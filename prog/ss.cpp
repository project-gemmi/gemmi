// Copyright 2017 Global Phasing Ltd.

#include <cstdio>     // for printf, fprintf
#include <exception>  // for exception
#include <iostream>   // for cout, cerr
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/modify.hpp>    // for remove_hydrogens
#include <gemmi/to_cif.hpp>    // for write_cif_to_file
#include <gemmi/to_mmcif.hpp>  // for make_mmcif_document
#include <gemmi/to_pdb.hpp>    // for write_pdb
#include <gemmi/monlib.hpp>    // for MonLib
#include <gemmi/topo.hpp>      // for Topo, prepare_topology
#include <gemmi/fstream.hpp>   // for Ofstream
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/dssp.hpp>      // for
#include <gemmi/neighbor.hpp>  // for NeighborSearch
#include "monlib_opt.h"

#define GEMMI_PROG ss
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex {
  FormatIn=AfterMonLibOptions, Sort, Update, RemoveH, KeepH, Water, Unique, NoChange
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nDetermines secondary structure of proteins using DSSP method."
    "\nUses CCP4 monomer library to determine peptide bonds.\n"
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  MonLibUsage[0], // Monomers
  MonLibUsage[1], // Libin
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  // --angle, --existing-h, --dssp2 and other parameters
  { 0, 0, 0, 0, 0, 0 }
};

size_t count_hydrogens(const gemmi::Structure& st) {
  size_t count = 0;
  for (const gemmi::Model& model : st.models)
    for (gemmi::const_CRA cra : model.all())
      if (cra.atom->is_hydrogen())
        ++count;
  return count;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  p.check_exclusive_group({KeepH, RemoveH, Water, Unique});

  gemmi::HydrogenChange h_change = gemmi::HydrogenChange::ReAddButWater;
  if (p.options[Water])
    h_change = gemmi::HydrogenChange::ReAdd;
  else if (p.options[Unique])
    h_change = gemmi::HydrogenChange::ReAddKnown;
  else if (p.options[KeepH])
    h_change = gemmi::HydrogenChange::Shift;
  else if (p.options[NoChange])
    h_change = gemmi::HydrogenChange::NoChange;
  else if (p.options[RemoveH])
    h_change = gemmi::HydrogenChange::Remove;

  gemmi::CoorFormat input_format = coor_format_as_enum(p.options[FormatIn]);
  if (input_format == gemmi::CoorFormat::Unknown)
    input_format = gemmi::coor_format_from_ext_gz(input);
  gemmi::CoorFormat output_format = gemmi::coor_format_from_ext_gz(output);
  bool preserve_doc = (input_format == gemmi::CoorFormat::Mmcif &&
                       output_format == gemmi::CoorFormat::Mmcif);

  MonArguments mon_args;
  if (h_change != gemmi::HydrogenChange::Remove) {
    bool needs_monlib = input_format != gemmi::CoorFormat::ChemComp;
    mon_args = get_monomer_args(p.options, needs_monlib);
  }

   if (p.options[Verbose])
     std::printf("Reading coordinates from %s\n", input.c_str());
   try {
     gemmi::Structure st;
     std::unique_ptr<cif::Document> doc;
     if (preserve_doc || mon_args.has_plus())
       doc.reset(new cif::Document);
     st = gemmi::read_structure_gz(input, input_format, doc.get());
     if (st.models.empty() || st.models[0].chains.empty()) {
       std::fprintf(stderr, "No atoms in the input file. Wrong format?\n");
       return 1;
       gemmi::setup_entities(st);
       size_t initial_h = 0;
     if (p.options[Verbose])
       initial_h = count_hydrogens(st);
     if (h_change == gemmi::HydrogenChange::Remove)
       gemmi::remove_hydrogens(st);
     if (h_change != gemmi::HydrogenChange::Remove || p.options[Sort] || p.options[Update]) {
       gemmi::MonLib monlib;
       std::vector<std::string> wanted = st.models[0].get_all_residue_names();
       read_monomer_lib_and_user_files(monlib, wanted, mon_args, doc.get());
       if (!wanted.empty())
         gemmi::fail("Please create definitions for missing monomers.");
       if (p.options[Update])
         monlib.update_old_atom_names(st, {&gemmi::Logger::to_stdout});
       for (size_t i = 0; i != st.models.size(); ++i) {
         gemmi::Model& model = st.models[i];
         // preparing topology modifies hydrogens in the model
         std::unique_ptr<gemmi::Topo> topo =
           prepare_topology(st, monlib, i, h_change, p.options[Sort], {&gemmi::Logger::to_stderr});
         gemmi::NeighborSearch ns(model, st.cell, 9);
         //size_t idx = ns.grid.index_q(u, v, w);
         for (gemmi::Topo::ChainInfo& ci : topo->chain_infos)
           // TODO: add Calphas atoms to ns
           gemmi::dssp_determine_hydrogen_bonds(ns, ci);
           // TODO: find nearby Calphas for each peptide
           //ns.for_each(const Position &pos, char alt, double radius, const Func &func);
           // TODO: check if NH-CO makes hydrogen bond
     }
     if (p.options[Verbose]) {
       std::printf("Hydrogen site count: %zu in input, %zu in output.\n",
                   initial_h, count_hydrogens(st));
       std::printf("Writing coordinates to %s\n", output.c_str());
     }
     gemmi::Ofstream os(output, &std::cout);
     if (output_format == gemmi::CoorFormat::Pdb) {
       shorten_ccd_codes(st);
       gemmi::write_pdb(st, os.ref());
     } else {
       if (preserve_doc) {
         gemmi::MmcifOutputGroups groups(false);
         groups.atoms = true;
         gemmi::update_mmcif_block(st, doc->blocks[0], groups);
       } else {
         doc.reset(new cif::Document(gemmi::make_mmcif_document(st)));
       }
       cif::WriteOptions cif_options;
       cif_options.prefer_pairs = true;
       cif::write_cif_to_stream(os.ref(), *doc, cif_options);
       }
     }
     }
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

