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
#include "monlib_opt.h"

#define GEMMI_PROG h
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
    "\n\nAdd hydrogens in positions specified by the monomer library."
    "\nBy default, it removes and re-adds all hydrogens."
    "\nBy default, hydrogens are not added to water."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  MonLibUsage[0], // Monomers
  MonLibUsage[1], // Libin
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tOrder atoms in residues according to _chem_comp_atom." },
  { Update, 0, "", "update", Arg::None,
    "  --update  \tIf deprecated atom names (from _chem_comp_atom.alt_atom_id)"
    " are used in the model, change them." },
  { NoOp, 0, "", "", Arg::None,
    "Hydrogen options, mutually exclusive. Default: add hydrogens, but not to water." },
  { Water, 0, "", "water", Arg::None,
    "  --water  \tAdd hydrogens also to water." },
  { Unique, 0, "", "unique", Arg::None,
    "  --unique  \tAdd only hydrogens with uniquely determined positions." },
  { KeepH, 0, "", "keep", Arg::None,
    "  --keep  \tDo not add/remove hydrogens, only change positions." },
  { NoChange, 0, "", "no-change", Arg::None,
    "  --no-change  \tDo not change hydrogens, not even positions." },
  { RemoveH, 0, "", "remove", Arg::None,
    "  --remove  \tOnly remove hydrogens." },
  MonLibUsage[2], // details about Libin (--lib)
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

  MonArguments mon_args;
  if (h_change != gemmi::HydrogenChange::Remove)
    mon_args = get_monomer_args(p.options);

  if (p.options[Verbose])
    std::printf("Reading coordinates from %s\n", input.c_str());
  gemmi::CoorFormat input_format = coor_format_as_enum(p.options[FormatIn]);
  if (input_format == gemmi::CoorFormat::Unknown)
    input_format = gemmi::coor_format_from_ext_gz(input);
  gemmi::CoorFormat output_format = gemmi::coor_format_from_ext_gz(output);
  bool preserve_doc = (input_format == gemmi::CoorFormat::Mmcif &&
                       output_format == gemmi::CoorFormat::Mmcif);
  try {
    gemmi::Structure st;
    std::unique_ptr<cif::Document> doc;
    if (preserve_doc)
      doc.reset(new cif::Document);
    st = gemmi::read_structure_gz(input, input_format, doc.get());
    if (st.models.empty() || st.models[0].chains.empty()) {
      std::fprintf(stderr, "No atoms in the input file. Wrong format?\n");
      return 1;
    }
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
      if (p.options[Update]) {
        std::string msg = monlib.update_old_atom_names(st);
        std::printf("%s", msg.c_str());
      }
      for (size_t i = 0; i != st.models.size(); ++i) {
        // preparing topology modifies hydrogens in the model
        prepare_topology(st, monlib, i, h_change, p.options[Sort], &std::cerr);
      }
    }
    if (p.options[Verbose])
      std::printf("Hydrogen site count: %zu in input, %zu in output.\n",
                  initial_h, count_hydrogens(st));
    if (p.options[Verbose])
      std::printf("Writing coordinates to %s\n", output.c_str());
    gemmi::Ofstream os(output, &std::cout);
    if (gemmi::coor_format_from_ext_gz(output) == gemmi::CoorFormat::Pdb) {
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
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
