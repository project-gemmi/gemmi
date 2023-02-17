// Copyright 2017 Global Phasing Ltd.

#include <cstdio>
#include <cstdlib>   // for getenv
#include <stdexcept>
#include <iostream>  // for cout
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/modify.hpp>    // for remove_hydrogens
#include <gemmi/to_cif.hpp>    // for write_cif_to_file
#include <gemmi/to_mmcif.hpp>  // for make_mmcif_document
#include <gemmi/to_pdb.hpp>    // for write_pdb
#include <gemmi/monlib.hpp>    // for MonLib, read_monomer_lib
#include <gemmi/topo.hpp>      // for Topo
#include <gemmi/fstream.hpp>   // for Ofstream
#include <gemmi/riding_h.hpp>  // for prepare_topology
#include <gemmi/read_cif.hpp>  // for read_cif_gz
#include <gemmi/mmread_gz.hpp> // for read_structure_gz

#define GEMMI_PROG h
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex { Monomers=4, FormatIn, RemoveH, KeepH, Water, Sort };

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
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  { RemoveH, 0, "", "remove", Arg::None,
    "  --remove  \tOnly remove hydrogens." },
  { KeepH, 0, "", "keep", Arg::None,
    "  --keep  \tDo not add/remove hydrogens, only change positions." },
  { Water, 0, "", "water", Arg::None,
    "  --water  \tAdd hydrogens also to waters." },
  { Sort, 0, "", "sort", Arg::None,
    "  --sort  \tOrder atoms in residues according to _chem_comp_atom." },
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
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  p.check_exclusive_pair(KeepH, RemoveH);

  gemmi::HydrogenChange h_change = gemmi::HydrogenChange::ReAddButWater;
  if (p.options[RemoveH])
    h_change = gemmi::HydrogenChange::Remove;
  else if (p.options[KeepH])
    h_change = gemmi::HydrogenChange::Shift;
  else if (p.options[Water])
    h_change = gemmi::HydrogenChange::ReAdd;

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
    if (h_change == gemmi::HydrogenChange::Remove) {
      gemmi::remove_hydrogens(st);
    } else {
      if (monomer_dir == nullptr || *monomer_dir == '\0') {
        std::fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
        return 1;
      }
      std::vector<std::string> res_names = st.models[0].get_all_residue_names();
      if (p.options[Verbose])
        std::printf("Reading %zu monomers and all links from %s\n",
                    res_names.size(), input.c_str());
      std::string libin;
      gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir, res_names,
                                                     gemmi::read_cif_gz, libin, true);
      for (size_t i = 0; i != st.models.size(); ++i)
        // preparing topology modifies hydrogens in the model
        prepare_topology(st, monlib, i, h_change, p.options[Sort], &std::cerr);
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
      cif::write_cif_to_stream(os.ref(), *doc, cif::Style::PreferPairs);
    }
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
