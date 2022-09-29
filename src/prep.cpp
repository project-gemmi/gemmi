// Copyright 2017-2022 Global Phasing Ltd.

#include <stdio.h>             // for printf, fprintf
#include <cstdlib>             // for getenv
#include <iostream>            // for cerr
#include <stdexcept>           // for exception
#include "gemmi/crd.hpp"       // for prepare_crd, prepare_rst
#include "gemmi/to_cif.hpp"    // for write_cif_block_to_stream
#include "gemmi/fstream.hpp"   // for Ofstream
#include "gemmi/polyheur.hpp"  // for setup_entities
#include "gemmi/monlib.hpp"    // for MonLib, read_monomer_lib
#include "gemmi/read_cif.hpp"  // for read_cif_gz
#include "gemmi/read_coor.hpp" // for read_structure_gz
#include "gemmi/contact.hpp"   // for ContactSearch

#define GEMMI_PROG prep
#include "options.h"

namespace cif = gemmi::cif;

namespace {

enum OptionIndex {
  Split=4, Monomers, Libin, AutoCis, AutoLink,
  NoZeroOccRestr, NoHydrogens, KeepHydrogens
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n " EXE_NAME " --split [options] INPUT_FILE OUTPUT_BASENAME"
    "\n\nPrepare intermediate Refmac files."
    "\nINPUT_FILE can be in PDB, mmCIF or mmJSON format."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Split, 0, "", "split", Arg::None,
    "  --split  \tSplit output into two files: crd and rst." },
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { Libin, 0, "", "libin", Arg::Required,
    "  --libin=CIF  \tCustom additions to the monomer library." },
  { AutoCis, 0, "", "auto-cis", Arg::YesNo,
    "  --auto-cis=Y|N  \tAssign cis/trans ignoring CISPEP record (default: Y)." },
  { AutoLink, 0, "", "auto-link", Arg::YesNo,
    "  --auto-link=Y|N  \tFind links not included in LINK/SSBOND (default: N)." },
  //{ NoZeroOccRestr, 0, "", "no-zero-occ", Arg::None,
  //  "  --no-zero-occ  \tNo restraints for zero-occupancy atoms." },
  { NoOp, 0, "", "", Arg::None,
    "\nHydrogen options (default: remove and add on riding positions):" },
  { NoHydrogens, 0, "H", "no-hydrogens", Arg::None,
    "  -H, --no-hydrogens  \tRemove (and do not add) hydrogens." },
  { KeepHydrogens, 0, "", "keep-hydrogens", Arg::None,
    "  --keep-hydrogens  \tPreserve hydrogens from the input file." },
  { 0, 0, 0, 0, 0, 0 }
};

void assign_connections(gemmi::Model& model, gemmi::Structure& st) {
  using namespace gemmi;
  NeighborSearch ns(model, st.cell, 5.0);
  ns.populate();
  ContactSearch contacts(3.5f);
  contacts.ignore = ContactSearch::Ignore::AdjacentResidues;
  int counter = 0;
  contacts.for_each_contact(ns, [&](const CRA& cra1, const CRA& cra2,
                                    int image_idx, float dist_sq) {
    float r1 = cra1.atom->element.covalent_r();
    float r2 = cra2.atom->element.covalent_r();
    if (dist_sq > sq((r1 + r2) + 0.5))
      return;
    if (st.find_connection_by_cra(cra1, cra2))
      return;
    Connection conn;
    conn.name = "added" + std::to_string(++counter);
    conn.type = Connection::Covale;
    conn.asu = (image_idx == 0 ? Asu::Same : Asu::Different);
    conn.partner1 = make_address(*cra1.chain, *cra1.residue, *cra1.atom);
    conn.partner2 = make_address(*cra2.chain, *cra2.residue, *cra2.atom);
    conn.reported_distance = std::sqrt(dist_sq);
    printf("Added link %s - %s\n", atom_str(cra1).c_str(), atom_str(cra2).c_str());
    st.connections.push_back(conn);
  });
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  p.check_exclusive_pair(KeepHydrogens, NoHydrogens);
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  if (monomer_dir == nullptr || *monomer_dir == '\0') {
    fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    return 1;
  }
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  bool verbose = p.options[Verbose];

  try {
    if (verbose)
      printf("Reading %s ...\n", input.c_str());
    gemmi::Structure st = gemmi::read_structure_gz(input, gemmi::CoorFormat::Detect);
    gemmi::setup_entities(st);

    if (st.models.empty()) {
      fprintf(stderr, "No models found in the input file.\n");
      return 1;
    }
    gemmi::Model& model0 = st.models[0];

    std::string libin;
    if (p.options[Libin])
      libin = p.options[Libin].arg;
    if (verbose)
      printf("Reading monomer library...\n");
    gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir,
                                                   model0.get_all_residue_names(),
                                                   gemmi::read_cif_gz,
                                                   libin);
    if (p.is_yes(AutoCis, true))
      assign_cis_flags(model0);

    if (p.is_yes(AutoLink, false)) {
      assign_connections(model0, st);
    }

    if (verbose)
      printf("Preparing topology, hydrogens, restraints...\n");
    bool reorder = true;
    bool ignore_unknown_links = false;
    gemmi::HydrogenChange h_change;
    if (p.options[NoHydrogens])
      h_change = gemmi::HydrogenChange::Remove;
    else if (p.options[KeepHydrogens])
      h_change = gemmi::HydrogenChange::NoChange;
    else
      h_change = gemmi::HydrogenChange::ReAddButWater;
    auto topo = gemmi::prepare_topology(st, monlib, 0, h_change, reorder,
                                        &std::cerr, ignore_unknown_links);

    if (verbose)
      printf("Preparing data for Refmac...\n");
    cif::Document crd = prepare_refmac_crd(st, *topo, monlib, h_change);
    if (p.options[Split]) {
      output += ".crd";
      if (verbose)
        printf("Writing %s\n", output.c_str());
      gemmi::Ofstream os(output);
      write_cif_block_to_stream(os.ref(), crd.blocks.at(0), cif::Style::NoBlankLines);
      cif::Block rst = prepare_rst(*topo, monlib, st.cell);
      if (p.options[Split])
        output.replace(output.size()-3, 3, "rst");
      if (verbose)
        printf("Writing %s\n", output.c_str());
      if (p.options[Split])
        os = gemmi::Ofstream(output);
      else
        os->write("\n\n", 2);
      write_cif_block_to_stream(os.ref(), crd.blocks.at(1), cif::Style::NoBlankLines);
    } else {
      if (verbose)
        printf("Writing %s\n", output.c_str());
      gemmi::Ofstream os(output);
      write_cif_to_stream(os.ref(), crd, cif::Style::NoBlankLines);
    }
  } catch (std::exception& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
