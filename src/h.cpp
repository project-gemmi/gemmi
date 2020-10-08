// Copyright 2017 Global Phasing Ltd.

#include <cstdio>
#include <cstdlib>   // for getenv
#include <stdexcept>
#include <iostream>  // for cout
#include <gemmi/gzread.hpp>
#include "gemmi/chemcomp.hpp"  // for ChemComp
#include "gemmi/polyheur.hpp"  // for remove_hydrogens
#include "gemmi/to_cif.hpp"    // for write_cif_to_file
#include "gemmi/to_mmcif.hpp"  // for update_cif_block
#include "gemmi/to_pdb.hpp"    // for write_pdb
#include "gemmi/monlib.hpp"    // for MonLib, read_monomer_lib
#include "gemmi/topo.hpp"      // for Topo
#include "gemmi/fstream.hpp"   // for Ofstream
#include <gemmi/placeh.hpp>    // for place_hydrogens

#define GEMMI_PROG h
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Topo;
using gemmi::Restraints;

namespace {

enum OptionIndex { Monomers=4, RemoveH, KeepH, Water, Sort };

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

void remove_add_sort(Topo& topo, const std::vector<option::Option>& options) {
  int serial = 0;
  for (Topo::ChainInfo& chain_info : topo.chain_infos)
    for (Topo::ResInfo& ri : chain_info.res_infos) {
      const gemmi::ChemComp &cc = ri.chemcomp;
      gemmi::Residue &res = *ri.res;
      if (!options[KeepH]) {
        gemmi::remove_hydrogens(res);
        if (!options[RemoveH])
          if (options[Water] || !res.is_water())
            add_hydrogens(cc, res);
      }
      if (options[Sort]) {
        for (gemmi::Atom& atom : res.atoms) {
          auto it = cc.find_atom(atom.name);
          if (it == cc.atoms.end())
            gemmi::fail("No atom ", atom.name, " expected in ", res.name);
          atom.serial = int(it - cc.atoms.begin()); // temporary, for sorting only
        }
        std::sort(res.atoms.begin(), res.atoms.end(),
                  [](const gemmi::Atom& a, const gemmi::Atom& b) {
                    return a.serial != b.serial ? a.serial < b.serial
                                                : a.altloc < b.altloc;
                  });
      }
      if (!options[KeepH] || options[Sort])
        for (gemmi::Atom& atom : res.atoms)
          atom.serial = ++serial;
    }
}

int count_h(const gemmi::Structure& st) {
  int n = 0;
  for (const gemmi::Model& model : st.models)
    for (const gemmi::Chain& chain : model.chains)
      for (const gemmi::Residue& residue : chain.residues)
        for (const gemmi::Atom& atom : residue.atoms)
          if (atom.is_hydrogen())
            ++n;
  return n;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  if (monomer_dir == nullptr || *monomer_dir == '\0') {
    std::fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    return 1;
  }
  std::string input = p.coordinate_input_file(0);
  std::string output = p.nonOption(1);
  if (p.options[KeepH] && p.options[RemoveH])
    gemmi::fail("cannot use both --remove and --keep");
  if (p.options[Verbose])
    std::printf("Reading coordinates from %s\n", input.c_str());
  try {
    gemmi::Structure st = gemmi::read_structure_gz(input,
                                            gemmi::CoorFormat::UnknownAny);
    if (st.models.empty() || st.models[0].chains.empty()) {
      std::fprintf(stderr, "No atoms in the input file. Wrong format?\n");
      return 1;
    }
    int initial_h = 0;
    if (p.options[Verbose])
      initial_h = count_h(st);
    std::vector<std::string> res_names = st.models[0].get_all_residue_names();
    if (p.options[Verbose])
      std::printf("Reading %zu monomers and all links from %s\n",
                  res_names.size(), input.c_str());
    gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir, res_names,
                                                   gemmi::read_cif_gz);
    for (gemmi::Model& model : st.models) {
      Topo topo;
      topo.initialize_refmac_topology(st, model, monlib);
      remove_add_sort(topo, p.options);
      topo.finalize_refmac_topology(monlib);
      for (Topo::ChainInfo& chain_info : topo.chain_infos)
        for (Topo::ResInfo& ri : chain_info.res_infos)
          for (gemmi::Atom& atom : ri.res->atoms)
            if (!atom.is_hydrogen()) {
              try {
                place_hydrogens(atom, ri, topo);
              } catch (const std::runtime_error& e) {
                std::string loc = gemmi::atom_str(chain_info.name, *ri.res,
                                                  atom.name, atom.altloc);
                std::printf("Placing of hydrogen bonded to %s failed:\n  %s\n",
                            loc.c_str(), e.what());
              }
            }
    }
    if (p.options[Verbose])
      std::printf("Hydrogen site count: %d in input, %d in output.\n",
                  initial_h, count_h(st));
    if (p.options[Verbose])
      std::printf("Writing coordinates to %s\n", output.c_str());
    gemmi::Ofstream os(output, &std::cout);
    if (gemmi::coor_format_from_ext_gz(output) == gemmi::CoorFormat::Pdb)
      gemmi::write_pdb(st, os.ref());
    else
      cif::write_cif_to_stream(os.ref(), gemmi::make_mmcif_document(st),
                               cif::Style::PreferPairs);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (std::out_of_range& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include
