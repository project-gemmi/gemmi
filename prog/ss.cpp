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
  FormatIn=AfterMonLibOptions, Sort, Update, RemoveH, KeepH, Water, Unique, NoChange,
  CalculateH, HBondDefinition, Cutoff, PiHelixPreference,
  SearchPolyproline, HBondEnergyThreshold, MinCADistance, BendAngleMin
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
  { CalculateH, 0, "", "calculate-h", Arg::None,
    "  --calculate-h  \tCalculate hydrogen positions (original DSSP method). Default: use existing H atoms." },
  { HBondDefinition, 0, "", "hbond", Arg::Required,
    "  --hbond=TYPE  \tHydrogen bond definition: energy (DSSP) or geometry (distance+angle)." },
  { Cutoff, 0, "", "cutoff", Arg::Float,
    "  --cutoff=DIST  \tCutoff distance for neighbor search (default: 0.9 nm)." },
  { PiHelixPreference, 0, "", "pihelix", Arg::None,
    "  --pihelix  \tPrefer pi-helices over alpha-helices." },
  { SearchPolyproline, 0, "", "polypro", Arg::None,
    "  --polypro  \tSearch for polyproline helices." },
  { HBondEnergyThreshold, 0, "", "energy-cutoff", Arg::Float,
    "  --energy-cutoff=VAL  \tHydrogen bond energy cutoff (default: -0.5 kcal/mol)." },
  { MinCADistance, 0, "", "min-ca-dist", Arg::Float,
    "  --min-ca-dist=DIST  \tMinimum CA distance for hydrogen bonding (default: 9.0 A)." },
  { BendAngleMin, 0, "", "bend-angle", Arg::Float,
    "  --bend-angle=ANGLE  \tMinimum angle for bend assignment (default: 70.0 degrees)." },
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

  // Setup DSSP options
  gemmi::DsspOptions dssp_options;

  // Parse hydrogen mode
  if (p.options[CalculateH]) {
    dssp_options.hydrogen_mode = gemmi::HydrogenMode::Calculate;
  } else {
    dssp_options.hydrogen_mode = gemmi::HydrogenMode::Existing;
  }

  // Parse hydrogen bond definition
  if (p.options[HBondDefinition]) {
    const char* hbond = p.options[HBondDefinition].arg;
    if (gemmi::alpha_up(*hbond) == 'E') {
      dssp_options.hbond_definition = gemmi::HBondDefinition::Energy;
    } else if (gemmi::alpha_up(*hbond) == 'G') {
      dssp_options.hbond_definition = gemmi::HBondDefinition::Geometry;
    } else {
      std::fprintf(stderr, "Error: Invalid hydrogen bond definition '%s'. Use 'energy' or 'geometry'.\n", hbond);
      return 1;
    }
  }

  // Parse other DSSP options
  if (p.options[Cutoff])
    dssp_options.cutoff = std::stod(p.options[Cutoff].arg);
  if (p.options[PiHelixPreference])
    dssp_options.pi_helix_preference = true;
  if (p.options[SearchPolyproline])
    dssp_options.search_polyproline = true;
  if (p.options[HBondEnergyThreshold])
    dssp_options.hbond_energy_cutoff = std::stod(p.options[HBondEnergyThreshold].arg);
  if (p.options[MinCADistance])
    dssp_options.min_ca_distance = std::stod(p.options[MinCADistance].arg);
  if (p.options[BendAngleMin])
    dssp_options.bend_angle_min = std::stod(p.options[BendAngleMin].arg);

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
       if (p.options[Update])
         monlib.update_old_atom_names(st, {&gemmi::Logger::to_stdout});
       for (size_t i = 0; i != st.models.size(); ++i) {
         gemmi::Model& model = st.models[i];
         // preparing topology modifies hydrogens in the model
         std::unique_ptr<gemmi::Topo> topo =
           prepare_topology(st, monlib, i, h_change, p.options[Sort], {&gemmi::Logger::to_stderr});
         gemmi::NeighborSearch ns(model, st.cell, dssp_options.min_ca_distance);
         for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
           gemmi::Chain& chain = model.chains[n_ch];
           for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
             gemmi::Residue& res = chain.residues[n_res];
             if (gemmi::Atom* ca = res.find_atom("CA", '*', gemmi::El::C))
               ns.add_atom(*ca, n_ch, n_res, int(ca - &res.atoms[0]));
           }
         }

         // Calculate secondary structure using DSSP
         gemmi::DsspCalculator dssp_calc(dssp_options);

         for (gemmi::Topo::ChainInfo& ci : topo->chain_infos) {
           if (ci.polymer) {
             std::string ss_string = dssp_calc.calculate_secondary_structure(ns, *topo);
             if (p.options[Verbose]) {
               std::printf("Chain %s secondary structure: %s\n",
                          ci.chain_ref.name.c_str(), ss_string.c_str());
               if (p.options[Verbose].count() > 1) {
                 for (size_t j = 0; j < dssp_calc.res_infos.size(); ++j) {
                   gemmi::Topo::ResInfo* resinfo = dssp_calc.res_infos[j];
                   gemmi::Residue& res = *resinfo->res;
                   auto offset = [&](gemmi::Topo::ResInfo* ri) {
                     if (!ri)
                       return 0;
                     auto it = std::find(dssp_calc.res_infos.begin(), dssp_calc.res_infos.end(), ri);
                     if (it == dssp_calc.res_infos.end())
                       return 0;
                     return int(it - dssp_calc.res_infos.begin()) - int(j);
                   };
                   std::printf("# %zu %s %s  %d, %g   %d, %g    %d, %g    %d, %g\n",
                               j+1, res.seqid.str().c_str(), ci.chain_ref.name.c_str(),
                               offset(resinfo->acceptors[0]), resinfo->acceptor_energies[0],
                               offset(resinfo->donors[0]), resinfo->donor_energies[0],
                               offset(resinfo->acceptors[1]), resinfo->acceptor_energies[1],
                               offset(resinfo->donors[1]), resinfo->donor_energies[1]);
                 }
               }
             }

             // Store secondary structure in structure
             // TODO: Add secondary structure annotation to the structure
           }
         }
       }
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
  } catch (std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

