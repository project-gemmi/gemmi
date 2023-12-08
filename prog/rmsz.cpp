// Copyright 2018 Global Phasing Ltd.

#include <stdio.h>   // for printf, fprintf, stderr
#include <cmath>     // for sqrt
#include <iostream>  // for cerr
#include <exception>
#include "gemmi/model.hpp"     // for Structure, Atom, etc
#include "gemmi/chemcomp.hpp"  // for ChemComp
#include "gemmi/monlib.hpp"    // for MonLib
#include "gemmi/topo.hpp"      // for Topo
#include "gemmi/calculate.hpp" // for find_best_plane, get_distance_from_plane
#include "gemmi/polyheur.hpp"  // for setup_entities
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#include <gemmi/sprintf.hpp>   // for snprintf_z
#include "monlib_opt.h"

#define GEMMI_PROG rmsz
#include "options.h"

using gemmi::Topo;

namespace {

enum OptionIndex { Quiet=AfterMonLibOptions, FormatIn, Cutoff, Sort, Missing };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] [FILE]..."
    "\n\nValidate geometry of a coordinate file with (Refmac) monomer library."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Quiet, 0, "q", "quiet", Arg::None,
    "  -q, --quiet  \tShow only summary." },
  MonLibUsage[0], // Monomers
  MonLibUsage[1], // Libin
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  { Cutoff, 0, "", "cutoff", Arg::Float,
    "  --cutoff=ZC  \tList bonds and angles with Z score > ZC (default: 2)." },
  { Sort, 0, "s", "sort", Arg::None,
    "  -s, --sort  \tSort output according to |Z|." },
  { Missing, 0, "", "missing", Arg::None,
    "  --missing  \tList missing atoms." },
  MonLibUsage[2], // details about Libin (--lib)
  { 0, 0, 0, 0, 0, 0 }
};

struct RMS {
  int n = 0;
  double sum_sq = 0.;
  void put(double x) { ++n; sum_sq += x * x; }
  double get_value() const { return std::sqrt(sum_sq / n); }
};

struct RMSes {
  RMS d_bond;
  RMS d_angle;
  RMS d_torsion;
  RMS d_plane;
  RMS z_bond;
  RMS z_angle;
  RMS z_torsion;
  RMS z_plane;
  int wrong_chirality = 0;
  int all_chiralities = 0;
  int wrong_bond = 0;
  int wrong_angle = 0;
  int wrong_torsion = 0;
  int wrong_plane = 0;
};

void check_restraint(const Topo::Rule rule,
                     const Topo& topo,
                     const gemmi::UnitCell& cell,
                     double cutoff,
                     const char* tag,
                     RMSes* rmses,
                     int verbosity,
                     std::multimap<double, std::string>* lines) {
  double z = 0;
  char buf[200];
  size_t pos = 0;
  #define PRINT(...) do { \
    pos = std::min(pos, (size_t)200); \
    pos += gemmi::snprintf_z(buf + pos, 200 - pos, __VA_ARGS__); \
  } while(0)
  auto end_line = [&]() {
    if (lines)
      lines->emplace(-z, buf);
    else
      printf("%s\n", buf);
    pos = 0;
  };
  switch (rule.rkind) {
    case Topo::RKind::Bond: {
      const Topo::Bond& t = topo.bonds[rule.index];
      double dist = t.calculate();
      if (t.asu == gemmi::Asu::Different)
        dist = cell.find_nearest_image(t.atoms[0]->pos, t.atoms[1]->pos, t.asu).dist();
      z = t.calculate_z_(dist);
      if (z > cutoff) {
        rmses->wrong_bond++;
        if (verbosity >= 0) {
          PRINT("%s bond %s: |Z|=%.1f", tag, t.restr->str().c_str(), z);
          if (verbosity > 0) {
            PRINT(verbosity < 2 ? " %*.3f -> %.3f" : " %*g -> %g",
                  std::max(50 - (int)pos, 7), t.restr->value, t.calculate());
            if (verbosity > 2)
              PRINT("%*s%g", std::max(70 - (int)pos, 5), "esd=", t.restr->esd);
          }
          end_line();
        }
      }
      rmses->z_bond.put(z);
      rmses->d_bond.put(z * t.restr->esd);
      break;
    }
    case Topo::RKind::Angle: {
      const Topo::Angle& t = topo.angles[rule.index];
      z = t.calculate_z();
      if (z > cutoff) {
        rmses->wrong_angle++;
        if (verbosity >= 0) {
          PRINT("%s angle %s: |Z|=%.1f", tag, t.restr->str().c_str(), z);
          if (verbosity > 0) {
            PRINT(verbosity < 2 ? " %*.1f -> %.1f" : " %*g -> %g",
                  std::max(50 - (int)pos, 7), t.restr->value, gemmi::deg(t.calculate()));
            if (verbosity > 2)
              PRINT("%*s%g", std::max(70 - (int)pos, 5), "esd=", t.restr->esd);
          }
          end_line();
        }
      }
      rmses->z_angle.put(z);
      rmses->d_angle.put(z * t.restr->esd);
      break;
    }
    case Topo::RKind::Torsion: {
      const Topo::Torsion& t = topo.torsions[rule.index];
      // if _chem_comp_tor.value_angle_esd is 0 we skip it silently
      // (Refmac also ignores such items)
      if (t.restr->esd == 0)
        return;
      z = t.calculate_z();  // takes into account t.restr->period
      if (z > cutoff) {
        rmses->wrong_torsion++;
        if (verbosity >= 0) {
          PRINT("%s torsion %s: |Z|=%.1f", tag, t.restr->str().c_str(), z);
          if (verbosity > 0) {
            PRINT(verbosity < 2 ? " %*.1f -> %.1f" : " %*g -> %g",
                  std::max(50 - (int)pos, 7), t.restr->value, gemmi::deg(t.calculate()));
            if (verbosity > 2)
              PRINT("%*s%g p=%d", std::max(70 - (int)pos, 5), "esd=",
                    t.restr->esd, t.restr->period);
          }
          end_line();
        }
      }
      rmses->z_torsion.put(z);
      rmses->d_torsion.put(z * t.restr->esd);
      break;
    }
    case Topo::RKind::Chirality: {
      const Topo::Chirality& t = topo.chirs[rule.index];
      rmses->all_chiralities++;
      if (!t.check()) {
        if (verbosity >= 0) {
          PRINT("%s wrong chirality of %s\n", tag, t.restr->str().c_str());
          end_line();
        }
        rmses->wrong_chirality++;
      }
      break;
    }
    case Topo::RKind::Plane: {
      const Topo::Plane& t = topo.planes[rule.index];
      auto coeff = find_best_plane(t.atoms);
      double max_z = 0;
      for (const gemmi::Atom* atom : t.atoms) {
        double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
        z = dist / t.restr->esd;
        if (verbosity >= 0)
          if (z > cutoff) {
            PRINT("%s atom %s not in plane %s  |Z|=%.1f", tag,
                  atom->name.c_str(), t.restr->str().c_str(), z);
            if (verbosity > 0) {
              PRINT((verbosity < 2 ? "  d=%.1f" : "  d=%g"), dist);
              if (verbosity > 2)
                PRINT("%*s%g", std::max(70 - (int)pos, 5), "esd=", t.restr->esd);
            }
            end_line();
          }
        if (z > max_z)
            max_z = z;
      }
      rmses->z_plane.put(max_z);
      rmses->d_plane.put(max_z * t.restr->esd);
      if (max_z > cutoff)
        rmses->wrong_plane++;
      break;
    }
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  MonArguments mon_args = get_monomer_args(p.options);
  double cutoff = 2.0;
  if (p.options[Cutoff])
    cutoff = std::strtod(p.options[Cutoff].arg, nullptr);
  int verbosity = p.options[Verbose].count() - p.options[Quiet].count();
  for (int i = 0; i < p.nonOptionsCount(); ++i) {
    std::string input = p.coordinate_input_file(i);
    printf("File: %s\n", input.c_str());
    gemmi::cif::Document st_doc;
    try {
      gemmi::CoorFormat format = coor_format_as_enum(p.options[FormatIn]);
      gemmi::Structure st = gemmi::read_structure_gz(input, format, &st_doc);
      if (st.models.empty() || st.models[0].chains.empty()) {
        fprintf(stderr, "No atoms in the input file. Wrong format?\n");
        return 1;
      }
      gemmi::setup_entities(st);

      gemmi::MonLib monlib;
      std::vector<std::string> wanted = st.models[0].get_all_residue_names();
      read_monomer_lib_and_user_files(monlib, wanted, mon_args, &st_doc);
      if (!wanted.empty())
        gemmi::fail("Please create definitions for missing monomers.");

      for (gemmi::Model& model : st.models) {
        if (st.models.size() > 1)
          printf("### Model %s ###\n", model.name.c_str());
        Topo topo;
        topo.warnings = &std::cerr;
        topo.initialize_refmac_topology(st, model, monlib);
        topo.apply_all_restraints(monlib);

        if (p.options[Missing]) {
          std::vector<gemmi::AtomAddress> vec = find_missing_atoms(topo);
          printf("%zu missing atoms.\n", vec.size());
          for (const gemmi::AtomAddress& aa : vec)
            printf("    %s\n", aa.str().c_str());
          continue;
        }

        RMSes rmses;
        std::multimap<double, std::string> line_storage;
        std::multimap<double, std::string>* lines = nullptr;
        if (p.options[Sort])
          lines = &line_storage;

        // We could iterate directly over Topo::bonds, Topo::angles, etc,
        // but then we couldn't output the provenance (res or "link" below).
        for (const Topo::ChainInfo& chain_info : topo.chain_infos)
          for (const Topo::ResInfo& ri : chain_info.res_infos) {
            for (const Topo::Link& prev : ri.prev) {
              assert(ri.res == prev.res2);
              std::string rtag = gemmi::cat(chain_info.chain_ref.name, ' ',
                                            prev.res1->str(), '-', ri.res->str());
              for (const Topo::Rule& rule : prev.link_rules)
                check_restraint(rule, topo, st.cell, cutoff, rtag.c_str(),
                                &rmses, verbosity, lines);
            }
            std::string rtag = chain_info.chain_ref.name + " ";
            if (st.input_format != gemmi::CoorFormat::ChemComp)
              rtag += ri.res->seqid.str();
            gemmi::cat_to(rtag, '(', ri.res->name,  ')');
            for (const Topo::Rule& rule : ri.monomer_rules)
              check_restraint(rule, topo, st.cell, cutoff, rtag.c_str(),
                              &rmses, verbosity, lines);
          }
        for (const Topo::Link& link : topo.extras)
          for (const Topo::Rule& rule : link.link_rules)
            check_restraint(rule, topo, st.cell, cutoff, "link",
                            &rmses, verbosity, lines);

        if (lines) {
          for (const auto& t : *lines)
            printf("%s\n", t.second.c_str());
        }
        printf("Model rmsZ: "
               "bond: %.3f, angle: %.3f, torsion: %.3f, planarity %.3f\n"
               "Model rmsD: "
               "bond: %.3f, angle: %.3f, torsion: %.3f, planarity %.3f\n"
               "wrong chirality: %d of %d\n",
               rmses.z_bond.get_value(),
               rmses.z_angle.get_value(),
               rmses.z_torsion.get_value(),
               rmses.z_plane.get_value(),
               rmses.d_bond.get_value(),
               rmses.d_angle.get_value(),
               rmses.d_torsion.get_value(),
               rmses.d_plane.get_value(),
               rmses.wrong_chirality, rmses.all_chiralities);
        printf("rmsZ > %g for:\n"
               "  %d of %d bonds,\n"
               "  %d of %d angles,\n"
               "  %d of %d torsion angles,\n"
               "  %d of %d planes.\n",
               cutoff,
               rmses.wrong_bond, rmses.z_bond.n,
               rmses.wrong_angle, rmses.z_angle.n,
               rmses.wrong_torsion, rmses.z_torsion.n,
               rmses.wrong_plane, rmses.z_plane.n);
      }
    } catch (std::exception& e) {
      fprintf(stderr, "ERROR: %s\n", e.what());
      return 1;
    }
  }
  return 0;
}
