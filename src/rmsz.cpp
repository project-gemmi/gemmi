// Copyright 2018 Global Phasing Ltd.

#include <stdio.h>
#include <cstdlib> // for getenv
#include <stdexcept>
#include "gemmi/model.hpp"     // for Structure, Atom, etc
#include "gemmi/chemcomp.hpp"  // for ChemComp
#include "gemmi/monlib.hpp"    // for MonLib, read_monomer_lib
#include "gemmi/topo.hpp"      // for Topo
#include "gemmi/calculate.hpp" // for find_best_plane, get_distance_from_plane
#include "gemmi/polyheur.hpp"  // for setup_entities
#include <gemmi/read_cif.hpp>  // for read_cif_gz
#include <gemmi/read_coor.hpp> // for read_structure_gz

#define GEMMI_PROG rmsz
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::Topo;

namespace {

enum OptionIndex { Monomers=4, FormatIn, Cutoff };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE"
    "\n\nValidate geometry of a coordinate file with (Refmac) monomer library."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library dir (default: $CLIBD_MON)." },
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  { Cutoff, 0, "", "cutoff", Arg::Float,
    "  --cutoff=ZC  \tList bonds and angles with Z score > ZC (default: 2)." },
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
};

double check_restraint(const Topo::Force force,
                       const Topo& topo,
                       double cutoff,
                       const char* tag,
                       RMSes* rmses,
                       bool verbose) {
  switch (force.rkind) {
    case Topo::RKind::Bond: {
      const Topo::Bond& t = topo.bonds[force.index];
      double z = t.calculate_z();
      if (z > cutoff) {
        int n = printf("%s bond %s: |Z|=%.1f", tag, t.restr->str().c_str(), z);
        if (verbose)
          printf(" %*.3f -> %.3f", std::max(50 - n, 7),
                 t.restr->value, t.calculate());
        puts("");
      }
      rmses->z_bond.put(z);
      rmses->d_bond.put(z * t.restr->esd);
      return z;
    }
    case Topo::RKind::Angle: {
      const Topo::Angle& t = topo.angles[force.index];
      double z = t.calculate_z();
      if (z > cutoff) {
        int n = printf("%s angle %s: |Z|=%.1f", tag, t.restr->str().c_str(), z);
        if (verbose)
          printf(" %*.1f -> %.1f", std::max(50 - n, 7),
                 t.restr->value, gemmi::deg(t.calculate()));
        puts("");
      }
      rmses->z_angle.put(z);
      rmses->d_angle.put(z * t.restr->esd);
      return z;
    }
    case Topo::RKind::Torsion: {
      const Topo::Torsion& t = topo.torsions[force.index];
      double z = t.calculate_z();  // takes into account t.restr->period
      if (z > cutoff) {
        int n = printf("%s torsion %s: |Z|=%.1f",
                       tag, t.restr->str().c_str(), z);
        if (verbose)
          printf(" %*.1f -> %.1f", std::max(50 - n, 7),
                 t.restr->value, gemmi::deg(t.calculate()));
        puts("");
      }
      rmses->z_torsion.put(z);
      rmses->d_torsion.put(z * t.restr->esd);
      return z;
    }
    case Topo::RKind::Chirality: {
      const Topo::Chirality& t = topo.chirs[force.index];
      rmses->all_chiralities++;
      if (!t.check()) {
        printf("%s wrong chirality of %s\n", tag, t.restr->str().c_str());
        rmses->wrong_chirality++;
        return 1.0;
      }
      return 0.0;
    }
    case Topo::RKind::Plane: {
      const Topo::Plane& t = topo.planes[force.index];
      auto coeff = find_best_plane(t.atoms);
      double max_z = 0;
      for (const gemmi::Atom* atom : t.atoms) {
        double dist = gemmi::get_distance_from_plane(atom->pos, coeff);
        double z = dist / t.restr->esd;
        if (z > cutoff)
          printf("%s atom %s not in plane %s, |Z|=%.1f\n", tag,
                 atom->name.c_str(), t.restr->str().c_str(), z);
        if (z > max_z)
            max_z = z;
      }
      rmses->z_plane.put(max_z);
      rmses->d_plane.put(max_z * t.restr->esd);
      return max_z;
    }
  }
  gemmi::unreachable();
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(1);
  const char* monomer_dir = p.options[Monomers] ? p.options[Monomers].arg
                                                : std::getenv("CLIBD_MON");
  if (monomer_dir == nullptr || *monomer_dir == '\0') {
    fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    return 1;
  }
  double cutoff = 2.0;
  if (p.options[Cutoff])
    cutoff = std::strtod(p.options[Cutoff].arg, nullptr);
  std::string input = p.coordinate_input_file(0);
  try {
    gemmi::CoorFormat format = coor_format_as_enum(p.options[FormatIn]);
    gemmi::Structure st = gemmi::read_structure_gz(input, format);
    if (st.input_format == gemmi::CoorFormat::Pdb ||
        st.input_format == gemmi::CoorFormat::ChemComp)
      gemmi::setup_entities(st);
    for (gemmi::Model& model : st.models) {
      if (st.models.size() > 1)
        printf("### Model %s ###\n", model.name.c_str());
      gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir,
                                                  model.get_all_residue_names(),
                                                  gemmi::read_cif_gz);
      Topo topo;
      topo.initialize_refmac_topology(st, model, monlib);
      topo.finalize_refmac_topology(monlib);

      RMSes rmses;
      // We could iterate directly over Topo::bonds, Topo::angles, etc,
      // but then we couldn't output the provenance (res or "link" below).
      for (const Topo::ChainInfo& chain_info : topo.chain_infos)
        for (const Topo::ResInfo& ri : chain_info.res_infos) {
          std::string res = chain_info.name + " " + ri.res->str();
          for (const Topo::Force& force : ri.forces)
            if (force.provenance == Topo::Provenance::PrevLink ||
                force.provenance == Topo::Provenance::Monomer)
              check_restraint(force, topo, cutoff, res.c_str(), &rmses,
                              p.options[Verbose]);
        }
      for (const Topo::ExtraLink& link : topo.extras) {
        for (const Topo::Force& force : link.forces)
          check_restraint(force, topo, cutoff, "link", &rmses,
                          p.options[Verbose]);
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
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (std::out_of_range& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
