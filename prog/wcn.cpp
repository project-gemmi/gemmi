// Copyright 2018 Global Phasing Ltd.
//
// B-factor model testing

#include <gemmi/neighbor.hpp>
#include <gemmi/elem.hpp>  // for is_hydrogen
#include <gemmi/math.hpp>  // for Correlation
#include <gemmi/resinfo.hpp>  // for find_tabulated_residue
#include <gemmi/polyheur.hpp> // for assign_subchains
#include <gemmi/fileutil.hpp> // for file_open
#include <gemmi/pdb_id.hpp>   // for expand_if_pdb_code
#include <gemmi/mmread_gz.hpp> // for read_structure_gz
#define GEMMI_PROG wcn
#include "options.h"
#include <stdio.h>
#include <cstdlib>    // for strtof, atoi
#include <algorithm>  // for sort

using namespace gemmi;

namespace {

enum OptionIndex { FromFile=4, ListResidues, MinDist, MaxDist,
                   Exponent, Blur, Rom, ChainName, Sanity, SideChains,
                   NoCrystal, OmitEnds, PrintRes, XyOut };

struct WcnArg {
  static option::ArgStatus SideChains(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"include", "exclude", "only"});
  }
};
const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]\n"
    "Calculation of local density / contact numbers: WCN, CN, ACN, LDM, etc."},
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { FromFile, 0, "f", "file", Arg::Required,
    "  -f, --file=FILE  \tObtain paths or PDB IDs from FILE, one per line." },
  { ListResidues, 0, "l", "list", Arg::None,
    "  -l, --list  \tList per-residue values." },
  { MinDist, 0, "", "min-dist", Arg::Float,
    "  --min-dist=DIST  \tMinimum distance for \"contacts\" (default: 0.8)." },
  { MaxDist, 0, "", "cutoff", Arg::Float,
    "  --cutoff=DIST  \tMaximum distance for \"contacts\" (default: 15)." },
  { Exponent, 0, "", "pow", Arg::Float,
    "  --pow=P  \tExponent in the weighting (default: 2)." },
  { Blur, 0, "", "blur", Arg::Float,
    "  --blur=SIGMA  \tApply Gaussian smoothing of predicted B-factors." },
  { Rom, 0, "", "rom", Arg::None,
    "  --rom  \tRotation only model: |pos-ctr_of_chain|^P instead of WCN." },
  { ChainName, 0, "", "chain", Arg::Required,
    "  --chain=CHAIN  \tUse only one chain from the INPUT file." },
  { Sanity, 0, "", "sanity", Arg::None,
    "  --sanity  \tRun sanity checks first." },
  { SideChains, 0, "", "sidechains", WcnArg::SideChains,
    "  --sidechains=X  \tOne of: include, exclude, only (default: include)." },
  { NoCrystal, 0, "", "no-crystal", Arg::None,
    "  --no-crystal  \tIgnore crystal symmetry and intermolecular contacts." },
  { OmitEnds, 0, "", "omit-ends", Arg::Int,
    "  --omit-ends=N  \tIgnore N terminal residues from each chain end." },
  { PrintRes, 0, "", "print-res", Arg::None,
    "  --print-res  \tPrint also resolution and R-free." },
  { XyOut, 0, "", "xy-out", Arg::Required,
    "  --xy-out=DIR  \tWrite DIR/name.xy files with WCN and B(exper)." },
  { 0, 0, 0, 0, 0, 0 }
};

struct Params {
  double min_dist = 0.8;
  double max_dist = 15.0;
  double exponent = 2.0;
  double blur = 0.0;
  std::string chain_name;
  bool rotation_only = false;
  char sidechains = 'i';
  int omit_ends = 0;
  std::string xy_out;
};

Position calculate_center_of_mass(const ResidueSpan& residue_span) {
  double mass = 0;
  Vec3 sum;
  for (const Residue& res : residue_span)
    for (const Atom& atom : res.atoms) {
      double w = atom.element.weight() * atom.occ;
      sum += atom.pos * w;
      mass += w;
    }
  return Position(sum / mass);
}

bool check_sanity(const Model& model) {
  for (const Chain& chain : model.chains)
    for (const Residue& res : chain.residues)
      for (const Atom& atom : res.atoms) {
        if (atom.occ < 0 || atom.occ > 1) {
          fprintf(stderr, "WRONG: atom %s in %s has occupancy: %g\n",
                          atom.name.c_str(), res.str().c_str(), atom.occ);
          return false;
        }
        if (atom.b_iso < 0 ||
            (atom.b_iso == 0 && !atom.is_hydrogen() && atom.occ != 0)) {
          fprintf(stderr, "WRONG: atom %s in %s has B_iso: %g\n",
                          atom.name.c_str(), res.str().c_str(), atom.b_iso);
          return false;
        }
      }

  return true;
}

// ranks are from 1 to data.size()
std::vector<int> get_ranks(const std::vector<double>& data) {
  std::vector<int> indices(data.size());
  for (size_t i = 0; i != indices.size(); ++i)
    indices[i] = (int) i;
  std::sort(indices.begin(), indices.end(),
            [&data](int a, int b) { return data[a] < data[b]; });
  std::vector<int> result(data.size());
  for (size_t i = 0; i < indices.size(); ++i)
    result[indices[i]] = (int) i + 1;
  return result;
}

template<typename T>
Correlation calculate_correlation(const std::vector<T>& a,
                                  const std::vector<T>& b) {
  assert(a.size() == b.size());
  Correlation cc;
  for (size_t i = 0; i != a.size(); ++i)
    cc.add_point(a[i], b[i]);
  return cc;
}

struct Result {
  int n_residues;
  int n;
  double b_mean;
  double b_stddev;
  double cc;
  double rank_cc;
  double mean_abs_dev;
  double relative_mean_abs_dev;
};

double calculate_weight(double dist_sq, const Params& params) {
  if (params.exponent == 2.0)  // canonical WCN
    return 1.0 / dist_sq;
  if (params.exponent == 0.0) // CN (a.k.a ACN)
    return 1.0;
  return std::pow(dist_sq, -0.5f * params.exponent);
}

// CA,N,C,O
bool is_protein_backbone(const std::string& name) {
  switch (name.size()) {
    case 1: return name[0] == 'N' || name[0] == 'C' || name[0] == 'O';
    case 2: return name[0] == 'C' && name[1] == 'A';
    default: return false;
  }
}

Result test_bfactor_models(Structure& st, const Params& params) {
  Model& model = st.first_model();

  // prepare cell lists for neighbour search
  std::unique_ptr<NeighborSearch> ns;
  if (!params.rotation_only || params.blur != 0) {
    ns.reset(new NeighborSearch(model, st.cell, params.max_dist));
    for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
      const Chain& chain = model.chains[n_ch];
      for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
        const Residue& res = chain.residues[n_res];
        if (find_tabulated_residue(res.name).is_buffer_or_water())
          continue;
        for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
          const Atom& atom = res.atoms[n_atom];
          if (!atom.is_hydrogen())
            ns->add_atom(atom, n_ch, n_res, n_atom);
        }
      }
    }
  }

  // calculate B-factor predictor
  std::vector<double> b_exper;
  std::vector<double> b_predict;
  std::vector<const Atom*> atom_ptr;
  int n_residues = 0;
  for (Chain& chain : model.chains) {
    if (!params.chain_name.empty() && chain.name != params.chain_name)
      continue;
    ResidueSpan polymer = chain.get_polymer();
    if (polymer.size() <= 2 * (size_t) params.omit_ends)
      continue;
    Position com = calculate_center_of_mass(polymer);
    auto p_end = polymer.end() - params.omit_ends;
    for (auto res = polymer.begin() + params.omit_ends; res != p_end; ++res) {
      ++n_residues;
      for (Atom& atom : res->atoms) {
        if (is_hydrogen(atom.element) || atom.occ == 0.0f)
          continue;
        if ((params.sidechains == 'e' && !is_protein_backbone(atom.name)) ||
            (params.sidechains == 'o' && is_protein_backbone(atom.name)))
          continue;
        double r2 = atom.pos.dist_sq(com);
        double value;
        if (params.rotation_only) {
          if (params.exponent == 2)
            value = r2;
          else
            value = std::pow(r2, 0.5 * params.exponent);
        } else {
          double wcn = 0;
          ns->for_each(atom.pos, atom.altloc, params.max_dist,
                       [&](const NeighborSearch::Mark& m, double dist_sq) {
              if (dist_sq > sq(params.min_dist)) {
                CRA cra = m.to_cra(model);
                double weight = calculate_weight(dist_sq, params);
                // if an atom is one of multiple conformations we iterate here
                // only over other atoms of the same conformation (and atoms
                // with no altloc) so we don't weight by occupancy.
                if (atom.altloc == '\0')
                  weight *= cra.atom->occ;
                wcn += weight;
              }
          });
          if (wcn == 0.0) {
            fprintf(stderr, "Warning: lonely atom %s %s %s\n",
                    chain.name.c_str(), res->str().c_str(), atom.name.c_str());
            continue;
          }
          value = 1 / wcn;
        }
        b_exper.push_back(atom.b_iso);
        b_predict.push_back(value);
        // re-purposing u11, u22 and name
        atom.flag = 1;
        atom.aniso.u11 = (float) value;
        atom.aniso.u22 = (float) r2;
        std::string names = chain.name;
        names += '\t';
        names += res->seqid.str();
        names += '\t';
        names += atom.name;
        atom.name = std::move(names);
        atom_ptr.push_back(&atom);
      }
    }
  }

  // smoothing - average weighted by Gaussian(dist)
  if (params.blur > 0) {
    double mult = -0.5 / (params.blur * params.blur);
    for (size_t i = 0; i != atom_ptr.size(); ++i) {
      const Atom& atom = *atom_ptr[i];
      double b_sum = 0;
      double weight_sum = 0;
      ns->for_each(atom.pos, atom.altloc, 3 * params.blur,
                   [&](const NeighborSearch::Mark& m, double dist_sq) {
          const_CRA cra = m.to_cra(model);
          if (cra.atom->flag) {
            double weight = std::exp(mult * dist_sq);
            if (atom.altloc == '\0' && cra.atom != &atom)
              weight *= cra.atom->occ;
            b_sum += weight * cra.atom->aniso.u11;
            weight_sum += weight;
          }
      });
      b_predict[i] = b_sum / weight_sum;
    }
  }

  // optionally, write xy file
  if (!params.xy_out.empty()) {
    std::string path = params.xy_out + "/" + st.name;
    if (!params.chain_name.empty())
      path += "-" + params.chain_name;
    path += ".xy";
    auto f = gemmi::file_open(path.c_str(), "w");
    fprintf(f.get(), "WCN\tB\tr2\tchain\tseqid\tatom\telem\n");
    for (size_t i = 0; i != b_predict.size(); ++i)
      fprintf(f.get(), "%g\t%g\t%.2f\t%s\t%s\n",
              b_predict[i], b_exper[i],
              atom_ptr[i]->aniso.u22, // squared distance to the center of mass
              // Atom::name was used here to store also chain and seqid
              atom_ptr[i]->name.c_str(),
              atom_ptr[i]->element.name());
  }

  // calculate statistics
  Correlation cc = calculate_correlation(b_predict, b_exper);
  Correlation rank_cc = calculate_correlation(get_ranks(b_predict),
                                              get_ranks(b_exper));
  Result r;
  r.b_mean = cc.mean_y;
  r.b_stddev = std::sqrt(cc.y_variance());
  r.n = (int) b_exper.size();
  r.cc = cc.coefficient();
  r.rank_cc = rank_cc.coefficient();
  {
    double slope = cc.slope();
    double intercept = cc.intercept();
    double a = 0, b = 0;
    for (size_t i = 0; i != b_exper.size(); ++i) {
      double diff = b_predict[i] * slope + intercept - b_exper[i];
      a += std::abs(diff);
      b += std::abs(b_exper[i] - cc.mean_y);
      //b_predict[i] = diff; // re-purposing b_predict
    }
    r.mean_abs_dev = a / b_exper.size();
    r.relative_mean_abs_dev = a / b;
  }
  r.n_residues = n_residues;
  return r;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  std::vector<std::string> paths = p.paths_from_args_or_file(FromFile, 0);

  int verbose = p.options[Verbose].count();
  Params params;
  if (p.options[MinDist])
    params.min_dist = std::strtof(p.options[MinDist].arg, nullptr);
  if (p.options[MaxDist])
    params.max_dist = std::strtof(p.options[MaxDist].arg, nullptr);
  if (p.options[Exponent])
    params.exponent = std::strtof(p.options[Exponent].arg, nullptr);
  if (p.options[Blur])
    params.blur = std::strtof(p.options[Blur].arg, nullptr);
  if (p.options[Rom])
    params.rotation_only = true;
  if (p.options[ChainName])
    params.chain_name = p.options[ChainName].arg;
  if (p.options[SideChains])
    params.sidechains = p.options[SideChains].arg[0];
  if (p.options[OmitEnds])
    params.omit_ends = std::max(std::atoi(p.options[OmitEnds].arg), 0);
  if (p.options[XyOut])
    params.xy_out = p.options[XyOut].arg;
  double sum_cc = 0;
  double sum_rmad = 0;
  double sum_rank_cc = 0;
  printf("PDB\tChain\t");
  if (p.options[PrintRes])
    printf("Res[A]\tRFree\t");
  printf("#res\tN\t<B>\tstd(B)\tCC\t1-RMAD\trankCC\n");
  try {
    for (std::string& path : paths) {
      if (verbose > 0)
        fprintf(stderr, "File: %s\n", path.c_str());
      if (p.options[FromFile] && !p.options[ChainName]) {
        size_t sep = path.find_first_of(" \t");
        if (sep != std::string::npos) {
          params.chain_name = gemmi::trim_str(path.substr(sep));
          path.resize(sep);
        }
      }
      Structure st = read_structure_gz(gemmi::expand_if_pdb_code(path));
      st.merge_chain_parts();
      if (p.options[NoCrystal])
        st.cell = UnitCell();
      if (p.options[Sanity]) {
        if (!check_sanity(st.models.at(0))) {
          fprintf(stderr, "Skipping %s\n", path.c_str());
          continue;
        }
      }
      gemmi::assign_subchains(st, false);
      Result r = test_bfactor_models(st, params);
      printf("%s\t%s\t",
             st.name.c_str(),
             params.chain_name.empty() ? "*" : params.chain_name.c_str());
      if (p.options[PrintRes]) {
        double rfree = 0;
        if (st.meta.refinement.size() > 0)
          rfree = st.meta.refinement[0].r_free;
        printf("%.2f\t%.2f\t", st.resolution, rfree);
      }
      printf("%d\t%d\t%.2f\t%.1f\t%.4f\t%.4f\t%.4f\n",
             r.n_residues, r.n, r.b_mean, r.b_stddev,
             r.cc, 1.0 - r.relative_mean_abs_dev, r.rank_cc);
      sum_cc += r.cc;
      sum_rmad += r.relative_mean_abs_dev;
      sum_rank_cc += r.rank_cc;
    }
    int N = (int) paths.size();
    if (N > 1)
      fprintf(stderr,
              "average of %4d files    CC=%#.4g  1-RMAD=%#.4g  rankCC=%#.4g\n",
              N, sum_cc / N, 1.0 - sum_rmad / N, sum_rank_cc / N);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}
