// Copyright 2018 Global Phasing Ltd.
//
// B-factor model testing

// TODO: calculation in asu only

#include <gemmi/subcells.hpp>
#include <gemmi/elem.hpp>  // for is_hydrogen
#include <gemmi/math.hpp>  // for Correlation
#include <gemmi/resinfo.hpp>  // for find_tabulated_residue
#include <gemmi/polyheur.hpp> // for assign_subchains
#include <gemmi/gzread.hpp>
#include <gemmi/fileutil.hpp> // for expand_if_pdb_code
#define GEMMI_PROG bfit
#include "options.h"
#include <stdio.h>
#include <cstdlib>  // for strtod
#include <algorithm>  // for sort

using namespace gemmi;

enum OptionIndex { Verbose=3, FromFile, ListResidues, MinDist, MaxDist,
                   Exponent, ChainName, Sanity, SideChains, NoCrystal,
                   OmitEnds, Tsv };

struct BfitArg {
  static option::ArgStatus SideChains(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"include", "exclude", "only"});
  }
};
static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nB-factor model testing."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
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
  { ChainName, 0, "", "chain", Arg::Required,
    "  --chain=CHAIN  \tUse only one chain from the INPUT file." },
  { Sanity, 0, "", "sanity", Arg::None,
    "  --sanity  \tRun sanity checks first." },
  { SideChains, 0, "", "sidechains", BfitArg::SideChains,
    "  --sidechains=X  \tOne of: include, exclude, only (default: include)." },
  { NoCrystal, 0, "", "no-crystal", Arg::None,
    "  --no-crystal  \tIgnore crystal symmetry and intermolecular contacts." },
  { OmitEnds, 0, "", "omit-ends", Arg::Int,
    "  --omit-ends=N  \tIgnore N terminal residues from each chain end." },
  { Tsv, 0, "", "tsv", Arg::None,
    "  --tsv  \tOutput as tab-separated values." },
  { 0, 0, 0, 0, 0, 0 }
};

struct Params {
  float min_dist = 0.8f;
  float max_dist = 15.0f;
  float exponent = 2.0f;
  std::string chain_name;
  char sidechains = 'i';
  int omit_ends = 0;
};

static bool check_sanity(const Model& model) {
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
static std::vector<int> get_ranks(const std::vector<double>& data) {
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

void normalize(std::vector<double>& values) {
  Variance variance;
  for (double x : values)
    variance.add_point(x);
  double stddev = std::sqrt(variance.for_population());
  for (double& x : values)
    x = (x - variance.mean_x) / stddev;
}

struct Result {
  int n_residues;
  int n;
  double b_mean;
  double cc;
  double rank_cc;
};

static float calculate_weight(float dist_sq, const Params& params) {
  if (params.exponent == 2.0)  // canonical WCN
    return 1.0f / dist_sq;
  if (params.exponent == 0.0) // CN (a.k.a ACN)
    return 1.0f;
  return std::pow(dist_sq, -0.5f * params.exponent);
}

// CA,N,C,O
static bool is_protein_backbone(const std::string& name) {
  switch (name.size()) {
    case 1: return name[0] == 'N' || name[0] == 'C' || name[0] == 'O';
    case 2: return name[0] == 'C' && name[1] == 'A';
    default: return false;
  }
}

static Result test_bfactor_models(const Structure& st, const Params& params) {
  const Model& model = st.models.at(0);

  SubCells sc(model, st.cell, params.max_dist);
  // populate sc
  for (int n_ch = 0; n_ch != (int) model.chains.size(); ++n_ch) {
    const Chain& chain = model.chains[n_ch];
    for (int n_res = 0; n_res != (int) chain.residues.size(); ++n_res) {
      const Residue& res = chain.residues[n_res];
      if (find_tabulated_residue(res.name).is_buffer_or_water())
        continue;
      for (int n_atom = 0; n_atom != (int) res.atoms.size(); ++n_atom) {
        const Atom& atom = res.atoms[n_atom];
        if (!atom.is_hydrogen())
          sc.add_atom(atom, n_ch, n_res, n_atom);
      }
    }
  }

  std::vector<double> b_exper;
  std::vector<double> b_predict;
  int n_residues = 0;
  for (const Chain& chain : model.chains) {
    if (!params.chain_name.empty() && chain.name != params.chain_name)
      continue;
    for (const Residue& res : chain.get_polymer()) {
      // TODO: params.omit_ends
      ++n_residues;
      for (const Atom& atom : res.atoms) {
        if (is_hydrogen(atom.element))
          continue;
        if ((params.sidechains == 'e' && !is_protein_backbone(atom.name)) ||
            (params.sidechains == 'o' && is_protein_backbone(atom.name)))
          continue;
        double wcn = 0;
        sc.for_each(atom.pos, atom.altloc, params.max_dist,
                    [&](const SubCells::Mark& m, float dist_sq) {
            if (dist_sq > sq(params.min_dist)) {
              const_CRA cra = m.to_cra(model);
              wcn += calculate_weight(dist_sq, params) * cra.atom->occ;
            }
        });
        if (wcn == 0.0) {
          fprintf(stderr, "Warning: lonely atom %s %s %s\n",
                  chain.name.c_str(), res.str().c_str(), atom.name.c_str());
          continue;
        }
        b_exper.push_back(atom.b_iso);
        b_predict.push_back(1 / wcn);
      }
    }
  }

  Correlation cc = calculate_correlation(b_exper, b_predict);
  Correlation rank_cc = calculate_correlation(get_ranks(b_exper),
                                              get_ranks(b_predict));
  Result r;
  r.b_mean = cc.mean_x;
  r.n = b_exper.size();
  r.cc = cc.coefficient();
  r.rank_cc = rank_cc.coefficient();
  r.n_residues = n_residues;
  return r;
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  std::vector<std::string> paths = p.paths_from_args_or_file(FromFile, 0);

  bool verbose = p.options[Verbose].count();
  Params params;
  if (p.options[MinDist])
    params.min_dist = std::strtof(p.options[MinDist].arg, nullptr);
  if (p.options[MaxDist])
    params.max_dist = std::strtof(p.options[MaxDist].arg, nullptr);
  if (p.options[Exponent])
    params.exponent = std::strtof(p.options[Exponent].arg, nullptr);
  if (p.options[ChainName])
    params.chain_name = p.options[ChainName].arg;
  if (p.options[SideChains])
    params.sidechains = p.options[SideChains].arg[0];
  if (p.options[OmitEnds])
    params.omit_ends = std::atoi(p.options[OmitEnds].arg);
  double sum_cc = 0;
  double sum_rank_cc = 0;
  const char* format = p.options[Tsv]
               ? "%s\t%s\t%d\t%d\t%.2f\t%.4f\t%.4f\n"
               : "%s %s #res=%3d N=%5d  <B>=%5.2f   CC=%.4f  rankCC=%.4f\n";
  if (p.options[Tsv])
    printf("PDB\tChain\tN\t<B>\tCC\trankCC\n");
  try {
    for (std::string& path : paths) {
      if (verbose > 0)
        fprintf(stderr, "File: %s\n", path.c_str());
      if (p.options[FromFile] && !p.options[ChainName]) {
        size_t sep = path.find_first_of(" \t");
        params.chain_name = gemmi::trim_str(path.substr(sep));
      }
      if (p.options[FromFile] && starts_with_pdb_code(path))
        path.resize(4);
      Structure st = read_structure_gz(gemmi::expand_if_pdb_code(path));
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
      printf(format,
             st.name.c_str(),
             params.chain_name.empty() ? "*" : params.chain_name.c_str(),
             r.n_residues, r.n, r.b_mean,
             r.cc, r.rank_cc);
      sum_cc += r.cc;
      sum_rank_cc += r.rank_cc;
    }
    if (paths.size() > 1)
      printf("average of %4zu files             CC=%#.4g  rankCC=%#.4g\n",
             paths.size(), sum_cc / paths.size(), sum_rank_cc / paths.size());
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
