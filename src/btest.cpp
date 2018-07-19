// Copyright 2018 Global Phasing Ltd.
//
// B-factor model testing

#include <gemmi/subcells.hpp>
#include <gemmi/elem.hpp>  // for is_hydrogen
#include <gemmi/math.hpp>  // for Correlation
#include <gemmi/resinfo.hpp>  // for find_tabulated_residue
#include <gemmi/gzread.hpp>
#define GEMMI_PROG btest
#include "options.h"
#include <stdio.h>
#include <cstdlib>  // for strtod
#include <algorithm>  // for sort

using namespace gemmi;

enum OptionIndex { Verbose=3, FromFile };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nB-factor model testing."},
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { FromFile, 0, "f", "file", Arg::Required,
    "  -f, --file=FILE  \tobtain file (or PDB ID) list from FILE" },
  { 0, 0, 0, 0, 0, 0 }
};

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
  int n;
  double b_mean;
  double cc;
  double rank_cc;
};

// for now checks only LDM
// from B. Halle (2002) http://www.pnas.org/content/99/3/1274
static Result test_bfactor_models(const Structure& st) {
  const float min_dist = 0.8f;
  const float max_dist = 15.0f;
  SubCells sc(st.models.at(0), st.cell, max_dist);
  const Model& model = st.models.at(0);
  std::vector<double> b_exper;
  std::vector<double> b_predict;
	for (const Chain& chain : model.chains) {
    for (const Residue& res : chain.residues) {
      if (!find_tabulated_residue(res.name).is_amino_acid())
        continue;
      for (const Atom& atom : res.atoms)
        if (!is_hydrogen(atom.element)) {
          double density = 0;
          sc.for_each(atom.pos, atom.altloc, max_dist,
                      [&](const SubCells::Mark& m, float dist_sq) {
              if (dist_sq > sq(min_dist) && !is_hydrogen(m.element)) {
                const_CRA cra = m.to_cra(model);
                ResidueInfo res_inf = find_tabulated_residue(cra.residue->name);
                if (res_inf.is_amino_acid()) {
                  double occ = cra.atom->occ;
                  //if (dist_sq < sq(7.35))
                  //  density += occ;
                  //density += occ * std::erfc(0.5 * (std::sqrt(dist_sq) - 7.0f)) / 2;
                  //density += occ * std::erfc(0.2 * (std::sqrt(dist_sq) - 4.0f)) / 2;
                  density += occ / dist_sq;
                  //density += occ / pow(dist_sq, 2.3/2.0);
                  //density += occ * std::exp(-std::sqrt(dist_sq) / 3.0f);
                  //density += occ * (1 - std::sqrt(dist_sq) / max_dist);
                  //density += occ / std::sqrt(dist_sq);
                }
              }
          });
          if (density == 0.0)
            continue;
          b_exper.push_back(atom.b_iso);
          b_predict.push_back(1 / density);
        }
    }
  }
  //normalize(b_exper);
  //normalize(b_predict);
  Correlation cc = calculate_correlation(b_exper, b_predict);
  Correlation rank_cc = calculate_correlation(get_ranks(b_exper),
                                              get_ranks(b_predict));
  Result r;
  r.b_mean = cc.mean_x;
  r.n = b_exper.size();
  r.cc = cc.coefficient();
  r.rank_cc = rank_cc.coefficient();
  return r;
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);

  std::vector<std::string> paths = p.paths_from_args_or_file(FromFile, 0, true);
  bool verbose = p.options[Verbose].count();
  double sum_cc = 0;
  double sum_rank_cc = 0;
  try {
    for (const std::string& path : paths) {
      if (verbose > 0)
        std::printf("File: %s\n", path.c_str());
      Structure st = read_structure_gz(path);
      Result r = test_bfactor_models(st);
      printf("%s <B>=%#.4g for %5d atoms   CC=%#.4g  rankCC=%#.4g\n",
             st.name.c_str(), r.b_mean, r.n, r.cc, r.rank_cc);
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
