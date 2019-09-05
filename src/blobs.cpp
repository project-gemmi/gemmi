// Copyright 2019 Global Phasing Ltd.

#include <stdio.h>
#include <algorithm>  // count
#include <stdexcept>
#include "gemmi/gzread.hpp"
#include "gemmi/polyheur.hpp"  // for remove_hydrogens
#include "gemmi/math.hpp"      // for Variance
#include "mapcoef.h"

#define GEMMI_PROG blobs
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { Sigma=AfterMapOptions, MaskRadius,
                   MinVolume, MinScore, MinPeak };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] PDB_OR_MMCIF MTZ_OR_MMCIF"
    "\n\nSearch for umodelled blobs of electron density."
    "\n\nOptions:" },
  MapUsage[Help],
  MapUsage[Version],
  MapUsage[Verbose],
  { MaskRadius, 0, "", "mask-radius", Arg::Float,
    "  --mask-radius=NUMBER  \tModel mask radius (default: 2.0 A)." },
  { Sigma, 0, "", "sigma", Arg::Float,
    "  --sigma=NUMBER  \tSigma level (default: 1.0)." },
  { NoOp, 0, "", "", Arg::None, "\nBlob criteria:" },
  { MinVolume, 0, "", "min-volume", Arg::Float,
    "  --min-volume=NUMBER  \tMinimum volume (default: 10.0 A^3)." },
  { MinScore, 0, "", "min-score", Arg::Float,
    "  --min-score=NUMBER  \tMinimum electron count (default: 50.0)." },
  { MinPeak, 0, "", "min-peak", Arg::Float,
    "  --min-peak=NUMBER  \tMinimum peak density (default: 0.0)." },
  { NoOp, 0, "", "", Arg::None, "\nOptions for map calculation:" },
  MapUsage[Diff],
  MapUsage[Section],
  MapUsage[FLabel],
  MapUsage[PhLabel],
  MapUsage[GridDims],
  MapUsage[Sample],
  MapUsage[GridQuery],
  { 0, 0, 0, 0, 0, 0 }
};


struct GridPos {
  int u, v, w;
  int idx;
};

struct Blob {
  std::vector<GridPos> points;
  gemmi::Position pos;
  float score = 0.0f;
  float max_value = 0.0f;
};

struct BlobCriteria {
  double min_volume = 10.0;
  double min_score = 50.0;
  double min_peak = 0.0;
  double cutoff;
};

inline bool finalize_blob(Blob& blob, const gemmi::Grid<float>& grid,
                          const BlobCriteria& criteria) {
  size_t size = blob.points.size();
  if (size < 3 ||
      size * grid.unit_cell.volume < criteria.min_volume * grid.point_count())
    return false;
  for (const GridPos& point : blob.points) {
    float value = grid.data[point.idx];
    blob.score += value;
    if (value > blob.max_value)
      blob.max_value = value;
  }
  // TODO: multiply score by volume per point?
  if (blob.score < criteria.min_score)
    return false;
  size_t sum[3] = {0, 0, 0};
  for (const GridPos& point : blob.points) {
    sum[0] += point.u;
    sum[1] += point.v;
    sum[2] += point.w;
  }
  gemmi::Fractional fract((double) sum[0] / (grid.nu * blob.points.size()),
                          (double) sum[1] / (grid.nv * blob.points.size()),
                          (double) sum[2] / (grid.nw * blob.points.size()));
  blob.pos = grid.unit_cell.orthogonalize(fract);
  return true;
}

std::vector<Blob> find_blobs_by_flood_fill(const gemmi::Grid<float>& grid,
                                           const BlobCriteria& criteria) {
  std::vector<Blob> blobs;
  std::array<std::array<int, 3>, 6> moves = {{{-1, 0, 0}, {1, 0, 0},
                                              {0 ,-1, 0}, {0, 1, 0},
                                              {0, 0, -1}, {0, 0, 1}}};
  std::vector<bool> visited(grid.data.size(), false);
  // TODO: iterate only over ASU
  int idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        assert(idx == grid.index_q(u, v, w));
        if (visited[idx])
          continue;
        float value = grid.data[idx];
        if (value < criteria.cutoff)
          continue;
        Blob blob;
        blob.points.push_back({u, v, w, idx});
        visited[idx] = true;
        for (size_t j = 0; j < blob.points.size()/*increasing!*/; ++j)
          for (const std::array<int, 3>& mv : moves) {
            GridPos nabe = { blob.points[j].u + mv[0],
                             blob.points[j].v + mv[1],
                             blob.points[j].w + mv[2],
                             0 };
            nabe.idx = grid.index_s(nabe.u, nabe.v, nabe.w);
            if (!visited[nabe.idx] && grid.data[nabe.idx] > criteria.cutoff) {
              visited[nabe.idx] = true;
              blob.points.push_back(nabe);
            }
          }
        if (finalize_blob(blob, grid, criteria))
          blobs.push_back(blob);
      }
  return blobs;
}


static int run(OptParser& p) {
  std::string model_path = p.coordinate_input_file(0);
  std::string sf_path = p.nonOption(1);

  // read model, remove hydrogens if any
  if (p.options[Verbose])
    printf("Reading coordinates from %s ...\n", model_path.c_str());
  gemmi::Structure st = gemmi::read_structure_gz(model_path);
  if (st.models.empty() || st.models[0].chains.empty()) {
    fprintf(stderr, "Not a coordinate file: %s\n", model_path.c_str());
    return 1;
  }
  gemmi::remove_hydrogens(st);
  gemmi::remove_waters(st);

  // read map (includes FFT)
  if (p.options[Verbose])
    printf("Reading reflections from %s ...\n", sf_path.c_str());
  gemmi::Grid<float> grid =
    read_sf_and_transform_to_map(sf_path.c_str(), p.options, true);
  if (p.options[Verbose]) {
    printf("Unit cell: %g A^3, grid points: %d, volume/point: %g A^3\n",
           grid.unit_cell.volume, grid.point_count(),
           grid.unit_cell.volume / grid.point_count());
  }

  // calculate map RMSD and setup blob criteria
  BlobCriteria criteria;
  gemmi::Variance grid_variance;
  for (float d : grid.data)
    grid_variance.add_point(d);
  double rmsd = std::sqrt(grid_variance.for_population());
  double sigma_level = 1.0;
  if (p.options[Sigma])
    sigma_level = std::strtod(p.options[Sigma].arg, nullptr);
  criteria.cutoff = sigma_level * rmsd;
  // search for blobs
  if (p.options[MinVolume])
    criteria.min_volume = std::strtod(p.options[MinVolume].arg, nullptr);
  if (p.options[MinScore])
    criteria.min_score = std::strtod(p.options[MinScore].arg, nullptr);
  if (p.options[MinPeak])
    criteria.min_peak = std::strtod(p.options[MinPeak].arg, nullptr);
  if (p.options[Verbose])
    printf("Map RMSD: %.3f, cut-off: %.3f e/A^3 (%g sigma)\n",
           rmsd, criteria.cutoff, sigma_level);

  // mask model by zeroing map values
  double radius = 2.0;
  if (p.options[MaskRadius])
    radius = std::strtod(p.options[MaskRadius].arg, nullptr);
  if (st.models.size() > 1)
    std::fprintf(stderr, "Note: only the first model is used.\n");
  for (const gemmi::Chain& chain : st.models[0].chains)
    for (const gemmi::Residue& res : chain.residues)
      for (const gemmi::Atom& atom : res.atoms)
        grid.set_points_around(atom.pos, radius, -INFINITY);
  grid.symmetrize([](float a, float b) { return std::min(a,b); });
  if (p.options[Verbose]) {
    int n = std::count(grid.data.begin(), grid.data.end(), -INFINITY);
    std::fprintf(stderr, "Masked points: %d of %d\n", n, grid.point_count());
  }

  // find and sort blobs
  std::vector<Blob> blobs = find_blobs_by_flood_fill(grid, criteria);
  std::sort(blobs.begin(), blobs.end(),
            [](const Blob& a, const Blob& b) { return a.score > b.score; });
  int n = 0;
  for (const Blob& b : blobs)
    printf("#%-4d %-3zu grid points, score %-8.4g  (%7.2f,%7.2f,%7.2f)\n",
           n++, b.points.size(), b.score, b.pos.x, b.pos.y, b.pos.z);

  // TODO:move blobs to symmetry image nearest to the model
  return 0;
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  //if (p.options[Sigma] && p.options[Absolute])
  //  gemmi::fail("cannot use both --sigma and --abs");
  try {
    return run(p);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
