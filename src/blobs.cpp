// Copyright 2019 Global Phasing Ltd.

#include <stdio.h>
#include <algorithm>  // count
#include <stdexcept>
#include "gemmi/gzread.hpp"
#include "gemmi/polyheur.hpp"  // for remove_hydrogens
#include "gemmi/math.hpp"      // for Variance
#include "gemmi/subcells.hpp"  // for SubCells
#include "mapcoef.h"

#define GEMMI_PROG blobs
#include "options.h"

namespace cif = gemmi::cif;

enum OptionIndex { SigmaCutoff=AfterMapOptions, AbsCutoff,
                   MaskRadius, MaskWater,
                   MinVolume, MinScore, MinPeak };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] MTZ_OR_MMCIF PDB_OR_MMCIF"
    "\n\nSearch for umodelled blobs of electron density."
    "\n\nOptions:" },
  MapUsage[Help],
  MapUsage[Version],
  MapUsage[Verbose],

  { NoOp, 0, "", "", Arg::None,
    "\nThe area around model is masked to search only unmodelled density." },
  { MaskRadius, 0, "", "mask-radius", Arg::Float,
    "  --mask-radius=NUMBER  \tMask radius (default: 2.0 A)." },
  { MaskWater, 0, "", "mask-waters", Arg::None,
    "  --mask-water  \tMask water (water is not masked by default)." },

  { NoOp, 0, "", "", Arg::None,
    "\nThe searched blobs must density above:" },
  { SigmaCutoff, 0, "", "sigma", Arg::Float,
    "  --sigma=NUMBER  \tSigma (RMSD) level (default: 1.0)." },
  { AbsCutoff, 0, "", "abs", Arg::Float,
    "  --abs=NUMBER  \tAbsolute level in electrons/A^3." },

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
  // the mask will be used as follows:
  // 1=in blob,  0=in asu, not in blob (so far),  -1=in neither
  std::vector<signed char> mask = grid.get_asu_mask<signed char>(0, -1);
  std::vector<gemmi::Op> ops = grid.get_scaled_ops_except_id();
  int idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        assert(idx == grid.index_q(u, v, w));
        if (mask[idx] != 0)
          continue;
        float value = grid.data[idx];
        if (value < criteria.cutoff)
          continue;
        Blob blob;
        blob.points.push_back({u, v, w, idx});
        mask[idx] = 1;
        for (size_t j = 0; j < blob.points.size()/*increasing!*/; ++j)
          for (const std::array<int, 3>& mv : moves) {
            GridPos nabe = { blob.points[j].u + mv[0],
                             blob.points[j].v + mv[1],
                             blob.points[j].w + mv[2],
                             0 };
            nabe.idx = grid.index_s(nabe.u, nabe.v, nabe.w);
            if (mask[nabe.idx] != 1 && grid.data[nabe.idx] > criteria.cutoff) {
              if (mask[nabe.idx] != 0)
                for (const gemmi::Op& op : ops) {
                  auto t = grid.transformed_uvw(op, nabe.u, nabe.v, nabe.w);
                  int mate_idx = grid.index_s(t[0], t[1], t[2]);
                  if (mask[mate_idx] == 0)
                    mask[mate_idx] = -1;
                }
              mask[nabe.idx] = 1;
              blob.points.push_back(nabe);
            }
          }
        if (finalize_blob(blob, grid, criteria))
          blobs.push_back(blob);
      }
  return blobs;
}


static int run(OptParser& p) {
  std::string sf_path = p.nonOption(0);
  std::string model_path = p.coordinate_input_file(1);

  // read model, remove hydrogens if any
  if (p.options[Verbose])
    printf("Reading coordinates from %s ...\n", model_path.c_str());
  gemmi::Structure st = gemmi::read_structure_gz(model_path);
  if (st.models.empty() || st.models[0].chains.empty()) {
    fprintf(stderr, "Not a coordinate file: %s\n", model_path.c_str());
    return 1;
  }
  if (st.models.size() > 1)
    std::fprintf(stderr, "Note: only the first model is used.\n");
  gemmi::Model& model = st.models[0];
  gemmi::remove_hydrogens(model);
  if (!p.options[MaskWater])
    gemmi::remove_waters(model);

  // read map (includes FFT)
  gemmi::Grid<float> grid =
    read_sf_and_fft_to_map(sf_path.c_str(), p.options,
                           p.options[Verbose] ? stdout : nullptr,
                           true);
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
  if (p.options[AbsCutoff]) {
    criteria.cutoff = std::strtod(p.options[SigmaCutoff].arg, nullptr);
    sigma_level = criteria.cutoff / rmsd;
  } else {
    if (p.options[SigmaCutoff])
      sigma_level = std::strtod(p.options[SigmaCutoff].arg, nullptr);
    criteria.cutoff = sigma_level * rmsd;
  }
  // search for blobs
  if (p.options[MinVolume])
    criteria.min_volume = std::strtod(p.options[MinVolume].arg, nullptr);
  if (p.options[MinScore])
    criteria.min_score = std::strtod(p.options[MinScore].arg, nullptr);
  if (p.options[MinPeak])
    criteria.min_peak = std::strtod(p.options[MinPeak].arg, nullptr);
  if (p.options[Verbose])
    printf("Map RMSD: %.3f, cut-off: %.3f e/A^3 (%.3f sigma)\n",
           rmsd, criteria.cutoff, sigma_level);

  // mask model by zeroing map values
  double radius = 2.0;
  if (p.options[MaskRadius])
    radius = std::strtod(p.options[MaskRadius].arg, nullptr);
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& res : chain.residues)
      for (const gemmi::Atom& atom : res.atoms)
        grid.set_points_around(atom.pos, radius, -INFINITY);
  grid.symmetrize_min();
  if (p.options[Verbose]) {
    int n = std::count(grid.data.begin(), grid.data.end(), -INFINITY);
    printf("Masked points: %d of %d\n", n, grid.point_count());
  }

  // find and sort blobs
  std::vector<Blob> blobs = find_blobs_by_flood_fill(grid, criteria);
  if (p.options[Verbose])
    printf("%zu blob%s found.\n", blobs.size(), blobs.size() == 1 ? "" : "s");
  std::sort(blobs.begin(), blobs.end(),
            [](const Blob& a, const Blob& b) { return a.score > b.score; });

  // move blob position to the symmetry image nearest to the model
  // assert(st.cell == grid.unit_cell);
  gemmi::SubCells sc(model, st.cell, 10.0);
  sc.populate();
  for (Blob& blob : blobs)
    if (const gemmi::SubCells::Mark* mark = sc.find_nearest_atom(blob.pos)) {
      const gemmi::Position& ref = mark->to_cra(model).atom->pos;
      gemmi::Fractional fpos = st.cell.fractionalize(blob.pos);
      st.cell.apply_transform_inverse(fpos, mark->image_idx);
      blob.pos = st.cell.orthogonalize_in_pbc(ref, fpos);
    }

  // output results
  int n = 0;
  for (const Blob& b : blobs)
    printf("#%-4d %-3zu grid points, score %-8.4g  (%7.2f,%7.2f,%7.2f)\n",
           n++, b.points.size(), b.score, b.pos.x, b.pos.y, b.pos.z);
  return 0;
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  if (p.options[SigmaCutoff] && p.options[AbsCutoff])
    gemmi::fail("cannot use both --sigma and --abs");
  try {
    return run(p);
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
