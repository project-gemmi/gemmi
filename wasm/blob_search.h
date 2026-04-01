// Copyright Global Phasing Ltd.
#pragma once

#include "common.h"
#include <emscripten/val.h>
#include <gemmi/blob.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <limits>
#include <vector>

namespace blob_wasm {

inline em::val float_blob_view(const std::vector<float>& v) {
  return em::val(emscripten::typed_memory_view(v.size(), v.data()));
}

struct BlobSearchResult {
  std::vector<float> centroids_;
  std::vector<float> peak_positions_;
  std::vector<float> scores_;
  std::vector<float> volumes_;
  std::vector<float> peak_values_;

  BlobSearchResult() = default;

  explicit BlobSearchResult(const std::vector<gemmi::Blob>& blobs) {
    centroids_.reserve(3 * blobs.size());
    peak_positions_.reserve(3 * blobs.size());
    scores_.reserve(blobs.size());
    volumes_.reserve(blobs.size());
    peak_values_.reserve(blobs.size());
    for (const gemmi::Blob& blob : blobs) {
      centroids_.push_back(static_cast<float>(blob.centroid.x));
      centroids_.push_back(static_cast<float>(blob.centroid.y));
      centroids_.push_back(static_cast<float>(blob.centroid.z));
      peak_positions_.push_back(static_cast<float>(blob.peak_pos.x));
      peak_positions_.push_back(static_cast<float>(blob.peak_pos.y));
      peak_positions_.push_back(static_cast<float>(blob.peak_pos.z));
      scores_.push_back(static_cast<float>(blob.score));
      volumes_.push_back(static_cast<float>(blob.volume));
      peak_values_.push_back(static_cast<float>(blob.peak_value));
    }
  }

  size_t size() const { return scores_.size(); }
  em::val centroids() const { return float_blob_view(centroids_); }
  em::val peak_positions() const { return float_blob_view(peak_positions_); }
  em::val scores() const { return float_blob_view(scores_); }
  em::val volumes() const { return float_blob_view(volumes_); }
  em::val peak_values() const { return float_blob_view(peak_values_); }
};

inline void move_near_model(gemmi::NeighborSearch& ns, gemmi::Position& pos) {
  if (const auto* mark = ns.find_nearest_atom(pos)) {
    gemmi::const_CRA cra = mark->to_cra(*ns.model);
    pos = ns.grid.unit_cell.find_nearest_pbc_position(cra.atom->pos, pos,
                                                      mark->image_idx, true);
  }
}

inline bool valid_model_index(const gemmi::Structure* st, int model_index) {
  return st != nullptr && model_index >= 0 && model_index < (int) st->models.size();
}

inline void mask_model(gemmi::Grid<float>& grid, const gemmi::Model& model,
                       double radius, bool mask_waters, bool negate) {
  const float mask_value = negate ? std::numeric_limits<float>::infinity()
                                  : -std::numeric_limits<float>::infinity();
  for (const gemmi::Chain& chain : model.chains)
    for (const gemmi::Residue& residue : chain.residues) {
      if (!mask_waters && residue.is_water())
        continue;
      for (const gemmi::Atom& atom : residue.atoms) {
        if (atom.is_hydrogen())
          continue;
        grid.set_points_around(atom.pos, radius, mask_value);
      }
    }
  if (negate)
    grid.symmetrize_max();
  else
    grid.symmetrize_min();
}

inline void move_blobs_near_model(std::vector<gemmi::Blob>& blobs,
                                  gemmi::Model& model,
                                  const gemmi::UnitCell& cell) {
  if (model.chains.empty())
    return;
  gemmi::NeighborSearch ns(model, cell, 10.0);
  ns.populate();
  for (gemmi::Blob& blob : blobs) {
    move_near_model(ns, blob.centroid);
    move_near_model(ns, blob.peak_pos);
  }
}

inline BlobSearchResult* find_blobs(const gemmi::Grid<float>& source_grid,
                                    double cutoff, double min_volume,
                                    double min_score, double min_peak,
                                    bool negate, gemmi::Structure* st,
                                    int model_index, double mask_radius,
                                    bool mask_waters) {
  gemmi::Grid<float> grid = source_grid;
  if (valid_model_index(st, model_index))
    mask_model(grid, st->models[model_index], mask_radius, mask_waters, negate);
  gemmi::BlobCriteria criteria;
  criteria.cutoff = cutoff;
  criteria.min_volume = min_volume;
  criteria.min_score = min_score;
  criteria.min_peak = min_peak;
  std::vector<gemmi::Blob> blobs = gemmi::find_blobs_by_flood_fill(grid, criteria, negate);
  if (valid_model_index(st, model_index))
    move_blobs_near_model(blobs, st->models[model_index], source_grid.unit_cell);
  return new BlobSearchResult(blobs);
}

}  // namespace blob_wasm
