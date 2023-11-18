// Copyright 2019 Global Phasing Ltd.

#include <cstdio>
#include <algorithm>  // count
#include <stdexcept>
#include "gemmi/blob.hpp"
#include "gemmi/assembly.hpp"  // for expand_ncs
#include "gemmi/polyheur.hpp"  // for remove_waters
#include "gemmi/modify.hpp"    // for remove_hydrogens
#include "gemmi/stats.hpp"     // for Variance
#include "gemmi/neighbor.hpp"  // for NeighborSearch
#include "gemmi/mmread_gz.hpp" // for read_structure_gz
#include "gemmi/calculate.hpp" // for calculate_center_of_mass
#include "mapcoef.h"

#define GEMMI_PROG blobs
#include "options.h"

namespace {

using std::printf;

enum OptionIndex { SigmaCutoff=AfterMapOptions, AbsCutoff,
                   MaskRadius, MaskWater,
                   MinVolume, MinScore, MinSigma, MinDensity,
                   Dimple };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] MTZ_OR_MMCIF PDB_OR_MMCIF"
    "\n\nSearch for umodelled blobs of electron density."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { NoOp, 0, "", "", Arg::None,
    "\nThe area around model is masked to search only unmodelled density." },
  { MaskRadius, 0, "", "mask-radius", Arg::Float,
    "  --mask-radius=NUMBER  \tMask radius (default: 2.0 A)." },
  { MaskWater, 0, "", "mask-waters", Arg::None,
    "  --mask-water  \tMask water (water is not masked by default)." },

  { NoOp, 0, "", "", Arg::None, "\nSearching blobs of density above:" },
  { SigmaCutoff, 0, "", "sigma", Arg::Float,
    "  --sigma=NUMBER  \tSigma (RMSD) level (default: 1.0)." },
  { AbsCutoff, 0, "", "abs", Arg::Float,
    "  --abs=NUMBER  \tAbsolute level in electrons/A^3." },

  { NoOp, 0, "", "", Arg::None, "\nBlob criteria:" },
  { MinVolume, 0, "", "min-volume", Arg::Float,
    "  --min-volume=NUMBER  \tMinimal volume (default: 10.0 A^3)." },
  { MinScore, 0, "", "min-score", Arg::Float,
    "  --min-score=NUMBER  \tMin. this electrons in blob (default: 15.0)." },
  { MinSigma, 0, "", "min-sigma", Arg::Float,
    "  --min-sigma=NUMBER  \tMin. peak rmsd (default: 0.0)." },
  { MinDensity, 0, "", "min-peak", Arg::Float,
    "  --min-peak=NUMBER  \tMin. peak density (default: 0.0 el/A^3)." },

  { NoOp, 0, "", "", Arg::None, "\nOptions for map calculation:" },
  MapUsage[Diff],
  MapUsage[Section],
  MapUsage[FLabel],
  MapUsage[PhLabel],
  MapUsage[WeightLabel],
  MapUsage[GridDims],
  MapUsage[ExactDims],
  MapUsage[Sample],
  MapUsage[GridQuery],
  MapUsage[TimingFft],

  { Dimple, 0, "", "dimple", Arg::None, nullptr }, // output for Dimple
  { 0, 0, 0, 0, 0, 0 }
};

gemmi::const_CRA move_near_model(gemmi::NeighborSearch& ns, gemmi::Position& pos) {
  if (const auto* mark = ns.find_nearest_atom(pos)) {
    gemmi::const_CRA cra = mark->to_cra(*ns.model);
    pos = ns.grid.unit_cell.find_nearest_pbc_position(cra.atom->pos, pos,
                                                      mark->image_idx, true);
    return cra;
  }
  return {nullptr, nullptr, nullptr};
}

int run(OptParser& p) {
  std::string sf_path = p.nonOption(0);
  std::string model_path = p.coordinate_input_file(1);

  // read model, remove hydrogens if any
  if (p.options[Verbose])
    printf("Reading coordinates from %s ...\n", model_path.c_str());
  gemmi::Structure st = gemmi::read_structure_gz(model_path);
  if (st.models.empty() || st.models[0].chains.empty()) {
    std::fprintf(stderr, "Not a coordinate file: %s\n", model_path.c_str());
    return 1;
  }
  if (st.models.size() > 1)
    std::fprintf(stderr, "Note: only the first model is used.\n");
  gemmi::Model& model = st.models[0];
  gemmi::remove_hydrogens(model);
  if (!p.options[MaskWater])
    gemmi::remove_waters(model);

  // read map (includes FFT)
  FILE* verbose_output = p.options[Verbose] ? stdout : nullptr;
  gemmi::Grid<float> grid = read_sf_and_fft_to_map(sf_path.c_str(), p.options,
                                                   verbose_output, true);
  if (p.options[Verbose])
    printf("Unit cell: %g A^3, grid points: %zu, volume/point: %g A^3.\n",
           grid.unit_cell.volume, grid.point_count(),
           grid.unit_cell.volume / grid.point_count());
  // move blob position to the symmetry image nearest to the model
  if (st.find_spacegroup() != grid.spacegroup)
    std::fprintf(stderr, "Warning: different space groups in model and data.\n");
  if (!st.cell.approx(grid.unit_cell, 0.1))
    std::fprintf(stderr, "Warning: different unit cells in model and data.\n");

  if (st.ncs_not_expanded()) {
    std::fprintf(stderr, "Note: NCS is expanded, listed blobs may be redundant.\n");
    gemmi::expand_ncs(st, gemmi::HowToNameCopiedChain::AddNumber);
  }

  // calculate map RMSD and setup blob criteria
  gemmi::BlobCriteria criteria;
  gemmi::Variance grid_variance(grid.data.begin(), grid.data.end());
  double rmsd = std::sqrt(grid_variance.for_population());
  double sigma_level = 1.0;
  if (p.options[AbsCutoff]) {
    criteria.cutoff = std::strtod(p.options[AbsCutoff].arg, nullptr);
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
  if (p.options[MinSigma])
    criteria.min_peak = std::strtod(p.options[MinSigma].arg, nullptr) * rmsd;
  if (p.options[MinDensity])
    criteria.min_peak = std::strtod(p.options[MinDensity].arg, nullptr);
  printf("Map RMSD: %.3f. Searching blobs above %.3f e/A^3 (%.3f sigma).\n",
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
    size_t n = std::count(grid.data.begin(), grid.data.end(), -INFINITY);
    printf("Masked points: %zu of %zu.\n", n, grid.point_count());
  }
  if (p.options[Dimple]) {
    gemmi::Position com = gemmi::calculate_center_of_mass(model).get();
    printf("Center of mass: %.2f %.2f %.2f\n", com.x, com.y, com.z);
  }

  // find and sort blobs
  std::vector<gemmi::Blob> blobs = gemmi::find_blobs_by_flood_fill(grid, criteria);
  if (p.options[Verbose])
    printf("%zu blob%s found.\n", blobs.size(), blobs.size() == 1 ? "" : "s");

  gemmi::NeighborSearch ns(model, grid.unit_cell, 10.0);
  ns.populate();

  // output results
  for (size_t i = 0; i != blobs.size(); ++i) {
    gemmi::const_CRA cra = move_near_model(ns, blobs[i].centroid);
    // Blob::peak_pos is left not moved, but we don't use it below
    const gemmi::Blob& b = blobs[i];
    std::string residue_info = "none";
    if (cra.chain && cra.residue)
      residue_info = cra.chain->name + " " + cra.residue->str();
    printf("#%-2zu %5.1f el in %5.1f A^3, %4.1f rmsd,"
           " (%6.1f,%6.1f,%6.1f) near %s\n",
           i, b.score, b.volume, b.peak_value / rmsd,
           b.centroid.x, b.centroid.y, b.centroid.z, residue_info.c_str());
  }
  return 0;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);
  p.check_exclusive_pair(SigmaCutoff, AbsCutoff);
  try {
    return run(p);
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
