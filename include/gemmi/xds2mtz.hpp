// Copyright Global Phasing Ltd.
//
// Convert XDS_ASCII to MTZ.

#ifndef GEMMI_XDS2MTZ_HPP_
#define GEMMI_XDS2MTZ_HPP_

#include <map>
#include "mtz.hpp"        // for Mtz
#include "xds_ascii.hpp"  // for XdsAscii
#include "intensit.hpp"   // for Intensities

namespace gemmi {

inline Mtz xds_to_mtz(XdsAscii& xds) {
  if (xds.is_merged()) {
    Intensities intensities;
    intensities.import_xds(xds);
    return intensities.prepare_merged_mtz(/*with_nobs=*/false);
  }
  Mtz mtz;
  mtz.cell.set_from_array(xds.cell_constants);
  mtz.spacegroup = find_spacegroup_by_number(xds.spacegroup_number);
  mtz.add_base();
  const char* pxd[3] = {"XDSproject", "XDScrystal", "XDSdataset"};
  if (xds.isets.empty()) {
    mtz.datasets.push_back({1, pxd[0], pxd[1], pxd[2], mtz.cell, xds.wavelength});
  } else {
    for (XdsAscii::Iset& iset : xds.isets) {
      double wavelength = iset.wavelength != 0 ? iset.wavelength : xds.wavelength;
      UnitCell cell;
      cell.set_from_array(iset.cell_constants[0] != 0 ? iset.cell_constants
                                                      : xds.cell_constants);
      mtz.datasets.push_back({iset.id, pxd[0], pxd[1], pxd[2], cell, wavelength});
    }
  }
  mtz.add_column("M/ISYM", 'Y', 0, -1, false);
  mtz.add_column("BATCH", 'B', 0, -1, false);
  mtz.add_column("I", 'J', 0, -1, false);
  mtz.add_column("SIGI", 'Q', 0, -1, false);
  mtz.add_column("XDET", 'R', 0, -1, false);
  mtz.add_column("YDET", 'R', 0, -1, false);
  mtz.add_column("ROT", 'R', 0, -1, false);
  if (xds.read_columns >= 11) {
    mtz.add_column("FRACTIONCALC", 'R', 0, -1, false);
    mtz.add_column("LP", 'R', 0, -1, false);
    mtz.add_column("CORR", 'R', 0, -1, false);
    if (xds.read_columns > 11)
      mtz.add_column("MAXC", 'I', 0, -1, false);
  }
  mtz.add_column("FLAG", 'I', 0, -1, false);
  mtz.nreflections = (int) xds.data.size();
  mtz.data.resize(mtz.columns.size() * xds.data.size());
  UnmergedHklMover hkl_mover(mtz.spacegroup);
  int max_frame = 0;
  for (const XdsAscii::Refl& refl : xds.data)
    max_frame = std::max(max_frame, refl.frame());
  int iset_offset = (max_frame + 11000) / 10000 * 10000;
  // iset,frame -> batch
  std::map<std::pair<int,int>, int> frames;
  size_t k = 0;
  for (const XdsAscii::Refl& refl : xds.data) {
    auto hkl = refl.hkl;
    int isym = hkl_mover.move_to_asu(hkl);
    for (size_t j = 0; j != 3; ++j)
      mtz.data[k++] = (float) hkl[j];
    mtz.data[k++] = (float) isym;
    int frame = refl.frame();
    int batch = frame + iset_offset * std::max(refl.iset - 1, 0);
    frames.emplace(std::make_pair(refl.iset, frame), batch);
    mtz.data[k++] = (float) batch;
    mtz.data[k++] = (float) refl.iobs;  // I
    mtz.data[k++] = (float) std::fabs(refl.sigma);  // SIGI
    mtz.data[k++] = (float) refl.xd;
    mtz.data[k++] = (float) refl.yd;
    mtz.data[k++] = (float) xds.rot_angle(refl);  // ROT
    if (xds.read_columns >= 11) {
      mtz.data[k++] = float(0.01 * refl.peak);  // FRACTIONCALC
      mtz.data[k++] = (float) refl.rlp;
      mtz.data[k++] = float(0.01 * refl.corr);
      if (xds.read_columns > 11)
        mtz.data[k++] = (float) refl.maxc;
    }
    mtz.data[k++] = refl.sigma < 0 ? 64.f : 0.f;  // FLAG
  }
  // Prepare a similar batch header as Pointless.
  Mtz::Batch batch;
  batch.set_dataset_id(1);

  // We don't set lbcell (refinement flags for unit cell),
  // because it's probably not used by any program anyway.
  // batch.ints[4] to [9] = left unset

  // We also skip jumpax, which is defined as:
  // reciprocal axis closest to principle goniostat axis E1
  // batch.ints[11] = left unset

  // We assume one crystal, one goniostat axis, one detector.
  batch.ints[12] = 1;  // ncryst
  batch.ints[14] = 2;  // ldtype 3D
  batch.ints[15] = 1;  // jsaxs - goniostat scan axis number
  batch.ints[17] = 1;  // ngonax - number of goniostat axes
  batch.ints[19] = 1;  // ndet

  batch.set_cell(mtz.cell);  // batch.floats[0] to [5]
  Vec3 s0(-1, 0, 0); // will be re-set if we have geometry info
  try {
    Mat33 Q = xds.calculate_conversion_from_cambridge().inverse();
    s0 = -Q.multiply(xds.get_s0_direction());
    // Orientation matrix U. It is calculated differently in Pointless,
    // so the results are slightly different (due to limited precision
    // of numbers in XDS file).
    Mat33 U = Q.multiply(xds.get_orientation());
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        batch.floats[6 + 3*i + j] = (float) U[j][i];
  } catch (std::runtime_error& e) {
    // if some of the headers are absent, U is not set
    // and s0 is set to "idealised" s0
    std::fprintf(stderr, "Note: %s, orientation matrix U not set.\n", e.what());
  }
  batch.floats[21] = float(xds.reflecting_range_esd);  // crydat(0)
  // In the so-called "Cambridge" frame (as used by Mosflm),
  // the principal rotation axis is along z
  // and the incident beam S0 is along x.
  // Therefore, we set the rotation axis phi (scanax) to:
  batch.floats[38+2] = 1.f;  // scanax = [0, 0, 1]
  batch.floats[47] = float(xds.oscillation_range);  // phi range
  // E1,E2,E3 vectors are the goniostat axes in Cambridge laboratory frame.
  // E1 is set to rotation axis. E2 and E3 are not set for ngonax==1.
  batch.floats[59+2] = 1.f;  // e1 = scanax

  // Idealised source vector is -x ([-1 0 0]), antiparallel to beam.
  batch.floats[80+0] = -1.f;  // source[0]

  // s0 source vector (including tilts), CMtz::MTZBAT::so in libccp4
  batch.floats[83] = (float) s0.x;
  batch.floats[84] = (float) s0.y;
  batch.floats[85] = (float) s0.z;

  batch.set_wavelength((float)xds.wavelength);  // batch.floats[86]
  // or batch.set_wavelength(iset.wavelength);
  // Detector geometry.
  batch.floats[111] = (float) xds.detector_distance;  // dx[0]
  batch.floats[113] = 1.f;  // detlm[0][0][0]
  batch.floats[114] = (float) xds.nx;
  batch.floats[115] = 1.f;
  batch.floats[116] = (float) xds.ny;
  batch.axes.push_back("PHI");  // gonlab[0]

  for (auto& t : frames) {
    batch.set_dataset_id(t.first.first);
    batch.number = t.second;
    int frame = t.first.second;
    double phistt = xds.starting_angle +
                    xds.oscillation_range * (frame - xds.starting_frame);
    batch.floats[36] = float(phistt);
    batch.floats[37] = float(phistt + xds.oscillation_range);  // phiend
    mtz.batches.push_back(batch);
  }
  mtz.sort(5);
  return mtz;
}


} // namespace gemmi
#endif
