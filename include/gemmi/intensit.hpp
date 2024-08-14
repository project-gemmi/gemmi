// Copyright 2020 Global Phasing Ltd.
//
// Class Intensities that reads multi-record data from MTZ, mmCIF or XDS_ASCII
// and merges it into mean or anomalous intensities.
// It can also read merged data.

#ifndef GEMMI_INTENSIT_HPP_
#define GEMMI_INTENSIT_HPP_

#include <cassert>
#include "symmetry.hpp"
#include "unitcell.hpp"
#include "util.hpp"     // for vector_remove_if
#include "mtz.hpp"      // for Mtz
#include "refln.hpp"    // for ReflnBlock
#include "stats.hpp"    // for Correlation
#include "xds_ascii.hpp" // for XdsAscii

namespace gemmi {

inline std::string miller_str(const Miller& hkl) {
  return cat('(', hkl[0], ' ', hkl[1], ' ', hkl[2], ')');
}

/// Returns STARANISO version or empty string.
GEMMI_DLL std::string read_staraniso_b_from_mtz(const Mtz& mtz, SMat33<double>& output);

struct GEMMI_DLL Intensities {
  struct Refl {
    Miller hkl;
    short isign;  // 1 for I(+), -1 for I(-), 0 for mean
    short nobs;
    double value;
    double sigma;

    bool operator<(const Refl& o) const {
      return std::tie(hkl[0], hkl[1], hkl[2], isign) <
             std::tie(o.hkl[0], o.hkl[1], o.hkl[2], o.isign);
    }
    const char* intensity_label() const {
      if (isign == 0)
        return "<I>";
      return isign > 0 ? "I(+)" : "I(-)";
    }
    std::string hkl_label() const {
      return cat(intensity_label(), ' ', miller_str(hkl));
    }
  };

  struct AnisoScaling {
    SMat33<double> b = {0., 0., 0., 0., 0., 0.};

    bool ok() const { return !b.all_zero(); }
    double scale(const Miller& hkl, const UnitCell& cell) const {
      Vec3 s = cell.frac.mat.left_multiply(Vec3(hkl[0], hkl[1], hkl[2]));
      return std::exp(0.5 * b.r_u_r(s));
    }
  };

  std::vector<Refl> data;
  const SpaceGroup* spacegroup = nullptr;
  UnitCell unit_cell;
  double unit_cell_rmsd[6] = {0., 0., 0., 0., 0., 0.};
  double wavelength;
  DataType type = DataType::Unknown;
  AnisoScaling staraniso_b;

  static const char* type_str(DataType data_type) {
    switch (data_type) {
      case DataType::Unknown: return "n/a";
      case DataType::Unmerged: return "I";
      case DataType::Mean: return "<I>";
      case DataType::Anomalous: return "I+/I-";
    }
    unreachable();
  }

  const char* type_str() const { return Intensities::type_str(type); }

  std::string spacegroup_str() const { return spacegroup ? spacegroup->xhm() : "none"; }

  // returns (d_max, d_min)
  std::array<double,2> resolution_range() const;

  // pre: both are sorted
  Correlation calculate_correlation(const Intensities& other) const;

  void remove_systematic_absences() {
    if (!spacegroup)
      return;
    GroupOps gops = spacegroup->operations();
    vector_remove_if(data, [&](Refl& x) { return gops.is_systematically_absent(x.hkl); });
  }

  void sort() { std::sort(data.begin(), data.end()); }

  void merge_in_place(DataType data_type);

  // for unmerged centric reflections set isign=1.
  void switch_to_asu_indices(bool merged=false);

  void read_unmerged_intensities_from_mtz(const Mtz& mtz);
  void read_mean_intensities_from_mtz(const Mtz& mtz);
  void read_anomalous_intensities_from_mtz(const Mtz& mtz, bool check_complete=false);

  void read_merged_intensities_from_mtz(const Mtz& mtz) {
    if (mtz.column_with_label("I(+)"))
      read_anomalous_intensities_from_mtz(mtz, true);
    else
      read_mean_intensities_from_mtz(mtz);
  }

  void read_mtz(const Mtz& mtz, DataType data_type) {
    switch (data_type) {
      case DataType::Unmerged:  read_unmerged_intensities_from_mtz(mtz); break;
      case DataType::Mean:      read_mean_intensities_from_mtz(mtz); break;
      case DataType::Anomalous: read_anomalous_intensities_from_mtz(mtz); break;
      case DataType::Unknown: assert(0);
    }
  }

  void read_unmerged_intensities_from_mmcif(const ReflnBlock& rb);
  void read_mean_intensities_from_mmcif(const ReflnBlock& rb);
  void read_anomalous_intensities_from_mmcif(const ReflnBlock& rb, bool check_complete=false);

  void read_merged_intensities_from_mmcif(const ReflnBlock& rb) {
    if (rb.find_column_index("pdbx_I_plus") != -1)
      read_anomalous_intensities_from_mmcif(rb, true);
    else if (rb.find_column_index("intensity_meas") != -1)
      read_mean_intensities_from_mmcif(rb);
    else
      fail("Intensities not found in the mmCIF file, block ", rb.block.name,
           " has neither intensity_meas nor pdbx_I_plus/minus");
  }

  void read_f_squared_from_mmcif(const ReflnBlock& rb);

  void read_mmcif(const ReflnBlock& rb, DataType data_type) {
    switch (data_type) {
      case DataType::Unmerged:  read_unmerged_intensities_from_mmcif(rb); break;
      case DataType::Mean:      read_mean_intensities_from_mmcif(rb); break;
      case DataType::Anomalous: read_anomalous_intensities_from_mmcif(rb); break;
      case DataType::Unknown:   break;
    }
  }

  void read_unmerged_intensities_from_xds(const XdsAscii& xds);

  // returns STARANISO version or empty string
  std::string take_staraniso_b_from_mtz(const Mtz& mtz) {
    return read_staraniso_b_from_mtz(mtz, staraniso_b.b);
  }
  bool take_staraniso_b_from_mmcif(const cif::Block& block);

  Mtz prepare_merged_mtz(bool with_nobs);
};

} // namespace gemmi
#endif
