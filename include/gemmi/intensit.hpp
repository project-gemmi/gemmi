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
    std::int8_t isign;  // 1 for I(+), -1 for I(-), 0 for mean or unmerged
    std::int8_t isym;   // for unmerged data: encodes symmetry op like M/ISYM in MTZ
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
  std::vector<Op> isym_ops;
  AnisoScaling staraniso_b;

  static const char* type_str(DataType data_type) {
    switch (data_type) {
      case DataType::Unmerged: return "I";
      case DataType::Mean: return "<I>";
      case DataType::Anomalous: return "I+/I-";
      case DataType::MergedAM:
      case DataType::MergedMA:
      case DataType::Unknown: return "n/a";
    }
    unreachable();
  }

  const char* type_str() const { return Intensities::type_str(type); }

  std::string spacegroup_str() const { return spacegroup ? spacegroup->xhm() : "none"; }

  // returns (d_max, d_min)
  std::array<double,2> resolution_range() const;

  // pre: both are sorted
  Correlation calculate_correlation(const Intensities& other) const;

  void add_if_valid(const Miller& hkl, int8_t isign, int8_t isym, double value, double sigma) {
    // XDS marks rejected reflections with negative sigma.
    // Sigma 0.0 rarely happens (e.g. 5tkn), but is also problematic.
    if (!std::isnan(value) && sigma > 0)
      data.push_back({hkl, isign, isym, /*nobs=*/0, value, sigma});
  }

  void remove_systematic_absences() {
    if (!spacegroup)
      return;
    GroupOps gops = spacegroup->operations();
    vector_remove_if(data, [&](Refl& x) { return gops.is_systematically_absent(x.hkl); });
  }

  void sort() { std::sort(data.begin(), data.end()); }

  void merge_in_place(DataType data_type);

  void switch_to_asu_indices();

  void read_unmerged_intensities_from_mtz(const Mtz& mtz);
  void read_mean_intensities_from_mtz(const Mtz& mtz);
  // with check_complete=true, throw if anomalous data is null where it shouldn't be
  void read_anomalous_intensities_from_mtz(const Mtz& mtz, bool check_complete=false);

  void read_mtz(const Mtz& mtz, DataType data_type) {
    bool check_anom_complete = false;
    if (data_type == DataType::Unknown)
      data_type = mtz.batches.empty() ? DataType::MergedMA : DataType::Unmerged;
    if (data_type == DataType::MergedAM) {
      data_type = mtz.iplus_column() ? DataType::Anomalous : DataType::Mean;
      // if I(+) and I(-) is empty where IMEAN is not, throw error
      check_anom_complete = true;
    }
    if (data_type == DataType::MergedMA)
      data_type = mtz.imean_column() ? DataType::Mean : DataType::Anomalous;

    if (data_type == DataType::Unmerged)
      read_unmerged_intensities_from_mtz(mtz);
    else if (data_type == DataType::Mean)
      read_mean_intensities_from_mtz(mtz);
    else // (data_type == DataType::Anomalous)
      read_anomalous_intensities_from_mtz(mtz, check_anom_complete);
  }

  void read_unmerged_intensities_from_mmcif(const ReflnBlock& rb);
  void read_mean_intensities_from_mmcif(const ReflnBlock& rb);
  void read_anomalous_intensities_from_mmcif(const ReflnBlock& rb, bool check_complete=false);

  void read_f_squared_from_mmcif(const ReflnBlock& rb);

  void read_mmcif(const ReflnBlock& rb, DataType data_type) {
    bool check_anom_complete = false;
    if (data_type == DataType::Unknown)
      data_type = rb.is_unmerged() ? DataType::Unmerged : DataType::MergedMA;

    if (data_type == DataType::MergedAM || data_type == DataType::MergedMA) {
      bool has_anom = rb.find_column_index("pdbx_I_plus") != -1;
      bool has_mean = rb.find_column_index("intensity_meas") != -1;
      if (!has_anom && !has_mean)
        fail("Intensities not found in the mmCIF file, block ", rb.block.name,
             " has neither intensity_meas nor pdbx_I_plus/minus");
      if (data_type == DataType::MergedAM) {
        data_type = has_anom ? DataType::Anomalous : DataType::Mean;
        // if both I(+) and I(-) is empty where IMEAN has value, throw error
        check_anom_complete = true;
      }
      if (data_type == DataType::MergedMA)
        data_type =  has_mean ? DataType::Mean : DataType::Anomalous;
    }
    if (data_type == DataType::Unmerged)
      read_unmerged_intensities_from_mmcif(rb);
    else if (data_type == DataType::Mean)
      read_mean_intensities_from_mmcif(rb);
    else // (data_type == DataType::Anomalous)
      read_anomalous_intensities_from_mmcif(rb, check_anom_complete);
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
