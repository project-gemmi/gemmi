// Copyright 2019 Global Phasing Ltd.
//
// Metadata from coordinate files.

#ifndef GEMMI_METADATA_HPP_
#define GEMMI_METADATA_HPP_

#include <array>
#include "math.hpp"

namespace gemmi {

// corresponds to the mmCIF _software category
struct SoftwareItem {
  enum Classification {
    DataCollection, DataExtraction, DataProcessing, DataReduction,
    DataScaling, ModelBuilding, Phasing, Refinement, Unspecified
  };
  std::string name;
  std::string version;
  Classification classification = Unspecified;
  int pdbx_ordinal = -1;
};

struct TlsGroup {
  std::string id;           // _pdbx_refine_tls.id
  std::string selection;    // _pdbx_refine_tls_group.selection_details
  Position origin;          // _pdbx_refine_tls.origin_x/y/z
  Mat33 T;                  // _pdbx_refine_tls.T[][]
  Mat33 L;                  // _pdbx_refine_tls.L[][]
  Mat33 S;                  // _pdbx_refine_tls.S[][]
};

// RefinementInfo corresponds to REMARK 3.
// BasicRefinementInfo is used for both total and per-bin statistics.
// For per-bin data, each values corresponds to one _refine_ls_shell.* tag.
struct BasicRefinementInfo {
  double resolution_high = NAN;      // _refine.ls_d_res_high
  double resolution_low = NAN;       // _refine.ls_d_res_low
  double completeness = NAN;         // _refine.ls_percent_reflns_obs
  int reflection_count = -1;         // _refine.ls_number_reflns_obs
  int rfree_set_count = -1;          // _refine.ls_number_reflns_R_free
  double r_all = NAN;                // _refine.ls_R_factor_obs
  double r_work = NAN;               // _refine.ls_R_factor_R_work
  double r_free = NAN;               // _refine.ls_R_factor_R_free
};

struct RefinementInfo : BasicRefinementInfo {
  struct Restr {
    std::string type;
    int count = -1;
    double weight = NAN;
    std::string function;
  };
  std::string id;
  std::string cross_validation_method; // _refine.pdbx_ls_cross_valid_method
  std::string rfree_selection_method;  // _refine.pdbx_R_Free_selection_details
  int bin_count = -1;        // _refine_ls_shell.pdbx_total_number_of_bins_used
  std::vector<BasicRefinementInfo> bins;
  double b_wilson = NAN;              // _reflns.B_iso_Wilson_estimate
  double mean_b = NAN;                // _refine.B_iso_mean
  Mat33 aniso_b;                      // _refine.aniso_B[][]
  double luzzati_error = NAN; // _refine_analyze.Luzzati_coordinate_error_obs
  double dpi_blow_r = NAN;            // _refine.pdbx_overall_SU_R_Blow_DPI
  double dpi_blow_rfree = NAN;        // _refine.pdbx_overall_SU_R_free_Blow_DPI
  double dpi_cruickshank_r = NAN;     // _refine.overall_SU_R_Cruickshank_DPI
  double dpi_cruickshank_rfree = NAN; // _refine.pdbx_overall_SU_R_free_Cruickshank_DPI
  double cc_fo_fc = NAN;              // _refine.correlation_coeff_Fo_to_Fc
  double cc_fo_fc_free = NAN;         // _refine.correlation_coeff_Fo_to_Fc_free
  std::vector<Restr> restr;           // _refine_ls_restr
  std::vector<TlsGroup> tls_groups;   // _pdbx_refine_tls
  std::string remarks;
};

// Corresponds to REMARK 200.
struct XRDInfo {
  std::string collection_date;        // _diffrn_detector.pdbx_collection_date
  double temperature;                 // _diffrn.ambient_temp
};


struct Metadata {
  std::vector<RefinementInfo> refinement;
  std::vector<SoftwareItem> software;

  bool has(double RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !std::isnan(r.*field); });
  }
  bool has(int RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return r.*field != -1; });
  }
  bool has(std::string RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !(r.*field).empty(); });
  }
  bool has(Mat33 RefinementInfo::*field) const {
    return std::any_of(refinement.begin(), refinement.end(),
        [&](const RefinementInfo& r) { return !std::isnan((r.*field)[0][0]); });
  }
  bool has_restr() const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !r.restr.empty(); });
  }
  bool has_tls() const {
    return std::any_of(refinement.begin(), refinement.end(),
            [&](const RefinementInfo& r) { return !r.tls_groups.empty(); });
  }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
