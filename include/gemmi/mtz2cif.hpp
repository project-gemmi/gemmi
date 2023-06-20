// Copyright 2020 Global Phasing Ltd.
//
// A class for converting MTZ (merged or unmerged) to SF-mmCIF

// TODO:
//  - cell parameters may be different in CELL and DCELL records, check for it
//  - check that the FP column is not from Refmac
//  - should we allow for repeated column names in MTZ?

#ifndef GEMMI_MTZ2CIF_HPP_
#define GEMMI_MTZ2CIF_HPP_

#include <ostream>
#include "mtz.hpp"       // for Mtz
#include "xds_ascii.hpp" // for XdsAscii
#include "merge.hpp"     // for Intensities

namespace gemmi {

class GEMMI_DLL MtzToCif {
public:
  // options that can be set directly
  std::vector<std::string> spec_lines; // conversion specification (cf. default_spec)
  const char* block_name = nullptr;  // NAME in data_NAME
  std::string entry_id = "xxxx";     // _entry.id
  bool with_comments = true;         // write comments
  bool with_history = true;          // write MTZ history in comments
  bool skip_empty = false;           // skip reflections with no values
  bool skip_negative_sigi = false;   // skip refl. with sigma(I) < 0 in unmerged
  bool enable_UB = false;            // write _diffrn_orient_matrix.UB
  bool write_staraniso_tensor = true; // write _reflns.pdbx_aniso_B_tensor_*
  bool write_special_marker_for_pdb = false;
  int less_anomalous = 0;            // skip (+)/(-) columns even if in spec
  std::string skip_empty_cols;       // columns used to determine "emptiness"
  double wavelength = NAN;           // user-specified wavelength
  int trim = 0;                      // output only reflections -N<=h,k,l<=N
  int free_flag_value = -1;          // -1 = auto: 0 or (if we have >50% of 0's) 1
  std::string staraniso_version;     // for _software.version in "special_marker"
  std::string gemmi_run_from;        // added to gemmi as _software.description

  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "H                          H index_h",
      "K                          H index_k",
      "L                          H index_l",
      "? FREE|RFREE|FREER|FreeR_flag|R-free-flags|FreeRflag I status S",
      "? IMEAN|I|IOBS|I-obs       J intensity_meas",
      "& SIG{prev}                Q intensity_sigma",
      "? I(+)|IOBS(+)|I-obs(+)    K pdbx_I_plus",
      "& SIG{prev}                M pdbx_I_plus_sigma",
      "? I(-)|IOBS(-)|I-obs(-)    K pdbx_I_minus",
      "& SIG{prev}                M pdbx_I_minus_sigma",
      // TODO: FP from Refmac should show warning or error
      "? F|FP|FOBS|F-obs          F F_meas_au",
      "& SIG{prev}                Q F_meas_sigma_au",
      "? F(+)|FOBS(+)|F-obs(+)    G pdbx_F_plus",
      "& SIG{prev}                L pdbx_F_plus_sigma",
      "? F(-)|FOBS(-)|F-obs(-)    G pdbx_F_minus",
      "& SIG{prev}                L pdbx_F_minus_sigma",
      "? DP                       D pdbx_anom_difference",
      "& SIGDP                    Q pdbx_anom_difference_sigma",
      "? FC                       F F_calc",
      "? PHIC                     P phase_calc",
      "? FOM                      W fom",
      "? HLA                      A pdbx_HL_A_iso",
      "& HLB                      A pdbx_HL_B_iso",
      "& HLC                      A pdbx_HL_C_iso",
      "& HLD                      A pdbx_HL_D_iso",
      "? FWT|2FOFCWT              F pdbx_FWT",
      "& PHWT|PH2FOFCWT           P pdbx_PHWT .3f",
      "? DELFWT|FOFCWT            F pdbx_DELFWT",
      "& DELPHWT|PHDELWT|PHFOFCWT P pdbx_DELPHWT .3f",
      nullptr
    };
    static const char* unmerged[] = {
      "$dataset    diffrn_id",  // diffrn_id - sweep id deduced from BATCH
      "$counter    id",         // reflection counter (1, 2, ...)
      "H         H index_h",
      "K         H index_k",
      "L         H index_l",
      "? I       J intensity_net",
      "& SIGI    Q intensity_sigma .5g",
      // new! https://github.com/wwpdb-dictionaries/mmcif_pdbx/pull/33
      "?ROT      R pdbx_scan_angle",
      "$image      pdbx_image_id",
      nullptr
    };
    return for_merged ? merged : unmerged;
  }

  void write_cif(const Mtz& mtz, const Mtz* mtz2,
                 SMat33<double>* staraniso_b, std::ostream& os);
  void write_cif_from_xds(const XdsAscii& xds, std::ostream& os);
};

GEMMI_DLL void write_staraniso_b_in_mmcif(const SMat33<double>& b,
                                          const std::string& entry_id,
                                          char* buf, std::ostream& os);

/// remove '_dataset_name' that can be appended to column names in ccp4i
GEMMI_DLL void remove_appendix_from_column_names(Mtz& mtz, std::ostream& out);

GEMMI_DLL bool validate_merged_mtz_deposition_columns(const Mtz& mtz, std::ostream& out);

// note: both mi and ui get modified
GEMMI_DLL bool validate_merged_intensities(Intensities& mi, Intensities& ui,
                                           bool relaxed_check, std::ostream& out);

} // namespace gemmi
#endif
