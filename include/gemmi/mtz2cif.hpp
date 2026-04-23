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
#include "intensit.hpp"  // for Intensities
#include "logger.hpp"    // for Logger
#include "metadata.hpp"  // for SoftwareItem

namespace gemmi {

/// @file
/// @brief Converter for MTZ reflection data to SF-mmCIF format.

/// @brief Converts MTZ files (merged or unmerged) to SF-mmCIF reflection tables.
///
/// This class provides configuration options for column selection, naming,
/// filtering, and metadata handling when converting MTZ format to mmCIF.
class GEMMI_DLL MtzToCif {
public:
  /// Column conversion specification lines (see default_spec for format).
  std::vector<std::string> spec_lines;
  /// CIF data block name (NAME in data_NAME).
  const char* block_name = nullptr;
  /// Entry identifier (_entry.id tag).
  std::string entry_id = "xxxx";
  /// Whether to write comments describing the conversion.
  bool with_comments = true;
  /// Whether to write MTZ history records in comments.
  bool with_history = true;
  /// Skip reflections where all selected numeric columns are missing.
  bool skip_empty = false;
  /// Skip unmerged reflections with sigma(I) < 0.
  bool skip_negative_sigi = false;
  /// Write _diffrn_orient_matrix.UB orientation matrix.
  bool enable_UB = false;
  /// Write _reflns.pdbx_aniso_B_tensor_* Starraniso B-tensor (if available).
  bool write_staraniso_tensor = true;
  /// Write PDB-specific special marker for validation.
  bool write_special_marker_for_pdb = false;
  /// If non-zero, skip anomalous (+/-) column pairs.
  int less_anomalous = 0;
  /// Columns used to determine if reflection is "empty" (when skip_empty=true).
  std::string skip_empty_cols;
  /// User-specified wavelength (NAN means use MTZ value).
  double wavelength = NAN;
  /// Trim reflections: output only those with -N<=h,k,l<=N (0 = no trim).
  int trim = 0;
  /// Free flag value: -1=auto, 0 or 1=explicit.
  int free_flag_value = -1;
  /// Starraniso version string for metadata.
  std::string staraniso_version;
  /// Description string appended to gemmi software entry.
  std::string gemmi_run_from;

  /// @brief Get default column specification for merged or unmerged data.
  ///
  /// The returned spec_lines describe MTZ-to-mmCIF column mapping.
  /// Format: [?|&|$|H][COLUMN_NAME] [TYPE] [mmCIF_TAG] [FORMAT]
  /// - ? = optional column (try alternatives separated by |)
  /// - & = required, uses previous column's result
  /// - $ = internal (dataset_id, counter)
  /// - H = required by IUCR standard
  /// @param for_merged If true, return spec for merged data; else unmerged.
  /// @return Null-terminated array of spec strings.
  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "H                              H index_h",
      "K                              H index_k",
      "L                              H index_l",
      "? FREE|RFREE|FREER|FreeR_flag|R-free-flags|FreeRflag I status S",
      "? IMEAN|I|IOBS|I-obs           J intensity_meas",
      "& SIG{prev}                    Q intensity_sigma",
      "? I(+)|IOBS(+)|I-obs(+)|Iplus  K pdbx_I_plus",
      "& SIG{prev}                    M pdbx_I_plus_sigma",
      "? I(-)|IOBS(-)|I-obs(-)|Iminus K pdbx_I_minus",
      "& SIG{prev}                    M pdbx_I_minus_sigma",
      // TODO: FP from Refmac should show warning or error
      "? F|FP|FOBS|F-obs              F F_meas_au",
      "& SIG{prev}                    Q F_meas_sigma_au",
      "? F(+)|FOBS(+)|F-obs(+)|Fplus  G pdbx_F_plus",
      "& SIG{prev}                    L pdbx_F_plus_sigma",
      "? F(-)|FOBS(-)|F-obs(-)|Fminus G pdbx_F_minus",
      "& SIG{prev}                    L pdbx_F_minus_sigma",
      "? DP                           D pdbx_anom_difference",
      "& SIGDP                        Q pdbx_anom_difference_sigma",
      "? FC                           F F_calc",
      "? PHIC                         P phase_calc",
      "? FOM                          W fom",
      "? HLA                          A pdbx_HL_A_iso",
      "& HLB                          A pdbx_HL_B_iso",
      "& HLC                          A pdbx_HL_C_iso",
      "& HLD                          A pdbx_HL_D_iso",
      "? FWT|2FOFCWT                  F pdbx_FWT",
      "& PHWT|PH2FOFCWT               P pdbx_PHWT .3f",
      "? DELFWT|FOFCWT                F pdbx_DELFWT",
      "& DELPHWT|PHDELWT|PHFOFCWT     P pdbx_DELPHWT .3f",
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
      "? ROT     R pdbx_scan_angle",
      "$image      pdbx_image_id",
      nullptr
    };
    return for_merged ? merged : unmerged;
  }

  /// @brief Write MTZ reflection data to CIF format.
  /// @param mtz First MTZ dataset (required).
  /// @param mtz2 Optional second MTZ dataset for anomalous comparison.
  /// @param staraniso_b Optional Starraniso B-tensor to include.
  /// @param os Output stream for CIF file.
  void write_cif(const Mtz& mtz, const Mtz* mtz2,
                 SMat33<double>* staraniso_b, std::ostream& os);
  /// @brief Write XDS reflection data to CIF format.
  /// @param xds XDS_ASCII data to convert.
  /// @param os Output stream for CIF file.
  void write_cif_from_xds(const XdsAscii& xds, std::ostream& os) const;
};

/// @brief Write Starraniso B-tensor to mmCIF format.
/// @param b 3x3 symmetric B-tensor matrix.
/// @param entry_id Entry identifier for tags.
/// @param buf Temporary buffer for formatting.
/// @param os Output stream.
GEMMI_DLL void write_staraniso_b_in_mmcif(const SMat33<double>& b,
                                          const std::string& entry_id,
                                          char* buf, std::ostream& os);

/// @brief Remove '_dataset_name' appendix from MTZ column labels.
///
/// This suffix is sometimes added by CCP4i and needs removal for proper conversion.
/// @param mtz MTZ file to modify.
/// @param logger For reporting changes.
GEMMI_DLL void remove_appendix_from_column_names(Mtz& mtz, const Logger& logger);

/// @brief Validate merged MTZ has required columns for PDB deposition.
/// @param mtz MTZ to check.
/// @param logger For reporting results.
/// @return True if all required columns present.
GEMMI_DLL bool validate_merged_mtz_deposition_columns(const Mtz& mtz, const Logger& logger);

/// @brief Validate merged intensity data for consistency and quality.
///
/// Compares merged and unmerged intensity columns for anomalous differences and completeness.
/// Modifies both Intensities objects.
/// @param mi Merged intensities (modified).
/// @param ui Unmerged intensities (modified).
/// @param relaxed_check If true, apply looser validation criteria.
/// @param logger For reporting issues.
/// @return True if validation passes.
GEMMI_DLL bool validate_merged_intensities(Intensities& mi, Intensities& ui,
                                           bool relaxed_check, const Logger& logger);

/// @brief Extract software information from MTZ history records.
/// @param history Vector of history strings from MTZ file.
/// @return Vector of SoftwareItem objects describing processing steps.
GEMMI_DLL std::vector<SoftwareItem>
get_software_from_mtz_history(const std::vector<std::string>& history);

} // namespace gemmi
#endif
