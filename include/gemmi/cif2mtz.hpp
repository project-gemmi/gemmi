// Copyright 2021 Global Phasing Ltd.
//
// A class for converting SF-mmCIF to MTZ (merged or unmerged).

#ifndef GEMMI_CIF2MTZ_HPP_
#define GEMMI_CIF2MTZ_HPP_

#include <ostream>
#include "cifdoc.hpp"   // for Loop, as_int, ...
#include "fail.hpp"     // for fail
#include "mtz.hpp"      // for Mtz
#include "numb.hpp"     // for as_number
#include "refln.hpp"    // for ReflnBlock

namespace gemmi {

struct CifToMtz {
  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "pdbx_r_free_flag FreeR_flag I 0",
      "status FreeR_flag s 0", // s is a special flag
      "intensity_meas IMEAN J 1",
      "intensity_sigma SIGIMEAN Q 1",
      "pdbx_I_plus I(+) K 1",
      "pdbx_I_plus_sigma SIGI(+) M 1",
      "pdbx_I_minus I(-) K 1",
      "pdbx_I_minus_sigma SIGI(-) M 1",
      "F_meas_au FP F 1",
      "F_meas_sigma_au SIGFP Q 1",
      "pdbx_F_plus F(+) G 1",
      "pdbx_F_plus_sigma SIGF(+) L 1",
      "pdbx_F_minus F(-) G 1",
      "pdbx_F_minus_sigma SIGF(-) L 1",
      "pdbx_anom_difference DP D 1",
      "pdbx_anom_difference_sigma SIGDP Q 1",
      "F_calc FC F 1",
      "phase_calc PHIC P 1",
      "fom FOM W 1",
      "weight FOM W 1",
      "pdbx_HL_A_iso HLA A 1",
      "pdbx_HL_B_iso HLB A 1",
      "pdbx_HL_C_iso HLC A 1",
      "pdbx_HL_D_iso HLD A 1",
      "pdbx_FWT FWT F 1",
      "pdbx_PHWT PHWT P 1",
      "pdbx_DELFWT DELFWT F 1",
      "pdbx_DELPHWT DELPHWT P 1",
      nullptr
    };
    static const char* unmerged[] = {
      "intensity_meas I J 1",  // for unmerged data is category refln
      "intensity_net I J 1",
      "intensity_sigma SIGI Q 1",
      "pdbx_detector_x XDET R 1",
      "pdbx_detector_y YDET R 1",
      "pdbx_scan_angle ROT R 1",
      nullptr
    };
    return for_merged ? merged : unmerged;
  };

  struct Entry {
    std::string refln_tag;
    std::string col_label;
    char col_type;
    int dataset_id;

    Entry(const std::string& line) {
      std::vector<std::string> tokens;
      tokens.reserve(4);
      split_str_into_multi(line, " \t\r\n", tokens);
      if (tokens.size() != 4)
        fail("line should have 4 words: " + line);
      if (tokens[2].size() != 1 || tokens[3].size() != 1 ||
          (tokens[3][0] != '0' && tokens[3][0] != '1'))
        fail("incorrect line: " + line);
      refln_tag = tokens[0];
      col_label = tokens[1];
      col_type = tokens[2][0];
      dataset_id = tokens[3][0] - '0';
    }
  };

  // Alternative mmCIF tags for the same MTZ label should be consecutive
  std::vector<Entry> spec_entries;
  bool verbose = false;
  bool force_unmerged = false;
  const char* title = nullptr;
  std::vector<std::string> history;
  std::vector<std::string> spec_lines;

  Mtz convert_block_to_mtz(const ReflnBlock& rb, std::ostream& out) {
    Mtz mtz;
    if (title)
      mtz.title = title;
    if (!history.empty()) {
      mtz.history.reserve(mtz.history.size() + history.size());
      mtz.history.insert(mtz.history.end(), history.begin(), history.end());
    }
    mtz.cell = rb.cell;
    mtz.spacegroup = rb.spacegroup;
    mtz.add_dataset("HKL_base");
    mtz.add_dataset("unknown").wavelength = rb.wavelength;
    const cif::Loop* loop = rb.refln_loop ? rb.refln_loop : rb.diffrn_refln_loop;
    if (!loop)
      fail("_refln category not found in mmCIF block: " + rb.block.name);
    if (verbose)
      out << "Searching tags with known MTZ equivalents ...\n";
    bool uses_status = false;
    std::vector<int> indices;
    std::string tag = loop->tags[0];
    const size_t len = tag.find('.') + 1;

    bool unmerged = force_unmerged || !rb.refln_loop;

    if (!spec_lines.empty()) {
      spec_entries.reserve(spec_lines.size());
      for (const std::string& line : spec_lines)
        spec_entries.emplace_back(line);
    } else {
      const char** line = default_spec(!unmerged);
      for (; *line != nullptr; ++line)
        spec_entries.emplace_back(*line);
    }

    // always start with H, K, L
    tag.replace(len, std::string::npos, "index_h");
    for (char c : {'h', 'k', 'l'}) {
      tag.back() = c;
      int index = loop->find_tag(tag);
      if (index == -1)
        fail("Miller index tag not found: " + tag);
      indices.push_back(index);
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'H';
      col->label = alpha_up(c);
    }

    std::unique_ptr<UnmergedHklMover> hkl_mover;
    // M/ISYM and BATCH
    if (unmerged) {
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'Y';
      col->label = "M/ISYM";

      col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'B';
      col->label = "BATCH";

      mtz.batches.emplace_back();
      mtz.batches.back().set_cell(mtz.cell);
      hkl_mover.reset(new UnmergedHklMover(mtz.spacegroup));
    }

    // other columns according to the spec
    bool column_added = false;
    for (const Entry& entry : spec_entries) {
      tag.replace(len, std::string::npos, entry.refln_tag);
      int index = loop->find_tag(tag);
      if (index == -1)
        continue;
      if (mtz.column_with_label(entry.col_label))
        continue;
      column_added = true;
      indices.push_back(index);
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = entry.dataset_id;
      col->type = entry.col_type;
      if (col->type == 's') {
        col->type = 'I';
        uses_status = true;
      }
      col->label = entry.col_label;
      if (verbose)
        out << "  " << tag << " -> " << col->label << '\n';
    }
    if (!column_added)
      fail(force_unmerged ? "Unmerged d" : "D", "ata not found in block ", rb.block.name);

    for (size_t i = 0; i != mtz.columns.size(); ++i) {
      mtz.columns[i].parent = &mtz;
      mtz.columns[i].idx = i;
    }
    mtz.nreflections = (int) loop->length();

    if (unmerged) {
      if (verbose)
        out << "The BATCH columns is set to dummy value 1.\n";
    }

    // fill in the data
    mtz.data.resize(mtz.columns.size() * mtz.nreflections);
    int k = 0;
    for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
      size_t j = 0;
      if (unmerged) {
        std::array<int, 3> hkl;
        for (int ii = 0; ii != 3; ++ii)
          hkl[ii] = cif::as_int(loop->values[i + indices[ii]]);
        int isym = hkl_mover->move_to_asu(hkl);
        for (; j != 3; ++j)
          mtz.data[k++] = (float) hkl[j];
        mtz.data[k++] = (float) isym;
        mtz.data[k++] = 1.0f; // batch number TODO "pdbx_image_id",
      } else {
        for (; j != 3; ++j)
          mtz.data[k++] = (float) cif::as_int(loop->values[i + indices[j]]);
      }
      if (uses_status)
        mtz.data[k++] = status_to_freeflag(loop->values[i + indices[j++]]);
      for (; j != indices.size(); ++j) {
        const std::string& v = loop->values[i + indices[j]];
        if (cif::is_null(v)) {
          mtz.data[k] = (float) NAN;
        } else {
          mtz.data[k] = (float) cif::as_number(v);
          if (std::isnan(mtz.data[k]))
            out << "Value #" << i + indices[j] << " in the loop is not a number: "
                << v << '\n';
        }
        ++k;
      }
    }
    return mtz;
  }

private:
  static float status_to_freeflag(const std::string& str) {
    char c = str[0];
    if (c == '\'' || c == '"')
      c = str[1];
    if (c == 'o')
      return 1.f;
    if (c == 'f')
      return 0.f;
    return NAN;
  }
};

} // namespace gemmi
#endif
