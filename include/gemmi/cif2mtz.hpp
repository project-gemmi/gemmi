// Copyright 2021 Global Phasing Ltd.
//
// A class for converting SF-mmCIF to MTZ (merged or unmerged).

#ifndef GEMMI_CIF2MTZ_HPP_
#define GEMMI_CIF2MTZ_HPP_

#include <ostream>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include "cifdoc.hpp"   // for Loop, as_int, ...
#include "fail.hpp"     // for fail
#include "mtz.hpp"      // for Mtz
#include "numb.hpp"     // for as_number
#include "refln.hpp"    // for ReflnBlock
#include "version.hpp"  // for GEMMI_VERSION

namespace gemmi {

template<typename DataProxy>
DataType check_data_type_under_symmetry(const DataProxy& proxy) {
  const SpaceGroup* sg = proxy.spacegroup();
  if (!sg)
    return DataType::Unknown;
  std::unordered_map<Op::Miller, int, MillerHash> seen;
  ReciprocalAsu asu(sg);
  GroupOps gops = sg->operations();
  bool centric = gops.is_centrosymmetric();
  DataType data_type = DataType::Mean;
  for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
    auto hkl_isym = asu.to_asu(proxy.get_hkl(i), gops);
    int sign = hkl_isym.second % 2 + 1;  // 2=positive, 1=negative
    auto r = seen.emplace(hkl_isym.first, sign);
    if (!r.second) {
      if ((r.first->second & sign) != 0 || centric)
        return DataType::Unmerged;
      r.first->second |= sign;
      data_type = DataType::Anomalous;
    }
  }
  return data_type;
}


struct CifToMtz {
  // Alternative mmCIF tags for the same MTZ label should be consecutive
  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "pdbx_r_free_flag FreeR_flag I 0",
      "status FreeR_flag I 0 o=1,f=0",
      "intensity_meas IMEAN J 1",
      "intensity_sigma SIGIMEAN Q 1",
      "pdbx_I_plus I(+) K 1",
      "pdbx_I_plus_sigma SIGI(+) M 1",
      "pdbx_I_minus I(-) K 1",
      "pdbx_I_minus_sigma SIGI(-) M 1",
      "F_meas FP F 1",
      "F_meas_au FP F 1",
      "F_meas_sigma SIGFP Q 1",
      "F_meas_sigma_au SIGFP Q 1",
      "pdbx_F_plus F(+) G 1",
      "pdbx_F_plus_sigma SIGF(+) L 1",
      "pdbx_F_minus F(-) G 1",
      "pdbx_F_minus_sigma SIGF(-) L 1",
      "pdbx_anom_difference DP D 1",
      "pdbx_anom_difference_sigma SIGDP Q 1",
      "F_calc FC F 1",
      "F_calc_au FC F 1",
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
      "intensity_meas I J 0",  // for unmerged data is category refln
      "intensity_net I J 0",
      "intensity_sigma SIGI Q 0",
      "pdbx_detector_x XDET R 0",
      "pdbx_detector_y YDET R 0",
      "pdbx_scan_angle ROT R 0",
      nullptr
    };
    return for_merged ? merged : unmerged;
  };

  struct Entry {
    std::string refln_tag;
    std::string col_label;
    char col_type;
    int dataset_id;
    std::vector<std::pair<std::string, float>> code_to_number;

    Entry(const std::string& line) {
      std::vector<std::string> tokens;
      tokens.reserve(4);
      split_str_into_multi(line, " \t\r\n", tokens);
      if (tokens.size() != 4 && tokens.size() != 5)
        fail("line should have 4  or 5 words: " + line);
      if (tokens[2].size() != 1 || tokens[3].size() != 1 ||
          (tokens[3][0] != '0' && tokens[3][0] != '1'))
        fail("incorrect line: " + line);
      refln_tag = tokens[0];
      col_label = tokens[1];
      col_type = tokens[2][0];
      dataset_id = tokens[3][0] - '0';
      // for compatibility with older spec
      if (col_type == 's' && tokens.size() == 4) {
        col_type = 'I';
        tokens.push_back("o=1,f=0");
      }
      if (tokens.size() == 5) {
        std::vector<std::string> items = split_str(tokens[4], ',');
        code_to_number.reserve(items.size());
        for (const std::string& item : items) {
          size_t pos = item.find('=');
          if (pos == std::string::npos)
            fail("wrong mapping (", item, ") in: ", line);
          float f;
          auto result = fast_float::from_chars(item.c_str() + pos + 1,
                                               item.c_str() + item.size(), f);
          if (result.ec != std::errc())
            fail("failed to parse value in ", item, " in: ", line);
          code_to_number.emplace_back(item.substr(0, pos), f);
        }
      }
    }

    float translate_code_to_number(const std::string& v) const {
      if (v.size() == 1) {
        for (const auto& c2n : code_to_number)
          if (c2n.first.size() == 1 && c2n.first[0] == v[0])
            return c2n.second;
      } else {
        std::string s = cif::as_string(v);
        for (const auto& c2n : code_to_number)
          if (c2n.first == s)
            return c2n.second;
      }
      return NAN;
    }
  };

  bool verbose = false;
  bool force_unmerged = false;
  std::string title;
  std::vector<std::string> history = { "From gemmi-cif2mtz " GEMMI_VERSION };
  std::vector<std::string> spec_lines;

  Mtz convert_block_to_mtz(const ReflnBlock& rb, std::ostream& out) const {
    Mtz mtz;
    mtz.title = title.empty() ? "Converted from mmCIF block " + rb.block.name : title;
    if (!history.empty()) {
      mtz.history.reserve(mtz.history.size() + history.size());
      mtz.history.insert(mtz.history.end(), history.begin(), history.end());
    }
    mtz.cell = rb.cell;
    mtz.spacegroup = rb.spacegroup;
    mtz.add_dataset("HKL_base");
    const cif::Loop* loop = rb.refln_loop ? rb.refln_loop : rb.diffrn_refln_loop;
    if (!loop)
      fail("_refln category not found in mmCIF block: " + rb.block.name);
    if (verbose)
      out << "Searching tags with known MTZ equivalents ...\n";
    std::vector<int> indices;
    std::vector<const Entry*> entries;  // used for code_to_number only
    std::string tag = loop->tags[0];
    const size_t tag_offset = rb.tag_offset();

    bool unmerged = force_unmerged || !rb.refln_loop;

    std::vector<Entry> spec_entries;
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
    tag.replace(tag_offset, std::string::npos, "index_h");
    for (char c : {'h', 'k', 'l'}) {
      tag.back() = c;
      int index = loop->find_tag(tag);
      if (index == -1)
        fail("Miller index tag not found: " + tag);
      indices.push_back(index);
      entries.push_back(nullptr);
      auto col = mtz.columns.emplace(mtz.columns.end());
      col->dataset_id = 0;
      col->type = 'H';
      col->label = alpha_up(c);
    }

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
    }

    if (!unmerged)
      mtz.add_dataset("unknown").wavelength = rb.wavelength;

    // other columns according to the spec
    bool column_added = false;
    for (const Entry& entry : spec_entries) {
      tag.replace(tag_offset, std::string::npos, entry.refln_tag);
      int index = loop->find_tag(tag);
      if (index == -1)
        continue;
      if (mtz.column_with_label(entry.col_label))
        continue;
      column_added = true;
      indices.push_back(index);
      entries.push_back(entry.code_to_number.empty() ? nullptr : &entry);
      auto col = mtz.columns.emplace(mtz.columns.end());
      // dataset_id is meaningless in unmerged MTZ files
      col->dataset_id = unmerged ? 0 : entry.dataset_id;
      col->type = entry.col_type;
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

    std::unique_ptr<UnmergedHklMover> hkl_mover;
    struct BatchInfo {
      int sweep_id;
      int frame_id;
    };
    std::vector<BatchInfo> batch_nums;
    if (unmerged) {
      hkl_mover.reset(new UnmergedHklMover(mtz.spacegroup));
      tag.replace(tag_offset, std::string::npos, "diffrn_id");
      int sweep_id_index = loop->find_tag(tag);
      if (sweep_id_index == -1 && verbose)
        out << "No diffrn_id. Assuming a single sweep.\n";
      tag.replace(tag_offset, std::string::npos, "pdbx_image_id");
      int image_id_index = loop->find_tag(tag);
      if (verbose) {
        if (image_id_index == -1)
          out << "No pdbx_image_id, setting BATCH to a dummy value.\n";
        else
          out << "  " << tag << " & diffrn_id -> BATCH\n";
      }
      struct SweepInfo {
        std::set<int> frame_ids;
        int offset;
        float wavelength = 0.f;
        std::string crystal_id = "unknown";
        int dataset_id;
      };
      std::map<int, SweepInfo> sweeps;
      cif::Block& block = const_cast<cif::Block&>(rb.block);
      cif::Table tab_w0 = block.find("_diffrn.", {"id", "crystal_id"});
      cif::Table tab_w1 = block.find("_diffrn_radiation.",
                                     {"diffrn_id", "wavelength_id"});
      cif::Table tab_w2 = block.find("_diffrn_radiation_wavelength.",
                                     {"id", "wavelength"});
      // store sweep and frame numbers corresponding to reflections
      batch_nums.reserve(loop->length());
      for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
        int sweep_id = 1;
        if (sweep_id_index >= 0)
          sweep_id = cif::as_int(loop->values[i + sweep_id_index], 1);
        int frame = 1;
        if (image_id_index >= 0) {
          double d = cif::as_number(loop->values[i + image_id_index], 1.);
          frame = (int) std::ceil(d);
        }
        batch_nums.push_back({sweep_id, frame});
        if (frame >= 0) {
          SweepInfo& sweep = sweeps[sweep_id];
          // if new sweep was added - try to set crystal_id and wavelength
          if (sweep.frame_ids.empty() && sweep_id_index >= 0) {
            const std::string& sweep_str = loop->values[i + sweep_id_index];
            try {
              sweep.crystal_id = tab_w0.find_row(sweep_str).str(1);
            } catch(std::exception&) {}
            try {
              const std::string& wave_id = tab_w1.find_row(sweep_str)[1];
              const std::string& wavelen = tab_w2.find_row(wave_id)[1];
              sweep.wavelength = (float) cif::as_number(wavelen, 0.);
            } catch(std::exception&) {}
          }
          sweep.frame_ids.insert(frame);
        }
      }

      // add datasets and set SweepInfo::dataset_id
      for (auto& sweep_pair : sweeps) {
        SweepInfo& si = sweep_pair.second;
        Mtz::Dataset* ds = mtz.dataset_with_name(si.crystal_id);
        if (ds == nullptr) {
          ds = &mtz.add_dataset(si.crystal_id);
          // both crystal name and dataset name are set to crystal_id
          ds->project_name = "unknown";
          ds->wavelength = si.wavelength;
        }
        si.dataset_id = ds->id;
      }

      // determine offset that makes frame numbers unique
      int cap = 0;
      for (auto& it : sweeps) {
        it.second.offset = cap;
        cap += *--it.second.frame_ids.end() + 1100;
        cap -= cap % 1000;
      }
      // add offset to BatchInfo::frame_id
      for (BatchInfo& p : batch_nums)
        if (p.sweep_id >= 0 && p.frame_id >= 0)
          p.frame_id += sweeps.at(p.sweep_id).offset;

      // add MTZ batches
      for (const auto& sweep_pair : sweeps) {
        const SweepInfo& sweep_info = sweep_pair.second;
        Mtz::Batch batch;
        batch.set_dataset_id(sweep_info.dataset_id);
        batch.set_cell(mtz.cell);
        batch.set_wavelength(sweep_info.wavelength);
        // FIXME should we set more properties in BATCH header?
        int min_frame = *sweep_info.frame_ids.begin();
        int max_frame = *--sweep_info.frame_ids.end();
        if (2 * sweep_info.frame_ids.size() > size_t(max_frame - min_frame)) {
          // probably consecutive range, even if some frames are missing
          for (int n = min_frame; n <= max_frame; ++n) {
            batch.number = n + sweep_info.offset;
            mtz.batches.push_back(batch);
          }
        } else {
          for (int n : sweep_info.frame_ids) {
            batch.number = n + sweep_info.offset;
            mtz.batches.push_back(batch);
          }
        }
      }
    }  // - if (unmerged)

    // fill in the data
    mtz.data.resize(mtz.columns.size() * mtz.nreflections);
    size_t k = 0, row = 0;
    for (size_t i = 0; i < loop->values.size(); i += loop->tags.size()) {
      if (unmerged) {
        std::array<int, 3> hkl;
        for (size_t ii = 0; ii != 3; ++ii)
          hkl[ii] = cif::as_int(loop->values[i + indices[ii]]);
        int isym = hkl_mover->move_to_asu(hkl);
        for (size_t j = 0; j != 3; ++j)
          mtz.data[k++] = (float) hkl[j];
        mtz.data[k++] = (float) isym;
        mtz.data[k++] = batch_nums.empty() ? 1.f : (float) batch_nums[row++].frame_id;
      } else {
        for (size_t j = 0; j != 3; ++j)
          mtz.data[k++] = (float) cif::as_int(loop->values[i + indices[j]]);
      }
      for (size_t j = 3; j != indices.size(); ++j) {
        const std::string& v = loop->values[i + indices[j]];
        if (cif::is_null(v)) {
          mtz.data[k] = (float) NAN;
        } else if (entries[j] != nullptr) {
          mtz.data[k] = entries[j]->translate_code_to_number(v);
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
