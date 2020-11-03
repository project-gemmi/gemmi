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
#include <map>
#include <algorithm>     // for all_of
#include "mtz.hpp"       // for Mtz
#include "atox.hpp"      // for read_word
#include "merge.hpp"     // for Intensities, read_unmerged_intensities_from_mtz
#include "sprintf.hpp"   // for gf_snprintf, to_str
#include "version.hpp"   // for GEMMI_VERSION

namespace gemmi {

struct MtzToCif {
  // options that can be set directly
  std::vector<std::string> spec_lines; // conversion specification (cf. default_spec)
  const char* block_name = nullptr;  // NAME in data_NAME
  std::string mtz_path;              // path written in a comment
  bool with_comments = true;         // write comments
  bool skip_empty = false;           // skip reflections with no values
  bool enable_UB = false;            // write _diffrn_orient_matrix.UB
  bool enable_angle_phi = false;     // write phistt as _diffrn_refln.angle_phi
  std::string skip_empty_cols;       // columns used to determine "emptiness"
  double wavelength = NAN;           // user-specified wavelength
  int trim = 0;                      // output only reflections -N<=h,k,l<=N

  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "H H index_h",
      "K H index_k",
      "L H index_l",
      "? IMEAN|I J intensity_meas",
      "& SIGIMEAN|SIGI Q intensity_sigma",
      "? I(+)    K pdbx_I_plus",
      "& SIGI(+) M pdbx_I_plus_sigma",
      "? I(-)    K pdbx_I_minus",
      "& SIGI(-) M pdbx_I_minus_sigma",
      "? FP      F F_meas_au",   // TODO: FP from Refmac should show warning or error
      "& SIGFP   Q F_meas_sigma_au",
      "? F(+)    G pdbx_F_plus",
      "& SIGF(+) L pdbx_F_plus_sigma",
      "? F(-)    G pdbx_F_minus",
      "& SIGF(-) L pdbx_F_minus_sigma",
      "? FREE|RFREE|FreeR_flag I status",
      "? FWT|2FOFCWT      F pdbx_FWT",
      "& PHWT|PH2FOFCWT   P pdbx_PHWT",
      "? DELFWT|FOFCWT    F pdbx_DELFWT",
      "& DELPHWT|PHDELWT|PHFOFCWT P pdbx_DELPHWT",
      nullptr
    };
    static const char* unmerged[] = {
      // The first few columns of _diffrn_refln are set automatically:
      //   diffrn_id - based on BATCH
      //   standard_code - dummy mandatory column (set to 1)
      //   scale_group_code - dummy mandatory column (set to 1)
      //   id - reflection counter (1, 2, ...)
      //   angle_phi - (optional) phistt from batch header, written if phistt is changing
      "H H index_h",
      "K H index_k",
      "L H index_l",
      "? I       J intensity_net",
      "& SIGI    Q intensity_sigma .5g",
      nullptr
    };
    return for_merged ? merged : unmerged;
  }

  void write_cif(const Mtz& mtz, const Mtz* mtz2, std::ostream& os);

private:
  // describes which MTZ column is to be translated to what mmCIF column
  struct Trans {
    int col_idx;
    bool is_status = false;
    std::string tag;  // excluding category
    std::string format = "%g";
    int min_width = 0;
  };

  // data corresponding to one sweep (dataset) in unmerged MTZ file
  struct SweepData {
    int id;
    int batch_count;
    const Mtz::Batch* first_batch;
    const Mtz::Dataset* dataset;
  };

  std::vector<Trans> recipe;

  static std::vector<SweepData> gather_sweep_data(const Mtz& mtz) {
    std::vector<SweepData> data;
    for (const Mtz::Batch& batch : mtz.batches) {
      int sweep_id = batch.dataset_id();
      auto it = std::find_if(data.begin(), data.end(),
                 [&](const SweepData& sweep) { return sweep.id == sweep_id; });
      if (it == data.end())
        data.push_back({sweep_id, 1, &batch, nullptr});
      else
        it->batch_count++;
    }
    for (SweepData& sweep : data) {
      for (const Mtz::Dataset& d : mtz.datasets)
        if (d.id == sweep.id) {
          sweep.dataset = &d;
          break;
        }
      if (!sweep.dataset) {
        mtz.warn("Reference to absent dataset: " + std::to_string(sweep.id));
      }
    }
    return data;
  }

  // Get the first (non-zero) DWAVEL corresponding to an intensity, amplitude
  // or sigma column from the template.
  static double get_wavelength(const Mtz& mtz, const std::vector<Trans>& spec) {
    for (const Trans& tr : spec) {
      const Mtz::Column& col = mtz.columns.at(tr.col_idx);
      if (col.type == 'F' || col.type == 'J' ||
          col.type == 'K' || col.type == 'G' ||
          col.type == 'Q' || col.type == 'M' || col.type == 'L') { // sigma
        double wavelength = mtz.dataset(col.dataset_id).wavelength;
        if (wavelength != 0.)
          return wavelength;
      }
    }
    return 0.;
  }

  static int find_column_index(const std::string& column, const Mtz& mtz) {
    int idx = -1;
    for (const std::string& label : split_str(column, '|')) {
      for (size_t i = 0; i != mtz.columns.size(); ++i) {
        if (mtz.columns[i].label == label) {
          if (idx == -1)
            idx = (int) i;
          else
            mtz.warn("Column label duplicated: " + label);
        }
      }
      if (idx != -1)
        break;
    }
    return idx;
  }

  static int check_format(const std::string& fmt) {
    // expected format: [#_+-]?\d*(\.\d+)?[fFgGeEc]
    int min_width = 0;
    if (fmt.find('%') != std::string::npos)
      fail("Specify format without %. Got: " + fmt);
    const char* p = fmt.c_str();
    if (*p == '_' || *p == '+' || *p == '-' || *p == '#')
     ++p;
    if (is_digit(*p)) {
      min_width = *p++ - '0';
      if (is_digit(*p)) // two digits of width number max
        min_width = min_width * 10 + (*p++ - '0');
    }
    if (*p == '.' && is_digit(*(p+1))) {
      p += 2;
      if (is_digit(*p)) // two digits of precision numbers max
        ++p;
    }
    if (!std::isalpha(*p) || *(p+1) != '\0')
      fail("wrong format : " + fmt + "\nCorrect examples: g, .4f, 12.5e");
    char c = alpha_up(*p);
    if (c != 'F' && c != 'G' && c != 'E')
      fail("expected floating-point format, got: " + fmt);
    return min_width;
  }

  // state for parse_spec_line
  struct SpecParserState {
    size_t verified_spec_size = 0;
    bool discard_next_line = false;
  };

  void prepare_recipe(const Mtz& mtz) {
    recipe.clear();
    SpecParserState state;
    if (!spec_lines.empty()) {
      for (const std::string& line : spec_lines)
        parse_spec_line(line.c_str(), mtz, state);
    } else {
      const char** lines = default_spec(/*for_merged=*/mtz.batches.empty());
      for (; *lines != nullptr; ++lines)
        parse_spec_line(*lines, mtz, state);
    }
    if (recipe.empty())
      fail("empty translation recipe");
    for (size_t i = 0; i != recipe.size(); ++i)
      for (size_t j = i + 1; j != recipe.size(); ++j)
        if (recipe[i].tag == recipe[j].tag)
          fail("duplicated output tag: " + recipe[i].tag);
    // H, K, L must be the first columns in MTZ and are required in _refln
    for (int i = 2; i != -1; --i)
      if (!in_vector_f([&](const Trans& t) { return t.col_idx == i; }, recipe)) {
        Trans tr;
        tr.col_idx = i;
        tr.tag = "index_";
        tr.tag += ('h' + i); // h, k or l
        recipe.insert(recipe.begin(), tr);
      }
  }

  // adds results to recipe
  void parse_spec_line(const char* line, const Mtz& mtz, SpecParserState& state) {
    Trans tr;
    const char* p = line;
    if (*p == '&') {
      if (state.discard_next_line)
        return;
    } else {
      state.verified_spec_size = recipe.size();
      state.discard_next_line = false;
    }
    bool optional = (*p == '?' || *p == '&');
    if (optional)
      ++p;
    std::string column = read_word(p, &p);
    std::string type = read_word(p, &p);
    if (type.size() != 1)
      fail("Spec error: MTZ type '" + type + "' is not one character,"
           "\nin line: " + line);
    tr.tag = read_word(p, &p);
    if (tr.tag[0] == '_' || tr.tag.find('.') != std::string::npos)
      fail("Spec error: expected tag part after _refln., got: " +
           tr.tag + "\nin line: " + line);
    tr.col_idx = find_column_index(column, mtz);
    if (tr.col_idx == -1) {
      if (!optional)
        fail("Column not found: " + column);
      recipe.resize(state.verified_spec_size);
      state.discard_next_line = true;
      return;
    }
    const Mtz::Column& col = mtz.columns[tr.col_idx];
    if (type[0] != '*' && col.type != type[0])
      fail("Column " + col.label + " has type " +
           std::string(1, col.type) + " not " + type);
    tr.is_status = iequal(tr.tag, "status");
    std::string fmt = read_word(p, &p);
    if (!fmt.empty() && !tr.is_status) {
      tr.min_width = check_format(fmt);
      tr.format = "%" + fmt;
      if (tr.format[1] == '_')
        tr.format[1] = ' ';
    }
    recipe.push_back(tr);
  }

  void write_main_loop(const Mtz& mtz, char* buf, std::ostream& os);
};

// pre: unmerged MTZ is after switch_to_original_hkl()
inline void validate_merged_intensities(const Mtz& mtz1, const Mtz& mtz2, std::ostream& out) {
  if (mtz1.is_merged() == mtz2.is_merged())
    fail("both files are ", mtz1.is_merged() ? "merged" : "unmerged");
  const Mtz& merged = mtz1.is_merged() ? mtz1 : mtz2;
  const Mtz& unmerged = mtz1.is_merged() ? mtz2 : mtz1;
  Intensities in1 = read_unmerged_intensities_from_mtz(unmerged);
  size_t before = in1.data.size();
  in1.merge_in_place(/*output_plus_minus=*/false);  // it also sorts
  size_t after = in1.data.size();
  in1.remove_systematic_absences();
  Intensities in2 = read_mean_intensities_from_mtz(merged);
  in2.sort();
  if (in1.spacegroup != in2.spacegroup)
    out << "Warning: different space groups in two MTZ files:\n"
        << in1.spacegroup_str() << " and " << in2.spacegroup_str() << '\n';
  out << "Reflections: " << before << " -> " << after << " (merged) -> "
      << in1.data.size() << " (no sysabs)  vs  " << in2.data.size() << '\n';
  gemmi::Correlation corr;
  auto r1 = in1.data.begin();
  auto r2 = in2.data.begin();
  while (r1 != in1.data.end() && r2 != in2.data.end()) {
    if (r1->hkl == r2->hkl) {
      corr.add_point(r1->value, r2->value);
      ++r1;
      ++r2;
    } else if (std::tie(r1->hkl[0], r1->hkl[1], r1->hkl[2]) <
               std::tie(r2->hkl[0], r2->hkl[1], r2->hkl[2])) {
      ++r1;
    } else {
      ++r2;
    }
  }
  out << "IMEAN CC of " << corr.n << " values: " << 100 * corr.coefficient()
      << "% (mean ratio: " << corr.mean_ratio() << ")\n";
}

inline void MtzToCif::write_cif(const Mtz& mtz, const Mtz* mtz2, std::ostream& os) {
  if (mtz2 && mtz.is_merged() == mtz2->is_merged())
    fail("If two MTZ files are given, one must be merged and one unmerged,\n"
         "got two ", mtz.is_merged() ? "merged" : "unmerged");
  const Mtz* merged = mtz.is_merged() ? &mtz : mtz2;
  const Mtz* unmerged = mtz.is_merged() ? mtz2 : &mtz;

  char buf[256];
#define WRITE(...) os.write(buf, gf_snprintf(buf, 255, __VA_ARGS__))
  std::string id = ".";
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    if (!mtz_path.empty())
      os << "# from: " << mtz_path << '\n';
    for (const Mtz* m : {merged, unmerged})
      if (m) {
        os << "# title of " << (m->is_merged() ? "" : "un") << "merged MTZ: "
           << m->title << '\n';
        for (size_t i = 0; i != m->history.size(); ++i)
          os << "# MTZ history #" << i << ": " << m->history[i] << '\n';
      }
  }
  os << "data_" << (block_name ? block_name : "mtz");

  os << "\n\n_entry.id " << id << "\n\n";

  if (unmerged) {
    os << "_exptl_crystal.id 1\n\n";
    bool scaled = (unmerged->column_with_label("SCALEUSED") != nullptr);
    std::vector<SweepData> sweep_data = gather_sweep_data(*unmerged);

    os << "loop_\n_diffrn.id\n_diffrn.crystal_id\n_diffrn.details\n";
    for (const SweepData& sweep : sweep_data) {
      os << sweep.id << " 1 '";
      if (scaled)
        os << "scaled ";
      os << "unmerged data'\n";
    }
    os << '\n';

    os << "loop_\n"
          "_diffrn_measurement.diffrn_id\n"
          "_diffrn_measurement.details\n";
    for (const SweepData& sweep : sweep_data)
      os << sweep.id << " '" << sweep.batch_count << " frames'\n";
    os << '\n';

    os << "loop_\n"
          "_diffrn_radiation.diffrn_id\n"
          "_diffrn_radiation.wavelength_id\n";
    for (const SweepData& sweep : sweep_data)
      os << sweep.id << ' ' << sweep.id << '\n';
    os << '\n';

    os << "loop_\n"
          "_diffrn_radiation_wavelength.id\n"
          "_diffrn_radiation_wavelength.wavelength\n";
    for (const SweepData& sweep : sweep_data)
      os << sweep.id << ' '
         << (sweep.dataset ? to_str(sweep.dataset->wavelength) : " ")
         << '\n';
    os << '\n';

    if (enable_UB) {
      os << "loop_\n"
            "_diffrn_orient_matrix.diffrn_id\n"
            "_diffrn_orient_matrix.UB[1][1]\n"
            "_diffrn_orient_matrix.UB[1][2]\n"
            "_diffrn_orient_matrix.UB[1][3]\n"
            "_diffrn_orient_matrix.UB[2][1]\n"
            "_diffrn_orient_matrix.UB[2][2]\n"
            "_diffrn_orient_matrix.UB[2][3]\n"
            "_diffrn_orient_matrix.UB[3][1]\n"
            "_diffrn_orient_matrix.UB[3][2]\n"
            "_diffrn_orient_matrix.UB[3][3]\n";
      for (const SweepData& sweep : sweep_data) {
        Mat33 u = sweep.first_batch->matrix_U();
        Mat33 b = sweep.first_batch->get_cell().calculate_matrix_B();
        Mat33 ub = u.multiply(b);
        WRITE("%d  %#g %#g %#g  %#g %#g %#g  %#g %#g %#g\n", sweep.id,
              ub.a[0][0], ub.a[0][1], ub.a[0][2],
              ub.a[1][0], ub.a[1][1], ub.a[1][2],
              ub.a[2][0], ub.a[2][1], ub.a[2][2]);
      }
      os << '\n';
    }
  } else {
    double w = std::isnan(wavelength) ? get_wavelength(mtz, recipe) : wavelength;
    if (w != 0.)
      os << "_diffrn_radiation.diffrn_id 1\n"
         << "_diffrn_radiation.wavelength_id 1\n"
         << "_diffrn_radiation_wavelength.id 1\n"
         << "_diffrn_radiation_wavelength.wavelength " << to_str(w) << "\n\n";
  }

  const UnitCell& cell = mtz.get_cell();
  os << "_cell.entry_id " << id << '\n';
  WRITE("_cell.length_a    %8.3f\n", cell.a);
  WRITE("_cell.length_b    %8.3f\n", cell.b);
  WRITE("_cell.length_c    %8.3f\n", cell.c);
  WRITE("_cell.angle_alpha %8.3f\n", cell.alpha);
  WRITE("_cell.angle_beta  %8.3f\n", cell.beta);
  WRITE("_cell.angle_gamma %8.3f\n\n", cell.gamma);

  if (const SpaceGroup* sg = mtz.spacegroup) {
    os << "_symmetry.entry_id " << id << "\n"
          "_symmetry.space_group_name_H-M '" << sg->hm << "'\n"
          "_symmetry.Int_Tables_number " << sg->number << '\n';
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }

  if (merged)
    write_main_loop(*merged, buf, os);
  if (unmerged)
    write_main_loop(*unmerged, buf, os);
}

inline void MtzToCif::write_main_loop(const Mtz& mtz, char* buf, std::ostream& os) {
  prepare_recipe(mtz);
  // prepare indices
  std::vector<int> value_indices;  // used for --skip_empty
  std::vector<int> sigma_indices;  // used for status 'x'
  for (const Trans& tr : recipe) {
    const Mtz::Column& col = mtz.columns[tr.col_idx];
    if (skip_empty) {
      if (skip_empty_cols.empty() ? col.type != 'H' && col.type != 'I'
                                  : is_in_list(col.label, skip_empty_cols))
        value_indices.push_back(tr.col_idx);
    }
    if (col.type != 'Q' && col.type != 'L' && col.type != 'M')
      sigma_indices.push_back(tr.col_idx);
  }

  bool write_angle_phi = false;
  if (enable_angle_phi) {
    // don't write angle_phi if it has the same value in all batches
    for (size_t i = 1; i < mtz.batches.size(); ++i)
      if (mtz.batches[i].phi_start() != mtz.batches[0].phi_start()) {
        write_angle_phi = true;
        break;
      }
    }

  bool unmerged = !mtz.is_merged();
  os << "\nloop_\n";
  if (unmerged) {  // prepended tags that are not in the recipe
    os << "_diffrn_refln.diffrn_id\n"
          "_diffrn_refln.standard_code\n"
          "_diffrn_refln.scale_group_code\n"
          "_diffrn_refln.id\n";
    if (write_angle_phi) {
      os << "_diffrn_refln.angle_phi";
      if (with_comments)
        os << "                  # phistt from batch header";
      os << '\n';
    }
  }
  for (const Trans& tr : recipe) {
    const Mtz::Column& col = mtz.columns.at(tr.col_idx);
    const Mtz::Dataset& ds = mtz.dataset(col.dataset_id);
    os << (unmerged ? "_diffrn_refln." : "_refln.");
    if (with_comments) {
      // dataset is assigned to column only in merged MTZ
      if (unmerged)
        WRITE("%-26s # %s\n", tr.tag.c_str(), col.label.c_str());
      else
        WRITE("%-26s # %-14s from dataset %s\n",
              tr.tag.c_str(), col.label.c_str(), ds.dataset_name.c_str());
    } else {
      os << tr.tag << '\n';
    }
  }
  if (unmerged) // appended tag that is not in the recipe
    // ccp4_centroid_of_image_numbers is a proposed real number corresponding
    // to ITEM_ZD in XDS, BATCH (integer, with subtracted offset) in MTZ, etc.
    os << "_diffrn_refln.ccp4_centroid_of_image_numbers\n";
  int batch_idx = find_column_index("BATCH", mtz);
  std::map<int, const Mtz::Batch*> batch_by_number;
  for (const Mtz::Batch& b : mtz.batches)
    batch_by_number.emplace(b.number, &b);

  // prepare offsets
  std::map<int, int> batch_offset;
  for (const Mtz::Batch& b : mtz.batches) {
    auto result = batch_offset.emplace(b.dataset_id(), b.number);
    if (!result.second && result.first->second < b.number)
      result.first->second = b.number;
  }
  for (auto& it : batch_offset)
    it.second -= it.second % 1000;

  for (int i = 0, idx = 0; i != mtz.nreflections; ++i) {
    const float* row = &mtz.data[i * mtz.columns.size()];
    if (trim > 0) {
      if (row[0] < -trim || row[0] > trim ||
          row[1] < -trim || row[1] > trim ||
          row[2] < -trim || row[2] > trim) {
        continue;
      }
    }
    if (!value_indices.empty())
      if (std::all_of(value_indices.begin(), value_indices.end(),
                      [&](int n) { return std::isnan(row[n]); }))
        continue;
    int batch_number = 0;
    if (unmerged) {
      if (batch_idx == -1)
        fail("BATCH column not found");
      batch_number = (int) row[batch_idx];
      auto it = batch_by_number.find(batch_number);
      if (it == batch_by_number.end())
        fail("unexpected values in column BATCH");
      const Mtz::Batch& batch = *it->second;
      // values for prepended tags
      os << batch.dataset_id() << " . . " << ++idx << ' ';
      if (write_angle_phi)
        WRITE("%.2f ", batch.phi_start());
      batch_number -= batch_offset.at(batch.dataset_id());
    }
    bool first = true;
    for (const Trans& tr : recipe) {
      if (first)
        first = false;
      else
        os << ' ';
      float v = row[tr.col_idx];
      if (tr.is_status) {
        char status = 'x';
        if (sigma_indices.empty() ||
            !std::all_of(sigma_indices.begin(), sigma_indices.end(),
                         [&](int n) { return std::isnan(row[n]); }))
          status = v == 0. ? 'f' : 'o';
        os << status;
      } else if (std::isnan(v)) {
        for (int j = 1; j < tr.min_width; ++j)
          os << ' ';
        os << '?';
      } else {
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
        WRITE(tr.format.c_str(), v);
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif
      }
    }
    if (unmerged) // value for appended tag
      os << ' ' << batch_number;
    os << '\n';
  }
}
#undef WRITE

} // namespace gemmi
#endif
