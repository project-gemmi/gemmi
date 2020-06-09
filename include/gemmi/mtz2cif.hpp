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
#include "sprintf.hpp"   // for gf_snprintf, to_str
#include "version.hpp"   // for GEMMI_VERSION

namespace gemmi {

struct MtzToCif {
  // options that can be set directly
  std::vector<std::string> spec_lines;  // conversion specification (cf. default_spec)
  const char* block_name = nullptr;     // NAME in data_NAME
  std::string mtz_path;                 // path written in a comment
  bool with_comments = true;            // write comments
  bool skip_empty = false;              // skip reflections with no values
  std::string skip_empty_cols;          // columns used to determine "emptiness"
  double wavelength = NAN;              // user-specified wavelength
  int trim = 0;                         // output only reflections -N<=h,k,l<=N

  static const char** default_spec(bool for_merged) {
    static const char* merged[] = {
      "H H index_h",
      "K H index_k",
      "L H index_l",
      "? I       J intensity_meas",
      "& SIGI    Q intensity_sigma",
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

  void write_cif(const Mtz& mtz, std::ostream& os);

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
};

inline void MtzToCif::write_cif(const Mtz& mtz, std::ostream& os) {
  recipe.clear();
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

  char buf[256];
#define WRITE(...) os.write(buf, gf_snprintf(buf, 255, __VA_ARGS__))
  std::string id = ".";
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    if (!mtz_path.empty())
      os << "# from: " << mtz_path << '\n';
    os << "# MTZ title: " << mtz.title << '\n';
    for (size_t i = 0; i != mtz.history.size(); ++i)
      os << "# MTZ history #" << i << ": " << mtz.history[i] << '\n';
  }
  const bool unmerged = !mtz.batches.empty();
  os << "data_"
     << (block_name ? block_name : unmerged ? "unmerged" : "merged")
     << "\n\n_entry.id " << id << "\n\n";

  bool write_angle_phi = false;
  if (unmerged) {
    // don't write angle_phi if it has the same value in all batches
    for (size_t i = 1; i != mtz.batches.size(); ++i)
      if (mtz.batches[i].phi_start() != mtz.batches[0].phi_start()) {
        write_angle_phi = true;
        break;
      }
    os << "_exptl_crystal.id 1\n\n";
    bool scaled = (mtz.column_with_label("SCALEUSED") != nullptr);
    std::vector<SweepData> sweep_data = gather_sweep_data(mtz);

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

  if (!unmerged) {
    double w = std::isnan(wavelength) ? get_wavelength(mtz, recipe) : wavelength;
    if (w != 0.)
      os << "_diffrn_radiation.diffrn_id 1\n"
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
          "_symmetry.Int_Tables_number " << sg->number << "\n\n";
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }
  os << "loop_\n";
  if (unmerged) {
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
  int batch_idx = find_column_index("BATCH", mtz);
  std::map<int, const Mtz::Batch*> batch_by_number;
  for (const Mtz::Batch& b : mtz.batches)
    batch_by_number.emplace(b.number, &b);
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
    if (unmerged) {
      if (batch_idx == -1)
        fail("BATCH column not found");
      int batch_number = (int) row[batch_idx];
      auto it = batch_by_number.find(batch_number);
      if (it == batch_by_number.end())
        fail("unexpected values in column BATCH");
      const Mtz::Batch& batch = *it->second;
      os << batch.dataset_id() << " . . " << ++idx << ' ';
      if (write_angle_phi)
        WRITE("%.2f ", batch.phi_start());
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
    os << '\n';
  }
#undef WRITE
}

} // namespace gemmi
#endif
