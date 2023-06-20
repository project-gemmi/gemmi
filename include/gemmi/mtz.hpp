// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cassert>
#include <cstdint>       // for int32_t
#include <cstring>       // for memcpy
#include <cmath>         // for isnan
#include <algorithm>     // for sort, any_of
#include <array>
#include <initializer_list>
#include <ostream>
#include <string>
#include <vector>
#include "atox.hpp"      // for simple_atoi, read_word
#include "atof.hpp"      // for fast_atof
#include "input.hpp"     // for FileStream, CharArray
#include "iterator.hpp"  // for StrideIter
#include "fail.hpp"      // for fail
#include "fileutil.hpp"  // for file_open, is_little_endian, fileptr_t, ...
#include "math.hpp"      // for rad, Mat33
#include "symmetry.hpp"  // for find_spacegroup_by_name, SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for ialpha4_id, rtrim_str, ialpha3_id, ...

namespace gemmi {

template <typename T, typename FP=typename std::iterator_traits<T>::value_type>
std::array<FP,2> calculate_min_max_disregarding_nans(T begin, T end) {
  std::array<FP,2> minmax = {{NAN, NAN}};
  T i = begin;
  while (i != end && std::isnan(*i))
    ++i;
  if (i != end) {
    minmax[0] = minmax[1] = *i;
    while (++i != end) {
      if (*i < minmax[0])
        minmax[0] = *i;
      else if (*i > minmax[1])
        minmax[1] = *i;
    }
  }
  return minmax;
}

// Unmerged MTZ files always store in-asu hkl indices and symmetry operation
// encoded in the M/ISYM column. Here is a helper for writing such files.
struct UnmergedHklMover {
  UnmergedHklMover(const SpaceGroup* spacegroup) : asu_(spacegroup) {
    if (spacegroup)
      group_ops_ = spacegroup->operations();
  }

  // Modifies hkl and returns ISYM value for M/ISYM
  int move_to_asu(std::array<int, 3>& hkl) {
    std::pair<Miller, int> hkl_isym = asu_.to_asu(hkl, group_ops_);
    hkl = hkl_isym.first;
    return hkl_isym.second;
  }

private:
  ReciprocalAsu asu_;
  GroupOps group_ops_;
};


struct GEMMI_DLL Mtz {
  struct Dataset {
    int id;
    std::string project_name;
    std::string crystal_name;
    std::string dataset_name;
    UnitCell cell;
    double wavelength;  // 0 means not set
  };

  struct Column {
    int dataset_id;
    char type;
    std::string label;
    float min_value = NAN;
    float max_value = NAN;
    std::string source;  // from COLSRC
    Mtz* parent;
    std::size_t idx;

    Dataset& dataset() { return parent->dataset(dataset_id); }
    const Dataset& dataset() const { return parent->dataset(dataset_id); }
    bool has_data() const { return parent->has_data(); }
    int size() const { return has_data() ? parent->nreflections : 0; }
    size_t stride() const { return parent->columns.size(); }
    float& operator[](std::size_t n) { return parent->data[idx + n * stride()]; }
    float operator[](std::size_t n) const { return parent->data[idx + n * stride()]; }
    float& at(std::size_t n) { return parent->data.at(idx + n * stride()); }
    float at(std::size_t n) const { return parent->data.at(idx + n * stride()); }
    bool is_integer() const {
      return type == 'H' || type == 'B' || type == 'Y' || type == 'I';
    }
    using iterator = StrideIter<float>;
    iterator begin() {
      assert(parent);
      assert(&parent->columns[idx] == this);
      return iterator({parent->data.data(), idx, stride()});
    }
    iterator end() {
      return iterator({parent->data.data() + parent->data.size(), idx,
                       stride()});
    }
    using const_iterator = StrideIter<const float>;
    const_iterator begin() const { return const_cast<Column*>(this)->begin(); }
    const_iterator end() const { return const_cast<Column*>(this)->end(); }
  };

  struct Batch {
    Batch() {
      ints.resize(29, 0);
      floats.resize(156, 0.);
      // write the same values that are written by CCP4 progs such as COMBAT
      ints[0] = 29 + 156;
      ints[1] = 29;
      ints[2] = 156;
      // COMBAT sets BSCALE=1, but Pointless sets it to 0.
      //floats[43] = 1.f; // batch scale
    }
    int number;
    std::string title;
    std::vector<int> ints;
    std::vector<float> floats;
    std::vector<std::string> axes;

    UnitCell get_cell() const {
      return UnitCell(floats[0], floats[1], floats[2],
                      floats[3], floats[4], floats[5]);
    }
    void set_cell(const UnitCell& uc) {
      floats[0] = (float) uc.a;
      floats[1] = (float) uc.b;
      floats[2] = (float) uc.c;
      floats[3] = (float) uc.alpha;
      floats[4] = (float) uc.beta;
      floats[5] = (float) uc.gamma;
    }

    int dataset_id() const { return ints[20]; }
    void set_dataset_id(int id) { ints[20] = id; }
    float wavelength() const { return floats[86]; }
    void set_wavelength(float lambda) { floats[86] = lambda; }
    float phi_start() const { return floats[36]; }
    float phi_end() const { return floats[37]; }
    Mat33 matrix_U() const {
      return Mat33(floats[6], floats[9],  floats[12],
                   floats[7], floats[10], floats[13],
                   floats[8], floats[11], floats[14]);
    }
  };

  std::string source_path;  // input file path, if known
  bool same_byte_order = true;
  bool indices_switched_to_original = false;
  std::int64_t header_offset = 0;
  std::string version_stamp;
  std::string title;
  int nreflections = 0;
  std::array<int, 5> sort_order = {};
  double min_1_d2 = NAN;
  double max_1_d2 = NAN;
  float valm = NAN;
  int nsymop = 0;
  UnitCell cell;
  int spacegroup_number = 0;
  std::string spacegroup_name;
  std::vector<Op> symops;
  const SpaceGroup* spacegroup = nullptr;
  std::vector<Dataset> datasets;
  std::vector<Column> columns;
  std::vector<Batch> batches;
  std::vector<std::string> history;
  std::string appended_text;
  std::vector<float> data;

  // stream used for warnings when reading mtz file (and also in mtz2cif)
  std::ostream* warnings = nullptr;

  explicit Mtz(bool with_base=false) {
    if (with_base)
      add_base();
  }
  Mtz(Mtz&& o) noexcept { *this = std::move(o); }
  Mtz& operator=(Mtz&& o) noexcept {
    same_byte_order = o.same_byte_order;
    header_offset = o.header_offset;
    version_stamp = std::move(o.version_stamp);
    title = std::move(o.title);
    nreflections = o.nreflections;
    sort_order = o.sort_order;
    min_1_d2 = o.min_1_d2;
    max_1_d2 = o.max_1_d2;
    valm = o.valm;
    nsymop = o.nsymop;
    cell = std::move(o.cell);
    spacegroup_number = o.spacegroup_number;
    spacegroup_name = std::move(o.spacegroup_name);
    symops = std::move(o.symops);
    spacegroup = o.spacegroup;
    datasets = std::move(o.datasets);
    columns = std::move(o.columns);
    batches = std::move(o.batches);
    history = std::move(o.history);
    appended_text = std::move(o.appended_text);
    data = std::move(o.data);
    warnings = o.warnings;
    for (Mtz::Column& col : columns)
      col.parent = this;
    return *this;
  }
  Mtz(Mtz const&) = delete;
  Mtz& operator=(Mtz const&) = delete;

  void add_base() {
    datasets.push_back({0, "HKL_base", "HKL_base", "HKL_base", cell, 0.});
    for (int i = 0; i != 3; ++i)
      add_column(std::string(1, "HKL"[i]), 'H', 0, i, false);
  }

  // Functions to use after MTZ headers (and data) is read.

  double resolution_high() const { return std::sqrt(1.0 / max_1_d2); }
  double resolution_low() const  { return std::sqrt(1.0 / min_1_d2); }

  UnitCell& get_cell(int dataset=-1) {
    for (Dataset& ds : datasets)
      if (ds.id == dataset && ds.cell.is_crystal() && ds.cell.a > 0)
        return ds.cell;
    return cell;
  }

  const UnitCell& get_cell(int dataset=-1) const {
    return const_cast<Mtz*>(this)->get_cell(dataset);
  }

  void set_cell_for_all(const UnitCell& new_cell) {
    cell = new_cell;
    cell.set_cell_images_from_spacegroup(spacegroup);  // probably not needed
    for (Dataset& ds : datasets)
      ds.cell = cell;
  }

  UnitCell get_average_cell_from_batch_headers(double* rmsd) const {
    if (rmsd)
      for (int i = 0; i < 6; ++i)
        rmsd[i] = 0.;
    double avg[6] = {0., 0., 0., 0., 0., 0.};
    for (const Batch& batch : batches)
      for (int i = 0; i < 6; ++i) {
        // if batch headers are not set correctly, return global cell
        if (batch.floats[i] <= 0)
          return cell;
        avg[i] += batch.floats[i];
      }
    if (avg[0] <= 0 || avg[1] <= 0 || avg[2] <= 0 ||
        avg[3] <= 0 || avg[4] <= 0 || avg[5] <= 0)
      return UnitCell();
    size_t n = batches.size();
    for (int i = 0; i < 6; ++i)
      avg[i] /= n;
    if (rmsd) {
      for (const Batch& batch : batches)
        for (int i = 0; i < 6; ++i)
          rmsd[i] += sq(avg[i] - batch.floats[i]);
      for (int i = 0; i < 6; ++i)
        rmsd[i] = std::sqrt(rmsd[i] / n);
    }
    return UnitCell(avg[0], avg[1], avg[2], avg[3], avg[4], avg[5]);
  }

  void set_spacegroup(const SpaceGroup* new_sg) {
    spacegroup = new_sg;
    spacegroup_number = new_sg ? spacegroup->ccp4 : 0;
    spacegroup_name = new_sg ? spacegroup->hm : "";
  }

  Dataset& last_dataset() {
    if (datasets.empty())
      fail("MTZ dataset not found (missing DATASET header line?).");
    return datasets.back();
  }
  Dataset& dataset(int id) {
    if ((size_t)id < datasets.size() && datasets[id].id == id)
      return datasets[id];
    for (Dataset& d : datasets)
      if (d.id == id)
        return d;
    fail("MTZ file has no dataset with ID " + std::to_string(id));
  }
  const Dataset& dataset(int id) const {
    return const_cast<Mtz*>(this)->dataset(id);
  }
  Dataset* dataset_with_name(const std::string& name) {
    for (Dataset& d : datasets)
      if (d.dataset_name == name)
        return &d;
    return nullptr;
  }
  const Dataset* dataset_with_name(const std::string& label) const {
    return const_cast<Mtz*>(this)->dataset_with_name(label);
  }
  int count(const std::string& label) const {
    int n = 0;
    for (const Column& col : columns)
      if (col.label == label)
        ++n;
    return n;
  }
  int count_type(char type) const {
    int n = 0;
    for (const Column& col : columns)
      if (col.type == type)
        ++n;
    return n;
  }
  Column* column_with_label(const std::string& label,
                            const Dataset* ds=nullptr) {
    for (Column& col : columns)
      if (col.label == label && (!ds || ds->id == col.dataset_id))
        return &col;
    return nullptr;
  }
  const Column* column_with_label(const std::string& label,
                                  const Dataset* ds=nullptr) const {
    return const_cast<Mtz*>(this)->column_with_label(label, ds);
  }
  const Column& get_column_with_label(const std::string& label,
                                      const Dataset* ds=nullptr) const {
    if (const Column* col = column_with_label(label, ds))
      return *col;
    fail("Column label not found: " + label);
  }
  std::vector<const Column*> columns_with_type(char type) const {
    std::vector<const Column*> cols;
    for (const Column& col : columns)
      if (col.type == type)
        cols.push_back(&col);
    return cols;
  }

  std::vector<int> positions_of_columns_with_type(char col_type) const {
    std::vector<int> cols;
    for (int i = 0; i < (int) columns.size(); ++i)
      if (columns[i].type == col_type)
        cols.push_back(i);
    return cols;
  }

  // F(+)/(-) pairs should have type G (and L for sigma),
  // I(+)/(-) -- K (M for sigma), but E(+)/(-) has no special column type,
  // so here we use column labels not types.
  std::vector<std::pair<int,int>> positions_of_plus_minus_columns() const {
    std::vector<std::pair<int,int>> r;
    for (int i = 0; i < (int) columns.size(); ++i) {
      const Column& col = columns[i];
      size_t sign_pos = col.label.find("(+)");
      if (sign_pos != std::string::npos) {
        std::string minus_label = columns[i].label;
        minus_label[sign_pos+1] = '-';
        for (int j = 0; j < (int) columns.size(); ++j)
          if (columns[j].label == minus_label &&
              columns[j].type == col.type &&
              columns[j].dataset_id == col.dataset_id) {
            r.emplace_back(i, j);
            break;
          }
      }
    }
    return r;
  }

  const Column* column_with_one_of_labels(std::initializer_list<const char*> labels) const {
    for (const char* label : labels) {
      if (const Column* col = column_with_label(label))
        return col;
    }
    return nullptr;
  }

  Column* column_with_type_and_one_of_labels(char type,
                                             std::initializer_list<const char*> labels) {
    for (Column& col : columns)
      if (col.type == type) {
        for (const char* label : labels)
          if (col.label == label)
            return &col;
      }
    return nullptr;
  }

  Column* rfree_column() {
    // cf. MtzToCif::default_spec in mtz2cif.hpp
    return column_with_type_and_one_of_labels('I',
        {"FREE", "RFREE", "FREER", "FreeR_flag", "R-free-flags", "FreeRflag"});
  }
  const Column* rfree_column() const {
    return const_cast<Mtz*>(this)->rfree_column();
  }

  Column* imean_column() {
    return column_with_type_and_one_of_labels('J', {"IMEAN", "I", "IOBS", "I-obs"});
  }
  const Column* imean_column() const {
    return const_cast<Mtz*>(this)->imean_column();
  }

  Column* iplus_column() {
    return column_with_type_and_one_of_labels('K', {"I(+)", "IOBS(+)", "I-obs(+)"});
  }
  const Column* iplus_column() const {
    return const_cast<Mtz*>(this)->iplus_column();
  }

  Column* iminus_column() {
    return column_with_type_and_one_of_labels('K', {"I(-)", "IOBS(-)", "I-obs(-)"});
  }
  const Column* iminus_column() const {
    return const_cast<Mtz*>(this)->iminus_column();
  }

  bool has_data() const {
    return data.size() == columns.size() * nreflections;
  }

  bool is_merged() const { return batches.empty(); }

  void extend_min_max_1_d2(const UnitCell& uc, double& min, double& max) const {
    for (size_t i = 0; i < data.size(); i += columns.size()) {
      double res = uc.calculate_1_d2_double(data[i+0], data[i+1], data[i+2]);
      if (res < min)
        min = res;
      if (res > max)
        max = res;
    }
  }

  std::array<double,2> calculate_min_max_1_d2() const {
    if (!has_data() || columns.size() < 3)
      fail("No data.");
    double min_value = INFINITY;
    double max_value = 0.;
    if (cell.is_crystal() && cell.a > 0)
      extend_min_max_1_d2(cell, min_value, max_value);
    const UnitCell* prev_cell = nullptr;
    for (const Dataset& ds : datasets)
      if (ds.cell.is_crystal() && ds.cell.a > 0 && ds.cell != cell &&
          (!prev_cell || ds.cell != *prev_cell)) {
        extend_min_max_1_d2(ds.cell, min_value, max_value);
        prev_cell = &ds.cell;
      }
    if (min_value == INFINITY)
      min_value = 0;
    return {{min_value, max_value}};
  }

  void update_reso() {
    std::array<double,2> reso = calculate_min_max_1_d2();
    min_1_d2 = reso[0];
    max_1_d2 = reso[1];
  }

  // Functions for reading MTZ headers and data.

  void toggle_endiannes() {
    same_byte_order = !same_byte_order;
    swap_eight_bytes(&header_offset);
  }

  template<typename Stream>
  void read_first_bytes(Stream& stream) {
    char buf[20] = {0};

    if (!stream.read(buf, 20))
      fail("Could not read the MTZ file (is it empty?)");
    if (buf[0] != 'M' || buf[1] != 'T' || buf[2] != 'Z' || buf[3] != ' ')
      fail("Not an MTZ file - it does not start with 'MTZ '");

    // Bytes 9-12 have so-called machine stamp:
    // "The first 4 half-bytes represent the real, complex, integer and
    // character formats".
    // We don't try to handle all the combinations here, only the two most
    // common: big endian (for all types) and little endian (for all types).
    // BE is denoted by 1 and LE by 4.
    // If we get a value different than 1 and 4 we assume the native byte order.
    if ((buf[9] & 0xf0) == (is_little_endian() ? 0x10 : 0x40))
      toggle_endiannes();

    std::int32_t tmp_header_offset;
    std::memcpy(&tmp_header_offset, buf + 4, 4);
    if (!same_byte_order)
      swap_four_bytes(&tmp_header_offset);

    if (tmp_header_offset == -1) {
      std::memcpy(&header_offset, buf + 12, 8);
      if (!same_byte_order) {
        swap_eight_bytes(&header_offset);
      }
    } else {
      header_offset = (int64_t) tmp_header_offset;
    }
  }

  static const char* skip_word(const char* line) {
    while (*line != '\0' && !std::isspace(*line))
      ++line;
    while (std::isspace(*line))
      ++line;
    return line;
  }

  static UnitCell read_cell_parameters(const char* line) {
    double a = fast_atof(line, &line);
    double b = fast_atof(line, &line);
    double c = fast_atof(line, &line);
    double alpha = fast_atof(line, &line);
    double beta = fast_atof(line, &line);
    double gamma = fast_atof(line, &line);
    return UnitCell(a, b, c, alpha, beta, gamma);
  }

  template<typename T> void warn(const T& text) const {
    if (warnings)
      *warnings << text << std::endl;
  }

  template<typename Stream>
  void seek_headers(Stream& stream) {
    std::ptrdiff_t pos = 4 * std::ptrdiff_t(header_offset - 1);
    if (!stream.seek(pos))
      fail("Cannot rewind to the MTZ header at byte " + std::to_string(pos));
  }

  // read headers until END
  template<typename Stream>
  void read_main_headers(Stream& stream) {
    char line[81] = {0};
    seek_headers(stream);
    int ncol = 0;
    bool has_batch = false;
    while (stream.read(line, 80) && ialpha3_id(line) != ialpha3_id("END")) {
      const char* args = skip_word(line);
      switch (ialpha4_id(line)) {
        case ialpha4_id("VERS"):
          version_stamp = rtrim_str(args);
          break;
        case ialpha4_id("TITL"):
          title = rtrim_str(args);
          break;
        case ialpha4_id("NCOL"): {
          ncol = simple_atoi(args, &args);
          nreflections = simple_atoi(args, &args);
          int nbatches = simple_atoi(args);
          if (nbatches < 0 || nbatches > 10000000)  // sanity check
            fail("Wrong NCOL header");
          batches.resize(nbatches);
          break;
        }
        case ialpha4_id("CELL"):
          cell = read_cell_parameters(args);
          break;
        case ialpha4_id("SORT"):
          for (int& n : sort_order)
            n = simple_atoi(args, &args);
          break;
        case ialpha4_id("SYMI"): {
          nsymop = simple_atoi(args, &args);
          symops.reserve(nsymop);
          simple_atoi(args, &args); // ignore number of primitive operations
          args = skip_word(skip_blank(args)); // ignore lattice type
          spacegroup_number = simple_atoi(args, &args);
          args = skip_blank(args);
          if (*args != '\'')
            spacegroup_name = read_word(args);
          else if (const char* end = std::strchr(++args, '\''))
            spacegroup_name.assign(args, end);
          // ignore point group which is at the end of args
          break;
        }
        case ialpha4_id("SYMM"):
          symops.push_back(parse_triplet(args));
          break;
        case ialpha4_id("RESO"):
          min_1_d2 = fast_atof(args, &args);
          max_1_d2 = fast_atof(args, &args);
          break;
        case ialpha4_id("VALM"):
          if (*args != 'N') {
            const char* endptr;
            float v = (float) fast_atof(args, &endptr);
            if (*endptr == '\0' || is_space(*endptr))
              valm = v;
            else
              warn("Unexpected VALM value: " + rtrim_str(args));
          }
          break;
        case ialpha4_id("COLU"): {
          columns.emplace_back();
          Column& col = columns.back();
          col.label = read_word(args, &args);
          col.type = read_word(args, &args)[0];
          col.min_value = (float) fast_atof(args, &args);
          col.max_value = (float) fast_atof(args, &args);
          col.dataset_id = simple_atoi(args);
          col.parent = this;
          col.idx = columns.size() - 1;
          break;
        }
        case ialpha4_id("COLS"):
          if (columns.empty())
            fail("MTZ: COLSRC before COLUMN?");
          args = skip_word(args);
          columns.back().source = read_word(args);
          break;
        case ialpha4_id("COLG"):
          // Column group - not used.
          break;
        case ialpha4_id("NDIF"):
          datasets.reserve(simple_atoi(args));
          break;
        case ialpha4_id("PROJ"):
          datasets.emplace_back();
          datasets.back().id = simple_atoi(args, &args);
          datasets.back().project_name = read_word(skip_word(args));
          datasets.back().wavelength = 0.0;
          break;
        case ialpha4_id("CRYS"):
          if (simple_atoi(args, &args) == last_dataset().id)
            datasets.back().crystal_name = read_word(args);
          else
            warn("MTZ CRYSTAL line: unusual numbering.");
          break;
        case ialpha4_id("DATA"):
          if (simple_atoi(args, &args) == last_dataset().id)
            datasets.back().dataset_name = read_word(args);
          else
            warn("MTZ DATASET line: unusual numbering.");
          break;
        case ialpha4_id("DCEL"):
          if (simple_atoi(args, &args) == last_dataset().id)
            datasets.back().cell = read_cell_parameters(args);
          else
            warn("MTZ DCELL line: unusual numbering.");
          break;
        // case("DRES"): not in use yet
        case ialpha4_id("DWAV"):
          if (simple_atoi(args, &args) == last_dataset().id)
            datasets.back().wavelength = fast_atof(args);
          else
            warn("MTZ DWAV line: unusual numbering.");
          break;
        case ialpha4_id("BATCH"):
          // We take number of batches from the NCOL record and serial numbers
          // from BH. This header could be used only to check consistency.
          has_batch = true;
          break;
        default:
          warn("Unknown header: " + rtrim_str(line));
      }
    }
    if (ncol != (int) columns.size())
      fail("Number of COLU records inconsistent with NCOL record.");
    if (has_batch != !batches.empty())
      fail("BATCH header inconsistent with NCOL record.");
  }

  // read the part between END and MTZENDOFHEADERS
  template<typename Stream>
  void read_history_and_batch_headers(Stream& stream) {
    char buf[81] = {0};
    int n_headers = 0;
    while (stream.read(buf, 80) && ialpha4_id(buf) != ialpha4_id("MTZE")) {
      if (n_headers != 0) {
        const char* start = skip_blank(buf);
        const char* end = rtrim_cstr(start, start+80);
        history.emplace_back(start, end);
        --n_headers;
      } else if (ialpha4_id(buf) == ialpha4_id("MTZH")) {
        n_headers = simple_atoi(skip_word(buf));
        if (n_headers < 0 || n_headers > 30) {
          warn("Wrong MTZ: number of headers should be between 0 and 30");
          return;
        }
        history.reserve(n_headers);
      } else if (ialpha4_id(buf) == ialpha4_id("MTZB")) {
        for (Batch& batch : batches) {
          stream.read(buf, 80);
          if (ialpha3_id(buf) != ialpha3_id("BH "))
            fail("Missing BH header");
          const char* args = skip_word(buf);
          batch.number = simple_atoi(args, &args);
          int total_words = simple_atoi(args, &args);
          int int_words = simple_atoi(args, &args);
          int float_words = simple_atoi(args);
          if (total_words != int_words + float_words || total_words > 1000)
            fail("Wrong BH header");
          stream.read(buf, 80); // TITLE
          const char* end = rtrim_cstr(buf + 6, buf+76);
          batch.title.assign(buf, end - buf);
          batch.ints.resize(int_words);
          stream.read(batch.ints.data(), int_words * 4);
          batch.floats.resize(float_words);
          stream.read(batch.floats.data(), float_words * 4);
          stream.read(buf, 80);
          if (ialpha4_id(buf) != ialpha4_id("BHCH"))
            fail("Missing BHCH header");
          split_str_into_multi(buf + 5, " \t", batch.axes);
        }
      }
    }
    appended_text = stream.read_rest();
  }

  void setup_spacegroup() {
    spacegroup = find_spacegroup_by_name(spacegroup_name,
                                         cell.alpha, cell.gamma);
    if (!spacegroup) {
      warn("MTZ: unrecognized spacegroup name: " + spacegroup_name);
      return;
    }
    if (spacegroup->ccp4 != spacegroup_number)
      warn("MTZ: inconsistent spacegroup name and number");
    cell.set_cell_images_from_spacegroup(spacegroup);
    for (Dataset& d : datasets)
      d.cell.set_cell_images_from_spacegroup(spacegroup);
  }

  template<typename Stream>
  void read_raw_data(Stream& stream) {
    size_t n = columns.size() * nreflections;
    data.resize(n);
    if (!stream.seek(80))
      fail("Cannot rewind to the MTZ data.");
    if (!stream.read(data.data(), 4 * n))
      fail("Error when reading MTZ data");
    if (!same_byte_order)
      for (float& f : data)
        swap_four_bytes(&f);
  }

  template<typename Stream>
  void read_all_headers(Stream& stream) {
    read_first_bytes(stream);
    read_main_headers(stream);
    read_history_and_batch_headers(stream);
    setup_spacegroup();
    if (datasets.empty())
      datasets.push_back({0, "HKL_base", "HKL_base", "HKL_base", cell, 0.});
  }

  template<typename Stream>
  void read_stream(Stream&& stream, bool with_data) {
    read_all_headers(stream);
    if (with_data)
      read_raw_data(stream);
  }

  void read_file(const std::string& path) {
    fileptr_t f = file_open(path.c_str(), "rb");
    try {
      source_path = path;
      read_stream(FileStream{f.get()}, true);
    } catch (std::runtime_error& e) {
      fail(std::string(e.what()) + ": " + path);
    }
  }

  template<typename Input>
  void read_input(Input&& input, bool with_data) {
    source_path = input.path();
    if (input.is_stdin()) {
      read_stream(FileStream{stdin}, with_data);
    } else if (CharArray mem = input.uncompress_into_buffer()) {
      read_stream(mem.stream(), with_data);
    } else {
      fileptr_t f = file_open(input.path().c_str(), "rb");
      read_stream(FileStream{f.get()}, true);
    }
  }

  std::vector<int> sorted_row_indices(int use_first=3) const {
    if (!has_data())
      fail("No data.");
    if (use_first <= 0 || use_first >= (int) columns.size())
      fail("Wrong use_first arg in Mtz::sort.");
    std::vector<int> indices(nreflections);
    for (int i = 0; i != nreflections; ++i)
      indices[i] = i;
    std::stable_sort(indices.begin(), indices.end(), [&](int i, int j) {
      int a = i * (int) columns.size();
      int b = j * (int) columns.size();
      for (int n = 0; n < use_first; ++n)
        if (data[a+n] != data[b+n])
          return data[a+n] < data[b+n];
      return false;
    });
    return indices;
  }

  bool sort(int use_first=3) {
    std::vector<int> indices = sorted_row_indices(use_first);
    sort_order = {{0, 0, 0, 0, 0}};
    for (int i = 0; i < use_first; ++i)
      sort_order[i] = i + 1;
    if (std::is_sorted(indices.begin(), indices.end()))
      return false;
    std::vector<float> new_data(data.size());
    size_t w = columns.size();
    for (size_t i = 0; i != indices.size(); ++i)
      std::memcpy(&new_data[i * w], &data[indices[i] * w], w * sizeof(float));
    data.swap(new_data);
    return true;
  }

  Miller get_hkl(size_t offset) const {
    return {{(int)data[offset], (int)data[offset+1], (int)data[offset+2]}};
  }
  void set_hkl(size_t offset, const Miller& hkl) {
    for (int i = 0; i != 3; ++i)
      data[offset + i] = static_cast<float>(hkl[i]);
  }

  /// Returns offset of the first hkl or (size_t)-1. Can be slow.
  size_t find_offset_of_hkl(const Miller& hkl, size_t start=0) const;

  /// (for merged MTZ only) change HKL to ASU equivalent, adjust phases, etc
  void ensure_asu(bool tnt_asu=false);

  /// reindex data, usually followed by ensure_asu()
  void reindex(const Op& op, std::ostream* out);

  /// Change symmetry to P1 and expand reflections. Does not sort.
  /// Similar to command EXPAND in SFTOOLS.
  void expand_to_p1();

  // (for unmerged MTZ only) change HKL according to M/ISYM
  bool switch_to_original_hkl() {
    if (indices_switched_to_original)
      return false;
    if (!has_data())
      fail("switch_to_original_hkl(): data not read yet");
    const Column* col = column_with_label("M/ISYM");
    if (col == nullptr || col->type != 'Y' || col->idx < 3)
      return false;
    std::vector<Op> inv_symops;
    inv_symops.reserve(symops.size());
    for (const Op& op : symops)
      inv_symops.push_back(op.inverse());
    for (size_t n = 0; n + col->idx < data.size(); n += columns.size()) {
      int isym = static_cast<int>(data[n + col->idx]) & 0xFF;
      const Op& op = inv_symops.at((isym - 1) / 2);
      Miller hkl = op.apply_to_hkl(get_hkl(n));
      int sign = (isym & 1) ? 1 : -1;
      for (int i = 0; i < 3; ++i)
        data[n+i] = static_cast<float>(sign * hkl[i]);
    }
    indices_switched_to_original = true;
    return true;
  }

  // (for unmerged MTZ only) change HKL to ASU equivalent and set ISYM
  bool switch_to_asu_hkl() {
    if (!indices_switched_to_original)
      return false;
    if (!has_data())
      fail("switch_to_asu_hkl(): data not read yet");
    const Column* col = column_with_label("M/ISYM");
    if (col == nullptr || col->type != 'Y' || col->idx < 3 || !spacegroup)
      return false;
    size_t misym_idx = col->idx;
    UnmergedHklMover hkl_mover(spacegroup);
    for (size_t n = 0; n + col->idx < data.size(); n += columns.size()) {
      Miller hkl = get_hkl(n);
      int isym = hkl_mover.move_to_asu(hkl);  // modifies hkl
      set_hkl(n, hkl);
      float& misym = data[n + misym_idx];
      misym = float(((int)misym & ~0xff) | isym);
    }
    indices_switched_to_original = false;
    return true;
  }

  Dataset& add_dataset(const std::string& name) {
    int id = 0;
    for (const Dataset& d : datasets)
      if (d.id >= id)
        id = d.id + 1;
    datasets.push_back({id, name, name, name, cell, 0.0});
    return datasets.back();
  }

  Column& add_column(const std::string& label, char type,
                     int dataset_id, int pos, bool expand_data) {
    if (datasets.empty())
      fail("No datasets.");
    if (dataset_id < 0)
      dataset_id = datasets.back().id;
    else
      dataset(dataset_id); // check if such dataset exist
    if (pos > (int) columns.size())
      fail("Requested column position after the end.");
    if (pos < 0)
      pos = (int) columns.size();
    auto col = columns.emplace(columns.begin() + pos);
    for (auto i = col + 1; i != columns.end(); ++i)
      i->idx++;
    col->dataset_id = dataset_id;
    col->type = type;
    col->label = label;
    col->parent = this;
    col->idx = pos;
    if (expand_data)
      expand_data_rows(1, pos);
    return *col;
  }

  // helper_functions
  void check_column(size_t idx, const char* msg) const {
    if (!has_data())
      fail(msg, ": data not read yet");
    if (idx >= columns.size())
      fail(msg, ": no column with 0-based index ", std::to_string(idx));
  }
  void check_trailing_cols(const Column& src_col,
                           const std::vector<std::string>& trailing_cols) const {
    assert(src_col.parent == this);
    if (!has_data())
      fail("data in source mtz not read yet");
    if (src_col.idx + trailing_cols.size() >= columns.size())
      fail("Not enough columns after " + src_col.label);
    for (size_t i = 0; i < trailing_cols.size(); ++i)
      if (!trailing_cols[i].empty() &&
          trailing_cols[i] != columns[src_col.idx + i + 1].label)
        fail("expected trailing column ", trailing_cols[i], ", found ", src_col.label);
  }
  void do_replace_column(size_t dest_idx, const Column& src_col,
                         const std::vector<std::string>& trailing_cols) {
    const Mtz* src_mtz = src_col.parent;
    for (size_t i = 0; i <= trailing_cols.size(); ++i) {
      Column& dst = columns[dest_idx + i];
      const Column& src = src_mtz->columns[src_col.idx + i];
      dst.type = src.type;
      dst.label = src.label;
      dst.min_value = src.min_value;
      dst.max_value = src.max_value;
      dst.source = src.source;
      dst.dataset_id = src.dataset_id;
    }
    if (src_mtz == this) {
      // internal copying
      for (size_t n = 0; n < data.size(); n += columns.size())
        for (size_t i = 0; i <= trailing_cols.size(); ++i)
          data[n + dest_idx + i] = data[n + src_col.idx + i];
    } else {
      // external copying - need to match indices
      std::vector<int> dst_indices = sorted_row_indices();
      std::vector<int> src_indices = src_mtz->sorted_row_indices();
      // cf. for_matching_reflections()
      size_t dst_stride = columns.size();
      size_t src_stride = src_mtz->columns.size();
      auto dst = dst_indices.begin();
      auto src = src_indices.begin();
      while (dst != dst_indices.end() && src != src_indices.end()) {
        Miller dst_hkl = get_hkl(*dst * dst_stride);
        Miller src_hkl = src_mtz->get_hkl(*src * src_stride);
        if (dst_hkl == src_hkl) {
          // copy values
          for (size_t i = 0; i <= trailing_cols.size(); ++i)
            data[*dst * dst_stride + dest_idx + i] =
              src_mtz->data[*src * src_stride + src_col.idx + i];
          ++dst;
          ++src;
        } else if (dst_hkl < src_hkl) {
          ++dst;
        } else {
          ++src;
        }
      }
    }
  }

  // extra_col are columns right after src_col that are also copied.
  Column& replace_column(size_t dest_idx, const Column& src_col,
                         const std::vector<std::string>& trailing_cols={}) {
    src_col.parent->check_trailing_cols(src_col, trailing_cols);
    check_column(dest_idx + trailing_cols.size(), "replace_column()");
    do_replace_column(dest_idx, src_col, trailing_cols);
    return columns[dest_idx];
  }

  // If dest_idx < 0 - columns are appended at the end
  // append new column(s), otherwise overwrite existing ones.
  Column& copy_column(int dest_idx, const Column& src_col,
                      const std::vector<std::string>& trailing_cols={}) {
    // check input consistency
    if (!has_data())
      fail("copy_column(): data not read yet");
    src_col.parent->check_trailing_cols(src_col, trailing_cols);
    // add new columns
    if (dest_idx < 0)
      dest_idx = (int) columns.size();
    // if src_col is from this Mtz it may get invalidated when adding columns
    int col_idx = -1;
    if (src_col.parent == this) {
      col_idx = (int) src_col.idx;
      if (col_idx >= dest_idx)
        col_idx += 1 + (int)trailing_cols.size();
    }
    for (int i = 0; i <= (int) trailing_cols.size(); ++i)
      add_column("", ' ', -1, dest_idx + i, false);
    expand_data_rows(1 + trailing_cols.size(), dest_idx);
    // copy the data
    const Column& src_col_now = col_idx < 0 ? src_col : columns[col_idx];
    // most of the work (hkl-based row matching and data copying) is done here:
    do_replace_column(dest_idx, src_col_now, trailing_cols);
    return columns[dest_idx];
  }

  void remove_column(size_t idx) {
    check_column(idx, "remove_column()");
    columns.erase(columns.begin() + idx);
    for (size_t i = idx; i < columns.size(); ++i)
      --columns[i].idx;
    for (size_t dest = idx, source = idx + 1; source < data.size(); ++source)
      for (size_t i = 0; i < columns.size() && source < data.size(); ++i)
        data[dest++] = data[source++];
    data.resize(columns.size() * nreflections);
  }

  template <typename Func>
  void remove_rows_if(Func condition) {
    if (!has_data())
      fail("No data.");
    auto out = data.begin();
    size_t width = columns.size();
    for (auto r = data.begin(); r < data.end(); r += width)
      if (!condition(&*r)) {
        if (r != out)
          std::copy(r, r + width, out);
        out += width;
      }
    data.erase(out, data.end());
    nreflections = int(data.size() / width);
  }

  void expand_data_rows(size_t added, int pos_=-1) {
    size_t old_row_size = columns.size() - added;
    if (data.size() != old_row_size * nreflections)
      fail("Internal error");
    data.resize(columns.size() * nreflections);
    size_t pos = pos_ == -1 ? old_row_size : (size_t) pos_;
    if (pos > old_row_size)
      fail("expand_data_rows(): pos out of range");
    std::vector<float>::iterator dst = data.end();
    for (int i = nreflections; i-- != 0; ) {
      for (size_t j = old_row_size; j-- != pos; )
        *--dst = data[i * old_row_size + j];
      for (size_t j = added; j-- != 0; )
        *--dst = NAN;
      for (size_t j = pos; j-- != 0; )
        *--dst = data[i * old_row_size + j];
    }
    assert(dst == data.begin());
  }

  void set_data(const float* new_data, size_t n) {
    size_t ncols = columns.size();
    if (n % ncols != 0)
      fail("Mtz.set_data(): expected " + std::to_string(ncols) + " columns.");
    nreflections = int(n / ncols);
    data.assign(new_data, new_data + n);
  }

  // Function for writing MTZ file
  void write_to_cstream(std::FILE* stream) const;
  void write_to_string(std::string& str) const;
  void write_to_file(const std::string& path) const;

private:
  template<typename Write> void write_to_stream(Write write) const;
};


inline Mtz read_mtz_file(const std::string& path) {
  Mtz mtz;
  mtz.read_file(path);
  return mtz;
}

template<typename Input>
Mtz read_mtz(Input&& input, bool with_data) {
  Mtz mtz;
  mtz.read_input(std::forward<Input>(input), with_data);
  return mtz;
}


// Abstraction of data source, cf. ReflnDataProxy.
struct MtzDataProxy {
  const Mtz& mtz_;
  size_t stride() const { return mtz_.columns.size(); }
  size_t size() const { return mtz_.data.size(); }
  using num_type = float;
  float get_num(size_t n) const { return mtz_.data[n]; }
  const UnitCell& unit_cell() const { return mtz_.cell; }
  const SpaceGroup* spacegroup() const { return mtz_.spacegroup; }
  Miller get_hkl(size_t offset) const { return mtz_.get_hkl(offset); }

  size_t column_index(const std::string& label) const {
    if (const Mtz::Column* col = mtz_.column_with_label(label))
      return col->idx;
    fail("MTZ file has no column with label: " + label);
  }
};

// Like above, but here the data is stored outside of the Mtz class
struct MtzExternalDataProxy : MtzDataProxy {
  const float* data_;
  MtzExternalDataProxy(const Mtz& mtz, const float* data)
    : MtzDataProxy{mtz}, data_(data) {}
  size_t size() const { return mtz_.columns.size() * mtz_.nreflections; }
  float get_num(size_t n) const { return data_[n]; }
  Miller get_hkl(size_t offset) const {
    return {{(int)data_[offset + 0],
             (int)data_[offset + 1],
             (int)data_[offset + 2]}};
  }
};

inline MtzDataProxy data_proxy(const Mtz& mtz) { return {mtz}; }

} // namespace gemmi

#endif
