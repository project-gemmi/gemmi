// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cassert>
#include <cmath>         // for isnan
#include <cstdint>       // for int32_t
#include <algorithm>     // for copy
#include <array>
#include <initializer_list>
#include <string>
#include <vector>
#include "fail.hpp"      // for fail
#include "input.hpp"     // for AnyStream, FileStream, CharArray
#include "iterator.hpp"  // for StrideIter
#include "logger.hpp"    // for Logger
#include "math.hpp"      // for rad, Mat33
#include "symmetry.hpp"  // for find_spacegroup_by_name, SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for ialpha4_id, rtrim_str, ialpha3_id, ...

namespace gemmi {

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

struct MtzMetadata {
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
  std::vector<std::string> history;
  std::string appended_text;
  // used to report non-critical problems when reading a file (also used in mtz2cif)
  Logger logger;
};

struct GEMMI_DLL Mtz : public MtzMetadata {
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

    const Column* get_next_column_if_type(char next_type) const {
      if (idx + 1 < parent->columns.size()) {
        const Column& next_col = parent->columns[idx + 1];
        if (next_col.dataset_id == dataset_id && next_col.type == next_type)
          return &next_col;
      }
      return nullptr;
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
    int number = 0;
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

  std::vector<Dataset> datasets;
  std::vector<Column> columns;
  std::vector<Batch> batches;
  std::vector<float> data;

  explicit Mtz(bool with_base=false) {
    if (with_base)
      add_base();
  }
  Mtz(Mtz&& o) noexcept { *this = std::move(o); }
  Mtz& operator=(Mtz&& o) noexcept {
    MtzMetadata::operator=(std::move(o));
    datasets = std::move(o.datasets);
    columns = std::move(o.columns);
    batches = std::move(o.batches);
    data = std::move(o.data);
    for (Mtz::Column& col : columns)
      col.parent = this;
    return *this;
  }

  // explicit to be aware where we make copies
  explicit Mtz(const Mtz& o) : MtzMetadata(o) {
    datasets = o.datasets;
    columns = o.columns;
    batches = o.batches;
    data = o.data;
    for (Mtz::Column& col : columns)
      col.parent = this;
  }

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

  UnitCellParameters get_average_cell_from_batch_headers(double* rmsd) const;

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

  Column* column_with_label(const std::string& label, const Dataset* ds=nullptr, char type='*') {
    for (Column& col : columns)
      if (col.label == label && (!ds || ds->id == col.dataset_id)
                             && (type == '*' || type == col.type))
        return &col;
    return nullptr;
  }
  const Column* column_with_label(const std::string& label, const Dataset* ds=nullptr,
                                  char type='*') const {
    return const_cast<Mtz*>(this)->column_with_label(label, ds, type);
  }

  const Column& get_column_with_label(const std::string& label, const Dataset* ds=nullptr) const {
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

  /// the order of labels matters
  const Column* column_with_one_of_labels(std::initializer_list<const char*> labels,
                                          char type='*') const {
    for (const char* label : labels)
      if (const Column* col = column_with_label(label, nullptr, type))
        return col;
    return nullptr;
  }

  /// the order of labels doesn't matter
  Column* column_with_type_and_any_of_labels(char type, std::initializer_list<const char*> labels) {
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
    return column_with_type_and_any_of_labels('I',
        {"FREE", "RFREE", "FREER", "FreeR_flag", "R-free-flags", "FreeRflag", "R_FREE_FLAGS"});
  }
  const Column* rfree_column() const {
    return const_cast<Mtz*>(this)->rfree_column();
  }

  Column* imean_column() {
    return column_with_type_and_any_of_labels('J', {"IMEAN", "I", "IOBS", "I-obs"});
  }
  const Column* imean_column() const {
    return const_cast<Mtz*>(this)->imean_column();
  }

  Column* iplus_column() {
    return column_with_type_and_any_of_labels('K', {"I(+)", "IOBS(+)", "I-obs(+)", "Iplus"});
  }
  const Column* iplus_column() const {
    return const_cast<Mtz*>(this)->iplus_column();
  }

  Column* iminus_column() {
    return column_with_type_and_any_of_labels('K', {"I(-)", "IOBS(-)", "I-obs(-)", "Iminus"});
  }
  const Column* iminus_column() const {
    return const_cast<Mtz*>(this)->iminus_column();
  }

  bool has_data() const {
    return data.size() == columns.size() * nreflections;
  }

  bool is_merged() const { return batches.empty(); }

  /// Calculates min/max for all combinations of reflections and unit cells,
  /// where unit cells are a global CELL and per-dataset DCELL.
  std::array<double,2> calculate_min_max_1_d2() const;

  void update_reso() {
    std::array<double,2> reso = calculate_min_max_1_d2();
    min_1_d2 = reso[0];
    max_1_d2 = reso[1];
  }

  // Functions for reading MTZ headers and data.

  void toggle_endianness() {
    same_byte_order = !same_byte_order;
    swap_eight_bytes(&header_offset);
  }

  void read_first_bytes(AnyStream& stream);

  /// read headers until END
  void read_main_headers(AnyStream& stream, std::vector<std::string>* save_headers);

  /// read the part between END and MTZENDOFHEADERS
  void read_history_and_batch_headers(AnyStream& stream);

  void setup_spacegroup();

  void read_raw_data(AnyStream& stream);

  void read_all_headers(AnyStream& stream);

  void read_stream(AnyStream&& stream, bool with_data) {
    read_all_headers(stream);
    if (with_data)
      read_raw_data(stream);
  }

  void read_file(const std::string& path) {
    try {
      source_path = path;
      read_stream(FileStream(path.c_str(), "rb"), true);
    } catch (std::system_error&) {
      throw;  // system_error::what() includes path, don't add anything
    } catch (std::runtime_error& e) {
      fail(std::string(e.what()) + ": " + path);
    }
  }

  template<typename Input>
  void read_input(Input&& input, bool with_data) {
    source_path = input.path();
    if (CharArray mem = input.uncompress_into_buffer())
      read_stream(MemoryStream(mem.data(), mem.size()), with_data);
    else
      read_stream(FileStream(input.path().c_str(), "rb"), with_data);
  }

  /// the same as read_input(MaybeGzipped(path), with_data)
  void read_file_gz(const std::string& path, bool with_data=true);

  std::vector<int> sorted_row_indices(int use_first=3) const;
  bool sort(int use_first=3);

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

  /// Reindex data, usually followed by ensure_asu(). Outputs messages through logger.
  void reindex(const Op& op);

  /// Change symmetry to P1 and expand reflections. Does not sort.
  /// Similar to command EXPAND in SFTOOLS.
  void expand_to_p1();

  /// (for unmerged MTZ only) change HKL according to M/ISYM
  bool switch_to_original_hkl();

  /// (for unmerged MTZ only) change HKL to ASU equivalent and set ISYM
  bool switch_to_asu_hkl();

  Dataset& add_dataset(const std::string& name) {
    int id = 0;
    for (const Dataset& d : datasets)
      if (d.id >= id)
        id = d.id + 1;
    datasets.push_back({id, name, name, name, cell, 0.0});
    return datasets.back();
  }

  Column& add_column(const std::string& label, char type,
                     int dataset_id, int pos, bool expand_data);

  // extra_col are columns right after src_col that are also copied.
  Column& replace_column(size_t dest_idx, const Column& src_col,
                         const std::vector<std::string>& trailing_cols={});

  // If dest_idx < 0 - columns are appended at the end
  // append new column(s), otherwise overwrite existing ones.
  Column& copy_column(int dest_idx, const Column& src_col,
                      const std::vector<std::string>& trailing_cols={});

  void remove_column(size_t idx);

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
    size_t pos = pos_ == -1 ? old_row_size : (size_t) pos_;
    if (pos > old_row_size)
      fail("expand_data_rows(): pos out of range");
    vector_insert_columns(data, old_row_size, (size_t)nreflections, added, pos, NAN);
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
