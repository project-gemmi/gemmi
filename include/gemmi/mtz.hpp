// Copyright 2019 Global Phasing Ltd.
//
// MTZ reflection file format.

#ifndef GEMMI_MTZ_HPP_
#define GEMMI_MTZ_HPP_

#include <cassert>
#include <cstdint>   // for int32_t
#include <cstdio>    // for FILE, fread
#include <cstring>   // for memcpy
#include <cmath>     // for isnan
#include <algorithm> // for sort
#include <complex>   // for complex
#include <array>
#include <string>
#include <vector>
#include "atox.hpp"      // for simple_atof, simple_atoi, read_word
#include "grid.hpp"      // for Grid
#include "iterator.hpp"  // for StrideIter
#include "fileutil.hpp"  // for file_open, is_little_endian, fileptr_t, ...
#include "symmetry.hpp"  // for find_spacegroup_by_name, SpaceGroup
#include "unitcell.hpp"  // for UnitCell
#include "util.hpp"      // for fail

#ifdef  __INTEL_COMPILER
# pragma warning push
# pragma warning disable 597  // for StrideIter, the same as in cifdoc.hpp.
#endif

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

struct Mtz {
  struct Dataset {
    int id;
    std::string project_name;
    std::string crystal_name;
    std::string dataset_name;
    UnitCell cell;
    double wavelength;
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
    unsigned stride() const { return (unsigned) parent->columns.size(); }
    float& operator[](int n) { return parent->data[idx + n * stride()]; }
    float operator[](int n) const { return parent->data[idx + n * stride()]; }
    float& at(int n) { return parent->data.at(idx + n * stride()); }
    float at(int n) const { return parent->data.at(idx + n * stride()); }
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

  bool same_byte_order = true;
  std::int32_t header_offset = 0;
  std::string version_stamp;
  std::string title;
  int nreflections = 0;
  int nbatches = 0;
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
  std::vector<std::string> history;
  std::vector<float> data;

  FILE* warnings = nullptr;

  Mtz() = default;
  Mtz(Mtz&& o) noexcept { *this = std::move(o); }
  Mtz& operator=(Mtz&& o) noexcept {
    same_byte_order = o.same_byte_order;
    header_offset = o.header_offset;
    version_stamp = std::move(o.version_stamp);
    title = std::move(o.title);
    nreflections = o.nreflections;
    nbatches = o.nbatches;
    sort_order = std::move(o.sort_order);
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
    history = std::move(o.history);
    data = std::move(o.data);
    warnings = o.warnings;
    for (Mtz::Column& col : columns)
      col.parent = this;
    return *this;
  }
  Mtz(Mtz const&) = delete;
  Mtz& operator=(Mtz const&) = delete;

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
  std::vector<const Column*> columns_with_type(char type) const {
    std::vector<const Column*> cols;
    for (const Column& col : columns)
      if (col.type == type)
        cols.push_back(&col);
    return cols;
  }

  bool has_data() const {
    return data.size() == columns.size() * nreflections;
  }

  void extend_min_max_1_d2(const UnitCell& uc, double& min, double& max) const {
    for (size_t i = 0; i < data.size(); i += columns.size()) {
      double res = uc.calculate_1_d2(data[i+0], data[i+1], data[i+2]);
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

  // Functions for reading MTZ headers and data.

  void toggle_endiannes() {
    same_byte_order = !same_byte_order;
    swap_four_bytes(&header_offset);
  }

  void read_first_bytes(std::FILE* stream) {
    char buf[12] = {0};

    if (std::fread(buf, 1, 12, stream) != 12)
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

    std::memcpy(&header_offset, buf + 4, 4);
    if (!same_byte_order)
      swap_four_bytes(&header_offset);
  }

  static const char* skip_word(const char* line) {
    while (*line != '\0' && !std::isspace(*line))
      ++line;
    while (std::isspace(*line))
      ++line;
    return line;
  }

  static UnitCell read_cell_parameters(const char* line) {
    double a = simple_atof(line, &line);
    double b = simple_atof(line, &line);
    double c = simple_atof(line, &line);
    double alpha = simple_atof(line, &line);
    double beta = simple_atof(line, &line);
    double gamma = simple_atof(line, &line);
    return UnitCell(a, b, c, alpha, beta, gamma);
  }

  void warn(const std::string& text) {
    if (warnings)
      std::fprintf(warnings, "%s\n", text.c_str());
  }

  void seek_headers(std::FILE* stream) {
    if (std::fseek(stream, 4 * (header_offset - 1), SEEK_SET) != 0)
      fail("Cannot rewind to the MTZ header at byte "
           + std::to_string(header_offset));
  }

  // read headers until END
  void read_main_headers(std::FILE* stream) {
    char line[81] = {0};
    seek_headers(stream);
    int ncol = 0;
    while (std::fread(line, 1, 80, stream) == 80 &&
           ialpha3_id(line) != ialpha3_id("END")) {
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
          nbatches = simple_atoi(args, &args);
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
          min_1_d2 = simple_atof(args, &args);
          max_1_d2 = simple_atof(args, &args);
          break;
        case ialpha4_id("VALM"):
          if (*args != 'N') {
            const char* endptr;
            float v = (float) simple_atof(args, &endptr);
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
          col.min_value = (float) simple_atof(args, &args);
          col.max_value = (float) simple_atof(args, &args);
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
            datasets.back().wavelength = simple_atof(args);
          else
            warn("MTZ DWAV line: unusual numbering.");
          break;
        case ialpha4_id("BATCH"):
          // this header is not used for anything?
          break;
        default:
          warn("Unknown header: " + rtrim_str(line));
      }
    }
    if (ncol != (int) columns.size())
      fail("Number of COLU records inconsistent with NCOL record.");
  }

  // read the part between END and MTZENDOFHEADERS
  void read_history_and_batch_headers(std::FILE* stream) {
    char buf[81] = {0};
    int n_headers = 0;
    while (std::fread(buf, 1, 80, stream) == 80 &&
           ialpha4_id(buf) != ialpha4_id("MTZE")) {
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
        for (int i = 0; i < nbatches; ++i) {
          // TODO: BH, etc
        }
      }
    }
  }

  void setup_spacegroup() {
    spacegroup = find_spacegroup_by_name(spacegroup_name);
    if (!spacegroup)
      warn("MTZ: unrecognized spacegroup name: " + spacegroup_name);
    if (spacegroup->ccp4 != spacegroup_number)
      warn("MTZ: inconsistent spacegroup name and number");
  }

  void read_raw_data(std::FILE* stream) {
    size_t n = columns.size() * nreflections;
    data.resize(n);
    if (std::fseek(stream, 80, SEEK_SET) != 0)
      fail("Cannot rewind to the MTZ data.");
    if (std::fread(data.data(), 4, n, stream) != n)
      fail("Error when reading MTZ data");
    if (!same_byte_order)
      for (float& f : data)
        swap_four_bytes(&f);
  }

  void read_all_headers(std::FILE* stream) {
    read_first_bytes(stream);
    read_main_headers(stream);
    read_history_and_batch_headers(stream);
    setup_spacegroup();
  }

  std::vector<int> sorted_row_indices() const {
    if (!has_data())
      fail("No data.");
    std::vector<int> indices(nreflections);
    for (int i = 0; i != nreflections; ++i)
      indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
      int a = i * (int) columns.size();
      int b = j * (int) columns.size();
      return data[a] < data[b] || (data[a] == data[b] && (
               data[a+1] < data[b+1] || (data[a+1] == data[b+1] && (
                 data[a+2] < data[b+2]))));
    });
    return indices;
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
                     int dataset_id=-1, int pos=-1) {
    if (datasets.empty())
      fail("No datasets.");
    if (dataset_id < 0)
      dataset_id = datasets.back().id;
    else
      dataset(dataset_id); // check if such dataset exist
    if (pos >= (int) columns.size())
      fail("Requested column position after the end.");
    if (pos < 0)
      pos = columns.size();
    auto col = columns.emplace(columns.begin() + pos);
    for (auto i = col + 1; i != columns.end(); ++i)
      i->idx++;
    col->dataset_id = dataset_id;
    col->type = type;
    col->label = label;
    col->parent = this;
    col->idx = pos;
    return *col;
  }

  void expand_data_rows(int added) {
    int old_row_size = columns.size() - added;
    if ((int) data.size() != old_row_size * nreflections)
      fail("Internal error");
    data.resize(columns.size() * nreflections);
    for (int i = nreflections; i-- != 0; ) {
      for (int j = added; j-- != 0; )
        data[i * columns.size() + old_row_size + j] = NAN;
      for (int j = old_row_size; j-- != 0; )
        data[i * columns.size() + j] = data[i * old_row_size + j];
    }
  }

  void set_data(const float* new_data, size_t n) {
    if (n % columns.size() != 0)
      fail("Mtz.set_data(): expected " +
           std::to_string(columns.size()) + " columns.");
    nreflections = n / columns.size();
    data.assign(new_data, new_data + n);
  }

  // Function for writing MTZ file
  void write_to_stream(std::FILE* stream) const;
  void write_to_file(const std::string& path) const;
};

inline Mtz read_mtz_stream(std::FILE* stream, bool with_data) {
  Mtz mtz;
  mtz.read_all_headers(stream);
  if (with_data)
    mtz.read_raw_data(stream);
  return mtz;
}

inline Mtz read_mtz_file(const std::string& path) {
  fileptr_t f = file_open(path.c_str(), "rb");
  try {
    return read_mtz_stream(f.get(), true);
  } catch (std::runtime_error& e) {
    fail(std::string(e.what()) + ": " + path);
  }
}

// Abstraction of data source, cf. ReflnDataProxy.
struct MtzDataProxy {
  const Mtz& mtz_;
  bool ok() const { return mtz_.has_data(); }
  constexpr std::array<size_t,3> hkl_col() const { return {{0, 1, 2}}; }
  size_t stride() const { return mtz_.columns.size(); }
  size_t size() const { return mtz_.data.size(); }
  int get_int(size_t n) const { return (int) mtz_.data[n]; }
  float get_num(size_t n) const { return mtz_.data[n]; }
  const UnitCell& unit_cell() const { return mtz_.cell; }
  const SpaceGroup* spacegroup() const { return mtz_.spacegroup; }
};


} // namespace gemmi

#ifdef GEMMI_WRITE_IMPLEMENTATION

#include "sprintf.hpp"

namespace gemmi {

#define WRITE(...) do { \
    int len = stbsp_snprintf(buf, 81, __VA_ARGS__); \
    std::memset(buf + len, ' ', 80 - len); \
    if (std::fwrite(buf, 80, 1, stream) != 1) \
      fail("Writing MTZ file failed"); \
  } while(0)

void Mtz::write_to_stream(std::FILE* stream) const {
  // uses: data, spacegroup, nreflections, nbatches, cell, sort_order,
  //       valm, columns, datasets, history
  if (!has_data())
    fail("Cannot write Mtz which has no data");
  if (!spacegroup)
    fail("Cannot write Mtz which has no space group");
  char buf[81] = {'M', 'T', 'Z', ' ', '\0'};
  std::int32_t header_start = (int) columns.size() * nreflections + 21;
  std::memcpy(buf + 4, &header_start, 4);
  std::int32_t machst = is_little_endian() ? 0x00004144 : 0x11110000;
  std::memcpy(buf + 8, &machst, 4);
  if (std::fwrite(buf, 80, 1, stream) != 1 ||
      std::fwrite(data.data(), 4, data.size(), stream) != data.size())
    fail("Writing MTZ file failed");
  WRITE("VERS MTZ:V1.1");
  WRITE("TITLE %s", title.c_str());
  WRITE("NCOL %8zu %12d %8d", columns.size(), nreflections, nbatches);
  if (cell.is_crystal())
    WRITE("CELL  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",
          cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
  WRITE("SORT  %3d %3d %3d %3d %3d", sort_order[0], sort_order[1],
        sort_order[2], sort_order[3], sort_order[4]);
  GroupOps ops = spacegroup->operations();
  WRITE("SYMINF %3d %2d %c %5d %*s'%s' PG%s",
        ops.order(),               // number of symmetry operations
        (int) ops.sym_ops.size(),  // number of primitive operations
        spacegroup->hm[0],         // lattice type
        spacegroup->ccp4,          // space group number 
        20 - (int) std::strlen(spacegroup->hm), "",
        spacegroup->hm,            // space group name
        spacegroup->point_group_hm()); // point group name
  for (Op op : ops)
    WRITE("SYMM %s", to_upper(op.triplet()).c_str());
  auto reso = calculate_min_max_1_d2();
  WRITE("RESO %-20.12f %-20.12f", reso[0], reso[1]);
  if (std::isnan(valm))
    WRITE("VALM NAN");
  else
    WRITE("VALM %f", valm);
  for (const Column& col : columns) {
    auto minmax = calculate_min_max_disregarding_nans(col.begin(), col.end());
    WRITE("COLUMN %-30s %c %17.9g %17.9g %4d",
          col.label.c_str(), col.type, minmax[0], minmax[1], col.dataset_id);
    if (!col.source.empty())
      WRITE("COLSRC %-30s %-36s  %4d",
            col.label.c_str(), col.source.c_str(), col.dataset_id);
  }
  WRITE("NDIF %8zu", datasets.size());
  for (const Dataset& ds : datasets) {
    WRITE("PROJECT %7d %s", ds.id, ds.project_name.c_str());
    WRITE("CRYSTAL %7d %s", ds.id, ds.crystal_name.c_str());
    WRITE("DATASET %7d %s", ds.id, ds.dataset_name.c_str());
    WRITE("DCELL %9d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",
          ds.id, ds.cell.a, ds.cell.b, ds.cell.c,
          ds.cell.alpha, ds.cell.beta, ds.cell.gamma);
    WRITE("DWAVEL %8d %10.5f", ds.id, ds.wavelength);
    //WRITE("BATCH");
  }
  WRITE("END");
  if (!history.empty()) {
    // According to mtzformat.html the file can have only up to 30 history
    // lines, but we don't enforce it here.
    WRITE("MTZHIST %3zu", history.size());
    for (const std::string& line : history)
      WRITE("%s", line.c_str());
  }
  WRITE("MTZENDOFHEADERS");
}

#undef WRITE

void Mtz::write_to_file(const std::string& path) const {
  fileptr_t f = file_open(path.c_str(), "wb");
  try {
    return write_to_stream(f.get());
  } catch (std::runtime_error& e) {
    fail(std::string(e.what()) + ": " + path);
  }
}

} // namespace gemmi

#endif // GEMMI_WRITE_IMPLEMENTATION

#ifdef  __INTEL_COMPILER
# pragma warning pop
#endif

#endif
