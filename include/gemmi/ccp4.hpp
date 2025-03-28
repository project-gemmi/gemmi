// Copyright 2018 Global Phasing Ltd.
//
// CCP4 format for maps and masks.

#ifndef GEMMI_CCP4_HPP_
#define GEMMI_CCP4_HPP_

#include <cassert>
#include <cmath>     // for ceil, fabs, floor, round
#include <cstdint>   // for uint16_t, int32_t
#include <cstdio>    // for FILE
#include <cstring>   // for memcpy
#include <array>
#include <string>
#include <type_traits>  // for is_same
#include <vector>
#include "symmetry.hpp"
#include "fail.hpp"      // for fail
#include "fileutil.hpp"  // for file_open, is_little_endian, ...
#include "input.hpp"     // for AnyStream, FileStream
#include "grid.hpp"

namespace gemmi {

using std::int32_t;

// options for Ccp4<>::setup
enum class MapSetup {
  Full,         // reorder and expand to the whole unit cell
  NoSymmetry,   // reorder and resize to the whole cell, but no symmetry ops
  ReorderOnly   // reorder axes to X, Y, Z
};

struct Ccp4Base {
  DataStats hstats;  // data statistics read from / written to ccp4 map
  // stores raw headers if the grid was read from ccp4 map
  std::vector<int32_t> ccp4_header;
  bool same_byte_order = true;

  // methods to access info from ccp4 headers, w is word number from the spec
  void* header_word(int w) { return &ccp4_header.at(w - 1); }
  const void* header_word(int w) const { return &ccp4_header.at(w - 1); }

  int32_t header_i32(int w) const {
    int32_t value = ccp4_header.at(w - 1);
    if (!same_byte_order)
      swap_four_bytes(&value);
    return value;
  }
  std::array<int, 3> header_3i32(int w) const {
    return {{ header_i32(w), header_i32(w+1), header_i32(w+2) }};
  }
  float header_float(int w) const {
    int32_t int_value = header_i32(w);
    float f;
    std::memcpy(&f, &int_value, 4);
    return f;
  }
  std::string header_str(int w, size_t len=80) const {
    if (4 * ccp4_header.size() < 4 * (w - 1) + len)
      fail("invalid end of string");
    return std::string(reinterpret_cast<const char*>(header_word(w)), len);
  }
  void set_header_i32(int w, int32_t value) {
    if (!same_byte_order)
      swap_four_bytes(&value);
    ccp4_header.at(w - 1) = value;
  }
  void set_header_3i32(int w, int32_t x, int32_t y, int32_t z) {
    set_header_i32(w, x);
    set_header_i32(w+1, y);
    set_header_i32(w+2, z);
  }
  void set_header_float(int w, float value) {
    int32_t int32_value;
    std::memcpy(&int32_value, &value, 4);
    set_header_i32(w, int32_value);
  }
  void set_header_str(int w, const std::string& str) {
    std::memcpy(header_word(w), str.c_str(), str.size());
  }

  std::array<int, 3> axis_positions() const {
    if (ccp4_header.empty())
      return {{0, 1, 2}}; // assuming it's X,Y,Z
    std::array<int, 3> pos{{-1, -1, -1}};
    for (int i = 0; i != 3; ++i) {
      int mapi = header_i32(17 + i);
      if (mapi <= 0 || mapi > 3 || pos[mapi - 1] != -1)
        fail("Incorrect MAPC/MAPR/MAPS records");
      pos[mapi - 1] = i;
    }
    return pos;
  }

  double header_rfloat(int w) const { // rounded to 5 digits
    return std::round(1e5 * header_float(w)) / 1e5;
  }

  Box<Fractional> get_extent() const {
    Box<Fractional> box;
    // cf. setup()
    auto pos = axis_positions();
    std::array<int, 3> start = header_3i32(5);
    std::array<int, 3> size = header_3i32(1);
    std::array<int, 3> sampl = header_3i32(8);
    for (int i = 0; i < 3; ++i) {
      double scale = 1. / sampl[i];
      int p = pos[i];
      box.minimum.at(i) = scale * start[p] - 1e-9;
      box.maximum.at(i) = scale * (start[p] + size[p] - 1) + 1e-9;
    }
    return box;
  }

  // Skew transformation (words 25-37) is supported by CCP4 maplib and PyMOL,
  // but it's not in the MRC format and is not supported by most programs.
  // From maplib.html: Skew transformation is from standard orthogonal
  // coordinate frame (as used for atoms) to orthogonal map frame, as
  //                            Xo(map) = S * (Xo(atoms) - t)
  bool has_skew_transformation() const {
    return header_i32(25) != 0;  // LSKFLG should be 0 or 1
  }
  Transform get_skew_transformation() const {
    return {
      // 26-34 SKWMAT
      { header_float(26), header_float(27), header_float(28),
        header_float(29), header_float(30), header_float(31),
        header_float(32), header_float(33), header_float(34) },
      // 35-37 SKWTRN
      { header_float(35), header_float(36), header_float(37) }
    };
  }

  // ORIGIN (words 50-52), used in MRC format, zeros in CCP4 format
  Position get_origin() const {
    return Position(header_float(50), header_float(51), header_float(52));
  }

  // this function assumes that the whole unit cell is covered with offset 0
  void prepare_ccp4_header_except_mode_and_stats(GridMeta& grid) {
    GroupOps ops;
    if (grid.spacegroup)
      ops = grid.spacegroup->operations();
    ccp4_header.clear();
    ccp4_header.resize(256 + ops.order() * 20, 0);
    set_header_3i32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
    set_header_3i32(5, 0, 0, 0); // NXSTART, NYSTART, NZSTART
    if (grid.axis_order == AxisOrder::XYZ)
      set_header_3i32(8, grid.nu, grid.nv, grid.nw);  // MX, MY, MZ
    else // grid.axis_order == AxisOrder::ZYX
      set_header_3i32(8, grid.nw, grid.nv, grid.nu);
    set_header_float(11, (float) grid.unit_cell.a);
    set_header_float(12, (float) grid.unit_cell.b);
    set_header_float(13, (float) grid.unit_cell.c);
    set_header_float(14, (float) grid.unit_cell.alpha);
    set_header_float(15, (float) grid.unit_cell.beta);
    set_header_float(16, (float) grid.unit_cell.gamma);
    if (grid.axis_order == AxisOrder::XYZ)
      set_header_3i32(17, 1, 2, 3); // MAPC, MAPR, MAPS
    else // grid.axis_order == AxisOrder::ZYX
      set_header_3i32(17, 3, 2, 1);
    set_header_i32(23, grid.spacegroup ? grid.spacegroup->ccp4 : 1); // ISPG
    set_header_i32(24, ops.order() * 80);  // NSYMBT
    set_header_str(27, "CCP4"); // EXTTYP
    set_header_i32(28, 20140);  // NVERSION
    set_header_str(53, "MAP ");
    set_header_i32(54, is_little_endian() ? 0x00004144 : 0x11110000); // MACHST
    set_header_i32(56, 1); // labels
    std::memset(header_word(57), ' ', 800 + ops.order() * 80);
    set_header_str(57, "written by GEMMI");
    int n = 257;
    for (Op op : ops) {
      set_header_str(n, op.triplet());
      n += 20;
    }
  }

  void update_header_mode_and_stats(int mode) {
    set_header_i32(4, mode);
    set_header_float(20, (float) hstats.dmin);
    set_header_float(21, (float) hstats.dmax);
    set_header_float(22, (float) hstats.dmean);
    set_header_float(55, (float) hstats.rms);
    // labels could be modified but it's not important
  }

  bool full_cell_(const GridMeta& grid) const {
    if (ccp4_header.empty())
      return true; // assuming it's full cell
    return
      // NXSTART et al. must be 0
      header_i32(5) == 0 && header_i32(6) == 0 && header_i32(7) == 0 &&
      // MX == NX
      header_i32(8) == grid.nu && header_i32(9) == grid.nv && header_i32(10) == grid.nw;
  }

  void read_ccp4_header_(GridMeta* grid, AnyStream& f, const std::string& path) {
    const size_t hsize = 256;
    ccp4_header.resize(hsize);
    if (!f.read(ccp4_header.data(), 4 * hsize))
      fail("Failed to read map header: " + path);
    if (header_str(53, 4) != "MAP ")
      fail("Not a CCP4 map: " + path);
    std::string machst = header_str(54, 4);
    if (machst[0] != 0x44 && machst[0] != 0x11)
      fail("Unsupported machine stamp (endianness) in the file?");
    same_byte_order = machst[0] == (is_little_endian() ? 0x44 : 0x11);
    size_t ext_w = header_i32(24) / 4;  // NSYMBT in words
    if (ext_w != 0) {
      if (ext_w > 1000000)
        fail("Unexpectedly long extended header: " + path);
      ccp4_header.resize(hsize + ext_w);
      if (!f.read(ccp4_header.data() + hsize, 4 * ext_w))
        fail("Failed to read extended header: " + path);
    }
    for (int i = 0; i < 3; ++i) {
      int axis = header_i32(17 + i);
      if (axis < 1 || axis > 3)
        fail("Unexpected axis value in word " + std::to_string(17 + i)
             + ": " + std::to_string(axis));
    }
    hstats.dmin = header_float(20);
    hstats.dmax = header_float(21);
    hstats.dmean = header_float(22);
    hstats.rms = header_float(55);
    if (grid) {
      grid->unit_cell.set(header_rfloat(11), header_rfloat(12), header_rfloat(13),
                          header_rfloat(14), header_rfloat(15), header_rfloat(16));
      grid->nu = header_i32(1);
      grid->nv = header_i32(2);
      grid->nw = header_i32(3);
      grid->spacegroup = find_spacegroup_by_number(header_i32(23));
      auto pos = axis_positions();
      grid->axis_order = AxisOrder::Unknown;
      if (pos[0] == 0 && pos[1] == 1 && pos[2] == 2 && full_cell_(*grid))
        grid->axis_order = AxisOrder::XYZ;
    }
  }
};

template<typename T=float>
struct Ccp4 : public Ccp4Base {
  Grid<T> grid;

  /// If the header is empty, prepare it; otherwise, update only MODE
  /// and, if update_stats==true, also DMIN, DMAX, DMEAN and RMS.
  void update_ccp4_header(int mode=-1, bool update_stats=true) {
    if (mode > 2 && mode != 6)
      fail("Only modes 0, 1, 2 and 6 are supported.");
    if (grid.point_count() == 0)
      fail("update_ccp4_header(): set the grid first (it has size 0)");
    if (grid.axis_order == AxisOrder::Unknown)
      fail("update_ccp4_header(): run setup() first");
    if (update_stats)
      hstats = calculate_data_statistics(grid.data);
    if (ccp4_header.empty())
      prepare_ccp4_header_except_mode_and_stats(grid);
    assert(ccp4_header.size() >= 256);
    if (mode < 0) {
      mode = mode_for_data();
      if (mode < 0)
        fail("update_ccp4_header: specify map mode explicitly (usually 2)");
    }
    update_header_mode_and_stats(mode);
  }

  static int mode_for_data() {
    if (std::is_same<T, std::int8_t>::value)
      return 0;
    if (std::is_same<T, std::int16_t>::value)
      return 1;
    if (std::is_same<T, float>::value)
      return 2;
    if (std::is_same<T, std::uint16_t>::value)
      return 6;
    return -1;
  }

  bool full_cell() const { return full_cell_(grid); }

  void read_ccp4_header(AnyStream& f, const std::string& path) {
    read_ccp4_header_(&grid, f, path);
    if (grid.axis_order != AxisOrder::Unknown)
      grid.calculate_spacing();
  }

  void setup(T default_value, MapSetup mode=MapSetup::Full);
  void set_extent(const Box<Fractional>& box);

  void read_ccp4_stream(AnyStream& f, const std::string& path);

  void read_ccp4_file(const std::string& path) {
    FileStream stream(path.c_str(), "rb");
    read_ccp4_stream(stream, path);
  }

  void read_ccp4_from_memory(const char* data, size_t size, const std::string& name) {
    MemoryStream stream(data, size);
    read_ccp4_stream(stream, name);
  }

  template<typename Input>
  void read_ccp4(Input&& input) {
    std::unique_ptr<AnyStream> stream = input.create_stream();
    read_ccp4_stream(*stream, input.path());
  }

  void write_ccp4_map(const std::string& path) const;
};


namespace impl {

template<typename From, typename To>
To translate_map_point(From f) { return static_cast<To>(f); }
// We convert map 2 to 0 by translating non-zero values to 1.
template<> inline
std::int8_t translate_map_point<float,std::int8_t>(float f) { return f != 0; }

template<typename TFile, typename TMem>
void read_data(AnyStream& f, std::vector<TMem>& content) {
  if (std::is_same<TFile, TMem>::value) {
    size_t len = content.size();
    if (!f.read(content.data(), sizeof(TMem) * len))
      fail("Failed to read all the data from the map file.");
  } else {
    constexpr size_t chunk_size = 64 * 1024;
    std::vector<TFile> work(chunk_size);
    for (size_t i = 0; i < content.size(); i += chunk_size) {
      size_t len = std::min(chunk_size, content.size() - i);
      if (!f.read(work.data(), sizeof(TFile) * len))
        fail("Failed to read all the data from the map file.");
      for (size_t j = 0; j < len; ++j)
        content[i+j] = translate_map_point<TFile,TMem>(work[j]);
    }
  }
}

template<typename TFile, typename TMem>
void write_data(const std::vector<TMem>& content, FILE* f) {
  if (std::is_same<TMem, TFile>::value) {
    size_t len = content.size();
    if (std::fwrite(content.data(), sizeof(TFile), len, f) != len)
      sys_fail("Failed to write data to the map file");
  } else {
    constexpr size_t chunk_size = 64 * 1024;
    std::vector<TFile> work(chunk_size);
    for (size_t i = 0; i < content.size(); i += chunk_size) {
      size_t len = std::min(chunk_size, content.size() - i);
      for (size_t j = 0; j < len; ++j)
        work[j] = static_cast<TFile>(content[i+j]);
      if (std::fwrite(work.data(), sizeof(TFile), len, f) != len)
        sys_fail("Failed to write data to the map file");
    }
  }
}

} // namespace impl

// This function was tested only on little-endian machines,
// let us know if you need support for other architectures.
template<typename T>
void Ccp4<T>::read_ccp4_stream(AnyStream& f, const std::string& path) {
  read_ccp4_header(f, path);
  grid.data.resize(grid.point_count());
  int mode = header_i32(4);
  if (mode == 0)
    impl::read_data<std::int8_t>(f, grid.data);
  else if (mode == 1)
    impl::read_data<std::int16_t>(f, grid.data);
  else if (mode == 2)
    impl::read_data<float>(f, grid.data);
  else if (mode == 6)
    impl::read_data<std::uint16_t>(f, grid.data);
  else
    fail("Mode " + std::to_string(mode) + " is not supported "
         "(only 0, 1, 2 and 6 are supported).");
  //if (std::fgetc(f) != EOF)
  //  fail("The map file is longer then expected.");

  if (!same_byte_order) {
    if (sizeof(T) == 2)
      for (T& value : grid.data)
        swap_two_bytes(&value);
    else if (sizeof(T) == 4)
      for (T& value : grid.data)
        swap_four_bytes(&value);
  }
}

template<typename T>
void Ccp4<T>::setup(T default_value, MapSetup mode) {
  if (grid.axis_order == AxisOrder::XYZ || ccp4_header.empty())
    return;
  // cell sampling does not change
  const std::array<int, 3> sampl = header_3i32(8);
  // get old metadata
  const std::array<int, 3> pos = axis_positions();
  std::array<int, 3> start = header_3i32(5);
  int end[3] = { start[0] + grid.nu, start[1] + grid.nv, start[2] + grid.nw };
  // set new metadata
  if (mode == MapSetup::ReorderOnly) {
    set_header_3i32(5, start[pos[0]], start[pos[1]], start[pos[2]]);
    for (int i = 0; i < 3; ++i) {
      end[i] -= start[i];
      start[i] = 0;
    }
    int crs[3] = { grid.nu, grid.nv, grid.nw };
    grid.nu = crs[pos[0]];
    grid.nv = crs[pos[1]];
    grid.nw = crs[pos[2]];
  } else {
    grid.nu = sampl[0];
    grid.nv = sampl[1];
    grid.nw = sampl[2];
    set_header_3i32(5, 0, 0, 0); // start
  }
  set_header_3i32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
  set_header_3i32(17, 1, 2, 3); // axes (MAPC, MAPR, MAPS)
  grid.axis_order = full_cell() ? AxisOrder::XYZ : AxisOrder::Unknown;
  if (grid.axis_order == AxisOrder::XYZ)
    grid.calculate_spacing();

  // now set the data
  {
    std::vector<T> full(grid.point_count(), default_value);
    int it[3];
    int idx = 0;
    for (it[2] = start[2]; it[2] < end[2]; it[2]++) // sections
      for (it[1] = start[1]; it[1] < end[1]; it[1]++) // rows
        for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
          T val = grid.data[idx++];
          size_t new_index = grid.index_s(it[pos[0]], it[pos[1]], it[pos[2]]);
          full[new_index] = val;
        }
    grid.data = std::move(full);
  }

  if (mode == MapSetup::Full &&
      // no need to apply symmetry if we started with the whole cell
      (end[pos[0]] - start[pos[0]] < sampl[0] ||
       end[pos[1]] - start[pos[1]] < sampl[1] ||
       end[pos[2]] - start[pos[2]] < sampl[2]))
    grid.symmetrize_nondefault(default_value);
}

template<typename T>
void Ccp4<T>::set_extent(const Box<Fractional>& box) {
  if (ccp4_header.empty())
    fail("set_extent(): no header in the map. Call update_ccp4_header() first");
  if (!full_cell())
    fail("Ccp4::set_extent() works only after setup()");
  if (grid.axis_order != AxisOrder::XYZ)
    fail("Ccp4::set_extent() works only with XYZ order");
  int u0 = (int)std::ceil(box.minimum.x * grid.nu);
  int v0 = (int)std::ceil(box.minimum.y * grid.nv);
  int w0 = (int)std::ceil(box.minimum.z * grid.nw);
  int nu = (int)std::floor(box.maximum.x * grid.nu) - u0 + 1;
  int nv = (int)std::floor(box.maximum.y * grid.nv) - v0 + 1;
  int nw = (int)std::floor(box.maximum.z * grid.nw) - w0 + 1;
  // set the data
  std::vector<T> new_data((size_t)nu * nv * nw);
  grid.get_subarray(new_data.data(), {u0, v0, w0}, {nu, nv, nw});
  grid.data.swap(new_data);
  // and metadata
  grid.nu = nu;
  grid.nv = nv;
  grid.nw = nw;
  set_header_3i32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
  set_header_3i32(5, u0, v0, w0);
  // AxisOrder::XYZ is used only for grid covering full cell
  grid.axis_order = AxisOrder::Unknown;
}

template<typename T>
void Ccp4<T>::write_ccp4_map(const std::string& path) const {
  assert(ccp4_header.size() >= 256);
  fileptr_t f = file_open(path.c_str(), "wb");
  std::fwrite(ccp4_header.data(), 4, ccp4_header.size(), f.get());
  int mode = header_i32(4);
  if (mode == 0)
    impl::write_data<std::int8_t>(grid.data, f.get());
  else if (mode == 1)
    impl::write_data<std::int16_t>(grid.data, f.get());
  else if (mode == 2)
    impl::write_data<float>(grid.data, f.get());
  else if (mode == 6)
    impl::write_data<std::uint16_t>(grid.data, f.get());
}

GEMMI_DLL Ccp4<float> read_ccp4_map(const std::string& path, bool setup);
GEMMI_DLL Ccp4<int8_t> read_ccp4_mask(const std::string& path, bool setup);
GEMMI_DLL Ccp4Base read_ccp4_header(const std::string& path);

} // namespace gemmi
#endif
