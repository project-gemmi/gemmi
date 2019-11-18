// Copyright 2018 Global Phasing Ltd.
//
// CCP4 format for maps and masks.

#ifndef GEMMI_CCP4_HPP_
#define GEMMI_CCP4_HPP_

#include <cassert>
#include <cmath>     // for NAN, sqrt
#include <cstdint>   // for uint16_t, uint32_t
#include <cstdio>    // for FILE
#include <cstring>   // for memcpy
#include <array>
#include <string>
#include <typeinfo>  // for typeid
#include <vector>
#include "symmetry.hpp"
#include "fail.hpp"      // for fail
#include "fileutil.hpp"  // for file_open, is_little_endian, ...
#include "input.hpp"     // for FileStream
#include "grid.hpp"

namespace gemmi {

using std::int32_t;

// options for Ccp4<>::setup
enum class GridSetup {
  ReorderOnly,  // reorder axes to X, Y, Z
  ResizeOnly,   // reorder and resize to the whole cell, but no symmetry ops
  Full,         // reorder and expand to the whole unit cell
  FullCheck     // additionally consistency of redundant data
};

template<typename T=float>
struct Ccp4 {
  Grid<T> grid;
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
  float header_float(int w) const {
    int32_t int_value = header_i32(w);
    float f;
    std::memcpy(&f, &int_value, 4);
    return f;
  }
  // ccp4 map header has mostly 80-byte strings
  std::string header_str(int w, size_t len=80) const {
    if (4 * ccp4_header.size() < 4 * (w - 1) + len)
      fail("invalid end of string");
    return std::string(reinterpret_cast<const char*>(header_word(w)), len);
  }
  void set_header_i32(int w, int32_t value) {
    if (!same_byte_order)
      swap_four_bytes(&value);
    ccp4_header.at(w - 1) = value;
  };
  void set_header_3i32(int w, int32_t x, int32_t y, int32_t z) {
    set_header_i32(w, x);
    set_header_i32(w+1, y);
    set_header_i32(w+2, z);
  };
  void set_header_float(int w, float value) {
    int32_t int32_value;
    std::memcpy(&int32_value, &value, 4);
    set_header_i32(w, int32_value);
  };
  void set_header_str(int w, const std::string& str) {
    std::memcpy(header_word(w), str.c_str(), str.size());
  }

  void prepare_ccp4_header(int mode) {
    GroupOps ops;
    if (grid.spacegroup)
      ops = grid.spacegroup->operations();
    ccp4_header.clear();
    ccp4_header.resize(256 + ops.order() * 20, 0);
    set_header_3i32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
    set_header_3i32(5, 0, 0, 0); // NXSTART, NYSTART, NZSTART
    if (grid.hkl_orient == HklOrient::HKL)
      set_header_3i32(8, grid.nu, grid.nv, grid.nw);  // MX, MY, MZ
    else // grid.hkl_orient == HklOrient::LKH
      set_header_3i32(8, grid.nw, grid.nv, grid.nu);
    set_header_float(11, (float) grid.unit_cell.a);
    set_header_float(12, (float) grid.unit_cell.b);
    set_header_float(13, (float) grid.unit_cell.c);
    set_header_float(14, (float) grid.unit_cell.alpha);
    set_header_float(15, (float) grid.unit_cell.beta);
    set_header_float(16, (float) grid.unit_cell.gamma);
    if (grid.hkl_orient == HklOrient::HKL)
      set_header_3i32(17, 1, 2, 3); // MAPC, MAPR, MAPS
    else // grid.hkl_orient == HklOrient::LKH
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
    for (const Op& op : ops) {
      set_header_str(n, op.triplet());
      n += 20;
    }
    update_ccp4_header(mode);
  }

  void update_ccp4_header(int mode, bool update_stats=false) {
    if (update_stats)
      hstats = calculate_data_statistics(grid.data);
    if (mode != 0 && mode != 1 && mode != 2 && mode != 6)
      fail("Only modes 0, 1, 2 and 6 are supported.");
    if (ccp4_header.empty()) {
      prepare_ccp4_header(mode);
      return;
    }
    assert(ccp4_header.size() >= 256);
    set_header_i32(4, mode);
    set_header_float(20, (float) hstats.dmin);
    set_header_float(21, (float) hstats.dmax);
    set_header_float(22, (float) hstats.dmean);
    set_header_float(55, (float) hstats.rms);
    // labels could be modified but it's not important
  }

  bool full_cell() const {
    if (ccp4_header.empty())
      return true; // assuming it's full cell
    return
      // NXSTART et al. must be 0
      header_i32(5) == 0 && header_i32(6) == 0 && header_i32(7) == 0 &&
      // MX == NX
      header_i32(8) == grid.nu && header_i32(9) == grid.nv &&
                                  header_i32(10) == grid.nw &&
      // just in case, check ORIGIN
      header_i32(50) == 0 && header_i32(51) == 0 && header_i32(52) == 0;
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

  template<typename Stream>
  void read_ccp4_header(Stream f, const std::string& path) {
    const size_t hsize = 256;
    ccp4_header.resize(hsize);
    if (!f.read(ccp4_header.data(), 4 * hsize))
      fail("Failed to read map header: " + path);
    if (header_str(53, 4) != "MAP ")
      fail("Not a CCP4 map: " + path);
    std::string machst = header_str(54, 4);
    if (machst[0] != 0x44 && machst[0] != 0x11)
      fail("Unsupported machine stamp (endiannes) in the file?");
    same_byte_order = machst[0] == (is_little_endian() ? 0x44 : 0x11);
    grid.unit_cell.set(header_rfloat(11), header_rfloat(12), header_rfloat(13),
                       header_rfloat(14), header_rfloat(15), header_rfloat(16));
    size_t ext_w = header_i32(24) / 4;  // NSYMBT in words
    if (ext_w != 0) {
      if (ext_w > 1000000)
        fail("Unexpectedly long extended header: " + path);
      ccp4_header.resize(hsize + ext_w);
      if (!f.read(ccp4_header.data() + hsize, 4 * ext_w))
        fail("Failed to read extended header: " + path);
    }
    grid.nu = header_i32(1);
    grid.nv = header_i32(2);
    grid.nw = header_i32(3);
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
    grid.spacegroup = find_spacegroup_by_number(header_i32(23));
    auto pos = axis_positions();
    grid.full_canonical = pos[0] == 0 && pos[1] == 1 && pos[2] == 2 &&
                          full_cell();
  }

  double setup(GridSetup mode, T default_value);

  template<typename Stream>
  void read_ccp4_stream(Stream f, const std::string& path);

  void read_ccp4_file(const std::string& path) {
    fileptr_t f = file_open(path.c_str(), "rb");
    read_ccp4_stream(FileStream{f.get()}, path);
  }

  inline void read_ccp4_from_memory(const char* data, size_t size,
                                    const std::string& name) {
    return read_ccp4_stream(MemoryStream{data, data + size}, name);
  }

  template<typename Input>
  void read_ccp4(Input&& input) {
    if (input.is_stdin())
      return read_ccp4_stream(FileStream{stdin}, "stdin");
    if (input.is_compressed())
      return read_ccp4_stream(input.get_uncompressing_stream(), input.path());
    return read_ccp4_file(input.path());
  }

  void write_ccp4_map(const std::string& path) const;
};


namespace impl {

template<typename Stream, typename TFile, typename TMem>
void read_data(Stream& f, std::vector<TMem>& content) {
  if (typeid(TFile) == typeid(TMem)) {
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
        content[i+j] = static_cast<TMem>(work[j]);
    }
  }
}

template<typename TFile, typename TMem>
void write_data(const std::vector<TMem>& content, FILE* f) {
  if (typeid(TMem) == typeid(TFile)) {
    size_t len = content.size();
    if (std::fwrite(content.data(), sizeof(TFile), len, f) != len)
      fail("Failed to write data to the map file.");
  } else {
    constexpr size_t chunk_size = 64 * 1024;
    std::vector<TFile> work(chunk_size);
    for (size_t i = 0; i < content.size(); i += chunk_size) {
      size_t len = std::min(chunk_size, content.size() - i);
      for (size_t j = 0; j < len; ++j)
        work[j] = static_cast<TFile>(content[i+j]);
      if (std::fwrite(work.data(), sizeof(TFile), len, f) != len)
        fail("Failed to write data to the map file.");
    }
  }
}

} // namespace impl

// This function was tested only on little-endian machines,
// let us know if you need support for other architectures.
template<typename T> template<typename Stream>
void Ccp4<T>::read_ccp4_stream(Stream f, const std::string& path) {
  read_ccp4_header(f, path);
  grid.data.resize(grid.nu * grid.nv * grid.nw);
  int mode = header_i32(4);
  if (mode == 0)
    impl::read_data<Stream, std::int8_t>(f, grid.data);
  else if (mode == 1)
    impl::read_data<Stream, std::int16_t>(f, grid.data);
  else if (mode == 2)
    impl::read_data<Stream, float>(f, grid.data);
  else if (mode == 6)
    impl::read_data<Stream, std::uint16_t>(f, grid.data);
  else
    fail("Only modes 0, 1, 2 and 6 are supported.");
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

namespace impl {

template<typename T> bool is_same(T a, T b) { return a == b; }
template<> inline bool is_same(float a, float b) {
  return std::isnan(b) ? std::isnan(a) : a == b;
}
template<> inline bool is_same(double a, double b) {
  return std::isnan(b) ? std::isnan(a) : a == b;
}

}

template<typename T>
double Ccp4<T>::setup(GridSetup mode, T default_value) {
  double max_error = 0.0;
  if (grid.full_canonical || ccp4_header.empty())
    return max_error;
  // cell sampling does not change
  int sampl[3] = { header_i32(8), header_i32(9), header_i32(10) };
  // get old metadata
  auto pos = axis_positions();
  int start[3] = { header_i32(5), header_i32(6), header_i32(7) };
  int end[3] = { start[0] + grid.nu, start[1] + grid.nv, start[2] + grid.nw };
  // set new metadata
  if (mode == GridSetup::ReorderOnly) {
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
  // now set the data
  std::vector<T> full(grid.nu * grid.nv * grid.nw, default_value);
  int it[3];
  int idx = 0;
  for (it[2] = start[2]; it[2] < end[2]; it[2]++) // sections
    for (it[1] = start[1]; it[1] < end[1]; it[1]++) // rows
      for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
        T val = grid.data[idx++];
        int new_index = grid.index_s(it[pos[0]], it[pos[1]], it[pos[2]]);
        full[new_index] = val;
      }
  grid.data = std::move(full);
  if (mode == GridSetup::Full) {
    grid.full_canonical = true;
    grid.symmetrize([&default_value](T a, T b) {
        return impl::is_same(a, default_value) ? b : a;
    });
  } else if (mode == GridSetup::FullCheck) {
    grid.full_canonical = true;
    grid.symmetrize([&max_error, &default_value](T a, T b) {
        if (impl::is_same(a, default_value)) {
          return b;
        } else {
          if (!impl::is_same(b, default_value))
            max_error = std::max(max_error, std::fabs(double(a - b)));
          return a;
        }
    });
  } else {
    grid.full_canonical = pos[0] == 0 && pos[1] == 1 && pos[2] == 2 &&
                          full_cell();
  }
  return max_error;
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

} // namespace gemmi
#endif
