// Copyright 2017 Global Phasing Ltd.
//
// 3d grid for volume density. CCP4 format for maps and masks.

#ifndef GEMMI_GRID_HH
#define GEMMI_GRID_HH

#include <cassert>
#include <cmath>     // for NAN, sqrt
#include <cstdint>   // for uint16_t, uint32_t
#include <cstdio>    // for FILE, fread
#include <cstring>   // for memcpy
#include <string>
#include <typeinfo>  // for typeid
#include <vector>
#include "unitcell.hpp"
#include "util.hpp"  // for fail, file_open

namespace gemmi {

inline bool is_little_endian() {
  std::uint32_t x = 1;
  return *reinterpret_cast<char *>(&x) == 1;
}

struct GridStats {
  double dmin = NAN;
  double dmax = NAN;
  double dmean = NAN;
  double rms = NAN;
};

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
template<typename T=float>
struct Grid {
  int nu, nv, nw;
  mol::UnitCell unit_cell;
  std::vector<T> data;
  double spacing[3];
  GridStats stats;
  // stores raw headers if the grid was read from ccp4 map
  std::vector<char> ccp4_header;

  void set_size(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
    spacing[0] = 1.0 / (nu * unit_cell.ar);
    spacing[1] = 1.0 / (nv * unit_cell.br);
    spacing[2] = 1.0 / (nw * unit_cell.cr);
  }
  void set_spacing(double sp) {
    set_size(iround(unit_cell.a / sp),
             iround(unit_cell.b / sp),
             iround(unit_cell.c / sp));
  }

  T& node(int u, int v, int w) {
#if 1
    if (u >= nu)
      u -= nu;
    else if (u < 0)
      u += nu;
    if (v >= nv)
      v -= nv;
    else if (v < 0)
      v += nv;
    if (w >= nw)
      w -= nw;
    else if (w < 0)
      w += nw;
#endif
    int idx = w * nu * nv + v * nu + u;
    return data[idx];
  }

  void set_points_around(const mol::Position& ctr, double radius, T value) {
    mol::Position fctr = unit_cell.fractionalize(ctr);
    fctr.x -= std::floor(fctr.x);
    fctr.y -= std::floor(fctr.y);
    fctr.z -= std::floor(fctr.z);
#if 1
    int du = (int) std::ceil(radius / spacing[0]);
    int dv = (int) std::ceil(radius / spacing[1]);
    int dw = (int) std::ceil(radius / spacing[2]);
    int u0 = iround(fctr.x * nu);
    int v0 = iround(fctr.y * nv);
    int w0 = iround(fctr.z * nw);
    for (int w = w0-dw; w <= w0+dw; ++w)
      for (int v = v0-dv; v <= v0+dv; ++v)
        for (int u = u0-du; u < u0+du; ++u) {
#else
    for (int w = 0; w < nw; ++w)
      for (int v = 0; v < nv; ++v)
        for (int u = 0; u < nu; ++u) {
#endif
          mol::Position fdelta = {fctr.x - double(u) / nu,
                                  fctr.y - double(v) / nv,
                                  fctr.z - double(w) / nw};
          if (fdelta.x > 0.5)
            fdelta.x -= 1.0;
          else if (fdelta.x < -0.5)
            fdelta.x += 1.0;
          if (fdelta.y > 0.5)
            fdelta.y -= 1.0;
          else if (fdelta.y < -0.5)
            fdelta.y += 1.0;
          if (fdelta.z > 0.5)
            fdelta.z -= 1.0;
          else if (fdelta.z < -0.5)
            fdelta.z += 1.0;
          mol::Position d = unit_cell.orthogonalize(fdelta);
          if (d.x*d.x + d.y*d.y + d.z*d.z < radius*radius) {
            node(u, v, w) = value;
          }
        }
  }

  GridStats calculate_statistics() const;
  void read_ccp4(const std::string& path);
  void write_ccp4_map(const std::string& path, int mode=2) const;
  size_t write_ccp4_mask(const std::string& path, double threshold) const;

  // methods to access info from ccp4 headers, w is word number from the spec
  int32_t header_i32(int w) const {
    return *reinterpret_cast<const int32_t*>(ccp4_header.data() + 4 * (w - 1));
  }
  float header_float(int w) const {
    return *reinterpret_cast<const float*>(ccp4_header.data() + 4 * (w - 1));
  }
  // ccp4 map header has mostly 80-byte strings
  std::string header_str(int w, size_t len=80) const {
    return std::string(ccp4_header.data() + 4 * (w - 1),
                       std::min(len, ccp4_header.size() - w));
  }
};


namespace impl {

template<typename TFile, typename TMem>
void read_data(FILE* f, std::vector<TMem>& content) {
  if (typeid(TFile) == typeid(TMem)) {
    size_t len = content.size();
    if (std::fread(content.data(), sizeof(TMem), len, f) != len)
      fail("Failed to read all the data from the map file.");
  } else {
    constexpr size_t chunk_size = 64 * 1024;
    std::vector<TFile> work(chunk_size);
    for (size_t i = 0; i < content.size(); i += chunk_size) {
      size_t len = std::min(chunk_size, content.size() - i);
      if (std::fread(work.data(), sizeof(TFile), len, f) != len)
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

template<typename Out, typename In>
void write_arrays(const std::string& path, const std::vector<char>& header,
                  const std::vector<In>& content) {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "wb");
  std::fwrite(header.data(), sizeof(char), header.size(), f.get());
  write_data<Out>(content, f.get());
}

template<typename T>
std::vector<char> make_ccp4_header(const Grid<T>& grid, int mode) {
  int nsym = 1;
  std::vector<char> header(1024 + nsym * 80, 0);
  auto word_ptr = [&header](int w) { return header.data() + 4 * (w - 1); };
  auto set_i32 = [&word_ptr](int word, int32_t value) {
    *reinterpret_cast<int32_t*>(word_ptr(word)) = value;
  };
  auto set_tri_i32 = [&set_i32](int word, int32_t x, int32_t y, int32_t z) {
    set_i32(word, x);
    set_i32(word+1, y);
    set_i32(word+2, z);
  };
  auto set_float = [&word_ptr](int word, float value) {
    *reinterpret_cast<float*>(word_ptr(word)) = value;
  };
  set_tri_i32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
  set_i32(4, mode);
  set_tri_i32(5, 0, 0, 0); // NXSTART, NYSTART, NZSTART
  set_tri_i32(8, grid.nu, grid.nv, grid.nw);  // MX, MY, MZ
  set_float(11, (float) grid.unit_cell.a);
  set_float(12, (float) grid.unit_cell.b);
  set_float(13, (float) grid.unit_cell.c);
  set_float(14, (float) grid.unit_cell.alpha);
  set_float(15, (float) grid.unit_cell.beta);
  set_float(16, (float) grid.unit_cell.gamma);
  set_tri_i32(17, 1, 2, 3); // MAPC, MAPR, MAPS
  set_float(20, (float) grid.stats.dmin);
  set_float(21, (float) grid.stats.dmax);
  set_float(22, (float) grid.stats.dmean);
  set_i32(23, 1); // ISPG
  set_i32(24, nsym * 80);
  std::memcpy(word_ptr(27), "CCP4", 4); // EXTTYP
  //set_i32(28, nversion);
  std::memcpy(word_ptr(53), "MAP ", 4);
  set_i32(54, is_little_endian() ? 0x00004144 : 0x11110000); // MACHST
  set_float(55, (float) grid.stats.rms);
  set_i32(56, 1); // labels
  memset(word_ptr(57), ' ', 800 + nsym * 80);
  strcpy(word_ptr(57), "from unfinished GEMMI code");
  memcpy(word_ptr(257), "X,  Y,  Z", 9);
  return header;
}

template<typename T>
std::vector<char> update_ccp4_header(const Grid<T>& grid, int mode) {
  if (grid.ccp4_header.empty())
    return make_ccp4_header(grid, mode);
  std::vector<char> header = grid.ccp4_header;
  assert(header.size() >= 1024);
  // selectively copy-pasted from make_ccp4_header()
  auto word_ptr = [&header](int w) { return header.data() + 4 * (w - 1); };
  auto set_i32 = [&word_ptr](int word, int32_t value) {
    *reinterpret_cast<int32_t*>(word_ptr(word)) = value;
  };
  auto set_float = [&word_ptr](int word, float value) {
    *reinterpret_cast<float*>(word_ptr(word)) = value;
  };
  set_i32(4, mode);
  set_float(20, (float) grid.stats.dmin);
  set_float(21, (float) grid.stats.dmax);
  set_float(22, (float) grid.stats.dmean);
  set_float(55, (float) grid.stats.rms);
  // labels could be modified but it's not important
  return header;
}

} // namespace impl

template<typename T>
GridStats Grid<T>::calculate_statistics() const {
  GridStats st;
  if (data.empty())
    return st;
  double sum = 0;
  double sq_sum = 0;
  st.dmin = st.dmax = data[0];
  for (float d : data) {
    sum += d;
    sq_sum += d * d;
    if (d < st.dmin)
      st.dmin = d;
    if (d > st.dmax)
      st.dmax = d;
  }
  st.dmean = sum / data.size();
  st.rms = std::sqrt(sq_sum / data.size() - st.dmean * st.dmean);
  return st;
}

// This function was tested only on little-endian machines,
// let us know if you need support for other architectures.
template<typename T>
void Grid<T>::read_ccp4(const std::string& path) {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "rb");
  const size_t hsize = 1024;
  ccp4_header.resize(hsize);
  if (std::fread(ccp4_header.data(), 1, hsize, f.get()) != hsize)
    fail("Failed to read map header: " + path);
  if (ccp4_header[208] != 'M' || ccp4_header[209] != 'A' ||
      ccp4_header[210] != 'P' || ccp4_header[211] != ' ')
    fail("Not a CCP4 map: " + path);
  int mode = header_i32(4);
  unit_cell.set(header_float(11), header_float(12), header_float(13),
                header_float(14), header_float(15), header_float(16));
  size_t nsymbt = header_i32(24);
  if (nsymbt > 1000000)
    fail("Unexpectedly long extendended header: " + path);
  ccp4_header.resize(hsize + nsymbt);
  if (std::fread(ccp4_header.data() + hsize, 1, nsymbt, f.get()) != nsymbt)
    fail("Failed to read extended header: " + path);
  nu = header_i32(1);
  nv = header_i32(2);
  nw = header_i32(3);
  for (int i = 17; i < 20; ++i) {
    int axis = header_i32(17);
    if (axis < 1 || axis > 3)
      fail("Unexpected axis value in word " + std::to_string(i));
  }
  stats.dmin = header_float(20);
  stats.dmax = header_float(21);
  stats.dmean = header_float(22);
  stats.rms = header_float(55);

  data.resize(nu * nv * nw);
  if (mode == 0)
    impl::read_data<signed char>(f.get(), data);
  else if (mode == 1)
    impl::read_data<int16_t>(f.get(), data);
  else if (mode == 2)
    impl::read_data<float>(f.get(), data);
  else if (mode == 6)
    impl::read_data<std::uint16_t>(f.get(), data);
  else
    fail("Only modes 0, 1, 2 and 6 are supported.");
#if 0
  if (std::fgetc(f.get()) != EOF)
    fail("The map file is longer then expected.");
#endif
}

template<typename T>
void Grid<T>::write_ccp4_map(const std::string& path, int mode) const {
  std::vector<char> header = impl::update_ccp4_header(*this, mode);
  if (mode == 0)
    impl::write_arrays<signed char>(path, header, data);
  else if (mode == 1)
    impl::write_arrays<int16_t>(path, header, data);
  else if (mode == 2)
    impl::write_arrays<float>(path, header, data);
  else if (mode == 6)
    impl::write_arrays<std::uint16_t>(path, header, data);
  else
    fail("Only modes 0, 1, 2 and 6 are supported.");
}

template<typename T>
size_t Grid<T>::write_ccp4_mask(const std::string& path,
                                double threshold) const {
  size_t counter = 0;
  std::vector<signed char> v(data.size());
  for (size_t i = 0; i < data.size(); i++) {
    v[i] = data[i] > threshold ? 1 : 0;
    counter += v[i];
  }
  impl::write_arrays<signed char>(path, impl::update_ccp4_header(*this, 0), v);
  return counter;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
