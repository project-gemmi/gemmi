// Copyright 2017 Global Phasing Ltd.
//
// 3d grid for volume density. CCP4 format for maps and masks.

// TODO: handle big-endian files

#ifndef GEMMI_GRID_HH
#define GEMMI_GRID_HH

#include <cmath>     // for NAN, sqrt
#include <cstdio>    // for FILE
#include <cstring>   // for memcpy
#include <memory>    // for unique_ptr
#include <vector>
#include "unitcell.hpp"
#include "util.hpp"  // for fail

namespace gemmi {

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
template<typename T=float>
struct Grid {
  int nu, nv, nw;
  mol::UnitCell unit_cell;
  std::vector<T> data;
  // statistics
  double dmin = NAN, dmax = NAN, dmean = NAN, rms = NAN;
  // stores raw headers if the grid was read from ccp4 map
  std::vector<char> ccp4_header;

  void set_size(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
  }
  void set_spacing(double sp) {
    set_size(iround(unit_cell.a / sp),
             iround(unit_cell.b / sp),
             iround(unit_cell.c / sp));
  }

  void set_points_around(const mol::Position& pos, double radius, T value) {
    mol::Position corners[8];
    int i = 0;
    for (double x : {pos.x - radius, pos.x + radius})
      for (double y : {pos.y - radius, pos.y + radius})
        for (double z : {pos.z - radius, pos.z + radius})
          corners[i++] = {x, y, z};
    //TODO
    int cnt = 0;
    for (int u = 0; u < nu; ++u)
      for (int v = 0; v < nv; ++v)
        for (int w = 0; w < nw; ++w) {
          mol::Position fr = {double(u) / nu, double(v) / nv, double(w) / nw};
          mol::Position orth = unit_cell.orthogonalize(fr);
          double dx = pos.x - orth.x;
          double dy = pos.y - orth.y;
          double dz = pos.z - orth.z;
          if (dx*dx + dy*dy + dz*dz < radius*radius) {
            ++cnt;
            data[u*nv*nw + v*nw + w] = value;
          }
        }
    fprintf(stderr, "radius: %g, points set: %d\n", radius, cnt);
  }

  void calculate_statistics();
  void read_ccp4(const std::string& path);
  void write_ccp4_map(const std::string& path, int mode=2) const;
  void write_ccp4_mask(const std::string& path, double threshold) const;

private:
  uint32_t header_u32(int w) {
    return *reinterpret_cast<uint32_t*>(ccp4_header.data() + 4 * (w - 1));
  }
  uint32_t header_float(int w) {
    return *reinterpret_cast<float*>(ccp4_header.data() + 4 * (w - 1));
  }
};


namespace impl {

typedef std::unique_ptr<FILE, decltype(&std::fclose)> fileptr_t;

template<typename T>
void write_arrays(const std::string& path, const std::vector<char>& header,
                  const std::vector<T>& content) {
  fileptr_t f(std::fopen(path.c_str(), "wb"), &std::fclose);
  if (!f)
    fail("Failed to open file for writing: " + path);
  std::fwrite(header.data(), sizeof(char), header.size(), f.get());
  std::fwrite(content.data(), sizeof(T), content.size(), f.get());
}

template<typename T>
std::vector<char> make_ccp4_header(const Grid<T>& grid, int mode) {
  int nsymbt = 0;
  std::vector<char> header(1024 + nsymbt, 0);
  auto word_ptr = [&header](int w) { return header.data() + 4 * (w - 1); };
  auto set_u32 = [&word_ptr](int word, uint32_t value) {
    *reinterpret_cast<uint32_t*>(word_ptr(word)) = value;
  };
  auto set_tri_u32 = [&set_u32](int word, uint32_t x, uint32_t y, uint32_t z) {
    set_u32(word, x);
    set_u32(word+1, y);
    set_u32(word+2, z);
  };
  auto set_float = [&word_ptr](int word, float value) {
    *reinterpret_cast<float*>(word_ptr(word)) = value;
  };
  set_tri_u32(1, grid.nu, grid.nv, grid.nw); // NX, NY, NZ
  set_u32(4, mode);
  set_tri_u32(5, 0, 0, 0); // NXSTART, NYSTART, NZSTART
  set_tri_u32(8, grid.nu, grid.nv, grid.nw);  // MX, MY, MZ
  set_float(11, grid.unit_cell.a);
  set_float(12, grid.unit_cell.b);
  set_float(13, grid.unit_cell.c);
  set_float(14, grid.unit_cell.alpha);
  set_float(15, grid.unit_cell.beta);
  set_float(16, grid.unit_cell.gamma);
  set_tri_u32(17, 1, 2, 3); // MAPC, MAPR, MAPS
  set_float(20, grid.dmin);
  set_float(21, grid.dmax);
  set_float(22, grid.dmean);
  set_u32(23, 1); // ISPG
  set_u32(24, nsymbt);
  std::memcpy(word_ptr(27), "CCP4", 4); // EXTTYP
  //set_u32(28, nversion);
  std::memcpy(word_ptr(53), "MAP ", 4);
  set_u32(54, 0x00004444); // MACHST for little endian (0x11110000 for BE)
  set_float(55, grid.rms);
  set_u32(56, 0); // labels
  return header;
}

} // namespace impl

template<typename T>
void Grid<T>::calculate_statistics() {
  if (data.empty())
    return;
  double sum = 0;
  double sq_sum = 0;
  dmin = dmax = data[0];
  for (float d : data) {
    sum += d;
    sq_sum += d * d;
    if (d < dmin)
      dmin = d;
    if (d > dmax)
      dmax = d;
  }
  dmean = sum / data.size();
  rms = std::sqrt(sq_sum / data.size() - dmean * dmean);
}

template<typename T>
void Grid<T>::read_ccp4(const std::string& path) {
  impl::fileptr_t f(std::fopen(path.c_str(), "rb"), &std::fclose);
  if (!f)
    fail("Failed to open file: " + path);
  const size_t hsize = 1024;
  ccp4_header.resize(hsize);
  if (fread(ccp4_header.data(), 1, hsize, f.get()) != hsize)
    fail("Failed to read map header: " + path);
  if (ccp4_header[208] != 'M' || ccp4_header[209] != 'A' ||
      ccp4_header[210] != 'P' || ccp4_header[211] != ' ')
    fail("Not a CCP4 map: " + path);
  uint32_t mode = header_u32(4);
  if (mode != 2 && mode != 0)
    fail("Only Mode 2 and Mode 0 of CCP4 map is supported.");
  uint32_t nsymbt = header_u32(24);
  if (nsymbt > 1000000)
    fail("Unexpectedly long extendended header: " + path);
  ccp4_header.resize(hsize + nsymbt);
  if (fread(ccp4_header.data() + hsize, 1, nsymbt, f.get()) != nsymbt)
    fail("Failed to read extended header: " + path);
  if (mode == 2) {
    size_t len = 16 ; // n_crs[0]*n_crs[1]*n_crs[2]
    data.resize(len);
    //if (c-r-s in normal order)
      if (fread(data.data(), 4, len, f.get()) != len)
        fail("Failed to read all the data from: " + path);
  }
}

template<typename T>
void Grid<T>::write_ccp4_map(const std::string& path, int mode) const {
  if (mode == 2) {
    impl::write_arrays(path, impl::make_ccp4_header(*this, mode), data);
  } else if (mode == 0) {
    std::vector<signed char> v(data.size());
    // TODO: scale data
    impl::write_arrays(path, impl::make_ccp4_header(*this, mode), v);
  } else {
    fail("Only modes 0 and 2 are supported.");
  }
}

template<typename T>
void Grid<T>::write_ccp4_mask(const std::string& path, double threshold) const {
  std::vector<signed char> v(data.size());
  for (size_t i = 0; i < data.size(); i++)
    v[i] = data[i] < threshold ? 0 : 1;
  impl::write_arrays(path, impl::make_ccp4_header(*this, 0), v);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
