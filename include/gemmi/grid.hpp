// Copyright 2017 Global Phasing Ltd.
//
// 3d grid for volume density. CCP4 format for maps and masks.

// TODO: handle big-endian files

#ifndef GEMMI_GRID_HH
#define GEMMI_GRID_HH

#include <cmath>     // for NAN, sqrt
#include <cstdio>    // for FILE
#include <memory>    // for unique_ptr
#include <vector>
#include "unitcell.hpp"

namespace gemmi {
//namespace mol {

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
struct Grid {
  int nu, nv, nw;
  mol::UnitCell unit_cell;
  std::vector<float> data;
  // statistics
  double dmin = NAN, dmax = NAN, dmean = NAN, rms = NAN;
  // ccp4 map metadata
  std::vector<std::string> labels;

  Grid() noexcept : nu(0), nv(0), nw(0) {}
  Grid(int u, int v, int w) { set_size(u, v, w); }

  void set_size(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
  }

  void calculate_statistics() {
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

  void read_ccp4(const std::string& path);

  void write_ccp4_map(const std::string& path, int mode=2) const;

  void write_ccp4_mask(const std::string& path, double threshold) const {
    std::vector<signed char> v(data.size());
    for (size_t i = 0; i < data.size(); i++)
      v[i] = data[i] < threshold ? 0 : 1;
    write_arrays(path, make_ccp4_header(0), v);
  }

private:
  typedef std::unique_ptr<FILE, decltype(&std::fclose)> fileptr_t;
  std::vector<char> make_ccp4_header(int mode) const;

  template<typename T>
  void write_arrays(const std::string& path, const std::vector<char>& header,
                    const std::vector<T>& content) const {
    fileptr_t f(std::fopen(path.c_str(), "wb"), &std::fclose);
    if (!f)
      throw std::runtime_error("Failed to open file for writing: " + path);
    std::fwrite(header.data(), sizeof(char), header.size(), f.get());
    std::fwrite(content.data(), sizeof(T), content.size(), f.get());
  }
};

inline void Grid::read_ccp4(const std::string& path) {
  fileptr_t f(std::fopen(path.c_str(), "rb"), &std::fclose);
  if (!f)
    throw std::runtime_error("Failed to open file: " + path);
}

inline void Grid::write_ccp4_map(const std::string& path, int mode) const {
  if (mode == 2) {
    write_arrays(path, make_ccp4_header(mode), data);
  } else if (mode == 0) {
    std::vector<signed char> v(data.size());
    // TODO: scale data
    write_arrays(path, make_ccp4_header(mode), v);
  } else {
    throw std::runtime_error("Only modes 0 and 2 are supported.");
  }
}

inline std::vector<char> Grid::make_ccp4_header(int mode) const {
  std::vector<char> h(1024, 0);
  return h;
}

//} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
