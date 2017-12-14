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
#include <functional> // for function
#include <string>
#include <typeinfo>  // for typeid
#include <vector>
#include "unitcell.hpp"
#include "symmetry.hpp"
#include "util.hpp"  // for fail, file_open

namespace gemmi {

using std::int32_t;

// options for Grid<>::setup
enum class GridSetup {
  ReorderOnly,  // reorder axes to X, Y, Z
  Full,         // reorder and expand to the whole unit cell
  FullCheck     // additionally consistency of redundant data
};

inline bool is_little_endian() {
  std::uint32_t x = 1;
  return *reinterpret_cast<char *>(&x) == 1;
}

int modulo(int a, int n) {
  if (a >= n)
    a %= n;
  else if (a < 0)
    a = a % n + n;
  return a;
}

struct GridStats {
  double dmin = NAN;
  double dmax = NAN;
  double dmean = NAN;
  double rms = NAN;
};

template<typename T>
GridStats calculate_grid_statistics(const std::vector<T>& data) {
  GridStats st;
  if (data.empty())
    return st;
  double sum = 0;
  double sq_sum = 0;
  st.dmin = st.dmax = data[0];
  for (double d : data) {
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


inline bool has_small_factorization(int n) {
  for (int k : {2, 3, 5})
    while (n % k == 0)
      n /= k;
  return n == 1 || n == -1;
}

struct GridMeta {
  int nu, nv, nw;
  UnitCell unit_cell;
  bool full_canonical; // grid for the whole unit cell with X,Y,Z order
  const SpaceGroup* space_group;
  double spacing[3];
  GridStats hstats;  // data statistics read from / written to ccp4 map
  // stores raw headers if the grid was read from ccp4 map
  std::vector<int32_t> ccp4_header;

  // methods to access info from ccp4 headers, w is word number from the spec
  int32_t* header_word(int w) { return &ccp4_header.at(w - 1); }
  const int32_t* header_word(int w) const { return &ccp4_header.at(w - 1); }

  int32_t header_i32(int w) const { return *header_word(w); }
  float header_float(int w) const {
    float f;
    std::memcpy(&f, header_word(w), 4);
    return f;
  }
  // ccp4 map header has mostly 80-byte strings
  std::string header_str(int w, size_t len=80) const {
    if (4 * ccp4_header.size() < 4 * (w - 1) + len)
      fail("invalid end of string");
    return std::string(reinterpret_cast<const char*>(header_word(w)), len);
  }
  void set_header_i32(int w, int32_t value) {
    *header_word(w) = value;
  };
  void set_header_3i32(int w, int32_t x, int32_t y, int32_t z) {
    set_header_i32(w, x);
    set_header_i32(w+1, y);
    set_header_i32(w+2, z);
  };
  void set_header_float(int w, float value) {
    std::memcpy(header_word(w), &value, 4);
  };
  void set_header_str(int w, const std::string& str) {
    // TODO: should the written string be null-terminated? probably not
    std::memcpy(header_word(w), str.c_str(), str.size() + 1);
  }

  void make_ccp4_header(int mode) {
    int nsym = 1;
    ccp4_header.clear();
    ccp4_header.resize(256 + nsym * 20, 0);
    set_header_3i32(1, nu, nv, nw); // NX, NY, NZ
    set_header_3i32(5, 0, 0, 0); // NXSTART, NYSTART, NZSTART
    set_header_3i32(8, nu, nv, nw);  // MX, MY, MZ
    set_header_float(11, (float) unit_cell.a);
    set_header_float(12, (float) unit_cell.b);
    set_header_float(13, (float) unit_cell.c);
    set_header_float(14, (float) unit_cell.alpha);
    set_header_float(15, (float) unit_cell.beta);
    set_header_float(16, (float) unit_cell.gamma);
    set_header_3i32(17, 1, 2, 3); // MAPC, MAPR, MAPS
    set_header_i32(23, 1); // ISPG
    set_header_i32(24, nsym * 80);
    std::memcpy(header_word(27), "CCP4", 4); // EXTTYP
    //set_header_i32(28, nversion);
    std::memcpy(header_word(53), "MAP ", 4);
    // MACHST
    set_header_i32(54, is_little_endian() ? 0x00004144 : 0x11110000);
    set_header_i32(56, 1); // labels
    std::memset(header_word(57), ' ', 800 + nsym * 80);
    set_header_str(57, "written by GEMMI");
    std::memcpy(header_word(257), "X,  Y,  Z", 9);
    update_ccp4_header(mode);
  }

  void update_ccp4_header(int mode) {
    if (mode != 0 && mode != 1 && mode != 2 && mode != 6)
      fail("Only modes 0, 1, 2 and 6 are supported.");
    if (ccp4_header.empty())
      return make_ccp4_header(mode);
    assert(ccp4_header.size() >= 256);
    // selectively copy-pasted from make_ccp4_header()
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
      header_i32(5) == 0 && header_i32(6) == 0 && header_i32(7)  == 0 &&
      // MX == NX
      header_i32(8) == nu && header_i32(9) == nv && header_i32(10) == nw &&
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
        gemmi::fail("Incorrect MAPC/MAPR/MAPS records");
      pos[mapi - 1] = i;
    }
    return pos;
  }

  void read_ccp4_header(FILE* f, const std::string& path) {
    const size_t hsize = 256;
    ccp4_header.resize(hsize);
    if (std::fread(ccp4_header.data(), 4, hsize, f) != hsize)
      fail("Failed to read map header: " + path);
    if (header_str(53, 4) != "MAP ")
      fail("Not a CCP4 map: " + path);
    unit_cell.set(header_float(11), header_float(12), header_float(13),
                  header_float(14), header_float(15), header_float(16));
    size_t ext_w = header_i32(24) / 4;  // nsymbt in words
    if (ext_w > 1000000)
      fail("Unexpectedly long extendended header: " + path);
    ccp4_header.resize(hsize + ext_w);
    if (std::fread(ccp4_header.data() + hsize, 4, ext_w, f) != ext_w)
      fail("Failed to read extended header: " + path);
    nu = header_i32(1);
    nv = header_i32(2);
    nw = header_i32(3);
    for (int i = 0; i < 3; ++i) {
      int axis = header_i32(17 + i);
      if (axis < 1 || axis > 3)
        fail("Unexpected axis value in word " + std::to_string(17 + i));
    }
    hstats.dmin = header_float(20);
    hstats.dmax = header_float(21);
    hstats.dmean = header_float(22);
    hstats.rms = header_float(55);
    space_group = find_spacegroup_by_number(header_i32(23));
    auto pos = axis_positions();
    full_canonical = pos[0] == 0 && pos[1] == 1 && pos[2] == 2 && full_cell();
  }
};

// For now, for simplicity, the grid covers whole unit cell
// and space group is P1.
template<typename T=float>
struct Grid : GridMeta {
  std::vector<T> data;

  void set_size_without_checking(int u, int v, int w) {
    nu = u, nv = v, nw = w;
    data.resize(u * v * w);
    spacing[0] = 1.0 / (nu * unit_cell.ar);
    spacing[1] = 1.0 / (nv * unit_cell.br);
    spacing[2] = 1.0 / (nw * unit_cell.cr);
    full_canonical = true;
  }

  void set_size(int u, int v, int w) {
    if (space_group) {
      auto factors = space_group->operations().find_grid_factors();
      if (u % factors[0] != 0 || v % factors[1] != 0 || w % factors[2] != 0)
        fail("Grid not compatible with the space group " + space_group->xhm());
    }
    set_size_without_checking(u, v, w);
  }

  void set_size_from_max_spacing(double max_spacing) {
    const SpaceGroup& sg = space_group ? *space_group : get_spacegroup_p1();
    std::array<int, 3> sg_fac = sg.operations().find_grid_factors();
    int m[3];
    for (int i = 0; i != 3; ++i) {
      int f = std::max(2, sg_fac[i]);
      int n = int(std::ceil(unit_cell[i] / (max_spacing * f)));
      while (!has_small_factorization(n))
        ++n;
      m[i] = n * f;
    }
    set_size_without_checking(m[0], m[1], m[2]);
  }

  int index(int u, int v, int w) const { return w * nu * nv + v * nu + u; }

  T get_value(int u, int v, int w) const { return data[index(u, v, w)]; }

  // quick-wrap assumes (for efficiency) that the index is not far from [0,nu).
  int quick_wrapped_index(int u, int v, int w) const {
    if (u >= nu) u -= nu; else if (u < 0) u += nu;
    if (v >= nv) v -= nv; else if (v < 0) v += nv;
    if (w >= nw) w -= nw; else if (w < 0) w += nw;
    return index(u, v, w);
  }

  int wrapped_index(int u, int v, int w) const {
    return index(modulo(u, nu), modulo(v, nv), modulo(w, nw));
  }

  int symmetric_index(int u, int v, int w, const Op& op) const {
    //TODO apply symmetry to u v w
    return quick_wrapped_index(u, v, w);
  }

  void set_points_around(const Position& ctr, double radius, T value) {
    int du = (int) std::ceil(radius / spacing[0]);
    int dv = (int) std::ceil(radius / spacing[1]);
    int dw = (int) std::ceil(radius / spacing[2]);
    if (du > nu || dv > nv || dw > nw)
      fail("Masking radius bigger than the unit cell?");
    Position fctr = unit_cell.fractionalize(ctr).wrap_to_unit();
    int u0 = iround(fctr.x * nu);
    int v0 = iround(fctr.y * nv);
    int w0 = iround(fctr.z * nw);
    for (int w = w0-dw; w <= w0+dw; ++w)
      for (int v = v0-dv; v <= v0+dv; ++v)
        for (int u = u0-du; u < u0+du; ++u) {
          Position fdelta{fctr.x - u * (1.0 / nu),
                          fctr.y - v * (1.0 / nv),
                          fctr.z - w * (1.0 / nw)};
          for (int i = 0; i < 3; ++i)
            if (fdelta[i] > 0.5)
              fdelta[i] -= 1.0;
            else if (fdelta[i] < -0.5)
              fdelta[i] += 1.0;
          Position d = unit_cell.orthogonalize(fdelta);
          if (d.x*d.x + d.y*d.y + d.z*d.z < radius*radius) {
            data[quick_wrapped_index(u, v, w)] = value;
          }
        }
  }

  void make_zeros_and_ones(double threshold) {
    for (auto& d : data)
      d = d > threshold ? 1 : 0;
  }

  void symmetrize(std::function<T(T, T)> func) {
    if (!space_group || space_group->number == 1 || !full_canonical)
      return;
    std::vector<Op> ops = space_group->operations().all_ops_sorted();
    std::vector<int> mates(ops.size(), 0);
    std::vector<bool> visited(data.size(), false);
    int idx = -1;
    for (int w = 0; w != nw; ++w)
      for (int v = 0; v != nv; ++v)
        for (int u = 0; u != nu; ++u) {
          ++idx;
          assert(idx == index(u, v, w));
          if (visited[idx])
            continue;
          mates[0] = idx;
          for (size_t k = 1; k < ops.size(); ++k)
            mates[k] = symmetric_index(u, v, w, ops[k]);
          T value = data[idx];
          for (int k : mates) {
            assert(!visited[k]);
            visited[k] = true;
            value = func(value, data[k]);
          }
          for (int k : mates)
            data[k] = value;
        }
    assert(idx == (int) data.size());
  }

  void setup(GridSetup mode);

  void read_ccp4_map(const std::string& path);
  void write_ccp4_map(const std::string& path) const;
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

} // namespace impl

// This function was tested only on little-endian machines,
// let us know if you need support for other architectures.
template<typename T>
void Grid<T>::read_ccp4_map(const std::string& path) {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "rb");
  read_ccp4_header(f.get(), path);
  data.resize(nu * nv * nw);
  int mode = header_i32(4);
  if (mode == 0)
    impl::read_data<std::int8_t>(f.get(), data);
  else if (mode == 1)
    impl::read_data<std::int16_t>(f.get(), data);
  else if (mode == 2)
    impl::read_data<float>(f.get(), data);
  else if (mode == 6)
    impl::read_data<std::uint16_t>(f.get(), data);
  else
    fail("Only modes 0, 1, 2 and 6 are supported.");
  //if (std::fgetc(f.get()) != EOF)
  //  fail("The map file is longer then expected.");
}

template<typename T>
void Grid<T>::setup(GridSetup mode) {
  if (full_canonical || ccp4_header.empty())
    return;
  // cell sampling does not change
  int sampl[3] = { header_i32(8), header_i32(9), header_i32(10) };
  // get old metadata
  auto pos = axis_positions();
  int start[3] = { header_i32(5), header_i32(6), header_i32(7) };
  int end[3] = { start[0] + nu, start[1] + nv, start[2] + nw };
  // set new metadata
  if (mode != GridSetup::ReorderOnly) {
    nu = sampl[0];
    nv = sampl[1];
    nw = sampl[2];
    set_header_3i32(5, 0, 0, 0); // start
  } else {
    set_header_3i32(5, start[pos[0]], start[pos[1]], start[pos[2]]);
    for (int i = 0; i < 3; ++i) {
      end[i] -= start[i];
      start[i] = 0;
    }
    int crs[3] = { nu, nv, nw };
    nu = crs[pos[0]];
    nv = crs[pos[1]];
    nw = crs[pos[2]];
  }
  set_header_3i32(1, nu, nv, nw); // NX, NY, NZ
  set_header_3i32(17, 1, 2, 3); // axes (MAPC, MAPR, MAPS)
  // now set the data
  std::vector<T> full(nu * nv * nw, NAN);
  int it[3];
  int idx = 0;
  for (it[2] = start[2]; it[2] < end[2]; it[2]++) // sections
    for (it[1] = start[1]; it[1] < end[1]; it[1]++) // rows
      for (it[0] = start[0]; it[0] < end[0]; it[0]++) { // cols
        T val = data[idx++];
        int new_index = wrapped_index(it[pos[0]], it[pos[1]], it[pos[2]]);
        if (mode == GridSetup::FullCheck &&
            !std::isnan(full[new_index]) && full[new_index] != val)
          fail("inconsistent data for the same point: " + std::to_string(val) +
               " != " + std::to_string(full[new_index]));
        full[new_index] = val;
      }
  data = full;
  if (mode == GridSetup::Full)
    symmetrize([](T a, T b) { return std::isnan(a) ? b : a; });
  else if (mode == GridSetup::FullCheck)
    symmetrize([](T a, T b) {
        if (!std::isnan(a) && !std::isnan(b) && a != b)
          fail("map values do not preserve symmetry: " + std::to_string(a) +
               " != " + std::to_string(b));
        return std::isnan(a) ? b : a;
    });
}

template<typename T>
void Grid<T>::write_ccp4_map(const std::string& path) const {
  gemmi::fileptr_t f = gemmi::file_open(path.c_str(), "wb");
  std::fwrite(ccp4_header.data(), 4, ccp4_header.size(), f.get());
  int mode = header_i32(4);
  if (mode == 0)
    impl::write_data<std::int8_t>(data, f.get());
  else if (mode == 1)
    impl::write_data<std::int16_t>(data, f.get());
  else if (mode == 2)
    impl::write_data<float>(data, f.get());
  else if (mode == 6)
    impl::write_data<std::uint16_t>(data, f.get());
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
