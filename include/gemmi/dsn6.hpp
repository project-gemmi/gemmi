// Copyright Global Phasing Ltd.
//
// DSN6/BRIX density map format.

#ifndef GEMMI_DSN6_HPP_
#define GEMMI_DSN6_HPP_

#include <cmath>      // for ceil
#include <cstdint>    // for int16_t, uint16_t
#include <stdexcept>
#include "grid.hpp"
#include "stats.hpp"

namespace gemmi {

namespace impl {

inline int16_t read_dsn6_i16(const char* buf, size_t buf_size,
                             size_t index, bool little_endian) {
  const size_t offset = 2 * index;
  if (offset + 1 >= buf_size)
    throw std::runtime_error("DSN6 header is truncated");
  const unsigned char b0 = static_cast<unsigned char>(buf[offset]);
  const unsigned char b1 = static_cast<unsigned char>(buf[offset + 1]);
  const uint16_t value = little_endian ? static_cast<uint16_t>(b0 | (b1 << 8))
                                       : static_cast<uint16_t>((b0 << 8) | b1);
  return static_cast<int16_t>(value);
}

}  // namespace impl

/// Reads a DSN6/BRIX map from memory into a Grid<float> and returns statistics.
inline DataStats read_dsn6_from_memory(const char* buf, size_t size,
                                       Grid<float>& grid) {
  if (size < 512)
    throw std::runtime_error("DSN6: header is truncated");

  bool little_endian = false;
  if (impl::read_dsn6_i16(buf, size, 18, false) == 100) {
    little_endian = false;
  } else if (impl::read_dsn6_i16(buf, size, 18, true) == 100) {
    little_endian = true;
  } else {
    throw std::runtime_error("DSN6: endian detection failed");
  }

  auto header = [buf, size, little_endian](size_t index) {
    return impl::read_dsn6_i16(buf, size, index, little_endian);
  };

  const std::array<int, 3> origin = {{header(0), header(1), header(2)}};
  const std::array<int, 3> n_real = {{header(3), header(4), header(5)}};
  const std::array<int, 3> n_grid = {{header(6), header(7), header(8)}};
  const int cell_scale = header(17);
  const int prod_word = header(15);
  if (cell_scale == 0)
    throw std::runtime_error("DSN6: invalid cell scale in header");
  if (prod_word == 0)
    throw std::runtime_error("DSN6: invalid density scale in header");

  const double cell_mult = 1.0 / cell_scale;
  grid.set_unit_cell(cell_mult * header(9),
                     cell_mult * header(10),
                     cell_mult * header(11),
                     cell_mult * header(12),
                     cell_mult * header(13),
                     cell_mult * header(14));
  grid.set_size(n_grid[0], n_grid[1], n_grid[2]);

  const float prod = static_cast<float>(prod_word) / 100.f;
  const int plus = header(16);
  size_t offset = 512;
  const std::array<int, 3> n_blocks = {{
    static_cast<int>(std::ceil(n_real[0] / 8.0)),
    static_cast<int>(std::ceil(n_real[1] / 8.0)),
    static_cast<int>(std::ceil(n_real[2] / 8.0)),
  }};

  for (int zz = 0; zz < n_blocks[2]; ++zz) {
    for (int yy = 0; yy < n_blocks[1]; ++yy) {
      for (int xx = 0; xx < n_blocks[0]; ++xx) {
        for (int k = 0; k < 8; ++k) {
          const int z = 8 * zz + k;
          for (int j = 0; j < 8; ++j) {
            const int y = 8 * yy + j;
            for (int i = 0; i < 8; ++i) {
              const int x = 8 * xx + i;
              if (offset >= size)
                throw std::runtime_error("DSN6: data is truncated");
              if (x < n_real[0] && y < n_real[1] && z < n_real[2]) {
                const unsigned char density_byte =
                    static_cast<unsigned char>(buf[offset]);
                const float density = (density_byte - plus) / prod;
                ++offset;
                grid.set_value(origin[0] + x,
                               origin[1] + y,
                               origin[2] + z,
                               density);
              } else {
                offset += static_cast<size_t>(8 - i);
                if (offset > size)
                  throw std::runtime_error("DSN6: data is truncated");
                break;
              }
            }
          }
        }
      }
    }
  }

  return calculate_data_statistics(grid.data);
}

}  // namespace gemmi
#endif
