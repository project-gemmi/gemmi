// Copyright 2022 Global Phasing Ltd.
//
// MaskedGrid used primarily as direct-space asu mask.

#ifndef GEMMI_ASUMASK_HPP_
#define GEMMI_ASUMASK_HPP_

#include "grid.hpp"

namespace gemmi {

template<typename T, typename V=std::int8_t> struct MaskedGrid {
  std::vector<V> mask;
  Grid<T>* grid;

  struct iterator {
    MaskedGrid& parent;
    size_t index;
    int u = 0, v = 0, w = 0;
    iterator(MaskedGrid& parent_, size_t index_)
      : parent(parent_), index(index_) {}
    iterator& operator++() {
      do {
        ++index;
        if (++u == parent.grid->nu) {
          u = 0;
          if (++v == parent.grid->nv) {
            v = 0;
            ++w;
          }
        }
      } while (index != parent.mask.size() && parent.mask[index] != 0);
      return *this;
    }
    typename GridBase<T>::Point operator*() {
      return {u, v, w, &parent.grid->data[index]};
    }
    bool operator==(const iterator &o) const { return index == o.index; }
    bool operator!=(const iterator &o) const { return index != o.index; }
  };
  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, mask.size()}; }

  std::array<int, 3> asu_max() const {
    std::array<int, 3> m = {{0, 0, 0}};
    for (auto pt : *const_cast<MaskedGrid*>(this))
      if (*pt.value == 0) {
        m[0] = std::max(pt.u, m[0]);
        m[1] = std::max(pt.v, m[1]);
        m[2] = std::max(pt.w, m[2]);
      }
    return m;
  }
};

template<typename T, typename V=std::int8_t>
std::vector<V> get_asu_mask(const Grid<T>& grid) {
  std::vector<V> mask(grid.data.size(), 0);
  std::vector<GridOp> ops = grid.get_scaled_ops_except_id();
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx)
        if (mask[idx] == 0)
          for (const GridOp& op : ops) {
            std::array<int, 3> t = op.apply(u, v, w);
            size_t mate_idx = grid.index_n(t[0], t[1], t[2]);
            // grid point can be on special position
            if (mate_idx != idx)
              mask[mate_idx] = 1;
          }
  return mask;
}

template<typename T>
MaskedGrid<T> masked_asu(Grid<T>& grid) {
  return {get_asu_mask(grid), &grid};
}

} // namespace gemmi
#endif
