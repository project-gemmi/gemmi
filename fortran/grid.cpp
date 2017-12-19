
#include <cstring>
#include "gemmi/grid.hpp"

extern "C" {
#include "grid.h"
}

using gemmi::Grid;

using GridC = Grid<signed char>;
static geGridC* as_c(GridC* sg) { return reinterpret_cast<geGridC*>(sg); }
static GridC* as_cpp(geGridC* sg) { return reinterpret_cast<GridC*>(sg); }

extern "C" {

geGridC* GridC_init(int nx, int ny, int nz) {
  GridC* grid = new GridC;
  grid->set_size(nx, ny, nz);
  return as_c(grid);
}

void GridC_set_unit_cell(geGridC* grid, double a, double b, double c,
                         double alpha, double beta, double gamma) {
  as_cpp(grid)->set_unit_cell(a, b, c, alpha, beta, gamma);
}

void GridC_mask_atom(geGridC* grid, double x, double y, double z,
                     double radius) {
  as_cpp(grid)->mask_atom(x, y, z, radius);
}

void GridC_apply_space_group(geGridC* grid, int ccp4_num) {
  GridC* g = as_cpp(grid);
  g->space_group = gemmi::find_spacegroup_by_number(ccp4_num);
  g->symmetrize([](signed char a, signed char b) { return a > b ? a : b; });
}

signed char* GridC_data(geGridC* grid) {
  return as_cpp(grid)->data.data();
}

void GridC_free(geGridC* grid) {
  delete as_cpp(grid);
}

}
// vim:sw=2:ts=2:et
