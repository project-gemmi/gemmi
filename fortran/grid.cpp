
#include <cstring>
#include "gemmi/grid.hpp"

extern "C" {
#include "grid.h"
}

using gemmi::Grid;

using Grid0 = Grid<std::int8_t>;
static geGrid0* as_c(Grid0* sg) { return reinterpret_cast<geGrid0*>(sg); }
static Grid0* as_cpp(geGrid0* sg) { return reinterpret_cast<Grid0*>(sg); }

extern "C" {

geGrid0* geGrid0_init(int nx, int ny, int nz) {
  Grid0* grid = new Grid0;
  grid->set_size(nx, ny, nz);
  return as_c(grid);
}

void geGrid0_set_unit_cell(geGrid0* grid, double a, double b, double c,
                           double alpha, double beta, double gamma) {
  as_cpp(grid)->set_unit_cell(a, b, c, alpha, beta, gamma);
}

void geGrid0_mask_atom(geGrid0* grid, double x, double y, double z,
                       double radius) {
  as_cpp(grid)->mask_atom(x, y, z, radius);
}

void geGrid0_apply_space_group(geGrid0* grid, int ccp4_num) {
  Grid0* g = as_cpp(grid);
  g->space_group = gemmi::find_spacegroup_by_number(ccp4_num);
  g->symmetrize([](std::int8_t a, std::int8_t b) { return a > b ? a : b; });
}

int8_t* geGrid0_data(geGrid0* grid) {
  return as_cpp(grid)->data.data();
}

void geGrid0_free(geGrid0* grid) {
  delete as_cpp(grid);
}

void geGrid0_prepare_ccp4_header(geGrid0* grid, int n) {
  as_cpp(grid)->prepare_ccp4_header(n);
}

void geGrid0_write_ccp4_map(geGrid0* grid, const char* path) {
  as_cpp(grid)->write_ccp4_map(path);
}

}
// vim:sw=2:ts=2:et
