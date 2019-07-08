
#include <cstring>
#include "gemmi/ccp4.hpp"

extern "C" {
#include "grid.h"
}

using Mask = gemmi::Ccp4<std::int8_t>;
static geMask* as_c(Mask* sg) { return reinterpret_cast<geMask*>(sg); }
static Mask* as_cpp(geMask* sg) { return reinterpret_cast<Mask*>(sg); }

extern "C" {

geMask* geMask_init(int nx, int ny, int nz) {
  Mask* mask = new Mask;
  mask->grid.set_size(nx, ny, nz);
  return as_c(mask);
}

void geMask_set_unit_cell(geMask* mask, double a, double b, double c,
                          double alpha, double beta, double gamma) {
  as_cpp(mask)->grid.set_unit_cell(a, b, c, alpha, beta, gamma);
}

void geMask_mask_atom(geMask* mask, double x, double y, double z,
                      double radius) {
  as_cpp(mask)->grid.mask_atom(x, y, z, radius);
}

void geMask_apply_space_group(geMask* mask, int ccp4_num) {
  gemmi::Grid<std::int8_t>& g = as_cpp(mask)->grid;
  g.spacegroup = gemmi::find_spacegroup_by_number(ccp4_num);
  g.symmetrize([](std::int8_t a, std::int8_t b) { return a > b ? a : b; });
}

int8_t* geMask_data(geMask* mask) {
  return as_cpp(mask)->grid.data.data();
}

int8_t geMask_get_value(geMask* mask, int u, int v, int w) {
  return as_cpp(mask)->grid.get_value(u, v, w);
}

void geMask_free(geMask* mask) {
  delete as_cpp(mask);
}

void geMask_update_ccp4_header(geMask* mask, int mode, bool update_stats) {
  as_cpp(mask)->update_ccp4_header(mode, update_stats);
}

void geMask_write_ccp4_map(geMask* mask, const char* path) {
  as_cpp(mask)->write_ccp4_map(path);
}

}
