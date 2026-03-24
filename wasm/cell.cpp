// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/unitcell.hpp>

int nearest_image_shift_x(const gemmi::NearestImage& image) { return image.pbc_shift[0]; }
int nearest_image_shift_y(const gemmi::NearestImage& image) { return image.pbc_shift[1]; }
int nearest_image_shift_z(const gemmi::NearestImage& image) { return image.pbc_shift[2]; }

void add_cell() {
  em::value_array<gemmi::Fractional>("Fractional")
    .element(&gemmi::Fractional::x)
    .element(&gemmi::Fractional::y)
    .element(&gemmi::Fractional::z)
    ;
  em::class_<gemmi::NearestImage>("NearestImage")
    .constructor<>()
    .function("dist", &gemmi::NearestImage::dist)
    .function("same_asu", &gemmi::NearestImage::same_asu)
    .function("symmetry_code", &gemmi::NearestImage::symmetry_code)
    .property("sym_idx", &gemmi::NearestImage::sym_idx)
    .property("pbc_shift_x", &nearest_image_shift_x)
    .property("pbc_shift_y", &nearest_image_shift_y)
    .property("pbc_shift_z", &nearest_image_shift_z)
    ;
  em::register_vector<gemmi::NearestImage>("NearestImageVector");
  em::class_<gemmi::UnitCellParameters>("UnitCellParameters")
    .property("a", &gemmi::UnitCellParameters::a)
    .property("b", &gemmi::UnitCellParameters::b)
    .property("c", &gemmi::UnitCellParameters::c)
    .property("alpha", &gemmi::UnitCellParameters::alpha)
    .property("beta", &gemmi::UnitCellParameters::beta)
    .property("gamma", &gemmi::UnitCellParameters::gamma)
    ;
  em::class_<gemmi::UnitCell, em::base<gemmi::UnitCellParameters>>("UnitCell")
    .constructor<double, double, double, double, double, double>()
    .property("volume", &gemmi::UnitCell::volume)
    .function("is_crystal", &gemmi::UnitCell::is_crystal)
    .function("fractionalize", &gemmi::UnitCell::fractionalize)
    .function("orthogonalize", &gemmi::UnitCell::orthogonalize)
    ;
}
