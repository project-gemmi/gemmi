// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/unitcell.hpp>

void add_cell() {
  em::value_array<gemmi::Fractional>("Fractional")
    .element(&gemmi::Fractional::x)
    .element(&gemmi::Fractional::y)
    .element(&gemmi::Fractional::z)
    ;
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
