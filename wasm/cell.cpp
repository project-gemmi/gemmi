// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/unitcell.hpp>

void add_cell() {
  em::class_<gemmi::UnitCellParameters>("UnitCellParameters")
    .property("a", &gemmi::UnitCellParameters::a)
    .property("b", &gemmi::UnitCellParameters::b)
    .property("c", &gemmi::UnitCellParameters::c)
    .property("alpha", &gemmi::UnitCellParameters::alpha)
    .property("beta", &gemmi::UnitCellParameters::beta)
    .property("gamma", &gemmi::UnitCellParameters::gamma)
    ;
  em::class_<gemmi::UnitCell, em::base<gemmi::UnitCellParameters>>("UnitCell")
    .property("volume", &gemmi::UnitCell::volume)
    .function("is_crystal", &gemmi::UnitCell::is_crystal)
    ;
}
