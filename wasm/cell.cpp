// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/unitcell.hpp>

void add_cell() {
  em::class_<gemmi::UnitCell>("UnitCell")
    .property("a", &gemmi::UnitCell::a)
    .property("b", &gemmi::UnitCell::b)
    .property("c", &gemmi::UnitCell::c)
    .property("alpha", &gemmi::UnitCell::alpha)
    .property("beta", &gemmi::UnitCell::beta)
    .property("gamma", &gemmi::UnitCell::gamma)
    ;
}
