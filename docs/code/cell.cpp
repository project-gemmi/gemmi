#include <gemmi/unitcell.hpp>

gemmi::UnitCell cell(25.14, 39.50, 45.07, 90, 90, 90);
auto param = {cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma};
double volume = cell.volume;
gemmi::Position p = cell.orthogonalize(gemmi::Fractional(0.5, 0.5, 0.5));
gemmi::Fractional f = cell.fractionalize(p);

