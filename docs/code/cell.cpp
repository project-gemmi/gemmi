#include <gemmi/unitcell.hpp>

// UnitCell has:
// * directly set properties a, b, c, alpha, beta, gamma,
// * calculated properties such as ar (a*), br (b*), ..., volume,
// * fractionalization and orthogonalization matrices,
// * a list of images (symmetry or NCS mates) that is set externally
//   when reading a file.
// * and a few functions such as orthogonalize(), fractionalize(),
//   is_special_position(), find_nearest_image().

gemmi::UnitCell cell(25.14, 39.50, 45.07, 90, 90, 90);
gemmi::Position p = cell.orthogonalize(gemmi::Fractional(0.5, 0.5, 0.5));

