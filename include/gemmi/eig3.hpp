// Eigen decomposition code for symmetric 3x3 matrices.

#ifndef GEMMI_EIG3_HPP_
#define GEMMI_EIG3_HPP_

#include "math.hpp"  // for SMat33, Mat33
#include "fail.hpp"  // for GEMMI_DLL

namespace gemmi {

GEMMI_DLL Mat33 eigen_decomposition(const SMat33<double>& A, double (&d)[3]);

} // namespace gemmi

#endif
