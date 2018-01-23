
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <cstdlib>  // for rand
#include <gemmi/unitcell.hpp>
#include <linalg.h>

TEST_CASE("Matrix44::determinant") {
  std::srand(12345);
  gemmi::Transform tr;
  linalg::mat<double,4,4> m44 = linalg::identity;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      tr.mat[i][j] = m44[i][j] = 10.0 * std::rand() / RAND_MAX - 5;
    tr.vec.at(i) = m44[i][3] = 10.0 * std::rand() / RAND_MAX - 5;
  }
  linalg::mat<double,4,4> inv_m44 = linalg::inverse(m44);
  gemmi::Transform inv_tr = tr.inverse();
  CHECK_EQ(inv_m44[3][3], doctest::Approx(1.0));
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      CHECK_EQ(inv_tr.mat[i][j], doctest::Approx(inv_m44[i][j]));
    CHECK_EQ(inv_tr.vec.at(i), doctest::Approx(inv_m44[i][3]));
    CHECK_EQ(inv_m44[3][i], 0);
  }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
