
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <cstdlib>  // for rand
#include <gemmi/unitcell.hpp>
#include <linalg.h>

TEST_CASE("Matrix44::determinant") {
  std::srand(12345);
  gemmi::Matrix44 g_m44;
  linalg::mat<double,4,4> la_m44;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      g_m44.a[i][j] = la_m44[i][j] = 10.0 * std::rand() / RAND_MAX - 5;
  CHECK_EQ(g_m44.determinant(), doctest::Approx(linalg::determinant(la_m44)));
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
