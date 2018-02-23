
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <cstdlib>  // for rand
#include <gemmi/unitcell.hpp>
#include <linalg.h>

static double draw() { return 10.0 * std::rand() / RAND_MAX - 5; }
static gemmi::Transform random_transform() {
  gemmi::Transform a;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      a.mat[i][j] = draw();
    a.vec.at(i) = draw();
  }
  return a;
}

TEST_CASE("Transform::inverse") {
  std::srand(12345);
  gemmi::Transform tr = random_transform();
  linalg::mat<double,4,4> m44 = linalg::identity;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      m44[i][j] = tr.mat[i][j];
    m44[i][3] = tr.vec.at(i);
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

TEST_CASE("Transform::combine") {
  std::srand(12345);
  gemmi::Transform a = random_transform();
  gemmi::Transform b = random_transform();
  gemmi::Vec3 v;
  for (int i = 0; i < 3; ++i)
    v.at(i) = draw();
  gemmi::Vec3 result1 = a.combine(b).apply(v);
  gemmi::Vec3 result2 = a.apply(b.apply(v));
  for (int i = 0; i < 3; ++i)
    CHECK_EQ(result1.at(i), doctest::Approx(result2.at(i)));
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
