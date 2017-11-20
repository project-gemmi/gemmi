
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <algorithm>
#include <gemmi/cif.hpp>
namespace cif = gemmi::cif;

template<typename T> void check_with_two_elements(T duo) {
  typename T::const_iterator it;
  it = duo.begin();
  CHECK_EQ(it++, duo.begin());
  CHECK_EQ(++it, duo.end());
  CHECK_EQ(it--, duo.end());
  CHECK_EQ(--it, duo.begin());
  typename T::const_iterator cit;
  const T& cduo = duo;
  cit = cduo.end();
  cit = it;
  CHECK_EQ(cit == it, true);
  CHECK_EQ(cit != it, false);
  std::string first = duo[0];
  std::string second = duo[1];
  std::string rev[2];
  std::reverse_copy(duo.begin(), duo.end(), rev);
  CHECK_EQ(rev[0], second);
  CHECK_EQ(rev[1], first);
  std::reverse(duo.begin(), duo.end());
  CHECK_EQ(duo[0], second);
  CHECK_EQ(duo[1], first);
  std::fill(duo.begin(), duo.end(), "filler");
  CHECK_EQ(std::count(duo.begin(), duo.end(), "filler"), 2);
}

TEST_CASE("testing cif::Column iterators") {
  cif::Document doc = cif::read_string("data_1"
          " loop_ _a _b _c _d a1 b1 c1 d1 a2 b2 c2 d2"
          " _pair 1");
  cif::Block& block = doc.blocks[0];
  {
    cif::Column col_c = block.find_loop("_c");
    CHECK_EQ(col_c.length(), 2);
    CHECK_EQ(col_c.at(0), "c1");
    CHECK_EQ(col_c.at(1), "c2");
    CHECK_EQ(std::count(col_c.begin(), col_c.end(), "c2"), 1);
    check_with_two_elements(col_c);
  }
  {
    cif::Table tab = block.find({"_b", "_d"});
    CHECK_EQ(tab.length(), 2);
    CHECK_EQ(tab.width(), 2);
    CHECK_EQ(tab.at(0).at(0), "b1");
    CHECK_EQ(tab.at(1).at(1), "d2");
    cif::Table::Row row = tab.at(1);
    CHECK_EQ(row.at(0), "b2");
    CHECK_EQ(row.at(1), "d2");
    CHECK_EQ(std::count(row.begin(), row.end(), "d2"), 1);
    check_with_two_elements(row);
  }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
