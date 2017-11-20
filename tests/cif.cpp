
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <algorithm>
#include <gemmi/cif.hpp>
namespace cif = gemmi::cif;

TEST_CASE("testing cif::Column iterators") {
  cif::Document doc = cif::read_string("data_1"
          " loop_ _a _b _c _d a1 b1 c1 d1 a2 b2 c2 d2"
          " _pair 1");
  cif::Block& block = doc.blocks[0];
  {
    cif::Column col_c = block.find_loop("_c");
    cif::Column::const_iterator it;
    it = col_c.begin();
    CHECK_EQ(it++, col_c.begin());
    CHECK_EQ(++it, col_c.end());
    CHECK_EQ(it--, col_c.end());
    CHECK_EQ(--it, col_c.begin());
    CHECK_EQ(col_c.length(), 2);
    CHECK_EQ(col_c.at(0), "c1");
    CHECK_EQ(col_c.at(1), "c2");
    CHECK_EQ(std::count(col_c.begin(), col_c.end(), "c2"), 1);
    CHECK_EQ(std::count(col_c.begin(), col_c.end(), "a2"), 0);
    std::string rev[2];
    std::reverse_copy(col_c.begin(), col_c.end(), rev);
    CHECK_EQ(rev[0], "c2");
    CHECK_EQ(rev[1], "c1");
    //std::fill(col_c.begin(), col_c.end(), "filler");
    //std::reverse(col_c.begin(), col_c.end());
  }
  {
    cif::Table tab = block.find({"_b", "_d"});
    CHECK_EQ(tab.length(), 2);
    CHECK_EQ(tab.width(), 2);
    CHECK_EQ(tab.at(0).at(0), "b1");
    CHECK_EQ(tab.at(1).at(1), "d2");
    cif::Table::Row row = tab.at(1);
    cif::Table::Row::const_iterator it;
    it = row.begin();
    CHECK_EQ(it++, row.begin());
    CHECK_EQ(++it, row.end());
    CHECK_EQ(it--, row.end());
    CHECK_EQ(--it, row.begin());
    CHECK_EQ(row.at(0), "b2");
    CHECK_EQ(row.at(1), "d2");
    CHECK_EQ(std::count(row.begin(), row.end(), "d2"), 1);
    CHECK_EQ(std::count(row.begin(), row.end(), "c2"), 0);
    std::string rev[2];
    std::reverse_copy(row.begin(), row.end(), rev);
    CHECK_EQ(rev[0], "d2");
    CHECK_EQ(rev[1], "b2");
    //std::fill(row.begin(), row.end(), "filler");
    //std::reverse(row.begin(), row.end());
  }
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
