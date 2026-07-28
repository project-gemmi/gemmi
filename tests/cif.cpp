
#include "doctest.h"

#include <algorithm>
#include <gemmi/read_cif.hpp>
#include <gemmi/refln.hpp>     // for ReflnBlock
#include <gemmi/intensit.hpp>  // for Intensities

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

TEST_CASE("cif::Column iterators") {
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
    CHECK_EQ(tab.at(1).at(0), "b2");
    CHECK_EQ(tab.column(0).at(1), "b2");
    CHECK_EQ(tab.find_column("_d").at(0), "d1");
    cif::Table::Row row = tab.at(1);
    CHECK_EQ(row.at(0), "b2");
    CHECK_EQ(row.at(1), "d2");
    CHECK_EQ(std::count(row.begin(), row.end(), "d2"), 1);
    check_with_two_elements(row);
  }
  {
    cif::Column col_a = block.find({"_a"}).column(0);
    CHECK_EQ(col_a.at(0), "a1");
    CHECK_EQ(col_a.at(1), "a2");
    check_with_two_elements(col_a);
  }
}

TEST_CASE("cif::Block::init_loop") {
  cif::Document doc = cif::read_string("data_1 _m.a 1 _m.b 2 _m.c 3 "
                                       "loop_ _p.u _p.v _p.w 5 6 7");
  cif::Block& block = doc.blocks[0];
  block.init_loop("_x.", {"one", "two"}).add_row({"1", "2"});
  CHECK_EQ(block.find_loop("_x.one").at(0), "1");
  block.init_loop("_m.", {"b"}).add_row({"10"});
  block.find({"_m.b"}).erase();
  CHECK_EQ(block.find_values("_m.b").item(), nullptr);
  CHECK_EQ(*block.find_value("_m.a"), "1");
  block.init_mmcif_loop("_m.", {"c"}).add_row({"20"});
  CHECK_EQ(block.find_value("_m.a"), nullptr);
  block.init_loop("_p.", {"x"}).add_row({"10"});
  CHECK_EQ(block.find_values("_p.x").at(0), "10");
  CHECK_EQ(block.find_values("_p.u").at(0), "5");
  block.set_pair("_p.v", "30");
  CHECK_EQ(block.find_values("_p.u").item(), nullptr);
  CHECK_EQ(block.find_values("_p.v").at(0), "30");
}

TEST_CASE("Intensities::import_f_squared_from_mmcif") {
  cif::Document doc = cif::read_string(
      "data_r\n"
      "_cell.length_a 10 _cell.length_b 10 _cell.length_c 10\n"
      "_cell.angle_alpha 90 _cell.angle_beta 90 _cell.angle_gamma 90\n"
      "_symmetry.space_group_name_H-M 'P 1'\n"
      "loop_ _refln.index_h _refln.index_k _refln.index_l\n"
      "_refln.F_meas_au _refln.F_meas_sigma_au\n"
      "1 0 0 10 0.5\n"
      "0 1 0 4 0.25\n");
  gemmi::ReflnBlock rblock(std::move(doc.blocks.at(0)));
  gemmi::Intensities intensities;
  intensities.import_f_squared_from_mmcif(rblock);
  REQUIRE(intensities.data.size() == 2);
  // I = F^2 and sigma(I) = 2 F sigma(F)
  CHECK(intensities.data[0].value == doctest::Approx(100.));
  CHECK(intensities.data[0].sigma == doctest::Approx(10.));
  CHECK(intensities.data[1].value == doctest::Approx(16.));
  CHECK(intensities.data[1].sigma == doctest::Approx(2.));
}
