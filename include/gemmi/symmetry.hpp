// Copyright 2017 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdint>
#include <cstdlib>    // for strtol
#include <cstring>    // for memchr
#include <array>
#include <algorithm>  // for count
#include <stdexcept>  // for runtime_error
#include <string>

namespace gemmi {
namespace sym {

// TRIPLET <-> SYM OP

//using RotOp = std::array<std::array<std::int8_t, 3>, 3>;
//using TranOp = std::array<std::int8_t, 3>;

struct Op {
  typedef std::array<std::array<std::int8_t, 3>, 3> Rot;
  typedef std::array<std::int8_t, 3> Tran;
  Rot rot;
  Tran tran;

  std::string triplet() const;

  Op inv() const {
    // TODO
    return { rot, tran };
  }
};
inline bool operator==(const Op& a, const Op& b) {
  return a.rot == b.rot && a.tran == b.tran;
}
inline bool operator!=(const Op& a, const Op& b) { return !(a == b); }

inline Op combine(const Op& a, const Op& b) {
  Op r;
  for (int i = 0; i != 3; ++i) {
    r.tran[i] = a.tran[i];
    for (int j = 0; j != 3; ++j) {
      r.rot[i][j] = a.rot[i][0] * b.rot[0][j] +
                    a.rot[i][1] * b.rot[1][j] +
                    a.rot[i][2] * b.rot[2][j];
      r.tran[i] += a.rot[i][j] * b.tran[j];
    }
    r.tran[i] %= 12;
    if (r.tran[i] < 0)
      r.tran[i] += 12;
  }
  return r;
}

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

inline std::array<std::int8_t, 4> parse_triplet_part(const std::string& s) {
  std::array<std::int8_t, 4> r = { 0, 0, 0, 0 };
  int sign = 1;
  for (const char* c = s.c_str(); c != s.c_str() + s.length(); ++c) {
    if (*c == '+' || *c == '-') {
      sign = (*c == '+' ? 1 : -1);
      continue;
    }
    if (*c == ' ' || *c == '\t')
      continue;
    if (sign == 0)
      fail("wrong triplet format in: " + s);
    if (*c >= '0' && *c <= '9') {
      char* endptr;
      r[3] = sign * strtol(c, &endptr, 10);
      int den = 1;
      if (*endptr == '/') {
        den = strtol(endptr + 1, &endptr, 10);
        if (den != 1 && den != 2 && den != 3 && den != 4 && den != 6)
          fail("Unexpected denominator " + std::to_string(den) + " in: " + s);
      }
      r[3] *= 12 / den;
      c = endptr - 1;
    } else if (std::memchr("xXhHaA", *c, 6)) {
      r[0] += sign;
    } else if (std::memchr("yYkKbB", *c, 6)) {
      r[1] += sign;
    } else if (std::memchr("zZlLcC", *c, 6)) {
      r[2] += sign;
    } else {
      fail("unexpected character '" + std::string(1, *c) + "' in: " + s);
    }
    sign = 0;
  }
  if (sign != 0)
    fail("trailing sign in: " + s);
  return r;
}

inline Op parse_triplet(const std::string& s) {
  if (std::count(s.begin(), s.end(), ',') != 2)
    fail("expected exactly two commas in triplet");
  size_t comma1 = s.find(',');
  size_t comma2 = s.find(',', comma1 + 1);
  auto a = parse_triplet_part(s.substr(0, comma1));
  auto b = parse_triplet_part(s.substr(comma1 + 1, comma2 - (comma1 + 1)));
  auto c = parse_triplet_part(s.substr(comma2 + 1));
  Op::Rot rot = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
  Op::Tran tran = {a[3], b[3], c[3]};
  return { rot, tran };
}

inline std::string make_triplet_part(int x, int y, int z, int w) {
  std::string s;
  int xyz[] = { x, y, z };
  for (int i = 0; i != 3; ++i)
    if (xyz[i] != 0)
      s += (xyz[i] < 0 ? "-" : s.empty() ? "" : "+") + std::string(1, 'x' + i);
  if (w != 0) {  // simplify w/12
    int denom = 1;
    for (int factor : {2, 2, 3})
      if (w % factor == 0)
        w /= factor;
      else
        denom *= factor;
    s += (w > 0 ? "+" : "") + std::to_string(w);
    if (denom != 1)
      s += "/" + std::to_string(denom);
  }
  return s;
}

inline std::string Op::triplet() const {
  return make_triplet_part(rot[0][0], rot[0][1], rot[0][2], tran[0]) +
   "," + make_triplet_part(rot[1][0], rot[1][1], rot[1][2], tran[1]) +
   "," + make_triplet_part(rot[2][0], rot[2][1], rot[2][2], tran[2]);
}


// CRYSTALLOGRAPHIC SPACE GROUPS

struct SpaceGroup {
  std::uint8_t number;
  std::uint8_t ccp4;
  char xHM[20];
  char Hall[16];
};

// the template here is only to substitute C++17 inline variables
// https://stackoverflow.com/questions/38043442/how-do-inline-variables-work
template<class Dummy>
struct Data_
{
  static const SpaceGroup arr[530];
};

template<class Dummy>
const SpaceGroup Data_<Dummy>::arr[530] = {
  {1, 0, "P 1", "P 1"},
};
using data = Data_<void>;


} // namespace sym
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
