// Copyright 2017 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdint>
#include <cstdlib>    // for strtol
#include <cstring>    // for memchr, strlen
#include <array>
#include <algorithm>  // for count
#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>

#include <iostream> // debug

namespace gemmi {
namespace sym {

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

// TRIPLET <-> SYM OP

struct Op {
  typedef std::array<std::array<std::int8_t, 3>, 3> Rot;
  typedef std::array<std::int8_t, 3> Tran;
  Rot rot;
  Tran tran;

  std::string triplet() const;
  std::string rot_triplet() const { return Op{rot, {0, 0, 0}}.triplet(); };

  Op invert() const;

  void normalize_tran() { // wrap the elements into [0,12)
    for (int i = 0; i != 3; ++i) {
      tran[i] %= 12;
      if (tran[i] < 0)
        tran[i] += 12;
    }
  }

  void translate(const Tran& a) {
    for (int i = 0; i != 3; ++i)
      tran[i] += a[i];
  }

  void shift_origin(const Tran& a) {
    // TODO
  }

  int det_rot() const { // should be 1 (rotation) or -1 (with inversion)
    return rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
         - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
         + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0]);
  }

  static Op identity() { return {{1,0,0, 0,1,0, 0,0,1}, {0,0,0}}; };
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
  }
  r.normalize_tran();
  return r;
}

inline Op Op::invert() const {
  int detr = det_rot();
  if (detr != 1 && detr != -1)
    fail("not a rotation/inversion: det|" + rot_triplet() + "|="
         + std::to_string(detr));
  Op inv;
  inv.rot[0][0] = detr * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2]);
  inv.rot[0][1] = detr * (rot[0][2] * rot[2][1] - rot[0][1] * rot[2][2]);
  inv.rot[0][2] = detr * (rot[0][1] * rot[1][2] - rot[0][2] * rot[1][1]);
  inv.rot[1][0] = detr * (rot[1][2] * rot[2][0] - rot[1][0] * rot[2][2]);
  inv.rot[1][1] = detr * (rot[0][0] * rot[2][2] - rot[0][2] * rot[2][0]);
  inv.rot[1][2] = detr * (rot[1][0] * rot[0][2] - rot[0][0] * rot[1][2]);
  inv.rot[2][0] = detr * (rot[1][0] * rot[2][1] - rot[2][0] * rot[1][1]);
  inv.rot[2][1] = detr * (rot[2][0] * rot[0][1] - rot[0][0] * rot[2][1]);
  inv.rot[2][2] = detr * (rot[0][0] * rot[1][1] - rot[1][0] * rot[0][1]);
  for (int i = 0; i != 3; ++i)
    inv.tran[i] = - tran[0] * inv.rot[i][0]
                  - tran[1] * inv.rot[i][1]
                  - tran[2] * inv.rot[i][2];
  return inv;
}

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


// INTERPRETING HALL SYMBOLS
// based on http://cci.lbl.gov/sginfo/hall_symbols.html

struct SymOps {
  bool centrosym = false;
  std::vector<Op> sym_ops;
  std::vector<Op::Tran> cen_ops;

  struct Iter {
    const SymOps& parent;
    unsigned inv, n_symop, n_cenop;
    void operator++() {
      if (++n_cenop == parent.cen_ops.size()) {
        n_cenop = 0;
        if (++n_symop = parent.sym_ops.size()) {
          n_symop = 0;
          if (inv == 0 && parent.centrosym)
            inv = 1;
          else
            inv = 2;
        }
      }
    }
    Op operator*() const {
      Op op = parent.sym_ops[n_symop];
      if (inv)
        ; //TODO op.change_signs();
      op.translate(parent.cen_ops[n_cenop]);
      return op;
    }
    bool operator==(const Iter& o) const {
      return inv == o.inv && n_symop == o.n_symop && n_cenop == o.n_cenop;
    }
    bool operator!=(const Iter& o) const { return !(*this == o); }
  };

  Iter begin() const { return {*this, 0, 0, 0}; };
  Iter end() const { return {*this, 2, 0, 0}; };
};

inline std::vector<Op::Tran> lattice_translations(char lattice_symbol) {
  switch (lattice_symbol & ~0x20) {
    case 'P': return {{0, 0, 0}};
    case 'A': return {{0, 0, 0}, {0, 6, 6}};
    case 'B': return {{0, 0, 0}, {6, 0, 6}};
    case 'C': return {{0, 0, 0}, {6, 6, 0}};
    case 'I': return {{0, 0, 0}, {6, 6, 6}};
    case 'R': return {{0, 0, 0}, {8, 4, 4}, {4, 8, 8}};
    case 'S': return {{0, 0, 0}, {4, 4, 8}, {8, 4, 8}};
    case 'T': return {{0, 0, 0}, {4, 8, 4}, {8, 4, 8}};
    case 'F': return {{0, 0, 0}, {0, 6, 6}, {6, 0, 6}, {6, 6, 0}};
    default: fail("not a lattice symbol: " + std::string(1, lattice_symbol));
  }
}

inline Op::Tran translation_from_symbol(char symbol) {
  switch (symbol) {
    case 'a': return {6, 0, 0};
    case 'b': return {0, 6, 0};
    case 'c': return {0, 0, 6};
    case 'n': return {6, 6, 6};
    case 'u': return {3, 0, 0};
    case 'v': return {0, 3, 0};
    case 'w': return {0, 0, 3};
    case 'd': return {3, 3, 3};
    default: fail("not a translation symbol: " + std::string(1, symbol));
  }
}

inline const char* skip_blank(const char* p) {
  if (p)
    while (*p == ' ' || *p == '\t')
      ++p;
  return p;
}

inline Op hall_matrix_symbol(const char* start, const char* end) {
  Op op;
  // TODO
  return op;
}

inline Op::Tran parse_translation(const char* start, const char* end) {
  Op::Tran t;
  char* endptr;
  for (int i = 0; i != 3; ++i) {
    t[i] = std::strtol(start, &endptr, 10) % 12;
    start = endptr;
  }
  if (endptr != end)
    fail("wrong format of translation: " + std::string(start, end));
  return t;
}

inline const char* find_chr(const char* s, int c, const char* end) {
  return static_cast<const char*>(std::memchr(s, c, end - s));
}

inline SymOps symops_from_hall(const char* hall) {
  if (hall == nullptr)
    fail("null");
  hall = skip_blank(hall);
  SymOps ops;
  ops.centrosym = (hall[0] == '-');
  ops.sym_ops.emplace_back(Op::identity());
  const char* lat = skip_blank(ops.centrosym ? hall + 1 : hall);
  ops.cen_ops = lattice_translations(*lat);
  const char* part = skip_blank(lat + 1);
  const char* end = part + std::strlen(part);
  while (part && part != end) {
    if (*part == '(') {
      const char* rb = find_chr(part, ')', end);
      if (!rb)
        fail("missing ')': " + std::string(hall));
      if (ops.sym_ops.empty())
        fail("misplaced translation: " + std::string(hall));
      Op::Tran tr = parse_translation(part + 1, rb);
      ops.sym_ops.back().shift_origin(tr);
      part = skip_blank(find_chr(rb + 1, ' ', end));
    } else {
      const char* space = find_chr(part, ' ', end);
      Op op = hall_matrix_symbol(part, space ? space : end);
      if (op != Op::identity())
        ops.sym_ops.emplace_back(op);
      part = skip_blank(space);
    }
  }
  return ops;
}


} // namespace sym
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
