// Copyright 2017 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdint>
#include <cstdlib>    // for strtol
#include <cstring>    // for memchr, strchr, strlen
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
  enum { TDEN=12 };  // may be changed to 24 for change-of-basis ops
  typedef int int_t;
  typedef std::array<std::array<int_t, 3>, 3> Rot;
  typedef std::array<int_t, 3> Tran;

  Rot rot;
  Tran tran;

  std::string triplet() const;
  std::string rot_triplet() const { return Op{rot, {0, 0, 0}}.triplet(); };

  Op inverted() const;

  // if the translation points outside of the unit cell, wrap it.
  Op& unitize() {
    for (int i = 0; i != 3; ++i) {
      if (tran[i] >= TDEN) // elements need to be in [0,TDEN)
        tran[i] %= TDEN;
      else if (tran[i] < 0)
        tran[i] = (tran[i] % TDEN) + TDEN;
    }
    return *this;
  }

  Op& translate(const Tran& a) {
    for (int i = 0; i != 3; ++i)
      tran[i] += a[i];
    return *this;
  }

  Op translated(const Tran& a) const { return Op(*this).translate(a); }

  Rot negated_rot() const {
    return { -rot[0][0], -rot[0][1], -rot[0][2],
             -rot[1][0], -rot[1][1], -rot[1][2],
             -rot[2][0], -rot[2][1], -rot[2][2] };
  }

  Op negated() { return { negated_rot(), { -tran[0], -tran[1], -tran[2] } }; }

  Op shifted_origin(const Tran& a) const;
  Op changed_basis(const Op& cob) const;

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
  r.unitize();
  return r;
}

inline Op Op::shifted_origin(const Tran& a) const {
  Op::Rot rid = Op::identity().rot;
  return combine(combine(Op{rid, a}, *this), Op{rid, {-a[0], -a[1], -a[2]}});
}

inline Op Op::changed_basis(const Op& cob) const {
  return combine(combine(cob, *this), cob.inverted());
}


inline Op Op::inverted() const {
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

inline std::array<Op::int_t, 4> parse_triplet_part(const std::string& s) {
  std::array<Op::int_t, 4> r = { 0, 0, 0, 0 };
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
        if (den < 1 || Op::TDEN % den != 0)
          fail("Unexpected denominator " + std::to_string(den) + " in: " + s);
      }
      r[3] *= Op::TDEN / den;
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
    for (int factor : {2, 2, 3})  // for Op::TDEN == 12
      if (w % factor == 0)
        w /= factor;
      else
        denom *= factor;
    s += (w > 0 && !s.empty() ? "+" : "") + std::to_string(w);
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


// LIST OF CRYSTALLOGRAPHIC SPACE GROUPS

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
// based on both ITfC vol.B ch.1.4 (2010)
// and http://cci.lbl.gov/sginfo/hall_symbols.html

struct SymOps {
  std::vector<Op> sym_ops;
  std::vector<Op::Tran> cen_ops;

  void add_missing_group_elements();

  const Op* find_by_rotation(const Op::Rot& r) const {
    for (const Op& op : sym_ops)
      if (op.rot == r)
        return &op;
    return nullptr;
  }

  struct Iter {
    const SymOps& parent;
    unsigned n_sym, n_cen;
    void operator++() {
      if (++n_cen == parent.cen_ops.size()) {
        ++n_sym;
        n_cen = 0;
      }
    }
    Op operator*() const {
      return parent.sym_ops[n_sym].translated(parent.cen_ops[n_cen]).unitize();
    }
    bool operator==(const Iter& other) const {
      return n_sym == other.n_sym && n_cen == other.n_cen;
    }
    bool operator!=(const Iter& other) const { return !(*this == other); }
  };

  Iter begin() const { return {*this, 0, 0}; };
  Iter end() const { return {*this, (unsigned) sym_ops.size(), 0}; };
};

// corresponds to Table A1.4.2.2 in ITfC vol.B (edition 2010)
inline std::vector<Op::Tran> lattice_translations(char lattice_symbol) {
  switch (lattice_symbol & ~0x20) {
    case 'P': return {{0, 0, 0}};
    case 'A': return {{0, 0, 0}, {0, 6, 6}};
    case 'B': return {{0, 0, 0}, {6, 0, 6}};
    case 'C': return {{0, 0, 0}, {6, 6, 0}};
    case 'I': return {{0, 0, 0}, {6, 6, 6}};
    case 'R': return {{0, 0, 0}, {8, 4, 4}, {4, 8, 8}};
    // hall_symbols.html has no H, ITfC 2010 has no S and T
    case 'S': return {{0, 0, 0}, {4, 4, 8}, {8, 4, 8}};
    case 'T': return {{0, 0, 0}, {4, 8, 4}, {8, 4, 8}};
    case 'H': return {{0, 0, 0}, {8, 4, 0}, {4, 8, 0}};
    case 'F': return {{0, 0, 0}, {0, 6, 6}, {6, 0, 6}, {6, 6, 0}};
    default: fail("not a lattice symbol: " + std::string(1, lattice_symbol));
  }
}

// matrices for Nz from Table 3 and 4 from hall_symbols.html
inline Op::Rot rotation_z(int N) {
  switch (N) {
    case 1: return {1,0,0,  0,1,0,  0,0,1};
    case 2: return {-1,0,0, 0,-1,0, 0,0,1};
    case 3: return {0,-1,0, 1,-1,0, 0,0,1};
    case 4: return {0,-1,0, 1,0,0,  0,0,1};
    case 6: return {1,-1,0, 1,0,0,  0,0,1};
    case '\'': return {0,-1,0, -1,0,0, 0,0,-1};
    case '"':  return {0,1,0,   1,0,0, 0,0,-1};
    case '*':  return {0,0,1,   1,0,0, 0,1,0};
    default: fail("incorrect axis definition");
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
    default: fail("unknown symbol: " + std::string(1, symbol));
  }
}

inline Op::Rot alter_order(const Op::Rot& r, int i, int j, int k) {
    return { r[i][i], r[i][j], r[i][k],
             r[j][i], r[j][j], r[j][k],
             r[k][i], r[k][j], r[k][k] };
}

inline Op hall_matrix_symbol(const char* start, const char* end,
                             int pos, char first) {
  Op op = Op::identity();
  bool neg = (*start == '-');
  const char* p = (neg ? start + 1 : start);
  if (*p < '1' || *p == '5' || *p > '6')
    fail("wrong n-fold order notation: " + std::string(start, end));
  int N = *p++ - '0';
  int fractional_tran = 0;
  char principal_axis = '\0';
  char diagonal_axis = '\0';
  for (; p < end; ++p) {
    if (*p >= '1' && *p <= '5') {
      if (fractional_tran != '\0')
        fail("two numeric subscripts");
      fractional_tran = *p - '0';
    } else if (*p == '\'' || *p == '"' || *p == '*') {
      if (N != (*p == '*' ? 3 : 2))
        fail("wrong symbol: " + std::string(start, end));
      diagonal_axis = *p;
    } else if (*p == 'x' || *p == 'y' || *p == 'z') {
      principal_axis = *p;
    } else {
      op.translate(translation_from_symbol(*p));
    }
  }
  // fill in implicit values
  if (!principal_axis && !diagonal_axis) {
    if (pos == 1) {
      principal_axis = 'z';
    } else if (pos == 2 && N == 2) {
      if (first == '2' || first == '4')
        principal_axis = 'x';
      else if (first == '3' || first == '6')
        diagonal_axis = '\'';
    } else if (pos == 3 && N == 3) {
      diagonal_axis = '*';
    } else {
      fail("missing axis");
    }
  }
  // get the operation
  op.rot = rotation_z(diagonal_axis ? diagonal_axis : N);
  if (neg)
    op.rot = op.negated_rot();
  if (principal_axis == 'x')
    op.rot = alter_order(op.rot, 2, 0, 1);
  else if (principal_axis == 'y')
    op.rot = alter_order(op.rot, 1, 2, 0);
  if (fractional_tran)
    op.tran[principal_axis - 'x'] += 12 / N * fractional_tran;
  return op;
}

inline Op parse_change_of_basis(const char* start, const char* end) {
  if (memchr(start, ',', end - start) != nullptr) { // long symbol (x,y,z+1/12)
    // note: we do not parse multipliers such as 1/2x
    return parse_triplet(std::string(start, end));
  }
  // short symbol (0 0 1)
  Op cob = Op::identity();
  char* endptr;
  for (int i = 0; i != 3; ++i) {
    cob.tran[i] = std::strtol(start, &endptr, 10) % 12;
    start = endptr;
  }
  if (endptr != end)
    fail("unexpected change-of-basis format: " + std::string(start, end));
  return cob;
}

void SymOps::add_missing_group_elements() {
  // Brute force. To be replaced with Dimino's algorithm
  // see Luc Bourhis' answer https://physics.stackexchange.com/a/351400/95713
  if (sym_ops.empty() || sym_ops[0] != Op::identity())
    fail("oops");
  size_t generator_count = sym_ops.size();
  size_t prev_size = 0;
  while (prev_size != sym_ops.size()) {
    prev_size = sym_ops.size();
    for (size_t i = 1; i != prev_size; ++i)
      for (size_t j = 1; j != generator_count; ++j) {
        Op new_op = combine(sym_ops[i], sym_ops[j]);
        if (find_by_rotation(new_op.rot) == nullptr)
          sym_ops.push_back(new_op);
      }
    if (sym_ops.size() > 1000)
      fail("1000+ elements in the group should not happen");
  }
}

inline const char* find_blank(const char* p) {
  while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '_') // '_' == ' '
    ++p;
  return p;
}

inline const char* skip_blank(const char* p) {
  if (p)
    while (*p == ' ' || *p == '\t' || *p == '_') // '_' can be used as space
      ++p;
  return p;
}

inline SymOps generators_from_hall(const char* hall) {
  if (hall == nullptr)
    fail("null");
  hall = skip_blank(hall);
  SymOps ops;
  ops.sym_ops.emplace_back(Op::identity());
  bool centrosym = (hall[0] == '-');
  if (centrosym)
    ops.sym_ops.emplace_back(Op::identity().negated());
  const char* lat = skip_blank(centrosym ? hall + 1 : hall);
  if (!lat)
    fail("not a hall symbol: " + std::string(hall));
  ops.cen_ops = lattice_translations(*lat);
  int counter = 0;
  char first = '\0';
  const char* part = skip_blank(lat + 1);
  while (*part != '\0' && *part != '(') {
    const char* space = find_blank(part);
    if (!first)
      first = part[0];
    ++counter;
    if (part[0] != '1' || (part[1] != ' ' && part[1] != '\0')) {
      Op op = hall_matrix_symbol(part, space, counter, first);
      ops.sym_ops.emplace_back(op);
    }
    part = skip_blank(space);
  }
  if (*part == '(') {
    const char* rb = std::strchr(part, ')');
    if (!rb)
      fail("missing ')': " + std::string(hall));
    if (ops.sym_ops.empty())
      fail("misplaced translation: " + std::string(hall));
    Op cob = parse_change_of_basis(part + 1, rb);
    for (auto op = ops.sym_ops.begin() + 1; op != ops.sym_ops.end(); ++op)
      *op = op->changed_basis(cob);
    if (*skip_blank(find_blank(rb + 1)) != '\0')
      fail("unexpected characters after ')': " + std::string(hall));
  }
  return ops;
}

inline SymOps symops_from_hall(const char* hall) {
  SymOps symops = generators_from_hall(hall);
  symops.add_missing_group_elements();
  return symops;
}


} // namespace sym
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
