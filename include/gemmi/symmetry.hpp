// Copyright 2017 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdint>
#include <cstdlib>    // for strtol
#include <cstring>    // for memchr, strchr, strlen
#include <array>
#include <algorithm>  // for count, sort, remove
#include <functional> // for hash
#include <stdexcept>  // for runtime_error, invalid_argument
#include <string>
#include <tuple>      // for tie
#include <vector>

// we use brace elision with std:array's
#ifdef __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
#endif

namespace gemmi {

// UTILS

namespace impl {

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

inline const char* skip_blank(const char* p) {
  if (p)
    while (*p == ' ' || *p == '\t' || *p == '_') // '_' can be used as space
      ++p;
  return p;
}

// much faster than s += std::to_string(n) for n in 0 ... 99
inline void append_small_number(std::string& s, int n) {
  if (n < 0 || n >= 100) {
    s += std::to_string(n);
  } else if (n < 10) {
    s += char('0' + n);
  } else { // 10 ... 99
    int tens = n / 10;
    s += char('0' + tens);
    s += char('0' + n - 10 * tens);
  }
}

} // namespace impl

// TRIPLET <-> SYM OP

struct Op {
  enum { TDEN=24 };  // 24 to handle 1/8 in change-of-basis
  typedef int int_t;
  typedef std::array<std::array<int_t, 3>, 3> Rot;
  typedef std::array<int_t, 3> Tran;

  Rot rot;
  Tran tran;

  std::string triplet() const;
  std::string rot_triplet() const { return Op{rot, {0, 0, 0}}.triplet(); };

  Op inverse(int* denom=nullptr) const;

  // if the translation points outside of the unit cell, wrap it.
  Op& wrap() {
    for (int i = 0; i != 3; ++i) {
      if (tran[i] >= TDEN) // elements need to be in [0,TDEN)
        tran[i] %= TDEN;
      else if (tran[i] < 0)
        tran[i] = ((tran[i] + 1) % TDEN) + TDEN - 1;
    }
    return *this;
  }

  Op& translate(const Tran& a) {
    for (int i = 0; i != 3; ++i)
      tran[i] += a[i];
    return *this;
  }

  Op translated(const Tran& a) const { return Op(*this).translate(a); }

  Op add_centering(const Tran& a) const { return translated(a).wrap(); }

  Rot negated_rot() const {
    return { -rot[0][0], -rot[0][1], -rot[0][2],
             -rot[1][0], -rot[1][1], -rot[1][2],
             -rot[2][0], -rot[2][1], -rot[2][2] };
  }

  Op negated() { return { negated_rot(), { -tran[0], -tran[1], -tran[2] } }; }

  int det_rot() const { // should be 1 (rotation) or -1 (with inversion)
    return rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
         - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
         + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0]);
  }

  Op combine(const Op& b) const {
    Op r;
    for (int i = 0; i != 3; ++i) {
      r.tran[i] = tran[i];
      for (int j = 0; j != 3; ++j) {
        r.rot[i][j] = rot[i][0] * b.rot[0][j] +
                      rot[i][1] * b.rot[1][j] +
                      rot[i][2] * b.rot[2][j];
        r.tran[i] += rot[i][j] * b.tran[j];
      }
    }
    return r;
  }

  template<typename T>
  void apply_in_place_mult(T& x, T& y, T& z, T mult) const {
    T ret[3];
    for (int i = 0; i != 3; ++i)
      ret[i] = rot[i][0] * x + rot[i][1] * y + rot[i][2] * z + tran[i] * mult;
    x = ret[0];
    y = ret[1];
    z = ret[2];
  }

  template<typename T>
  void apply_in_place(T& x, T& y, T& z) const {
    apply_in_place_mult(x, y, z, 1.0 / Op::TDEN);
  }

  double phase_shift(int h, int k, int l) const {
    constexpr double mult = -2 * 3.1415926535897932384626433832795 / Op::TDEN;
    return mult * (h * tran[0] + k * tran[1] + l * tran[2]);
  }

  std::array<std::array<int_t, 4>, 4> seitz() const {
    std::array<std::array<int_t, 4>, 4> t;
    for (int i = 0; i < 3; ++i)
      t[i] = { rot[i][0], rot[i][1], rot[i][2], tran[i] };
    t[3] = { 0, 0, 0, 1 };
    return t;
  }

  std::array<std::array<double, 4>, 4> float_seitz() const {
    std::array<std::array<double, 4>, 4> t;
    for (int i = 0; i < 3; ++i)
      t[i] = { double(rot[i][0]), double(rot[i][1]), double(rot[i][2]),
               tran[i] / double(Op::TDEN) };
    t[3] = { 0., 0., 0., 1. };
    return t;
  }

  static constexpr Op identity() { return {{1,0,0, 0,1,0, 0,0,1}, {0,0,0}}; };
  bool operator<(const Op& rhs) const {
    return std::tie(rot, tran) < std::tie(rhs.rot, rhs.tran);
  }
};

inline bool operator==(const Op& a, const Op& b) {
  return a.rot == b.rot && a.tran == b.tran;
}
inline bool operator!=(const Op& a, const Op& b) { return !(a == b); }

inline Op operator*(const Op& a, const Op& b) { return a.combine(b).wrap(); }
inline Op& operator*=(Op& a, const Op& b) { a = a * b; return a; }

inline Op Op::inverse(int* denom) const {
  int detr = det_rot();
  int sign = detr >= 0 ? 1 : -1;
  if (detr != 1 && detr != -1) {
    if (detr == 0 || denom == nullptr)
      impl::fail("not a rotation/inversion: det|" + rot_triplet() + "|="
                 + std::to_string(detr));
    *denom = detr * sign;
  }
  Op inv;
  inv.rot[0][0] = sign * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2]);
  inv.rot[0][1] = sign * (rot[0][2] * rot[2][1] - rot[0][1] * rot[2][2]);
  inv.rot[0][2] = sign * (rot[0][1] * rot[1][2] - rot[0][2] * rot[1][1]);
  inv.rot[1][0] = sign * (rot[1][2] * rot[2][0] - rot[1][0] * rot[2][2]);
  inv.rot[1][1] = sign * (rot[0][0] * rot[2][2] - rot[0][2] * rot[2][0]);
  inv.rot[1][2] = sign * (rot[1][0] * rot[0][2] - rot[0][0] * rot[1][2]);
  inv.rot[2][0] = sign * (rot[1][0] * rot[2][1] - rot[2][0] * rot[1][1]);
  inv.rot[2][1] = sign * (rot[2][0] * rot[0][1] - rot[0][0] * rot[2][1]);
  inv.rot[2][2] = sign * (rot[0][0] * rot[1][1] - rot[1][0] * rot[0][1]);
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
    if (*c == '+') {
      sign = 1;
      continue;
    }
    if (*c == '-') {
      sign = -1;
      continue;
    }
    if (*c == ' ' || *c == '\t')
      continue;
    if (sign == 0)
      impl::fail("wrong or unsupported triplet format: " + s);
    if (*c >= '0' && *c <= '9') {
      char* endptr;
      r[3] = sign * strtol(c, &endptr, 10);
      int den = 1;
      if (*endptr == '/') {
        den = strtol(endptr + 1, &endptr, 10);
        if (den < 1 || Op::TDEN % den != 0)
          impl::fail("Unexpected denominator " + std::to_string(den) + " in: "
                     + s);
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
      impl::fail(std::string("unexpected character '") + *c + "' in: " + s);
    }
    sign = 0;
  }
  if (sign != 0)
    impl::fail("trailing sign in: " + s);
  return r;
}

inline Op parse_triplet(const std::string& s) {
  if (std::count(s.begin(), s.end(), ',') != 2)
    impl::fail("expected exactly two commas in triplet");
  size_t comma1 = s.find(',');
  size_t comma2 = s.find(',', comma1 + 1);
  auto a = parse_triplet_part(s.substr(0, comma1));
  auto b = parse_triplet_part(s.substr(comma1 + 1, comma2 - (comma1 + 1)));
  auto c = parse_triplet_part(s.substr(comma2 + 1));
  Op::Rot rot = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
  Op::Tran tran = {a[3], b[3], c[3]};
  return { rot, tran };
}

inline std::string make_triplet_part(int x, int y, int z, int w,
                                     char style='x') {
  std::string s;
  int xyz[] = { x, y, z };
  for (int i = 0; i != 3; ++i)
    if (xyz[i] != 0) {
      if (xyz[i] < 0)
        s += '-';
      else if (!s.empty())
        s += '+';
      if (xyz[i] != 1 && xyz[i] != -1)
        impl::append_small_number(s, xyz[i] < 0 ? -xyz[i] : xyz[i]);
      s += char(style + i);
    }
  if (w != 0) { // reduce the w/TDEN fraction to lowest terms
    int denom = 1;
    for (int factor : {2, 2, 2, 3})  // for Op::TDEN == 24
      if (w % factor == 0)
        w /= factor;
      else
        denom *= factor;
    if (w > 0) {
      if (!s.empty())
        s += '+';
      impl::append_small_number(s, w);
    } else {
      s += '-';
      impl::append_small_number(s, -w);
    }
    if (denom != 1) {
      s += '/';
      impl::append_small_number(s, denom);
    }
  }
  return s;
}

inline std::string Op::triplet() const {
  return make_triplet_part(rot[0][0], rot[0][1], rot[0][2], tran[0]) +
   "," + make_triplet_part(rot[1][0], rot[1][1], rot[1][2], tran[1]) +
   "," + make_triplet_part(rot[2][0], rot[2][1], rot[2][2], tran[2]);
}


// GROUPS OF OPERATIONS

// corresponds to Table A1.4.2.2 in ITfC vol.B (edition 2010)
inline std::vector<Op::Tran> centring_vectors(char lattice_symbol) {
  constexpr int h = Op::TDEN / 2;
  constexpr int t = Op::TDEN / 3;
  constexpr int d = 2 * t;
  switch (lattice_symbol & ~0x20) {
    case 'P': return {{0, 0, 0}};
    case 'A': return {{0, 0, 0}, {0, h, h}};
    case 'B': return {{0, 0, 0}, {h, 0, h}};
    case 'C': return {{0, 0, 0}, {h, h, 0}};
    case 'I': return {{0, 0, 0}, {h, h, h}};
    case 'R': return {{0, 0, 0}, {t, d, d}, {d, t, t}};
    // hall_symbols.html has no H, ITfC 2010 has no S and T
    case 'S': return {{0, 0, 0}, {t, t, d}, {d, t, d}};
    case 'T': return {{0, 0, 0}, {t, d, t}, {d, t, d}};
    case 'H': return {{0, 0, 0}, {t, d, 0}, {d, t, 0}};
    case 'F': return {{0, 0, 0}, {0, h, h}, {h, 0, h}, {h, h, 0}};
    default: impl::fail(std::string("not a lattice symbol: ") + lattice_symbol);
  }
}


struct GroupOps {
  std::vector<Op> sym_ops;
  std::vector<Op::Tran> cen_ops;

  int order() const { return static_cast<int>(sym_ops.size()*cen_ops.size()); }

  void add_missing_elements();

  char find_centering() const {
    if (cen_ops.size() == 1 && cen_ops[0] == Op::Tran{0, 0, 0})
      return 'P';
    std::vector<Op::Tran> trans = cen_ops;
    std::sort(trans.begin(), trans.end());
    for (char c : {'A', 'B', 'C', 'I', 'F', 'R', 'S', 'T', 'H'})
      if (trans == centring_vectors(c))
        return c;
    return 0;
  }

  const Op* find_by_rotation(const Op::Rot& r) const {
    for (const Op& op : sym_ops)
      if (op.rot == r)
        return &op;
    return nullptr;
  }

  bool is_centric() const { return find_by_rotation({-1,0,0, 0,-1,0, 0,0,-1}); }

  void change_basis(const Op& cob) {
    if (sym_ops.empty() || cen_ops.empty())
      return;
    int den = 1;
    Op cob_inv = cob.inverse(&den);
    // assuming the first items in sym_ops and cen_ops are identities
    for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
      *op = cob.combine(*op).combine(cob_inv).wrap();
    constexpr Op::Rot rid = Op::identity().rot;
    for (auto tr = cen_ops.begin() + 1; tr != cen_ops.end(); ++tr)
      *tr = cob.combine(Op{rid, *tr}).combine(cob_inv).wrap().tran;
    if (den > 1) {  // this part was tested only on R 3 :R -> :H
      for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
        for (auto& i : op->rot)
          for (auto& j : i)
            j /= den;
      // the centering must have changed
      for (int i = static_cast<int>(cen_ops.size()) - 1; i >= 0; --i)
        for (int j = i - 1; j >= 0; --j)
          if (cen_ops[i] == cen_ops[j]) {
            cen_ops.erase(cen_ops.begin() + i);
            break;
          }
    }
  }

  std::vector<Op> all_ops_sorted() const {
    std::vector<Op> ops;
    ops.reserve(sym_ops.size() * cen_ops.size());
    for (const Op& so : sym_ops)
      for (const Op::Tran& co : cen_ops)
        ops.push_back(so.add_centering(co));
    std::sort(ops.begin(), ops.end());
    return ops;
  }

  Op get_op(int n) const {
    int n_cen = n / sym_ops.size();
    int n_sym = n % sym_ops.size();
    return sym_ops.at(n_sym).add_centering(cen_ops.at(n_cen));
  }

  bool is_same_as(const GroupOps& other) const {
    if (cen_ops.size() != other.cen_ops.size() ||
        sym_ops.size() != other.sym_ops.size())
      return false;
    return all_ops_sorted() == other.all_ops_sorted();
  }

  // minimal multiplicity for real-space grid in each direction
  // examples: 1,2,1 for P21, 1,1,6 for P61
  std::array<int, 3> find_grid_factors() const {
    const int T = Op::TDEN;
    int r[3] = {T, T, T};
    for (const Op& op : *this)
      for (int i = 0; i != 3; ++i)
        if (op.tran[i] != 0 && op.tran[i] < r[i])
          r[i] = op.tran[i];
    return {T / r[0], T / r[1], T / r[2]};
  }

  bool are_directions_symmetry_related(int u, int v) const {
    for (const Op& op : sym_ops)
      if (op.rot[u][v] != 0)
        return true;
    return false;
  }

  struct Iter {
    const GroupOps& gops;
    int n_sym, n_cen;
    void operator++() {
      if (++n_sym == (int) gops.sym_ops.size()) {
        ++n_cen;
        n_sym = 0;
      }
    }
    Op operator*() const {
      return gops.sym_ops.at(n_sym).translated(gops.cen_ops.at(n_cen)).wrap();
    }
    bool operator==(const Iter& other) const {
      return n_sym == other.n_sym && n_cen == other.n_cen;
    }
    bool operator!=(const Iter& other) const { return !(*this == other); }
  };

  Iter begin() const { return {*this, 0, 0}; };
  Iter end() const { return {*this, 0, (int) cen_ops.size()}; };
};

inline void GroupOps::add_missing_elements() {
  // We always keep identity as sym_ops[0].
  if (sym_ops.empty() || sym_ops[0] != Op::identity())
    impl::fail("oops");
  if (sym_ops.size() == 1)
    return;
  auto check_size = [&]() {
    if (sym_ops.size() > 1023)
      impl::fail("1000+ elements in the group should not happen");
  };
  // Below we assume that all centring vectors are already known (in cen_ops)
  // so when checking for a new element we compare only the 3x3 matrix.
#if 1 // Dimino's algorithm
  // see Luc Bourhis' answer https://physics.stackexchange.com/a/351400/95713
  std::vector<Op> gen(sym_ops.begin() + 1, sym_ops.end());
  sym_ops.resize(2);
  const Op::Rot idrot = Op::identity().rot;
  for (Op g = sym_ops[1] * sym_ops[1]; g.rot != idrot; g *= sym_ops[1]) {
    sym_ops.push_back(g);
    check_size();
  }
  for (size_t i = 1; i < gen.size(); ++i) {
    std::vector<Op> coset_repr(1, Op::identity());
    size_t init_size = sym_ops.size();
    for (;;) {
      size_t len = coset_repr.size();
      for (size_t j = 0; j != len; ++j) {
        for (size_t n = 0; n != i + 1; ++n) {
          Op sg = gen[n] * coset_repr[j];
          if (find_by_rotation(sg.rot) == nullptr) {
            sym_ops.push_back(sg);
            for (size_t k = 1; k != init_size; ++k)
              sym_ops.push_back(sg * sym_ops[k]);
            coset_repr.push_back(sg);
          }
        }
      }
      if (len == coset_repr.size())
        break;
      check_size();
    }
  }
#else // brute force
  size_t generator_count = sym_ops.size();
  size_t prev_size = 0;
  while (prev_size != sym_ops.size()) {
    prev_size = sym_ops.size();
    for (size_t i = 1; i != prev_size; ++i)
      for (size_t j = 1; j != generator_count; ++j) {
        Op new_op = sym_ops[i].combine(sym_ops[j]);
        if (find_by_rotation(new_op.rot) == nullptr)
          sym_ops.push_back(new_op.wrap());
      }
    check_size();
  }
#endif
}


// INTERPRETING HALL SYMBOLS
// based on both ITfC vol.B ch.1.4 (2010)
// and http://cci.lbl.gov/sginfo/hall_symbols.html

// matrices for Nz from Table 3 and 4 from hall_symbols.html
inline Op::Rot hall_rotation_z(int N) {
  switch (N) {
    case 1: return {1,0,0,  0,1,0,  0,0,1};
    case 2: return {-1,0,0, 0,-1,0, 0,0,1};
    case 3: return {0,-1,0, 1,-1,0, 0,0,1};
    case 4: return {0,-1,0, 1,0,0,  0,0,1};
    case 6: return {1,-1,0, 1,0,0,  0,0,1};
    case '\'': return {0,-1,0, -1,0,0, 0,0,-1};
    case '"':  return {0,1,0,   1,0,0, 0,0,-1};
    case '*':  return {0,0,1,   1,0,0, 0,1,0};
    default: impl::fail("incorrect axis definition");
  }
}
inline Op::Tran hall_translation_from_symbol(char symbol) {
  constexpr int h = Op::TDEN / 2;
  constexpr int q = Op::TDEN / 4;
  switch (symbol) {
    case 'a': return {h, 0, 0};
    case 'b': return {0, h, 0};
    case 'c': return {0, 0, h};
    case 'n': return {h, h, h};
    case 'u': return {q, 0, 0};
    case 'v': return {0, q, 0};
    case 'w': return {0, 0, q};
    case 'd': return {q, q, q};
    default: impl::fail(std::string("unknown symbol: ") + symbol);
  }
}

inline Op hall_matrix_symbol(const char* start, const char* end,
                             int pos, int& prev) {
  Op op = Op::identity();
  bool neg = (*start == '-');
  const char* p = (neg ? start + 1 : start);
  if (*p < '1' || *p == '5' || *p > '6')
    impl::fail("wrong n-fold order notation: " + std::string(start, end));
  int N = *p++ - '0';
  int fractional_tran = 0;
  char principal_axis = '\0';
  char diagonal_axis = '\0';
  for (; p < end; ++p) {
    if (*p >= '1' && *p <= '5') {
      if (fractional_tran != '\0')
        impl::fail("two numeric subscripts");
      fractional_tran = *p - '0';
    } else if (*p == '\'' || *p == '"' || *p == '*') {
      if (N != (*p == '*' ? 3 : 2))
        impl::fail("wrong symbol: " + std::string(start, end));
      diagonal_axis = *p;
    } else if (*p == 'x' || *p == 'y' || *p == 'z') {
      principal_axis = *p;
    } else {
      op.translate(hall_translation_from_symbol(*p));
    }
  }
  // fill in implicit values
  if (!principal_axis && !diagonal_axis) {
    if (pos == 1) {
      principal_axis = 'z';
    } else if (pos == 2 && N == 2) {
      if (prev == 2 || prev == 4)
        principal_axis = 'x';
      else if (prev == 3 || prev == 6)
        diagonal_axis = '\'';
    } else if (pos == 3 && N == 3) {
      diagonal_axis = '*';
    } else if (N != 1) {
      impl::fail("missing axis");
    }
  }
  // get the operation
  op.rot = hall_rotation_z(diagonal_axis ? diagonal_axis : N);
  if (neg)
    op.rot = op.negated_rot();
  auto alter_order = [](const Op::Rot& r, int i, int j, int k) {
    return Op::Rot{ r[i][i], r[i][j], r[i][k],
                    r[j][i], r[j][j], r[j][k],
                    r[k][i], r[k][j], r[k][k] };
  };
  if (principal_axis == 'x')
    op.rot = alter_order(op.rot, 2, 0, 1);
  else if (principal_axis == 'y')
    op.rot = alter_order(op.rot, 1, 2, 0);
  if (fractional_tran)
    op.tran[principal_axis - 'x'] += Op::TDEN / N * fractional_tran;
  prev = N;
  return op;
}

// Parses either short (0 0 1) or long notation (x,y,z+1/12)
// but without multpliers (such as 1/2x) to keep things simple for now.
inline Op parse_hall_change_of_basis(const char* start, const char* end) {
  if (memchr(start, ',', end - start) != nullptr) // long symbol
    return parse_triplet(std::string(start, end));
  // short symbol (0 0 1)
  Op cob = Op::identity();
  char* endptr;
  for (int i = 0; i != 3; ++i) {
    cob.tran[i] = std::strtol(start, &endptr, 10) % 12 * (Op::TDEN / 12);
    start = endptr;
  }
  if (endptr != end)
    impl::fail("unexpected change-of-basis format: " + std::string(start, end));
  return cob;
}

inline GroupOps generators_from_hall(const char* hall) {
  auto find_blank = [](const char* p) {
    while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '_') // '_' == ' '
      ++p;
    return p;
  };
  if (hall == nullptr)
    impl::fail("null");
  hall = impl::skip_blank(hall);
  GroupOps ops;
  ops.sym_ops.emplace_back(Op::identity());
  bool centrosym = (hall[0] == '-');
  if (centrosym)
    ops.sym_ops.emplace_back(Op::identity().negated());
  const char* lat = impl::skip_blank(centrosym ? hall + 1 : hall);
  if (!lat)
    impl::fail("not a hall symbol: " + std::string(hall));
  ops.cen_ops = centring_vectors(*lat);
  int counter = 0;
  int prev = 0;
  const char* part = impl::skip_blank(lat + 1);
  while (*part != '\0' && *part != '(') {
    const char* space = find_blank(part);
    ++counter;
    if (part[0] != '1' || (part[1] != ' ' && part[1] != '\0')) {
      Op op = hall_matrix_symbol(part, space, counter, prev);
      ops.sym_ops.emplace_back(op);
    }
    part = impl::skip_blank(space);
  }
  if (*part == '(') {
    const char* rb = std::strchr(part, ')');
    if (!rb)
      impl::fail("missing ')': " + std::string(hall));
    if (ops.sym_ops.empty())
      impl::fail("misplaced translation: " + std::string(hall));
    ops.change_basis(parse_hall_change_of_basis(part + 1, rb));

    if (*impl::skip_blank(find_blank(rb + 1)) != '\0')
      impl::fail("unexpected characters after ')': " + std::string(hall));
  }
  return ops;
}

inline GroupOps symops_from_hall(const char* hall) {
  GroupOps ops = generators_from_hall(hall);
  ops.add_missing_elements();
  return ops;
}

// CRYSTAL SYSTEMS AND POINT GROUPS

enum class CrystalSystem : unsigned char {
  Triclinic=0, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic
};

inline const char* crystal_system_str(CrystalSystem system) {
  static const char* names[7] = {
    "triclinic", "monoclinic", "orthorhombic", "tetragonal",
    "trigonal", "hexagonal", "cubic"
  };
  return names[static_cast<int>(system)];
}

enum class PointGroup : unsigned char {
  C1=0, Ci, C2, Cs, C2h, D2, C2v, D2h, C4, S4, C4h, D4, C4v, D2d, D4h, C3,
  C3i, D3, C3v, D3d, C6, C3h, C6h, D6, C6v, D3h, D6h, T, Th, O, Td, Oh
};

inline const char* point_group_hm(PointGroup pg) {
  static const char hm_pointgroup_names[32][6] = {
    "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm",
    "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm", "3",
    "-3", "32", "3m", "-3m", "6", "-6", "6/m", "622",
    "6mm", "-62m", "6/mmm", "23", "m-3", "432", "-43m", "m-3m",
  };
  return hm_pointgroup_names[static_cast<int>(pg)];
}

inline CrystalSystem crystal_system(PointGroup pg) {
  static const CrystalSystem crystal_systems[32] = {
    CrystalSystem::Triclinic,    CrystalSystem::Triclinic,
    CrystalSystem::Monoclinic,   CrystalSystem::Monoclinic,
    CrystalSystem::Monoclinic,   CrystalSystem::Orthorhombic,
    CrystalSystem::Orthorhombic, CrystalSystem::Orthorhombic,
    CrystalSystem::Tetragonal,   CrystalSystem::Tetragonal,
    CrystalSystem::Tetragonal,   CrystalSystem::Tetragonal,
    CrystalSystem::Tetragonal,   CrystalSystem::Tetragonal,
    CrystalSystem::Tetragonal,   CrystalSystem::Trigonal,
    CrystalSystem::Trigonal,     CrystalSystem::Trigonal,
    CrystalSystem::Trigonal,     CrystalSystem::Trigonal,
    CrystalSystem::Hexagonal,    CrystalSystem::Hexagonal,
    CrystalSystem::Hexagonal,    CrystalSystem::Hexagonal,
    CrystalSystem::Hexagonal,    CrystalSystem::Hexagonal,
    CrystalSystem::Hexagonal,    CrystalSystem::Cubic,
    CrystalSystem::Cubic,        CrystalSystem::Cubic,
    CrystalSystem::Cubic,        CrystalSystem::Cubic
  };
  return crystal_systems[static_cast<int>(pg)];
}

inline PointGroup point_group(int space_group_number) {
  static char indices[230] = {
     0,  1,  2,  2,  2,  3,  3,  3,  3,  4,
     4,  4,  4,  4,  4,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  6,  6,  6,  6,  6,  6,
     6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
     6,  6,  6,  6,  6,  6,  7,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  8,  8,  8,  8,  8,  8,
     9,  9, 10, 10, 10, 10, 10, 10, 11, 11,
    11, 11, 11, 11, 11, 11, 11, 11, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 14, 14, 14, 14, 14, 14, 14, 14,
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
    14, 14, 15, 15, 15, 15, 16, 16, 17, 17,
    17, 17, 17, 17, 17, 18, 18, 18, 18, 18,
    18, 19, 19, 19, 19, 19, 19, 20, 20, 20,
    20, 20, 20, 21, 22, 22, 23, 23, 23, 23,
    23, 23, 24, 24, 24, 24, 25, 25, 25, 25,
    26, 26, 26, 26, 27, 27, 27, 27, 27, 28,
    28, 28, 28, 28, 28, 28, 29, 29, 29, 29,
    29, 29, 29, 29, 30, 30, 30, 30, 30, 30,
    31, 31, 31, 31, 31, 31, 31, 31, 31, 31
  };
  return static_cast<PointGroup>(indices[space_group_number-1]);
}

// LIST OF CRYSTALLOGRAPHIC SPACE GROUPS

struct SpaceGroup { // typically 40 bytes
  int number;
  int ccp4;
  char hm[11];  // Hermann–Mauguin (international) notation
  char ext;
  char qualifier[5];
  char hall[15];

  std::string colon_ext() const { return ext ? std::string(":") + ext : ""; }
  std::string xhm() const { return hm + colon_ext(); }

  // P 1 2 1 -> P2, but P 1 1 2 -> P112. R 3:H -> H3.
  std::string short_name() const {
    std::string s(hm);
    size_t len = s.size();
    if (len > 6 && s[2] == '1' && s[len - 2] == ' ' && s[len - 1] == '1')
      s = s[0] + s.substr(4, len - 4 - 2);
    if (ext == 'H')
      s[0] = 'H';
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return s;
  }

  PointGroup point_group() const { return gemmi::point_group(number); }
  const char* point_group_hm() const {
    return gemmi::point_group_hm(point_group());
  }
  CrystalSystem crystal_system() const {
    return gemmi::crystal_system(point_group());
  }
  const char* crystal_system_str() const {
    return gemmi::crystal_system_str(crystal_system());
  }

  GroupOps operations() const { return symops_from_hall(hall); }
};

struct SpaceGroupAltName {
  char hm[11];
  char ext;
  int pos;
};

// the template here is only to substitute C++17 inline variables
// https://stackoverflow.com/questions/38043442/how-do-inline-variables-work
namespace impl {

template<class Dummy>
struct Tables_
{
  static const SpaceGroup main[554];
  static const SpaceGroupAltName alt_names[27];
};

template<class Dummy>
const SpaceGroup Tables_<Dummy>::main[554] = {
  // This table was generated by tools/gen_sg_table.py.
  // First 530 entries in the same order as in SgInfo, sgtbx and ITB.
  // Note: spacegroup 68 has three duplicates with different H-M names.
  {  1,    1, "P 1"       ,   0,     "", "P 1"           }, //   0
  {  2,    2, "P -1"      ,   0,     "", "-P 1"          }, //   1
  {  3,    3, "P 1 2 1"   ,   0,    "b", "P 2y"          }, //   2
  {  3, 1003, "P 1 1 2"   ,   0,    "c", "P 2"           }, //   3
  {  3,    0, "P 2 1 1"   ,   0,    "a", "P 2x"          }, //   4
  {  4,    4, "P 1 21 1"  ,   0,    "b", "P 2yb"         }, //   5
  {  4, 1004, "P 1 1 21"  ,   0,    "c", "P 2c"          }, //   6
  {  4,    0, "P 21 1 1"  ,   0,    "a", "P 2xa"         }, //   7
  {  5,    5, "C 1 2 1"   ,   0,   "b1", "C 2y"          }, //   8
  {  5, 2005, "A 1 2 1"   ,   0,   "b2", "A 2y"          }, //   9
  {  5, 4005, "I 1 2 1"   ,   0,   "b3", "I 2y"          }, //  10
  {  5,    0, "A 1 1 2"   ,   0,   "c1", "A 2"           }, //  11
  {  5, 1005, "B 1 1 2"   ,   0,   "c2", "B 2"           }, //  12
  {  5,    0, "I 1 1 2"   ,   0,   "c3", "I 2"           }, //  13
  {  5,    0, "B 2 1 1"   ,   0,   "a1", "B 2x"          }, //  14
  {  5,    0, "C 2 1 1"   ,   0,   "a2", "C 2x"          }, //  15
  {  5,    0, "I 2 1 1"   ,   0,   "a3", "I 2x"          }, //  16
  {  6,    6, "P 1 m 1"   ,   0,    "b", "P -2y"         }, //  17
  {  6, 1006, "P 1 1 m"   ,   0,    "c", "P -2"          }, //  18
  {  6,    0, "P m 1 1"   ,   0,    "a", "P -2x"         }, //  19
  {  7,    7, "P 1 c 1"   ,   0,   "b1", "P -2yc"        }, //  20
  {  7,    0, "P 1 n 1"   ,   0,   "b2", "P -2yac"       }, //  21
  {  7,    0, "P 1 a 1"   ,   0,   "b3", "P -2ya"        }, //  22
  {  7,    0, "P 1 1 a"   ,   0,   "c1", "P -2a"         }, //  23
  {  7,    0, "P 1 1 n"   ,   0,   "c2", "P -2ab"        }, //  24
  {  7, 1007, "P 1 1 b"   ,   0,   "c3", "P -2b"         }, //  25
  {  7,    0, "P b 1 1"   ,   0,   "a1", "P -2xb"        }, //  26
  {  7,    0, "P n 1 1"   ,   0,   "a2", "P -2xbc"       }, //  27
  {  7,    0, "P c 1 1"   ,   0,   "a3", "P -2xc"        }, //  28
  {  8,    8, "C 1 m 1"   ,   0,   "b1", "C -2y"         }, //  29
  {  8,    0, "A 1 m 1"   ,   0,   "b2", "A -2y"         }, //  30
  {  8,    0, "I 1 m 1"   ,   0,   "b3", "I -2y"         }, //  31
  {  8,    0, "A 1 1 m"   ,   0,   "c1", "A -2"          }, //  32
  {  8, 1008, "B 1 1 m"   ,   0,   "c2", "B -2"          }, //  33
  {  8,    0, "I 1 1 m"   ,   0,   "c3", "I -2"          }, //  34
  {  8,    0, "B m 1 1"   ,   0,   "a1", "B -2x"         }, //  35
  {  8,    0, "C m 1 1"   ,   0,   "a2", "C -2x"         }, //  36
  {  8,    0, "I m 1 1"   ,   0,   "a3", "I -2x"         }, //  37
  {  9,    9, "C 1 c 1"   ,   0,   "b1", "C -2yc"        }, //  38
  {  9,    0, "A 1 n 1"   ,   0,   "b2", "A -2yab"       }, //  39
  {  9,    0, "I 1 a 1"   ,   0,   "b3", "I -2ya"        }, //  40
  {  9,    0, "A 1 a 1"   ,   0,  "-b1", "A -2ya"        }, //  41
  {  9,    0, "C 1 n 1"   ,   0,  "-b2", "C -2yac"       }, //  42
  {  9,    0, "I 1 c 1"   ,   0,  "-b3", "I -2yc"        }, //  43
  {  9,    0, "A 1 1 a"   ,   0,   "c1", "A -2a"         }, //  44
  {  9,    0, "B 1 1 n"   ,   0,   "c2", "B -2ab"        }, //  45
  {  9,    0, "I 1 1 b"   ,   0,   "c3", "I -2b"         }, //  46
  {  9, 1009, "B 1 1 b"   ,   0,  "-c1", "B -2b"         }, //  47
  {  9,    0, "A 1 1 n"   ,   0,  "-c2", "A -2ab"        }, //  48
  {  9,    0, "I 1 1 a"   ,   0,  "-c3", "I -2a"         }, //  49
  {  9,    0, "B b 1 1"   ,   0,   "a1", "B -2xb"        }, //  50
  {  9,    0, "C n 1 1"   ,   0,   "a2", "C -2xac"       }, //  51
  {  9,    0, "I c 1 1"   ,   0,   "a3", "I -2xc"        }, //  52
  {  9,    0, "C c 1 1"   ,   0,  "-a1", "C -2xc"        }, //  53
  {  9,    0, "B n 1 1"   ,   0,  "-a2", "B -2xab"       }, //  54
  {  9,    0, "I b 1 1"   ,   0,  "-a3", "I -2xb"        }, //  55
  { 10,   10, "P 1 2/m 1" ,   0,    "b", "-P 2y"         }, //  56
  { 10, 1010, "P 1 1 2/m" ,   0,    "c", "-P 2"          }, //  57
  { 10,    0, "P 2/m 1 1" ,   0,    "a", "-P 2x"         }, //  58
  { 11,   11, "P 1 21/m 1",   0,    "b", "-P 2yb"        }, //  59
  { 11, 1011, "P 1 1 21/m",   0,    "c", "-P 2c"         }, //  60
  { 11,    0, "P 21/m 1 1",   0,    "a", "-P 2xa"        }, //  61
  { 12,   12, "C 1 2/m 1" ,   0,   "b1", "-C 2y"         }, //  62
  { 12,    0, "A 1 2/m 1" ,   0,   "b2", "-A 2y"         }, //  63
  { 12,    0, "I 1 2/m 1" ,   0,   "b3", "-I 2y"         }, //  64
  { 12,    0, "A 1 1 2/m" ,   0,   "c1", "-A 2"          }, //  65
  { 12, 1012, "B 1 1 2/m" ,   0,   "c2", "-B 2"          }, //  66
  { 12,    0, "I 1 1 2/m" ,   0,   "c3", "-I 2"          }, //  67
  { 12,    0, "B 2/m 1 1" ,   0,   "a1", "-B 2x"         }, //  68
  { 12,    0, "C 2/m 1 1" ,   0,   "a2", "-C 2x"         }, //  69
  { 12,    0, "I 2/m 1 1" ,   0,   "a3", "-I 2x"         }, //  70
  { 13,   13, "P 1 2/c 1" ,   0,   "b1", "-P 2yc"        }, //  71
  { 13,    0, "P 1 2/n 1" ,   0,   "b2", "-P 2yac"       }, //  72
  { 13,    0, "P 1 2/a 1" ,   0,   "b3", "-P 2ya"        }, //  73
  { 13,    0, "P 1 1 2/a" ,   0,   "c1", "-P 2a"         }, //  74
  { 13,    0, "P 1 1 2/n" ,   0,   "c2", "-P 2ab"        }, //  75
  { 13, 1013, "P 1 1 2/b" ,   0,   "c3", "-P 2b"         }, //  76
  { 13,    0, "P 2/b 1 1" ,   0,   "a1", "-P 2xb"        }, //  77
  { 13,    0, "P 2/n 1 1" ,   0,   "a2", "-P 2xbc"       }, //  78
  { 13,    0, "P 2/c 1 1" ,   0,   "a3", "-P 2xc"        }, //  79
  { 14,   14, "P 1 21/c 1",   0,   "b1", "-P 2ybc"       }, //  80
  { 14, 2014, "P 1 21/n 1",   0,   "b2", "-P 2yn"        }, //  81
  { 14,    0, "P 1 21/a 1",   0,   "b3", "-P 2yab"       }, //  82
  { 14,    0, "P 1 1 21/a",   0,   "c1", "-P 2ac"        }, //  83
  { 14,    0, "P 1 1 21/n",   0,   "c2", "-P 2n"         }, //  84
  { 14, 1014, "P 1 1 21/b",   0,   "c3", "-P 2bc"        }, //  85
  { 14,    0, "P 21/b 1 1",   0,   "a1", "-P 2xab"       }, //  86
  { 14,    0, "P 21/n 1 1",   0,   "a2", "-P 2xn"        }, //  87
  { 14,    0, "P 21/c 1 1",   0,   "a3", "-P 2xac"       }, //  88
  { 15,   15, "C 1 2/c 1" ,   0,   "b1", "-C 2yc"        }, //  89
  { 15,    0, "A 1 2/n 1" ,   0,   "b2", "-A 2yab"       }, //  90
  { 15,    0, "I 1 2/a 1" ,   0,   "b3", "-I 2ya"        }, //  91
  { 15,    0, "A 1 2/a 1" ,   0,  "-b1", "-A 2ya"        }, //  92
  { 15,    0, "C 1 2/n 1" ,   0,  "-b2", "-C 2yac"       }, //  93
  { 15,    0, "I 1 2/c 1" ,   0,  "-b3", "-I 2yc"        }, //  94
  { 15,    0, "A 1 1 2/a" ,   0,   "c1", "-A 2a"         }, //  95
  { 15,    0, "B 1 1 2/n" ,   0,   "c2", "-B 2ab"        }, //  96
  { 15,    0, "I 1 1 2/b" ,   0,   "c3", "-I 2b"         }, //  97
  { 15, 1015, "B 1 1 2/b" ,   0,  "-c1", "-B 2b"         }, //  98
  { 15,    0, "A 1 1 2/n" ,   0,  "-c2", "-A 2ab"        }, //  99
  { 15,    0, "I 1 1 2/a" ,   0,  "-c3", "-I 2a"         }, // 100
  { 15,    0, "B 2/b 1 1" ,   0,   "a1", "-B 2xb"        }, // 101
  { 15,    0, "C 2/n 1 1" ,   0,   "a2", "-C 2xac"       }, // 102
  { 15,    0, "I 2/c 1 1" ,   0,   "a3", "-I 2xc"        }, // 103
  { 15,    0, "C 2/c 1 1" ,   0,  "-a1", "-C 2xc"        }, // 104
  { 15,    0, "B 2/n 1 1" ,   0,  "-a2", "-B 2xab"       }, // 105
  { 15,    0, "I 2/b 1 1" ,   0,  "-a3", "-I 2xb"        }, // 106
  { 16,   16, "P 2 2 2"   ,   0,     "", "P 2 2"         }, // 107
  { 17,   17, "P 2 2 21"  ,   0,     "", "P 2c 2"        }, // 108
  { 17, 1017, "P 21 2 2"  ,   0,  "cab", "P 2a 2a"       }, // 109
  { 17, 2017, "P 2 21 2"  ,   0,  "bca", "P 2 2b"        }, // 110
  { 18,   18, "P 21 21 2" ,   0,     "", "P 2 2ab"       }, // 111
  { 18, 3018, "P 2 21 21" ,   0,  "cab", "P 2bc 2"       }, // 112
  { 18, 2018, "P 21 2 21" ,   0,  "bca", "P 2ac 2ac"     }, // 113
  { 19,   19, "P 21 21 21",   0,     "", "P 2ac 2ab"     }, // 114
  { 20,   20, "C 2 2 21"  ,   0,     "", "C 2c 2"        }, // 115
  { 20,    0, "A 21 2 2"  ,   0,  "cab", "A 2a 2a"       }, // 116
  { 20,    0, "B 2 21 2"  ,   0,  "bca", "B 2 2b"        }, // 117
  { 21,   21, "C 2 2 2"   ,   0,     "", "C 2 2"         }, // 118
  { 21,    0, "A 2 2 2"   ,   0,  "cab", "A 2 2"         }, // 119
  { 21,    0, "B 2 2 2"   ,   0,  "bca", "B 2 2"         }, // 120
  { 22,   22, "F 2 2 2"   ,   0,     "", "F 2 2"         }, // 121
  { 23,   23, "I 2 2 2"   ,   0,     "", "I 2 2"         }, // 122
  { 24,   24, "I 21 21 21",   0,     "", "I 2b 2c"       }, // 123
  { 25,   25, "P m m 2"   ,   0,     "", "P 2 -2"        }, // 124
  { 25,    0, "P 2 m m"   ,   0,  "cab", "P -2 2"        }, // 125
  { 25,    0, "P m 2 m"   ,   0,  "bca", "P -2 -2"       }, // 126
  { 26,   26, "P m c 21"  ,   0,     "", "P 2c -2"       }, // 127
  { 26,    0, "P c m 21"  ,   0, "ba-c", "P 2c -2c"      }, // 128
  { 26,    0, "P 21 m a"  ,   0,  "cab", "P -2a 2a"      }, // 129
  { 26,    0, "P 21 a m"  ,   0, "-cba", "P -2 2a"       }, // 130
  { 26,    0, "P b 21 m"  ,   0,  "bca", "P -2 -2b"      }, // 131
  { 26,    0, "P m 21 b"  ,   0, "a-cb", "P -2b -2"      }, // 132
  { 27,   27, "P c c 2"   ,   0,     "", "P 2 -2c"       }, // 133
  { 27,    0, "P 2 a a"   ,   0,  "cab", "P -2a 2"       }, // 134
  { 27,    0, "P b 2 b"   ,   0,  "bca", "P -2b -2b"     }, // 135
  { 28,   28, "P m a 2"   ,   0,     "", "P 2 -2a"       }, // 136
  { 28,    0, "P b m 2"   ,   0, "ba-c", "P 2 -2b"       }, // 137
  { 28,    0, "P 2 m b"   ,   0,  "cab", "P -2b 2"       }, // 138
  { 28,    0, "P 2 c m"   ,   0, "-cba", "P -2c 2"       }, // 139
  { 28,    0, "P c 2 m"   ,   0,  "bca", "P -2c -2c"     }, // 140
  { 28,    0, "P m 2 a"   ,   0, "a-cb", "P -2a -2a"     }, // 141
  { 29,   29, "P c a 21"  ,   0,     "", "P 2c -2ac"     }, // 142
  { 29,    0, "P b c 21"  ,   0, "ba-c", "P 2c -2b"      }, // 143
  { 29,    0, "P 21 a b"  ,   0,  "cab", "P -2b 2a"      }, // 144
  { 29,    0, "P 21 c a"  ,   0, "-cba", "P -2ac 2a"     }, // 145
  { 29,    0, "P c 21 b"  ,   0,  "bca", "P -2bc -2c"    }, // 146
  { 29,    0, "P b 21 a"  ,   0, "a-cb", "P -2a -2ab"    }, // 147
  { 30,   30, "P n c 2"   ,   0,     "", "P 2 -2bc"      }, // 148
  { 30,    0, "P c n 2"   ,   0, "ba-c", "P 2 -2ac"      }, // 149
  { 30,    0, "P 2 n a"   ,   0,  "cab", "P -2ac 2"      }, // 150
  { 30,    0, "P 2 a n"   ,   0, "-cba", "P -2ab 2"      }, // 151
  { 30,    0, "P b 2 n"   ,   0,  "bca", "P -2ab -2ab"   }, // 152
  { 30,    0, "P n 2 b"   ,   0, "a-cb", "P -2bc -2bc"   }, // 153
  { 31,   31, "P m n 21"  ,   0,     "", "P 2ac -2"      }, // 154
  { 31,    0, "P n m 21"  ,   0, "ba-c", "P 2bc -2bc"    }, // 155
  { 31,    0, "P 21 m n"  ,   0,  "cab", "P -2ab 2ab"    }, // 156
  { 31,    0, "P 21 n m"  ,   0, "-cba", "P -2 2ac"      }, // 157
  { 31,    0, "P n 21 m"  ,   0,  "bca", "P -2 -2bc"     }, // 158
  { 31,    0, "P m 21 n"  ,   0, "a-cb", "P -2ab -2"     }, // 159
  { 32,   32, "P b a 2"   ,   0,     "", "P 2 -2ab"      }, // 160
  { 32,    0, "P 2 c b"   ,   0,  "cab", "P -2bc 2"      }, // 161
  { 32,    0, "P c 2 a"   ,   0,  "bca", "P -2ac -2ac"   }, // 162
  { 33,   33, "P n a 21"  ,   0,     "", "P 2c -2n"      }, // 163
  { 33,    0, "P b n 21"  ,   0, "ba-c", "P 2c -2ab"     }, // 164
  { 33,    0, "P 21 n b"  ,   0,  "cab", "P -2bc 2a"     }, // 165
  { 33,    0, "P 21 c n"  ,   0, "-cba", "P -2n 2a"      }, // 166
  { 33,    0, "P c 21 n"  ,   0,  "bca", "P -2n -2ac"    }, // 167
  { 33,    0, "P n 21 a"  ,   0, "a-cb", "P -2ac -2n"    }, // 168
  { 34,   34, "P n n 2"   ,   0,     "", "P 2 -2n"       }, // 169
  { 34,    0, "P 2 n n"   ,   0,  "cab", "P -2n 2"       }, // 170
  { 34,    0, "P n 2 n"   ,   0,  "bca", "P -2n -2n"     }, // 171
  { 35,   35, "C m m 2"   ,   0,     "", "C 2 -2"        }, // 172
  { 35,    0, "A 2 m m"   ,   0,  "cab", "A -2 2"        }, // 173
  { 35,    0, "B m 2 m"   ,   0,  "bca", "B -2 -2"       }, // 174
  { 36,   36, "C m c 21"  ,   0,     "", "C 2c -2"       }, // 175
  { 36,    0, "C c m 21"  ,   0, "ba-c", "C 2c -2c"      }, // 176
  { 36,    0, "A 21 m a"  ,   0,  "cab", "A -2a 2a"      }, // 177
  { 36,    0, "A 21 a m"  ,   0, "-cba", "A -2 2a"       }, // 178
  { 36,    0, "B b 21 m"  ,   0,  "bca", "B -2 -2b"      }, // 179
  { 36,    0, "B m 21 b"  ,   0, "a-cb", "B -2b -2"      }, // 180
  { 37,   37, "C c c 2"   ,   0,     "", "C 2 -2c"       }, // 181
  { 37,    0, "A 2 a a"   ,   0,  "cab", "A -2a 2"       }, // 182
  { 37,    0, "B b 2 b"   ,   0,  "bca", "B -2b -2b"     }, // 183
  { 38,   38, "A m m 2"   ,   0,     "", "A 2 -2"        }, // 184
  { 38,    0, "B m m 2"   ,   0, "ba-c", "B 2 -2"        }, // 185
  { 38,    0, "B 2 m m"   ,   0,  "cab", "B -2 2"        }, // 186
  { 38,    0, "C 2 m m"   ,   0, "-cba", "C -2 2"        }, // 187
  { 38,    0, "C m 2 m"   ,   0,  "bca", "C -2 -2"       }, // 188
  { 38,    0, "A m 2 m"   ,   0, "a-cb", "A -2 -2"       }, // 189
  { 39,   39, "A b m 2"   ,   0,     "", "A 2 -2b"       }, // 190
  { 39,    0, "B m a 2"   ,   0, "ba-c", "B 2 -2a"       }, // 191
  { 39,    0, "B 2 c m"   ,   0,  "cab", "B -2a 2"       }, // 192
  { 39,    0, "C 2 m b"   ,   0, "-cba", "C -2a 2"       }, // 193
  { 39,    0, "C m 2 a"   ,   0,  "bca", "C -2a -2a"     }, // 194
  { 39,    0, "A c 2 m"   ,   0, "a-cb", "A -2b -2b"     }, // 195
  { 40,   40, "A m a 2"   ,   0,     "", "A 2 -2a"       }, // 196
  { 40,    0, "B b m 2"   ,   0, "ba-c", "B 2 -2b"       }, // 197
  { 40,    0, "B 2 m b"   ,   0,  "cab", "B -2b 2"       }, // 198
  { 40,    0, "C 2 c m"   ,   0, "-cba", "C -2c 2"       }, // 199
  { 40,    0, "C c 2 m"   ,   0,  "bca", "C -2c -2c"     }, // 200
  { 40,    0, "A m 2 a"   ,   0, "a-cb", "A -2a -2a"     }, // 201
  { 41,   41, "A b a 2"   ,   0,     "", "A 2 -2ab"      }, // 202
  { 41,    0, "B b a 2"   ,   0, "ba-c", "B 2 -2ab"      }, // 203
  { 41,    0, "B 2 c b"   ,   0,  "cab", "B -2ab 2"      }, // 204
  { 41,    0, "C 2 c b"   ,   0, "-cba", "C -2ac 2"      }, // 205
  { 41,    0, "C c 2 a"   ,   0,  "bca", "C -2ac -2ac"   }, // 206
  { 41,    0, "A c 2 a"   ,   0, "a-cb", "A -2ab -2ab"   }, // 207
  { 42,   42, "F m m 2"   ,   0,     "", "F 2 -2"        }, // 208
  { 42,    0, "F 2 m m"   ,   0,  "cab", "F -2 2"        }, // 209
  { 42,    0, "F m 2 m"   ,   0,  "bca", "F -2 -2"       }, // 210
  { 43,   43, "F d d 2"   ,   0,     "", "F 2 -2d"       }, // 211
  { 43,    0, "F 2 d d"   ,   0,  "cab", "F -2d 2"       }, // 212
  { 43,    0, "F d 2 d"   ,   0,  "bca", "F -2d -2d"     }, // 213
  { 44,   44, "I m m 2"   ,   0,     "", "I 2 -2"        }, // 214
  { 44,    0, "I 2 m m"   ,   0,  "cab", "I -2 2"        }, // 215
  { 44,    0, "I m 2 m"   ,   0,  "bca", "I -2 -2"       }, // 216
  { 45,   45, "I b a 2"   ,   0,     "", "I 2 -2c"       }, // 217
  { 45,    0, "I 2 c b"   ,   0,  "cab", "I -2a 2"       }, // 218
  { 45,    0, "I c 2 a"   ,   0,  "bca", "I -2b -2b"     }, // 219
  { 46,   46, "I m a 2"   ,   0,     "", "I 2 -2a"       }, // 220
  { 46,    0, "I b m 2"   ,   0, "ba-c", "I 2 -2b"       }, // 221
  { 46,    0, "I 2 m b"   ,   0,  "cab", "I -2b 2"       }, // 222
  { 46,    0, "I 2 c m"   ,   0, "-cba", "I -2c 2"       }, // 223
  { 46,    0, "I c 2 m"   ,   0,  "bca", "I -2c -2c"     }, // 224
  { 46,    0, "I m 2 a"   ,   0, "a-cb", "I -2a -2a"     }, // 225
  { 47,   47, "P m m m"   ,   0,     "", "-P 2 2"        }, // 226
  { 48,   48, "P n n n"   , '1',     "", "P 2 2 -1n"     }, // 227
  { 48,    0, "P n n n"   , '2',     "", "-P 2ab 2bc"    }, // 228
  { 49,   49, "P c c m"   ,   0,     "", "-P 2 2c"       }, // 229
  { 49,    0, "P m a a"   ,   0,  "cab", "-P 2a 2"       }, // 230
  { 49,    0, "P b m b"   ,   0,  "bca", "-P 2b 2b"      }, // 231
  { 50,   50, "P b a n"   , '1',     "", "P 2 2 -1ab"    }, // 232
  { 50,    0, "P b a n"   , '2',     "", "-P 2ab 2b"     }, // 233
  { 50,    0, "P n c b"   , '1',  "cab", "P 2 2 -1bc"    }, // 234
  { 50,    0, "P n c b"   , '2',  "cab", "-P 2b 2bc"     }, // 235
  { 50,    0, "P c n a"   , '1',  "bca", "P 2 2 -1ac"    }, // 236
  { 50,    0, "P c n a"   , '2',  "bca", "-P 2a 2c"      }, // 237
  { 51,   51, "P m m a"   ,   0,     "", "-P 2a 2a"      }, // 238
  { 51,    0, "P m m b"   ,   0, "ba-c", "-P 2b 2"       }, // 239
  { 51,    0, "P b m m"   ,   0,  "cab", "-P 2 2b"       }, // 240
  { 51,    0, "P c m m"   ,   0, "-cba", "-P 2c 2c"      }, // 241
  { 51,    0, "P m c m"   ,   0,  "bca", "-P 2c 2"       }, // 242
  { 51,    0, "P m a m"   ,   0, "a-cb", "-P 2 2a"       }, // 243
  { 52,   52, "P n n a"   ,   0,     "", "-P 2a 2bc"     }, // 244
  { 52,    0, "P n n b"   ,   0, "ba-c", "-P 2b 2n"      }, // 245
  { 52,    0, "P b n n"   ,   0,  "cab", "-P 2n 2b"      }, // 246
  { 52,    0, "P c n n"   ,   0, "-cba", "-P 2ab 2c"     }, // 247
  { 52,    0, "P n c n"   ,   0,  "bca", "-P 2ab 2n"     }, // 248
  { 52,    0, "P n a n"   ,   0, "a-cb", "-P 2n 2bc"     }, // 249
  { 53,   53, "P m n a"   ,   0,     "", "-P 2ac 2"      }, // 250
  { 53,    0, "P n m b"   ,   0, "ba-c", "-P 2bc 2bc"    }, // 251
  { 53,    0, "P b m n"   ,   0,  "cab", "-P 2ab 2ab"    }, // 252
  { 53,    0, "P c n m"   ,   0, "-cba", "-P 2 2ac"      }, // 253
  { 53,    0, "P n c m"   ,   0,  "bca", "-P 2 2bc"      }, // 254
  { 53,    0, "P m a n"   ,   0, "a-cb", "-P 2ab 2"      }, // 255
  { 54,   54, "P c c a"   ,   0,     "", "-P 2a 2ac"     }, // 256
  { 54,    0, "P c c b"   ,   0, "ba-c", "-P 2b 2c"      }, // 257
  { 54,    0, "P b a a"   ,   0,  "cab", "-P 2a 2b"      }, // 258
  { 54,    0, "P c a a"   ,   0, "-cba", "-P 2ac 2c"     }, // 259
  { 54,    0, "P b c b"   ,   0,  "bca", "-P 2bc 2b"     }, // 260
  { 54,    0, "P b a b"   ,   0, "a-cb", "-P 2b 2ab"     }, // 261
  { 55,   55, "P b a m"   ,   0,     "", "-P 2 2ab"      }, // 262
  { 55,    0, "P m c b"   ,   0,  "cab", "-P 2bc 2"      }, // 263
  { 55,    0, "P c m a"   ,   0,  "bca", "-P 2ac 2ac"    }, // 264
  { 56,   56, "P c c n"   ,   0,     "", "-P 2ab 2ac"    }, // 265
  { 56,    0, "P n a a"   ,   0,  "cab", "-P 2ac 2bc"    }, // 266
  { 56,    0, "P b n b"   ,   0,  "bca", "-P 2bc 2ab"    }, // 267
  { 57,   57, "P b c m"   ,   0,     "", "-P 2c 2b"      }, // 268
  { 57,    0, "P c a m"   ,   0, "ba-c", "-P 2c 2ac"     }, // 269
  { 57,    0, "P m c a"   ,   0,  "cab", "-P 2ac 2a"     }, // 270
  { 57,    0, "P m a b"   ,   0, "-cba", "-P 2b 2a"      }, // 271
  { 57,    0, "P b m a"   ,   0,  "bca", "-P 2a 2ab"     }, // 272
  { 57,    0, "P c m b"   ,   0, "a-cb", "-P 2bc 2c"     }, // 273
  { 58,   58, "P n n m"   ,   0,     "", "-P 2 2n"       }, // 274
  { 58,    0, "P m n n"   ,   0,  "cab", "-P 2n 2"       }, // 275
  { 58,    0, "P n m n"   ,   0,  "bca", "-P 2n 2n"      }, // 276
  { 59,   59, "P m m n"   , '1',     "", "P 2 2ab -1ab"  }, // 277
  { 59, 1059, "P m m n"   , '2',     "", "-P 2ab 2a"     }, // 278
  { 59,    0, "P n m m"   , '1',  "cab", "P 2bc 2 -1bc"  }, // 279
  { 59,    0, "P n m m"   , '2',  "cab", "-P 2c 2bc"     }, // 280
  { 59,    0, "P m n m"   , '1',  "bca", "P 2ac 2ac -1ac"}, // 281
  { 59,    0, "P m n m"   , '2',  "bca", "-P 2c 2a"      }, // 282
  { 60,   60, "P b c n"   ,   0,     "", "-P 2n 2ab"     }, // 283
  { 60,    0, "P c a n"   ,   0, "ba-c", "-P 2n 2c"      }, // 284
  { 60,    0, "P n c a"   ,   0,  "cab", "-P 2a 2n"      }, // 285
  { 60,    0, "P n a b"   ,   0, "-cba", "-P 2bc 2n"     }, // 286
  { 60,    0, "P b n a"   ,   0,  "bca", "-P 2ac 2b"     }, // 287
  { 60,    0, "P c n b"   ,   0, "a-cb", "-P 2b 2ac"     }, // 288
  { 61,   61, "P b c a"   ,   0,     "", "-P 2ac 2ab"    }, // 289
  { 61,    0, "P c a b"   ,   0, "ba-c", "-P 2bc 2ac"    }, // 290
  { 62,   62, "P n m a"   ,   0,     "", "-P 2ac 2n"     }, // 291
  { 62,    0, "P m n b"   ,   0, "ba-c", "-P 2bc 2a"     }, // 292
  { 62,    0, "P b n m"   ,   0,  "cab", "-P 2c 2ab"     }, // 293
  { 62,    0, "P c m n"   ,   0, "-cba", "-P 2n 2ac"     }, // 294
  { 62,    0, "P m c n"   ,   0,  "bca", "-P 2n 2a"      }, // 295
  { 62,    0, "P n a m"   ,   0, "a-cb", "-P 2c 2n"      }, // 296
  { 63,   63, "C m c m"   ,   0,     "", "-C 2c 2"       }, // 297
  { 63,    0, "C c m m"   ,   0, "ba-c", "-C 2c 2c"      }, // 298
  { 63,    0, "A m m a"   ,   0,  "cab", "-A 2a 2a"      }, // 299
  { 63,    0, "A m a m"   ,   0, "-cba", "-A 2 2a"       }, // 300
  { 63,    0, "B b m m"   ,   0,  "bca", "-B 2 2b"       }, // 301
  { 63,    0, "B m m b"   ,   0, "a-cb", "-B 2b 2"       }, // 302
  { 64,   64, "C m c a"   ,   0,     "", "-C 2ac 2"      }, // 303
  { 64,    0, "C c m b"   ,   0, "ba-c", "-C 2ac 2ac"    }, // 304
  { 64,    0, "A b m a"   ,   0,  "cab", "-A 2ab 2ab"    }, // 305
  { 64,    0, "A c a m"   ,   0, "-cba", "-A 2 2ab"      }, // 306
  { 64,    0, "B b c m"   ,   0,  "bca", "-B 2 2ab"      }, // 307
  { 64,    0, "B m a b"   ,   0, "a-cb", "-B 2ab 2"      }, // 308
  { 65,   65, "C m m m"   ,   0,     "", "-C 2 2"        }, // 309
  { 65,    0, "A m m m"   ,   0,  "cab", "-A 2 2"        }, // 310
  { 65,    0, "B m m m"   ,   0,  "bca", "-B 2 2"        }, // 311
  { 66,   66, "C c c m"   ,   0,     "", "-C 2 2c"       }, // 312
  { 66,    0, "A m a a"   ,   0,  "cab", "-A 2a 2"       }, // 313
  { 66,    0, "B b m b"   ,   0,  "bca", "-B 2b 2b"      }, // 314
  { 67,   67, "C m m a"   ,   0,     "", "-C 2a 2"       }, // 315
  { 67,    0, "C m m b"   ,   0, "ba-c", "-C 2a 2a"      }, // 316
  { 67,    0, "A b m m"   ,   0,  "cab", "-A 2b 2b"      }, // 317
  { 67,    0, "A c m m"   ,   0, "-cba", "-A 2 2b"       }, // 318
  { 67,    0, "B m c m"   ,   0,  "bca", "-B 2 2a"       }, // 319
  { 67,    0, "B m a m"   ,   0, "a-cb", "-B 2a 2"       }, // 320
  { 68,   68, "C c c a"   , '1',     "", "C 2 2 -1ac"    }, // 321
  { 68,    0, "C c c a"   , '2',     "", "-C 2a 2ac"     }, // 322
  { 68,    0, "C c c b"   , '1', "ba-c", "C 2 2 -1ac"    }, // 323
  { 68,    0, "C c c b"   , '2', "ba-c", "-C 2a 2c"      }, // 324
  { 68,    0, "A b a a"   , '1',  "cab", "A 2 2 -1ab"    }, // 325
  { 68,    0, "A b a a"   , '2',  "cab", "-A 2a 2b"      }, // 326
  { 68,    0, "A c a a"   , '1', "-cba", "A 2 2 -1ab"    }, // 327
  { 68,    0, "A c a a"   , '2', "-cba", "-A 2ab 2b"     }, // 328
  { 68,    0, "B b c b"   , '1',  "bca", "B 2 2 -1ab"    }, // 329
  { 68,    0, "B b c b"   , '2',  "bca", "-B 2ab 2b"     }, // 330
  { 68,    0, "B b a b"   , '1', "a-cb", "B 2 2 -1ab"    }, // 331
  { 68,    0, "B b a b"   , '2', "a-cb", "-B 2b 2ab"     }, // 332
  { 69,   69, "F m m m"   ,   0,     "", "-F 2 2"        }, // 333
  { 70,   70, "F d d d"   , '1',     "", "F 2 2 -1d"     }, // 334
  { 70,    0, "F d d d"   , '2',     "", "-F 2uv 2vw"    }, // 335
  { 71,   71, "I m m m"   ,   0,     "", "-I 2 2"        }, // 336
  { 72,   72, "I b a m"   ,   0,     "", "-I 2 2c"       }, // 337
  { 72,    0, "I m c b"   ,   0,  "cab", "-I 2a 2"       }, // 338
  { 72,    0, "I c m a"   ,   0,  "bca", "-I 2b 2b"      }, // 339
  { 73,   73, "I b c a"   ,   0,     "", "-I 2b 2c"      }, // 340
  { 73,    0, "I c a b"   ,   0, "ba-c", "-I 2a 2b"      }, // 341
  { 74,   74, "I m m a"   ,   0,     "", "-I 2b 2"       }, // 342
  { 74,    0, "I m m b"   ,   0, "ba-c", "-I 2a 2a"      }, // 343
  { 74,    0, "I b m m"   ,   0,  "cab", "-I 2c 2c"      }, // 344
  { 74,    0, "I c m m"   ,   0, "-cba", "-I 2 2b"       }, // 345
  { 74,    0, "I m c m"   ,   0,  "bca", "-I 2 2a"       }, // 346
  { 74,    0, "I m a m"   ,   0, "a-cb", "-I 2c 2"       }, // 347
  { 75,   75, "P 4"       ,   0,     "", "P 4"           }, // 348
  { 76,   76, "P 41"      ,   0,     "", "P 4w"          }, // 349
  { 77,   77, "P 42"      ,   0,     "", "P 4c"          }, // 350
  { 78,   78, "P 43"      ,   0,     "", "P 4cw"         }, // 351
  { 79,   79, "I 4"       ,   0,     "", "I 4"           }, // 352
  { 80,   80, "I 41"      ,   0,     "", "I 4bw"         }, // 353
  { 81,   81, "P -4"      ,   0,     "", "P -4"          }, // 354
  { 82,   82, "I -4"      ,   0,     "", "I -4"          }, // 355
  { 83,   83, "P 4/m"     ,   0,     "", "-P 4"          }, // 356
  { 84,   84, "P 42/m"    ,   0,     "", "-P 4c"         }, // 357
  { 85,   85, "P 4/n"     , '1',     "", "P 4ab -1ab"    }, // 358
  { 85,    0, "P 4/n"     , '2',     "", "-P 4a"         }, // 359
  { 86,   86, "P 42/n"    , '1',     "", "P 4n -1n"      }, // 360
  { 86,    0, "P 42/n"    , '2',     "", "-P 4bc"        }, // 361
  { 87,   87, "I 4/m"     ,   0,     "", "-I 4"          }, // 362
  { 88,   88, "I 41/a"    , '1',     "", "I 4bw -1bw"    }, // 363
  { 88,    0, "I 41/a"    , '2',     "", "-I 4ad"        }, // 364
  { 89,   89, "P 4 2 2"   ,   0,     "", "P 4 2"         }, // 365
  { 90,   90, "P 4 21 2"  ,   0,     "", "P 4ab 2ab"     }, // 366
  { 91,   91, "P 41 2 2"  ,   0,     "", "P 4w 2c"       }, // 367
  { 92,   92, "P 41 21 2" ,   0,     "", "P 4abw 2nw"    }, // 368
  { 93,   93, "P 42 2 2"  ,   0,     "", "P 4c 2"        }, // 369
  { 94,   94, "P 42 21 2" ,   0,     "", "P 4n 2n"       }, // 370
  { 95,   95, "P 43 2 2"  ,   0,     "", "P 4cw 2c"      }, // 371
  { 96,   96, "P 43 21 2" ,   0,     "", "P 4nw 2abw"    }, // 372
  { 97,   97, "I 4 2 2"   ,   0,     "", "I 4 2"         }, // 373
  { 98,   98, "I 41 2 2"  ,   0,     "", "I 4bw 2bw"     }, // 374
  { 99,   99, "P 4 m m"   ,   0,     "", "P 4 -2"        }, // 375
  {100,  100, "P 4 b m"   ,   0,     "", "P 4 -2ab"      }, // 376
  {101,  101, "P 42 c m"  ,   0,     "", "P 4c -2c"      }, // 377
  {102,  102, "P 42 n m"  ,   0,     "", "P 4n -2n"      }, // 378
  {103,  103, "P 4 c c"   ,   0,     "", "P 4 -2c"       }, // 379
  {104,  104, "P 4 n c"   ,   0,     "", "P 4 -2n"       }, // 380
  {105,  105, "P 42 m c"  ,   0,     "", "P 4c -2"       }, // 381
  {106,  106, "P 42 b c"  ,   0,     "", "P 4c -2ab"     }, // 382
  {107,  107, "I 4 m m"   ,   0,     "", "I 4 -2"        }, // 383
  {108,  108, "I 4 c m"   ,   0,     "", "I 4 -2c"       }, // 384
  {109,  109, "I 41 m d"  ,   0,     "", "I 4bw -2"      }, // 385
  {110,  110, "I 41 c d"  ,   0,     "", "I 4bw -2c"     }, // 386
  {111,  111, "P -4 2 m"  ,   0,     "", "P -4 2"        }, // 387
  {112,  112, "P -4 2 c"  ,   0,     "", "P -4 2c"       }, // 388
  {113,  113, "P -4 21 m" ,   0,     "", "P -4 2ab"      }, // 389
  {114,  114, "P -4 21 c" ,   0,     "", "P -4 2n"       }, // 390
  {115,  115, "P -4 m 2"  ,   0,     "", "P -4 -2"       }, // 391
  {116,  116, "P -4 c 2"  ,   0,     "", "P -4 -2c"      }, // 392
  {117,  117, "P -4 b 2"  ,   0,     "", "P -4 -2ab"     }, // 393
  {118,  118, "P -4 n 2"  ,   0,     "", "P -4 -2n"      }, // 394
  {119,  119, "I -4 m 2"  ,   0,     "", "I -4 -2"       }, // 395
  {120,  120, "I -4 c 2"  ,   0,     "", "I -4 -2c"      }, // 396
  {121,  121, "I -4 2 m"  ,   0,     "", "I -4 2"        }, // 397
  {122,  122, "I -4 2 d"  ,   0,     "", "I -4 2bw"      }, // 398
  {123,  123, "P 4/m m m" ,   0,     "", "-P 4 2"        }, // 399
  {124,  124, "P 4/m c c" ,   0,     "", "-P 4 2c"       }, // 400
  {125,  125, "P 4/n b m" , '1',     "", "P 4 2 -1ab"    }, // 401
  {125,    0, "P 4/n b m" , '2',     "", "-P 4a 2b"      }, // 402
  {126,  126, "P 4/n n c" , '1',     "", "P 4 2 -1n"     }, // 403
  {126,    0, "P 4/n n c" , '2',     "", "-P 4a 2bc"     }, // 404
  {127,  127, "P 4/m b m" ,   0,     "", "-P 4 2ab"      }, // 405
  {128,  128, "P 4/m n c" ,   0,     "", "-P 4 2n"       }, // 406
  {129,  129, "P 4/n m m" , '1',     "", "P 4ab 2ab -1ab"}, // 407
  {129,    0, "P 4/n m m" , '2',     "", "-P 4a 2a"      }, // 408
  {130,  130, "P 4/n c c" , '1',     "", "P 4ab 2n -1ab" }, // 409
  {130,    0, "P 4/n c c" , '2',     "", "-P 4a 2ac"     }, // 410
  {131,  131, "P 42/m m c",   0,     "", "-P 4c 2"       }, // 411
  {132,  132, "P 42/m c m",   0,     "", "-P 4c 2c"      }, // 412
  {133,  133, "P 42/n b c", '1',     "", "P 4n 2c -1n"   }, // 413
  {133,    0, "P 42/n b c", '2',     "", "-P 4ac 2b"     }, // 414
  {134,  134, "P 42/n n m", '1',     "", "P 4n 2 -1n"    }, // 415
  {134,    0, "P 42/n n m", '2',     "", "-P 4ac 2bc"    }, // 416
  {135,  135, "P 42/m b c",   0,     "", "-P 4c 2ab"     }, // 417
  {136,  136, "P 42/m n m",   0,     "", "-P 4n 2n"      }, // 418
  {137,  137, "P 42/n m c", '1',     "", "P 4n 2n -1n"   }, // 419
  {137,    0, "P 42/n m c", '2',     "", "-P 4ac 2a"     }, // 420
  {138,  138, "P 42/n c m", '1',     "", "P 4n 2ab -1n"  }, // 421
  {138,    0, "P 42/n c m", '2',     "", "-P 4ac 2ac"    }, // 422
  {139,  139, "I 4/m m m" ,   0,     "", "-I 4 2"        }, // 423
  {140,  140, "I 4/m c m" ,   0,     "", "-I 4 2c"       }, // 424
  {141,  141, "I 41/a m d", '1',     "", "I 4bw 2bw -1bw"}, // 425
  {141,    0, "I 41/a m d", '2',     "", "-I 4bd 2"      }, // 426
  {142,  142, "I 41/a c d", '1',     "", "I 4bw 2aw -1bw"}, // 427
  {142,    0, "I 41/a c d", '2',     "", "-I 4bd 2c"     }, // 428
  {143,  143, "P 3"       ,   0,     "", "P 3"           }, // 429
  {144,  144, "P 31"      ,   0,     "", "P 31"          }, // 430
  {145,  145, "P 32"      ,   0,     "", "P 32"          }, // 431
  {146,  146, "R 3"       , 'H',     "", "R 3"           }, // 432
  {146, 1146, "R 3"       , 'R',     "", "P 3*"          }, // 433
  {147,  147, "P -3"      ,   0,     "", "-P 3"          }, // 434
  {148,  148, "R -3"      , 'H',     "", "-R 3"          }, // 435
  {148, 1148, "R -3"      , 'R',     "", "-P 3*"         }, // 436
  {149,  149, "P 3 1 2"   ,   0,     "", "P 3 2"         }, // 437
  {150,  150, "P 3 2 1"   ,   0,     "", "P 3 2\""       }, // 438
  {151,  151, "P 31 1 2"  ,   0,     "", "P 31 2 (0 0 4)"}, // 439
  {152,  152, "P 31 2 1"  ,   0,     "", "P 31 2\""      }, // 440
  {153,  153, "P 32 1 2"  ,   0,     "", "P 32 2 (0 0 2)"}, // 441
  {154,  154, "P 32 2 1"  ,   0,     "", "P 32 2\""      }, // 442
  {155,  155, "R 3 2"     , 'H',     "", "R 3 2\""       }, // 443
  {155, 1155, "R 3 2"     , 'R',     "", "P 3* 2"        }, // 444
  {156,  156, "P 3 m 1"   ,   0,     "", "P 3 -2\""      }, // 445
  {157,  157, "P 3 1 m"   ,   0,     "", "P 3 -2"        }, // 446
  {158,  158, "P 3 c 1"   ,   0,     "", "P 3 -2\"c"     }, // 447
  {159,  159, "P 3 1 c"   ,   0,     "", "P 3 -2c"       }, // 448
  {160,  160, "R 3 m"     , 'H',     "", "R 3 -2\""      }, // 449
  {160, 1160, "R 3 m"     , 'R',     "", "P 3* -2"       }, // 450
  {161,  161, "R 3 c"     , 'H',     "", "R 3 -2\"c"     }, // 451
  {161, 1161, "R 3 c"     , 'R',     "", "P 3* -2n"      }, // 452
  {162,  162, "P -3 1 m"  ,   0,     "", "-P 3 2"        }, // 453
  {163,  163, "P -3 1 c"  ,   0,     "", "-P 3 2c"       }, // 454
  {164,  164, "P -3 m 1"  ,   0,     "", "-P 3 2\""      }, // 455
  {165,  165, "P -3 c 1"  ,   0,     "", "-P 3 2\"c"     }, // 456
  {166,  166, "R -3 m"    , 'H',     "", "-R 3 2\""      }, // 457
  {166, 1166, "R -3 m"    , 'R',     "", "-P 3* 2"       }, // 458
  {167,  167, "R -3 c"    , 'H',     "", "-R 3 2\"c"     }, // 459
  {167, 1167, "R -3 c"    , 'R',     "", "-P 3* 2n"      }, // 460
  {168,  168, "P 6"       ,   0,     "", "P 6"           }, // 461
  {169,  169, "P 61"      ,   0,     "", "P 61"          }, // 462
  {170,  170, "P 65"      ,   0,     "", "P 65"          }, // 463
  {171,  171, "P 62"      ,   0,     "", "P 62"          }, // 464
  {172,  172, "P 64"      ,   0,     "", "P 64"          }, // 465
  {173,  173, "P 63"      ,   0,     "", "P 6c"          }, // 466
  {174,  174, "P -6"      ,   0,     "", "P -6"          }, // 467
  {175,  175, "P 6/m"     ,   0,     "", "-P 6"          }, // 468
  {176,  176, "P 63/m"    ,   0,     "", "-P 6c"         }, // 469
  {177,  177, "P 6 2 2"   ,   0,     "", "P 6 2"         }, // 470
  {178,  178, "P 61 2 2"  ,   0,     "", "P 61 2 (0 0 5)"}, // 471
  {179,  179, "P 65 2 2"  ,   0,     "", "P 65 2 (0 0 1)"}, // 472
  {180,  180, "P 62 2 2"  ,   0,     "", "P 62 2 (0 0 4)"}, // 473
  {181,  181, "P 64 2 2"  ,   0,     "", "P 64 2 (0 0 2)"}, // 474
  {182,  182, "P 63 2 2"  ,   0,     "", "P 6c 2c"       }, // 475
  {183,  183, "P 6 m m"   ,   0,     "", "P 6 -2"        }, // 476
  {184,  184, "P 6 c c"   ,   0,     "", "P 6 -2c"       }, // 477
  {185,  185, "P 63 c m"  ,   0,     "", "P 6c -2"       }, // 478
  {186,  186, "P 63 m c"  ,   0,     "", "P 6c -2c"      }, // 479
  {187,  187, "P -6 m 2"  ,   0,     "", "P -6 2"        }, // 480
  {188,  188, "P -6 c 2"  ,   0,     "", "P -6c 2"       }, // 481
  {189,  189, "P -6 2 m"  ,   0,     "", "P -6 -2"       }, // 482
  {190,  190, "P -6 2 c"  ,   0,     "", "P -6c -2c"     }, // 483
  {191,  191, "P 6/m m m" ,   0,     "", "-P 6 2"        }, // 484
  {192,  192, "P 6/m c c" ,   0,     "", "-P 6 2c"       }, // 485
  {193,  193, "P 63/m c m",   0,     "", "-P 6c 2"       }, // 486
  {194,  194, "P 63/m m c",   0,     "", "-P 6c 2c"      }, // 487
  {195,  195, "P 2 3"     ,   0,     "", "P 2 2 3"       }, // 488
  {196,  196, "F 2 3"     ,   0,     "", "F 2 2 3"       }, // 489
  {197,  197, "I 2 3"     ,   0,     "", "I 2 2 3"       }, // 490
  {198,  198, "P 21 3"    ,   0,     "", "P 2ac 2ab 3"   }, // 491
  {199,  199, "I 21 3"    ,   0,     "", "I 2b 2c 3"     }, // 492
  {200,  200, "P m -3"    ,   0,     "", "-P 2 2 3"      }, // 493
  {201,  201, "P n -3"    , '1',     "", "P 2 2 3 -1n"   }, // 494
  {201,    0, "P n -3"    , '2',     "", "-P 2ab 2bc 3"  }, // 495
  {202,  202, "F m -3"    ,   0,     "", "-F 2 2 3"      }, // 496
  {203,  203, "F d -3"    , '1',     "", "F 2 2 3 -1d"   }, // 497
  {203,    0, "F d -3"    , '2',     "", "-F 2uv 2vw 3"  }, // 498
  {204,  204, "I m -3"    ,   0,     "", "-I 2 2 3"      }, // 499
  {205,  205, "P a -3"    ,   0,     "", "-P 2ac 2ab 3"  }, // 500
  {206,  206, "I a -3"    ,   0,     "", "-I 2b 2c 3"    }, // 501
  {207,  207, "P 4 3 2"   ,   0,     "", "P 4 2 3"       }, // 502
  {208,  208, "P 42 3 2"  ,   0,     "", "P 4n 2 3"      }, // 503
  {209,  209, "F 4 3 2"   ,   0,     "", "F 4 2 3"       }, // 504
  {210,  210, "F 41 3 2"  ,   0,     "", "F 4d 2 3"      }, // 505
  {211,  211, "I 4 3 2"   ,   0,     "", "I 4 2 3"       }, // 506
  {212,  212, "P 43 3 2"  ,   0,     "", "P 4acd 2ab 3"  }, // 507
  {213,  213, "P 41 3 2"  ,   0,     "", "P 4bd 2ab 3"   }, // 508
  {214,  214, "I 41 3 2"  ,   0,     "", "I 4bd 2c 3"    }, // 509
  {215,  215, "P -4 3 m"  ,   0,     "", "P -4 2 3"      }, // 510
  {216,  216, "F -4 3 m"  ,   0,     "", "F -4 2 3"      }, // 511
  {217,  217, "I -4 3 m"  ,   0,     "", "I -4 2 3"      }, // 512
  {218,  218, "P -4 3 n"  ,   0,     "", "P -4n 2 3"     }, // 513
  {219,  219, "F -4 3 c"  ,   0,     "", "F -4a 2 3"     }, // 514
  {220,  220, "I -4 3 d"  ,   0,     "", "I -4bd 2c 3"   }, // 515
  {221,  221, "P m -3 m"  ,   0,     "", "-P 4 2 3"      }, // 516
  {222,  222, "P n -3 n"  , '1',     "", "P 4 2 3 -1n"   }, // 517
  {222,    0, "P n -3 n"  , '2',     "", "-P 4a 2bc 3"   }, // 518
  {223,  223, "P m -3 n"  ,   0,     "", "-P 4n 2 3"     }, // 519
  {224,  224, "P n -3 m"  , '1',     "", "P 4n 2 3 -1n"  }, // 520
  {224,    0, "P n -3 m"  , '2',     "", "-P 4bc 2bc 3"  }, // 521
  {225,  225, "F m -3 m"  ,   0,     "", "-F 4 2 3"      }, // 522
  {226,  226, "F m -3 c"  ,   0,     "", "-F 4a 2 3"     }, // 523
  {227,  227, "F d -3 m"  , '1',     "", "F 4d 2 3 -1d"  }, // 524
  {227,    0, "F d -3 m"  , '2',     "", "-F 4vw 2vw 3"  }, // 525
  {228,  228, "F d -3 c"  , '1',     "", "F 4d 2 3 -1ad" }, // 526
  {228,    0, "F d -3 c"  , '2',     "", "-F 4ud 2vw 3"  }, // 527
  {229,  229, "I m -3 m"  ,   0,     "", "-I 4 2 3"      }, // 528
  {230,  230, "I a -3 d"  ,   0,     "", "-I 4bd 2c 3"   }, // 529
  // And extra entries from syminfo.lib
  {  5, 5005, "I 1 21 1"  ,   0,     "", "I 2yb"         }, // 530
  {  5, 3005, "C 1 21 1"  ,   0,     "", "C 2yb"         }, // 531
  { 18, 1018, "P 21212(a)",   0,     "", "P 2ab 2a"      }, // 532
  { 20, 1020, "C 2 2 21a)",   0,     "", "C 2ac 2"       }, // 533
  { 21, 1021, "C 2 2 2a"  ,   0,     "", "C 2ab 2b"      }, // 534
  { 22, 1022, "F 2 2 2a"  ,   0,     "", "F 2 2c"        }, // 535
  { 23, 1023, "I 2 2 2a"  ,   0,     "", "I 2ab 2bc"     }, // 536
  { 94, 1094, "P 42 21 2a",   0,     "", "P 4bc 2a"      }, // 537
  {197, 1197, "I 2 3a"    ,   0,     "", "I 2ab 2bc 3"   }, // 538
  // And extra entries from Open Babel, double checked with
  // Crystallographic Space Group Diagrams and Tables,
  // http://img.chem.ucl.ac.uk/sgp/
  // triclinic - enlarged unit cells
  {  1,    0, "A 1"       ,   0,     "", "A 1"           },
  {  1,    0, "B 1"       ,   0,     "", "B 1"           },
  {  1,    0, "C 1"       ,   0,     "", "C 1"           },
  {  1,    0, "F 1"       ,   0,     "", "F 1"           },
  {  1,    0, "I 1"       ,   0,     "", "I 1"           },
  {  2,    0, "A -1"      ,   0,     "", "-A 1"          },
  {  2,    0, "B -1"      ,   0,     "", "-B 1"          },
  {  2,    0, "C -1"      ,   0,     "", "-C 1"          },
  {  2,    0, "F -1"      ,   0,     "", "-F 1"          },
  {  2,    0, "I -1"      ,   0,     "", "-I 1"          },
  // monoclinic
  {  4,    0, "C 1 1 21"  ,   0,     "", "C 2c"          },
  { 12,    0, "F 1 2/m 1" ,   0,     "", "-F 2y"         },
  { 64,    0, "A b a m"   ,   0,     "", "-A 2 2ab"      },
  // tetragonal - enlarged C- and F-centred unit cells
  {117,    0, "C -4 2 b"  ,   0,     "", "C -4 2ya"      },
  {139,    0, "F 4/m m m" ,   0,     "", "-F 4 2"        },
};

template<class Dummy>
const SpaceGroupAltName Tables_<Dummy>::alt_names[27] = {
  // In 1990's ITfC vol.A changed some of the standard names, introducing
  // symbols 'e' and 'g'. sgtbx interprets these new symbols with
  // option ad_hoc_1992. spglib uses only the new symbols.
  {"A e m 2",   0, 190}, // A b m 2
  {"B m e 2",   0, 191}, // B m a 2
  {"B 2 e m",   0, 192}, // B 2 c m
  {"C 2 m e",   0, 193}, // C 2 m b
  {"C m 2 e",   0, 194}, // C m 2 a
  {"A e 2 m",   0, 195}, // A c 2 m
  {"A e a 2",   0, 202}, // A b a 2
  {"B b e 2",   0, 203}, // B b a 2
  {"B 2 e b",   0, 204}, // B 2 c b
  {"C 2 c e",   0, 205}, // C 2 c b
  {"C c 2 e",   0, 206}, // C c 2 a
  {"A e 2 a",   0, 207}, // A c 2 a
  {"C m c e",   0, 303}, // C m c a
  {"C c m e",   0, 304}, // C c m b
  {"A e m a",   0, 305}, // A b m a
  {"A e a m",   0, 306}, // A c a m
  {"B b e m",   0, 307}, // B b c m
  {"B m e b",   0, 308}, // B m a b
  {"C m m e",   0, 315}, // C m m a
  {"A e m m",   0, 317}, // A b m m
  {"B m e m",   0, 319}, // B m c m
  {"C c c e", '1', 321}, // C c c a
  {"C c c e", '2', 322}, // C c c a
  {"A e a a", '1', 325}, // A b a a
  {"A e a a", '2', 326}, // A b a a
  {"B b e b", '1', 329}, // B b c b
  {"B b e b", '2', 330}, // B b c b
};

} // namespace impl

using spacegroup_tables = impl::Tables_<void>;

inline const SpaceGroup* find_spacegroup_by_number(int ccp4) noexcept {
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.ccp4 == ccp4)
      return &sg;
  return nullptr;
}

inline const SpaceGroup& get_spacegroup_by_number(int ccp4) {
  const SpaceGroup* sg = find_spacegroup_by_number(ccp4);
  if (sg == nullptr)
    throw std::invalid_argument("Invalid space-group number: "
                                + std::to_string(ccp4));
  return *sg;
}

inline const SpaceGroup* find_spacegroup_by_name(std::string name) noexcept {
  if (name[0] == 'H')
    name[0] = 'R';
  const char* p = impl::skip_blank(name.c_str());
  if (*p >= '0' && *p <= '9') { // handle numbers
    char *endptr;
    long n = std::strtol(p, &endptr, 10);
    return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
  }
  char first = *p & ~0x20; // to uppercase
  if (first == '\0')
    return nullptr;
  p = impl::skip_blank(p+1);
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.hm[0] == first && sg.hm[2] == *p) {
      const char* a = impl::skip_blank(p + 1);
      const char* b = impl::skip_blank(sg.hm + 3);
      while (*a == *b && *b != '\0') {
        a = impl::skip_blank(a+1);
        b = impl::skip_blank(b+1);
      }
      if (*b == '\0' &&
          (*a == '\0' || (*a == ':' && *impl::skip_blank(a+1) == sg.ext)))
        return &sg;
    } else if (sg.hm[0] == first && sg.hm[2] == '1' && sg.hm[3] == ' ' &&
               sg.hm[4] != '1') {
      // check monoclinic short names
      const char* a = impl::skip_blank(p);
      const char* b = sg.hm + 4;
      while (*a == *b && *b != ' ') {
        a = impl::skip_blank(a+1);
        ++b;
      }
      if (*a == '\0' && *b == ' ')
        return &sg;
    }
  for (const SpaceGroupAltName& sg : spacegroup_tables::alt_names)
    if (sg.hm[0] == first && sg.hm[2] == *p) {
      const char* a = impl::skip_blank(p + 1);
      const char* b = impl::skip_blank(sg.hm + 3);
      while (*a == *b && *b != '\0') {
        a = impl::skip_blank(a+1);
        b = impl::skip_blank(b+1);
      }
      if (*b == '\0' &&
          (*a == '\0' || (*a == ':' && *impl::skip_blank(a+1) == sg.ext)))
        return &spacegroup_tables::main[sg.pos];
    }
  return nullptr;
}

inline const SpaceGroup& get_spacegroup_by_name(const std::string& name) {
  const SpaceGroup* sg = find_spacegroup_by_name(name);
  if (sg == nullptr)
    throw std::invalid_argument("Unknown space-group name: " + name);
  return *sg;
}

inline const SpaceGroup& get_spacegroup_p1() {
  return spacegroup_tables::main[0];
}

inline const SpaceGroup* find_spacegroup_by_ops(const GroupOps& gops) {
  char c = gops.find_centering();
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if ((c == sg.hall[0] || c == sg.hall[1]) &&
        gops.is_same_as(sg.operations()))
      return &sg;
  return nullptr;
}

} // namespace gemmi

namespace std {
template<> struct hash<gemmi::Op> {
  size_t operator()(const gemmi::Op& op) const {
    size_t h = 0;
    for (int i = 0; i != 3; ++i)
      for (int j = 0; j != 3; ++j)
        h = (h << 2) ^ (op.rot[i][j] + 1);
    for (int i = 0; i != 3; ++i)
      h = (h << 5) ^ op.tran[i];
    return h;
  }
};
} // namespace std

#ifdef __clang__
# pragma clang diagnostic pop  // ignored -Wmissing-braces
#endif

#endif
// vim:sw=2:ts=2:et
