// Copyright 2017-2019 Global Phasing Ltd.
//
// Crystallographic Symmetry. Space Groups. Coordinate Triplets.
//
// If this is all that you need from Gemmi you can just copy this file,
// src/symmetry.cpp fail.hpp and LICENSE.txt to your project.

#ifndef GEMMI_SYMMETRY_HPP_
#define GEMMI_SYMMETRY_HPP_

#include <cstdlib>    // for strtol, abs
#include <cstring>    // for memchr, strchr
#include <cmath>      // for fabs
#include <array>
#include <algorithm>  // for count, sort, remove
#include <functional> // for hash
#include <stdexcept>  // for runtime_error, invalid_argument
#include <string>
#include <tuple>      // for tie
#include <vector>

#include "fail.hpp"   // for fail, unreachable

// we use brace elision for Op::Rot = std:array<std::array<int,3>,3>
#ifdef __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wmissing-braces"
#endif

namespace gemmi {

// UTILS

namespace impl {

// copied a helper function from atox.hpp to keep it a two-header lib
inline const char* skip_blank(const char* p) {
  if (p)
    while (*p == ' ' || *p == '\t' || *p == '_') // '_' can be used as space
      ++p;
  return p;
}

} // namespace impl


// OP

// Op is a symmetry operation, or a change-of-basic transformation,
// or a different operation of similar kind.
// Both "rotation" matrix and translation vector are fractional, with DEN
// used as the denominator.
struct Op {
  static constexpr int DEN = 24;  // 24 to handle 1/8 in change-of-basis
  typedef std::array<std::array<int, 3>, 3> Rot;
  typedef std::array<int, 3> Tran;

  Rot rot;
  Tran tran;

  std::string triplet(char style='x') const;

  Op inverse() const;

  Op::Tran wrapped_tran() const {
    Op::Tran t = tran;
    for (int i = 0; i != 3; ++i) {
      if (t[i] >= DEN) // elements need to be in [0,DEN)
        t[i] %= DEN;
      else if (t[i] < 0)
        t[i] = ((t[i] + 1) % DEN) + DEN - 1;
    }
    return t;
  }

  // If the translation points outside of the unit cell, wrap it.
  Op& wrap() {
    tran = wrapped_tran();
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

  Rot transposed_rot() const {
    return { rot[0][0], rot[1][0], rot[2][0],
             rot[0][1], rot[1][1], rot[2][1],
             rot[0][2], rot[1][2], rot[2][2] };
  }

  // DEN^3 for rotation, -DEN^3 for rotoinversion
  int det_rot() const {
    return rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
         - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
         + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0]);
  }

  // Rotation-part type based on Table 1 in RWGK, Acta Cryst. A55, 383 (1999)
  int rot_type() const {
    int det = det_rot();
    int tr_den = rot[0][0] + rot[1][1] + rot[2][2];
    int tr = tr_den / DEN;
    const int table[] = {0, 0, 2, 3, 4, 6, 1};
    if (std::abs(det) == DEN * DEN * DEN && tr * DEN == tr_den && std::abs(tr) <= 3)
      return det > 0 ? table[3 + tr] : -table[3 - tr];
    return 0;
  }

  Op combine(const Op& b) const {
    Op r;
    for (int i = 0; i != 3; ++i) {
      r.tran[i] = tran[i] * Op::DEN;
      for (int j = 0; j != 3; ++j) {
        r.rot[i][j] = (rot[i][0] * b.rot[0][j] +
                       rot[i][1] * b.rot[1][j] +
                       rot[i][2] * b.rot[2][j]) / Op::DEN;
        r.tran[i] += rot[i][j] * b.tran[j];
      }
      r.tran[i] /= Op::DEN;
    }
    return r;
  }

  std::array<double, 3> apply_to_xyz(const std::array<double, 3>& xyz) const {
    std::array<double, 3> out;
    for (int i = 0; i != 3; ++i)
      out[i] = (rot[i][0] * xyz[0] + rot[i][1] * xyz[1] + rot[i][2] * xyz[2] +
                tran[i]) / Op::DEN;
    return out;
  }

  // Miller is defined in the same way in namespace gemmi in unitcell.hpp
  using Miller = std::array<int, 3>;

  Miller apply_to_hkl_without_division(const Miller& hkl) const {
    Miller r;
    for (int i = 0; i != 3; ++i)
      r[i] = (rot[0][i] * hkl[0] + rot[1][i] * hkl[1] + rot[2][i] * hkl[2]);
    return r;
  }
  static Miller divide_hkl_by_DEN(const Miller& hkl) {
    return {{ hkl[0] / DEN, hkl[1] / DEN, hkl[2] / DEN }};
  }
  Miller apply_to_hkl(const Miller& hkl) const {
    return divide_hkl_by_DEN(apply_to_hkl_without_division(hkl));
  }

  double phase_shift(const Miller& hkl) const {
    constexpr double mult = -2 * 3.1415926535897932384626433832795 / Op::DEN;
    return mult * (hkl[0] * tran[0] + hkl[1] * tran[1] + hkl[2] * tran[2]);
  }

  std::array<std::array<int, 4>, 4> int_seitz() const {
    std::array<std::array<int, 4>, 4> t;
    for (int i = 0; i < 3; ++i)
      t[i] = { rot[i][0], rot[i][1], rot[i][2], tran[i] };
    t[3] = { 0, 0, 0, 1 };
    return t;
  }

  std::array<std::array<double, 4>, 4> float_seitz() const {
    std::array<std::array<double, 4>, 4> t;
    double m = 1.0 / Op::DEN;
    for (int i = 0; i < 3; ++i)
      t[i] = { m * rot[i][0], m * rot[i][1], m * rot[i][2], m * tran[i] };
    t[3] = { 0., 0., 0., 1. };
    return t;
  }

  static constexpr Op identity() {
    return {{DEN,0,0, 0,DEN,0, 0,0,DEN}, {0,0,0}};
  }
  static constexpr Op::Rot inversion_rot() {
    return {-DEN,0,0, 0,-DEN,0, 0,0,-DEN};
  }
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

inline Op Op::inverse() const {
  int detr = det_rot();
  if (detr == 0)
    fail("cannot invert matrix: " + Op{rot, {0,0,0}}.triplet());
  int d2 = Op::DEN * Op::DEN;
  Op inv;
  inv.rot[0][0] = d2 * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2]) / detr;
  inv.rot[0][1] = d2 * (rot[0][2] * rot[2][1] - rot[0][1] * rot[2][2]) / detr;
  inv.rot[0][2] = d2 * (rot[0][1] * rot[1][2] - rot[0][2] * rot[1][1]) / detr;
  inv.rot[1][0] = d2 * (rot[1][2] * rot[2][0] - rot[1][0] * rot[2][2]) / detr;
  inv.rot[1][1] = d2 * (rot[0][0] * rot[2][2] - rot[0][2] * rot[2][0]) / detr;
  inv.rot[1][2] = d2 * (rot[1][0] * rot[0][2] - rot[0][0] * rot[1][2]) / detr;
  inv.rot[2][0] = d2 * (rot[1][0] * rot[2][1] - rot[2][0] * rot[1][1]) / detr;
  inv.rot[2][1] = d2 * (rot[2][0] * rot[0][1] - rot[0][0] * rot[2][1]) / detr;
  inv.rot[2][2] = d2 * (rot[0][0] * rot[1][1] - rot[1][0] * rot[0][1]) / detr;
  for (int i = 0; i != 3; ++i)
    inv.tran[i] = (-tran[0] * inv.rot[i][0]
                   -tran[1] * inv.rot[i][1]
                   -tran[2] * inv.rot[i][2]) / Op::DEN;
  return inv;
}

// inverse of Op::float_seitz()
inline Op seitz_to_op(const std::array<std::array<double,4>, 4>& t) {
  static_assert(Op::DEN == 24, "");
  auto check_round = [](double d) {
    double r = std::round(d * Op::DEN);
    if (std::fabs(r - d * Op::DEN) > 0.05)
      fail("all numbers in Seitz matrix must be equal Z/24");
    return static_cast<int>(r);
  };
  Op op;
  if (std::fabs(t[3][0]) + std::fabs(t[3][1]) + std::fabs(t[3][2]) +
      std::fabs(t[3][3] - 1) > 1e-3)
    fail("the last row in Seitz matrix must be [0 0 0 1]");
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      op.rot[i][j] = check_round(t[i][j]);
    op.tran[i] = check_round(t[i][3]);
  }
  return op;
}


// TRIPLET -> OP

inline int interpret_miller_character(char c, const std::string& s) {
  static const signed char values[] =
    //a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
    { 1, 2, 3, 0, 0, 0, 0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3 };
  size_t idx = size_t((c | 0x20) - 'a');  // "|0x20" = to lower
  if (idx >= sizeof(values) || values[idx] == 0)
    fail("unexpected character '", c, "' in: ", s);
  return values[idx] - 1;
}

// decimal_fract is useful only for non-crystallographic ops (such as x+0.12)
inline std::array<int, 4> parse_triplet_part(const std::string& s,
                                             double* decimal_fract=nullptr) {
  std::array<int, 4> r = { 0, 0, 0, 0 };
  int num = Op::DEN;
  const char* c = s.c_str();
  while (*(c = impl::skip_blank(c))) {
    if (*c == '+' || *c == '-') {
      num = (*c == '+' ? Op::DEN : -Op::DEN);
      c = impl::skip_blank(++c);
    }
    if (num == 0)
      fail("wrong or unsupported triplet format: " + s);
    int r_idx;
    int den = 1;
    double fract = 0;
    if ((*c >= '0' && *c <= '9') || *c == '.') {
      // syntax examples in this branch: "1", "-1/2", "+2*x", "1/2 * b"
      char* endptr;
      int n = std::strtol(c, &endptr, 10);
      // some COD CIFs have decimal fractions ("-x+0.25", ".5+Y", "1.25000-y")
      if (*endptr == '.') {
        // avoiding strtod() etc which is locale-dependent
        fract = n;
        for (double denom = 0.1; *++endptr >= '0' && *endptr <= '9'; denom *= 0.1)
          fract += int(*endptr - '0') * denom;
        double rounded = std::round(fract * num);
        if (!decimal_fract) {
          if (std::fabs(rounded - fract * num) > 0.05)
            fail("unexpected number in a symmetry triplet part: " + s);
          num = int(rounded);
        }
      } else {
        num *= n;
      }
      if (*endptr == '/')
        den = std::strtol(endptr + 1, &endptr, 10);
      if (*endptr == '*') {
        c = impl::skip_blank(endptr + 1);
        r_idx = interpret_miller_character(*c, s);
        ++c;
      } else {
        c = endptr;
        r_idx = 3;
      }
    } else {
      // syntax examples in this branch: "x", "+a", "-k/3"
      r_idx = interpret_miller_character(*c, s);
      c = impl::skip_blank(++c);
      if (*c == '/') {
        char* endptr;
        den = std::strtol(c + 1, &endptr, 10);
        c = endptr;
      }
    }
    if (den != 1) {
      if (den <= 0 || Op::DEN % den != 0 || fract != 0)
        fail("Wrong denominator " + std::to_string(den) + " in: " + s);
      num /= den;
    }
    r[r_idx] += num;
    if (decimal_fract)
      decimal_fract[r_idx] = num > 0 ? fract : -fract;
    num = 0;
  }
  if (num != 0)
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


// OP -> TRIPLET

namespace impl {

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

inline void append_sign_of(std::string& s, int n) {
  if (n < 0)
    s += '-';
  else if (!s.empty())
    s += '+';
}

// append w/DEN fraction reduced to the lowest terms
inline std::pair<int,int> get_op_fraction(int w) {
  // Op::DEN == 24 == 2 * 2 * 2 * 3
  int denom = 1;
  for (int i = 0; i != 3; ++i)
    if (w % 2 == 0)  // 2, 2, 2
      w /= 2;
    else
      denom *= 2;
  if (w % 3 == 0)    // 3
    w /= 3;
  else
    denom *= 3;
  return {w, denom};
}

inline void append_fraction(std::string& s, std::pair<int,int> frac) {
  append_small_number(s, frac.first);
  if (frac.second != 1) {
    s += '/';
    append_small_number(s, frac.second);
  }
}

} // namespace impl

inline std::string make_triplet_part(const std::array<int, 3>& xyz, int w,
                                     char style='x') {
  std::string s;
  const char* letters = "xyz hkl abc XYZ HKL ABC";
  switch(style | 0x20) {  // |0x20 converts to lower case
    case 'x': break;
    case 'h': letters += 4; break;
    case 'a': letters += 8; break;
    default: fail("unexpected triplet style: ", style);
  }
  if (!(style & 0x20))  // not lower
    letters += 12;
  for (int i = 0; i != 3; ++i)
    if (xyz[i] != 0) {
      impl::append_sign_of(s, xyz[i]);
      int a = std::abs(xyz[i]);
      if (a != Op::DEN) {
        std::pair<int,int> frac = impl::get_op_fraction(a);
        if (frac.first == 1) {  // e.g. "x/3"
          s += letters[i];
          s += '/';
          impl::append_small_number(s, frac.second);
        } else {  // e.g. "2/3*x"
          impl::append_fraction(s, frac);
          s += '*';
          s += letters[i];
        }
      } else {
        s += letters[i];
      }
    }
  if (w != 0) {
    impl::append_sign_of(s, w);
    std::pair<int,int> frac = impl::get_op_fraction(std::abs(w));
    impl::append_fraction(s, frac);
  }
  return s;
}

inline std::string Op::triplet(char style) const {
  return make_triplet_part(rot[0], tran[0], style) +
   "," + make_triplet_part(rot[1], tran[1], style) +
   "," + make_triplet_part(rot[2], tran[2], style);
}


// GROUPS OF OPERATIONS

// corresponds to Table A1.4.2.2 in ITfC vol.B (edition 2010)
inline std::vector<Op::Tran> centring_vectors(char centring_type) {
  constexpr int h = Op::DEN / 2;
  constexpr int t = Op::DEN / 3;
  constexpr int d = 2 * t;
  // note: find_centering() depends on the order of operations in vector
  switch (centring_type & ~0x20) {
    case 'P': return {{0, 0, 0}};
    case 'A': return {{0, 0, 0}, {0, h, h}};
    case 'B': return {{0, 0, 0}, {h, 0, h}};
    case 'C': return {{0, 0, 0}, {h, h, 0}};
    case 'I': return {{0, 0, 0}, {h, h, h}};
    case 'R': return {{0, 0, 0}, {d, t, t}, {t, d, d}};
    // hall_symbols.html has no H, ITfC 2010 has no S and T
    case 'H': return {{0, 0, 0}, {d, t, 0}, {t, d, 0}};
    case 'S': return {{0, 0, 0}, {t, t, d}, {d, t, d}};
    case 'T': return {{0, 0, 0}, {t, d, t}, {d, t, d}};
    case 'F': return {{0, 0, 0}, {0, h, h}, {h, 0, h}, {h, h, 0}};
    default: fail("not a centring type: ", centring_type);
  }
}


struct GroupOps {
  std::vector<Op> sym_ops;
  std::vector<Op::Tran> cen_ops;

  int order() const { return static_cast<int>(sym_ops.size()*cen_ops.size()); }

  void add_missing_elements();
  void add_missing_elements_part2(const std::vector<Op>& gen,
                                  size_t max_size, bool ignore_bad_gen);

  bool add_inversion() {
    size_t init_size = sym_ops.size();
    sym_ops.reserve(2 * init_size);
    for (const Op& op : sym_ops) {
      Op::Rot neg = op.negated_rot();
      if (find_by_rotation(neg)) {
        sym_ops.resize(init_size);
        return false;
      }
      sym_ops.push_back({neg, op.tran});
    }
    return true;
  }

  char find_centering() const {
    if (cen_ops.size() == 1 && cen_ops[0] == Op::Tran{0, 0, 0})
      return 'P';
    std::vector<Op::Tran> trans = cen_ops;
    std::sort(trans.begin(), trans.end());
    for (char c : {'A', 'B', 'C', 'I', 'F', 'R', 'H', 'S', 'T'}) {
      std::vector<Op::Tran> c_vectors = centring_vectors(c);
      if (c == 'R' || c == 'H') // these two are returned not sorted
        std::swap(c_vectors[1], c_vectors[2]);
      if (trans == c_vectors)
        return c;
    }
    return 0;
  }

  Op* find_by_rotation(const Op::Rot& r) {
    for (Op& op : sym_ops)
      if (op.rot == r)
        return &op;
    return nullptr;
  }

  const Op* find_by_rotation(const Op::Rot& r) const {
    return const_cast<GroupOps*>(this)->find_by_rotation(r);
  }

  bool is_centrosymmetric() const {
    return find_by_rotation(Op::inversion_rot()) != nullptr;
  }

  bool is_reflection_centric(const Op::Miller& hkl) const {
    Op::Miller mhkl = {{-Op::DEN * hkl[0], -Op::DEN * hkl[1], -Op::DEN * hkl[2]}};
    for (const Op& op : sym_ops)
      if (op.apply_to_hkl_without_division(hkl) == mhkl)
        return true;
    return false;
  }

  int epsilon_factor_without_centering(const Op::Miller& hkl) const {
    Op::Miller denh = {{Op::DEN * hkl[0], Op::DEN * hkl[1], Op::DEN * hkl[2]}};
    int epsilon = 0;
    for (const Op& op : sym_ops)
      if (op.apply_to_hkl_without_division(hkl) == denh)
        ++epsilon;
    return epsilon;
  }
  int epsilon_factor(const Op::Miller& hkl) const {
    return epsilon_factor_without_centering(hkl) * (int) cen_ops.size();
  }

  static bool has_phase_shift(const Op::Tran& c, const Op::Miller& hkl) {
    return (hkl[0] * c[0] + hkl[1] * c[1] + hkl[2] * c[2]) % Op::DEN != 0;
  }

  bool is_systematically_absent(const Op::Miller& hkl) const {
    for (auto i = cen_ops.begin() + 1; i != cen_ops.end(); ++i)
      if (has_phase_shift(*i, hkl))
        return true;
    Op::Miller denh = {{Op::DEN * hkl[0], Op::DEN * hkl[1], Op::DEN * hkl[2]}};
    for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
      if (op->apply_to_hkl_without_division(hkl) == denh) {
        for (const Op::Tran& c : cen_ops)
          if (has_phase_shift({{op->tran[0] + c[0],
                                op->tran[1] + c[1],
                                op->tran[2] + c[2]}}, hkl))
            return true;
      }
    return false;
  }

  void change_basis_impl(const Op& cob, const Op& inv) {
    if (sym_ops.empty() || cen_ops.empty())
      return;

    // Apply change-of-basis to sym_ops.
    // Ignore the first item in sym_ops -- it's identity.
    for (auto op = sym_ops.begin() + 1; op != sym_ops.end(); ++op)
      *op = cob.combine(*op).combine(inv).wrap();

    // The number of centering vectors may be different.
    // As an ad-hoc method (not proved to be robust) add lattice points
    // from a super-cell.
    int idet = inv.det_rot() / (Op::DEN * Op::DEN * Op::DEN);
    if (idet > 1) {
      std::vector<Op::Tran> new_cen_ops;
      new_cen_ops.reserve(cen_ops.size() * idet * idet * idet);
      for (int i = 0; i < idet; ++i)
        for (int j = 0; j < idet; ++j)
          for (int k = 0; k < idet; ++k)
            for (Op::Tran& cen : cen_ops)
              new_cen_ops.push_back({i * Op::DEN + cen[0],
                                     j * Op::DEN + cen[1],
                                     k * Op::DEN + cen[2]});
      cen_ops.swap(new_cen_ops);
    }

    // Apply change-of-basis to centering vectors
    Op cvec = Op::identity();
    for (auto tr = cen_ops.begin() + 1; tr != cen_ops.end(); ++tr) {
      cvec.tran = *tr;
      *tr = cob.combine(cvec).combine(inv).wrap().tran;
    }

    // Remove redundant centering vectors.
    for (int i = static_cast<int>(cen_ops.size()) - 1; i > 0; --i)
      for (int j = i - 1; j >= 0; --j)
        if (cen_ops[i] == cen_ops[j]) {
          cen_ops.erase(cen_ops.begin() + i);
          break;
        }
  }

  void change_basis_forward(const Op& cob) { change_basis_impl(cob, cob.inverse()); }
  void change_basis_backward(const Op& inv) { change_basis_impl(inv.inverse(), inv); }

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
    int n_cen = n / (int) sym_ops.size();
    int n_sym = n % (int) sym_ops.size();
    return sym_ops.at(n_sym).add_centering(cen_ops.at(n_cen));
  }

  bool is_same_as(const GroupOps& other) const {
    if (cen_ops.size() != other.cen_ops.size() ||
        sym_ops.size() != other.sym_ops.size())
      return false;
    return all_ops_sorted() == other.all_ops_sorted();
  }

  bool has_same_centring(const GroupOps& other) const {
    if (cen_ops.size() != other.cen_ops.size())
      return false;
    if (std::is_sorted(cen_ops.begin(), cen_ops.end()) &&
        std::is_sorted(other.cen_ops.begin(), other.cen_ops.end()))
      return cen_ops == other.cen_ops;
    std::vector<Op::Tran> v1 = cen_ops;
    std::vector<Op::Tran> v2 = other.cen_ops;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    return v1 == v2;
  }

  bool has_same_rotations(const GroupOps& other) const {
    if (sym_ops.size() != other.sym_ops.size())
      return false;
    auto sorted_rotations = [](const GroupOps& g) {
      std::vector<Op::Rot> r(g.sym_ops.size());
      for (size_t i = 0; i != r.size(); ++i)
        r[i] = g.sym_ops[i].rot;
      std::sort(r.begin(), r.end());
      return r;
    };
    return sorted_rotations(*this) == sorted_rotations(other);
  }

  // minimal multiplicity for real-space grid in each direction
  // examples: 1,2,1 for P21, 1,1,6 for P61
  std::array<int, 3> find_grid_factors() const {
    const int T = Op::DEN;
    int r[3] = {T, T, T};
    for (Op op : *this)
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

  // remove translation part of sym_ops
  GroupOps derive_symmorphic() const {
    GroupOps r(*this);
    for (Op& op : r.sym_ops)
      op.tran[0] = op.tran[1] = op.tran[2] = 0;
    return r;
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

  Iter begin() const { return {*this, 0, 0}; }
  Iter end() const { return {*this, 0, (int) cen_ops.size()}; }
};

inline void GroupOps::add_missing_elements() {
  // We always keep identity as sym_ops[0].
  if (sym_ops.empty() || sym_ops[0] != Op::identity())
    fail("oops");
  if (sym_ops.size() == 1)
    return;
  constexpr size_t max_size = 1024;
  // Below we assume that all centring vectors are already known (in cen_ops)
  // so when checking for a new element we compare only the 3x3 matrix.
  // Dimino's algorithm. https://physics.stackexchange.com/a/351400/95713
  std::vector<Op> gen(sym_ops.begin() + 1, sym_ops.end());
  sym_ops.resize(2);
  const Op::Rot idrot = Op::identity().rot;
  for (Op g = sym_ops[1] * sym_ops[1]; g.rot != idrot; g *= sym_ops[1]) {
    sym_ops.push_back(g);
    if (sym_ops.size() > max_size)
      fail("Too many elements in the group - bad generators");
  }
  // the rest is in separate function b/c it's reused in twin.hpp
  add_missing_elements_part2(gen, max_size, false);
}

inline void GroupOps::add_missing_elements_part2(const std::vector<Op>& gen,
                                                 size_t max_size, bool ignore_bad_gen) {
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
      if (sym_ops.size() > max_size) {
        if (!ignore_bad_gen)
          fail("Too many elements in the group - bad generators");
        // ignore this generator and continue with the next one
        sym_ops.resize(init_size);
        break;
      }
    }
  }
}

// Create GroupOps from Ops by separating centering vectors
inline GroupOps split_centering_vectors(const std::vector<Op>& ops) {
  const Op identity = Op::identity();
  GroupOps go;
  go.sym_ops.push_back(identity);
  for (const Op& op : ops)
    if (Op* old_op = go.find_by_rotation(op.rot)) {
      Op::Tran tran = op.wrapped_tran();
      if (op.rot == identity.rot)  // pure shift
        go.cen_ops.push_back(tran);
      if (tran == identity.tran)  // or rather |tran| < |old_op->tran| ?
        old_op->tran = op.tran;
    } else {
      go.sym_ops.push_back(op);
    }
  return go;
}

// INTERPRETING HALL SYMBOLS
// based on both ITfC vol.B ch.1.4 (2010)
// and http://cci.lbl.gov/sginfo/hall_symbols.html

// matrices for Nz from Table 3 and 4 from hall_symbols.html
inline Op::Rot hall_rotation_z(int N) {
  constexpr int d = Op::DEN;
  switch (N) {
    case 1: return {d,0,0,  0,d,0,  0,0,d};
    case 2: return {-d,0,0, 0,-d,0, 0,0,d};
    case 3: return {0,-d,0, d,-d,0, 0,0,d};
    case 4: return {0,-d,0, d,0,0,  0,0,d};
    case 6: return {d,-d,0, d,0,0,  0,0,d};
    case '\'': return {0,-d,0, -d,0,0, 0,0,-d};
    case '"':  return {0,d,0,   d,0,0, 0,0,-d};
    case '*':  return {0,0,d,   d,0,0, 0,d,0};
    default: fail("incorrect axis definition");
  }
}
inline Op::Tran hall_translation_from_symbol(char symbol) {
  constexpr int h = Op::DEN / 2;
  constexpr int q = Op::DEN / 4;
  switch (symbol) {
    case 'a': return {h, 0, 0};
    case 'b': return {0, h, 0};
    case 'c': return {0, 0, h};
    case 'n': return {h, h, h};
    case 'u': return {q, 0, 0};
    case 'v': return {0, q, 0};
    case 'w': return {0, 0, q};
    case 'd': return {q, q, q};
    default: fail(std::string("unknown symbol: ") + symbol);
  }
}

inline Op hall_matrix_symbol(const char* start, const char* end,
                             int pos, int& prev) {
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
      fail("missing axis");
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
    op.tran[principal_axis - 'x'] += Op::DEN / N * fractional_tran;
  prev = N;
  return op;
}

// Parses either short (0 0 1) or long notation (x,y,z+1/12)
// but without multipliers (such as 1/2x) to keep things simple for now.
inline Op parse_hall_change_of_basis(const char* start, const char* end) {
  if (std::memchr(start, ',', end - start) != nullptr) // long symbol
    return parse_triplet(std::string(start, end));
  // short symbol (0 0 1)
  Op cob = Op::identity();
  char* endptr;
  for (int i = 0; i != 3; ++i) {
    cob.tran[i] = std::strtol(start, &endptr, 10) % 12 * (Op::DEN / 12);
    start = endptr;
  }
  if (endptr != end)
    fail("unexpected change-of-basis format: " + std::string(start, end));
  return cob;
}

inline GroupOps generators_from_hall(const char* hall) {
  auto find_blank = [](const char* p) {
    while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '_') // '_' == ' '
      ++p;
    return p;
  };
  if (hall == nullptr)
    fail("null");
  hall = impl::skip_blank(hall);
  GroupOps ops;
  ops.sym_ops.emplace_back(Op::identity());
  bool centrosym = (hall[0] == '-');
  const char* lat = impl::skip_blank(centrosym ? hall + 1 : hall);
  if (!lat)
    fail("not a hall symbol: " + std::string(hall));
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
  if (centrosym)
    ops.sym_ops.push_back({Op::identity().negated_rot(), {0,0,0}});
  if (*part == '(') {
    const char* rb = std::strchr(part, ')');
    if (!rb)
      fail("missing ')': " + std::string(hall));
    if (ops.sym_ops.empty())
      fail("misplaced translation: " + std::string(hall));
    ops.change_basis_forward(parse_hall_change_of_basis(part + 1, rb));

    if (*impl::skip_blank(find_blank(rb + 1)) != '\0')
      fail("unexpected characters after ')': " + std::string(hall));
  }
  return ops;
}

inline GroupOps symops_from_hall(const char* hall) {
  GroupOps ops = generators_from_hall(hall);
  ops.add_missing_elements();
  return ops;
}

// CRYSTAL SYSTEMS, POINT GROUPS AND LAUE CLASSES

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

// http://reference.iucr.org/dictionary/Laue_class
enum class Laue : unsigned char {
  L1=0, L2m, Lmmm, L4m, L4mmm, L3, L3m, L6m, L6mmm, Lm3, Lm3m
};

inline Laue pointgroup_to_laue(PointGroup pg) {
  static const Laue laue[32] = {
    Laue::L1, Laue::L1,
    Laue::L2m, Laue::L2m, Laue::L2m,
    Laue::Lmmm, Laue::Lmmm, Laue::Lmmm,
    Laue::L4m, Laue::L4m, Laue::L4m,
    Laue::L4mmm, Laue::L4mmm, Laue::L4mmm, Laue::L4mmm,
    Laue::L3, Laue::L3,
    Laue::L3m, Laue::L3m, Laue::L3m,
    Laue::L6m, Laue::L6m, Laue::L6m,
    Laue::L6mmm, Laue::L6mmm, Laue::L6mmm, Laue::L6mmm,
    Laue::Lm3, Laue::Lm3,
    Laue::Lm3m, Laue::Lm3m, Laue::Lm3m,
  };
  return laue[static_cast<int>(pg)];
}

// return centrosymmetric pointgroup from the Laue class
inline PointGroup laue_to_pointgroup(Laue laue) {
  static const PointGroup pg[11] = {
    PointGroup::Ci, PointGroup::C2h, PointGroup::D2h, PointGroup::C4h,
    PointGroup::D4h, PointGroup::C3i, PointGroup::D3d, PointGroup::C6h,
    PointGroup::D6h, PointGroup::Th, PointGroup::Oh
  };
  return pg[static_cast<int>(laue)];
}

inline const char* laue_class_str(Laue laue) {
  return point_group_hm(laue_to_pointgroup(laue));
}

inline CrystalSystem crystal_system(Laue laue) {
  static const CrystalSystem crystal_systems[11] = {
    CrystalSystem::Triclinic,
    CrystalSystem::Monoclinic,
    CrystalSystem::Orthorhombic,
    CrystalSystem::Tetragonal, CrystalSystem::Tetragonal,
    CrystalSystem::Trigonal,   CrystalSystem::Trigonal,
    CrystalSystem::Hexagonal,  CrystalSystem::Hexagonal,
    CrystalSystem::Cubic,      CrystalSystem::Cubic
  };
  return crystal_systems[static_cast<int>(laue)];
}

inline CrystalSystem crystal_system(PointGroup pg) {
  return crystal_system(pointgroup_to_laue(pg));
}

inline unsigned char point_group_index_and_category(int space_group_number) {
  // 0x20=Sohncke, 0x40=enantiomorphic, 0x80=symmorphic
  enum : unsigned char { S=0x20, E=(0x20|0x40), Y=0x80, Z=(0x20|0x80) };
  static const unsigned char indices[230] = {
     0|Z,  1|Y,  2|Z,  2|S,  2|Z,  3|Y,  3,    3|Y,  3,    4|Y,  // 1-10
     4,    4|Y,  4,    4,    4,    5|Z,  5|S,  5|S,  5|S,  5|S,  // 11-20
     5|Z,  5|Z,  5|Z,  5|S,  6|Y,  6,    6,    6,    6,    6,    // 21-30
     6,    6,    6,    6,    6|Y,  6,    6,    6|Y,  6,    6,    // 31-40
     6,    6|Y,  6,    6|Y,  6,    6,    7|Y,  7,    7,    7,    // 41-50
     7,    7,    7,    7,    7,    7,    7,    7,    7,    7,    // 51-60
     7,    7,    7,    7,    7|Y,  7,    7,    7,    7|Y,  7,    // 61-70
     7|Y,  7,    7,    7,    8|Z,  8|E,  8|S,  8|E,  8|Z,  8|S,  // 71-80
     9|Y,  9|Y, 10|Y, 10,   10,   10,   10|Y, 10,   11|Z, 11|S,  // 81-90
    11|E, 11|E, 11|S, 11|S, 11|E, 11|E, 11|Z, 11|S, 12|Y, 12,    // 91-100
    12,   12,   12,   12,   12,   12,   12|Y, 12,   12,   12,    // 101-110
    13|Y, 13,   13,   13,   13|Y, 13,   13,   13,   13|Y, 13,    // 111-120
    13|Y, 13,   14|Y, 14,   14,   14,   14,   14,   14,   14,    // 121-130
    14,   14,   14,   14,   14,   14,   14,   14,   14|Y, 14,    // 131-140
    14,   14,   15|Z, 15|E, 15|E, 15|Z, 16|Y, 16|Y, 17|Z, 17|Z,  // 141-150
    17|E, 17|E, 17|E, 17|E, 17|Z, 18|Y, 18|Y, 18,   18,   18|Y,  // 151-160
    18,   19|Y, 19,   19|Y, 19,   19|Y, 19,   20|Z, 20|E, 20|E,  // 161-170
    20|E, 20|E, 20|S, 21|Y, 22|Y, 22,   23|Z, 23|E, 23|E, 23|E,  // 171-180
    23|E, 23|S, 24|Y, 24,   24,   24,   25|Y, 25,   25|Y, 25,    // 181-190
    26|Y, 26,   26,   26,   27|Z, 27|Z, 27|Z, 27|S, 27|S, 28|Y,  // 191-200
    28,   28|Y, 28,   28|Y, 28,   28,   29|Z, 29|S, 29|Z, 29|S,  // 201-210
    29|Z, 29|E, 29|E, 29|S, 30|Y, 30|Y, 30|Y, 30,   30,   30,    // 211-220
    31|Y, 31,   31,   31,   31|Y, 31,   31,   31,   31|Y, 31     // 221-230
  };
  return indices[space_group_number-1];
}

inline PointGroup point_group(int space_group_number) {
  auto n = point_group_index_and_category(space_group_number);
  return static_cast<PointGroup>(n & 0x1f);
}

// true for 65 Sohncke (non-enantiogenic) space groups
inline bool is_sohncke(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x20) != 0;
}

// true for 22 space groups (11 enantiomorphic pairs)
inline bool is_enantiomorphic(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x40) != 0;
}

// true for 73 space groups
inline bool is_symmorphic(int space_group_number) {
  return (point_group_index_and_category(space_group_number) & 0x80) != 0;
}

/// Inversion center of the Euclidean normalizer that is not at the origin of
/// reference settings. Returns (0,0,0) if absent. Based on tables in ch. 3.5
/// of ITA (2016) doi:10.1107/97809553602060000933 (column "Inversion through
/// a centre at").
inline Op::Tran nonzero_inversion_center(int space_group_number) {
  constexpr int D = Op::DEN;
  switch (space_group_number) {
    case 43: return {D/8, D/8, 0};
    case 80: return {D/4, 0, 0};
    case 98: return {D/4, 0, D/8};
    case 109: return {D/4, 0, 0};
    case 110: return {D/4, 0, 0};
    case 122: return {D/4, 0, D/8};
    case 210: return {D/8, D/8, D/8};
    default: return {0, 0, 0};
  }
}

GEMMI_DLL const char* get_basisop(int basisop_idx);


// Returns a change-of-basis operator for centred -> primitive transformation.
// The same operator as inverse of z2p_op in sgtbx.
inline Op::Rot centred_to_primitive(char centring_type) {
  constexpr int D = Op::DEN;
  constexpr int H = Op::DEN / 2;
  constexpr int T = Op::DEN / 3;
  switch (centring_type) {
    case 'P': return {D,0,0, 0,D,0, 0,0,D};
    case 'A': return {-D,0,0, 0,-H,H, 0,H,H};
    case 'B': return {-H,0,H, 0,-D,0, H,0,H};
    case 'C': return {H,H,0, H,-H,0, 0,0,-D};
    case 'I': return {-H,H,H, H,-H,H, H,H,-H};
    case 'R': return {2*T,-T,-T, T,T,-2*T, T,T,T};
    case 'H': return {2*T,-T,0, T,T,0, 0,0,D};  // not used normally
    case 'F': return {0,H,H, H,0,H, H,H,0};
    default: fail("not a centring type: ", centring_type);
  }
}


// LIST OF CRYSTALLOGRAPHIC SPACE GROUPS

struct SpaceGroup { // typically 44 bytes
  int number;
  int ccp4;
  char hm[11];  // Hermann-Mauguin (international) notation
  char ext;
  char qualifier[5];
  char hall[15];
  int basisop_idx;

  std::string xhm() const {
    std::string ret = hm;
    if (ext) {
      ret += ':';
      ret += ext;
    }
    return ret;
  }

  char centring_type() const { return ext == 'R' ? 'P' : hm[0]; }

  // (old) CCP4 spacegroup names start with H for hexagonal setting
  char ccp4_lattice_type() const { return ext == 'H' ? 'H' : hm[0]; }

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

  // As explained in Phenix newsletter CCN_2011_01.pdf#page=12
  // the PDB uses own, non-standard symbols for rhombohedral space groups.
  std::string pdb_name() const {
    std::string s;
    s += ccp4_lattice_type();
    s += hm+1;
    return s;
  }

  bool is_sohncke() const { return gemmi::is_sohncke(number); }
  bool is_enantiomorphic() const { return gemmi::is_enantiomorphic(number); }
  bool is_symmorphic() const { return gemmi::is_symmorphic(number); }
  PointGroup point_group() const { return gemmi::point_group(number); }
  const char* point_group_hm() const {
    return gemmi::point_group_hm(point_group());
  }
  Laue laue_class() const { return pointgroup_to_laue(point_group()); }
  const char* laue_str() const { return laue_class_str(laue_class()); }
  CrystalSystem crystal_system() const {
    return gemmi::crystal_system(point_group());
  }
  const char* crystal_system_str() const {
    return gemmi::crystal_system_str(crystal_system());
  }
  bool is_centrosymmetric() const {
    return laue_to_pointgroup(laue_class()) == point_group();
  }

  /// returns 'a', 'b' or 'c' for monoclinic SG, '\0' otherwise
  char monoclinic_unique_axis() const {
    if (crystal_system() == CrystalSystem::Monoclinic)
      return qualifier[qualifier[0] == '-' ? 1 : 0];
    return '\0';
  }

  const char* basisop_str() const { return get_basisop(basisop_idx); }
  Op basisop() const { return parse_triplet(basisop_str()); }
  bool is_reference_setting() const { return basisop_idx == 0; }

  Op centred_to_primitive() const {
    return {gemmi::centred_to_primitive(centring_type()), {0,0,0}};
  }

  /// Returns change-of-hand operator. Compatible with similar sgtbx function.
  Op change_of_hand_op() const {
    if (is_centrosymmetric())
      return Op::identity();
    Op::Tran t = nonzero_inversion_center(number);
    Op op{Op::inversion_rot(), {2*t[0], 2*t[1], 2*t[2]}};
    if (!is_reference_setting()) {
      Op b = basisop();
      op = b.combine(op).combine(b.inverse());
    }
    return op;
  }

  GroupOps operations() const { return symops_from_hall(hall); }
};

struct SpaceGroupAltName {
  char hm[11];
  char ext;
  int pos;
};

struct GEMMI_DLL spacegroup_tables {
  static const SpaceGroup main[564];
  static const SpaceGroupAltName alt_names[28];
  static const unsigned char ccp4_hkl_asu[230];
};

inline const SpaceGroup* find_spacegroup_by_number(int ccp4) noexcept {
  if (ccp4 == 0)
    return &spacegroup_tables::main[0];
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

inline const SpaceGroup& get_spacegroup_reference_setting(int number) {
  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.number == number && sg.is_reference_setting())
      return sg;
  throw std::invalid_argument("Invalid space-group number: "
                              + std::to_string(number));
}

/// If angles alpha and gamma are provided, they are used to
/// distinguish hexagonal and rhombohedral settings (e.g. for "R 3").
/// \param prefer can specify preferred H/R settings and 1/2 origin choice.
/// For example, prefer="2H" means the origin choice 2 and hexagonal
/// settings. The default is "1H".
inline const SpaceGroup* find_spacegroup_by_name(std::string name,
                                  double alpha=0., double gamma=0.,
                                  const char* prefer=nullptr) {
  bool prefer_2 = false;
  bool prefer_R = false;
  if (prefer)
    for (const char* p = prefer; *p != '\0'; ++p) {
      if (*p == '2')
        prefer_2 = true;
      else if (*p == 'R')
        prefer_R = true;
      else if (*p != '1' && *p != 'H')
        throw std::invalid_argument("find_spacegroup_by_name(): invalid arg 'prefer'");
    }
  const char* p = impl::skip_blank(name.c_str());
  if (*p >= '0' && *p <= '9') { // handle numbers
    char *endptr;
    long n = std::strtol(p, &endptr, 10);
    return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
  }
  char first = *p & ~0x20; // to uppercase
  if (first == '\0')
    return nullptr;
  if (first == 'H')
    first = 'R';
  p = impl::skip_blank(p+1);
  size_t start = p - name.c_str();
  // change letters to lower case, except the letter after :
  for (size_t i = start; i < name.size(); ++i) {
    if (name[i] >= 'A' && name[i] <= 'Z')
      name[i] |= 0x20;  // to lowercase
    else if (name[i] == ':')
      while (++i < name.size())
        if (name[i] >= 'a' && name[i] <= 'z')
          name[i] &= ~0x20;  // to uppercase
  }
  // allow names ending with R or H, such as R3R instead of R3:R
  if (name.back() == 'h' || name.back() == 'r') {
    name.back() &= ~0x20;  // to uppercase
    name.insert(name.end() - 1, ':');
  }
  // The string that const char* p points to was just modified.
  // This confuses some compilers (GCC 4.8), so let's re-assign p.
  p = name.c_str() + start;

  for (const SpaceGroup& sg : spacegroup_tables::main)
    if (sg.hm[0] == first) {
      if (sg.hm[2] == *p) {
        const char* a = impl::skip_blank(p + 1);
        const char* b = impl::skip_blank(sg.hm + 3);
        // In IT 1935 and 1952, symbols of centrosymmetric, cubic space groups
        // 200-206 and 221-230 had symbol 3 (not -3), e.g. Pm3 instead of Pm-3,
        // as listed in Table 3.3.3.1 in ITfC (2016) vol. A, p.788.
        while ((*a == *b && *b != '\0') ||
               (*a == '3' && *b == '-' && b == sg.hm + 4 && *++b == '3')) {
          a = impl::skip_blank(a+1);
          b = impl::skip_blank(b+1);
        }
        if (*b == '\0') {
          if (*a == '\0') {
            // Change hexagonal settings to rhombohedral if the unit cell
            // angles are more consistent with the latter.
            // We have possible ambiguity in the hexagonal crystal family.
            // For instance, "R 3" may mean "R 3:H" (hexagonal setting) or
            // "R 3:R" (rhombohedral setting). The :H symbols come first
            // in the table and are used by default. The ratio gamma:alpha
            // is 120:90 in the hexagonal system and 1:1 in rhombohedral.
            // We assume that the 'R' entry follows directly the 'H' entry.
            if (sg.ext == 'H' && (alpha == 0. ? prefer_R : gamma < 1.125 * alpha))
              return &sg + 1;
            // Similarly, the origin choice #2 follows directly #1.
            if (sg.ext == '1' && prefer_2)
              return &sg + 1;
            return &sg;
          }
          if (*a == ':' && *impl::skip_blank(a+1) == sg.ext)
            return &sg;
        }
      } else if (sg.hm[2] == '1' && sg.hm[3] == ' ') {
        // check monoclinic short names, matching P2 to "P 1 2 1";
        // as an exception "B 2" == "B 1 1 2" (like in the PDB)
        const char* b = sg.hm + 4;
        if (*b != '1' || (first == 'B' && *++b == ' ' && *++b != '1')) {
          char end = (b == sg.hm + 4 ? ' ' : '\0');
          const char* a = impl::skip_blank(p);
          while (*a == *b && *b != end) {
            ++a;
            ++b;
          }
          if (*impl::skip_blank(a) == '\0' && *b == end)
            return &sg;
        }
      }
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

// Reciprocal space asu (asymmetric unit).
// The same 12 choices of ASU as in CCP4 symlib and cctbx.
struct ReciprocalAsu {
  int idx;
  Op::Rot rot{};  // value-initialized only to avoid -Wmaybe-uninitialized
  bool is_ref;

  ReciprocalAsu(const SpaceGroup* sg, bool tnt=false) {
    if (sg == nullptr)
      fail("Missing space group");
    idx = spacegroup_tables::ccp4_hkl_asu[sg->number - 1];
    if (tnt) {
      idx += 10;
      is_ref = true; // TNT ASU is given wrt current (not standard) settings
    } else {
      is_ref = sg->is_reference_setting();
      if (!is_ref)
        rot = sg->basisop().rot;
    }
  }

  bool is_in(const Op::Miller& hkl) const {
    if (is_ref)
      return is_in_reference_setting(hkl[0], hkl[1], hkl[2]);
    Op::Miller r;
    for (int i = 0; i != 3; ++i)
      r[i] = rot[0][i] * hkl[0] + rot[1][i] * hkl[1] + rot[2][i] * hkl[2];
    return is_in_reference_setting(r[0], r[1], r[2]);
  }

  bool is_in_reference_setting(int h, int k, int l) const {
    switch (idx) {
      // 0-9: CCP4 hkl asu,  10-19: TNT hkl asu
      case 0: return l>0 || (l==0 && (h>0 || (h==0 && k>=0)));
      case 1: return k>=0 && (l>0 || (l==0 && h>=0));
      case 12: // orthorhombic-D
      case 2: return h>=0 && k>=0 && l>=0;
      case 3: return l>=0 && ((h>=0 && k>0) || (h==0 && k==0));
      case 14: // tetragonal-D, hexagonal-D
      case 4: return h>=k && k>=0 && l>=0;
      case 5: return (h>=0 && k>0) || (h==0 && k==0 && l>=0);
      case 16: // trigonal-D P312
      case 6: return h>=k && k>=0 && (k>0 || l>=0);
      case 17: // trigonal-D P321
      case 7: return h>=k && k>=0 && (h>k || l>=0);
      case 8: return h>=0 && ((l>=h && k>h) || (l==h && k==h));
      case 9: return k>=l && l>=h && h>=0;
      case 10: return k>0 || (k==0 && (h>0 || (h==0 && l>=0))); // triclinic
      case 11: return k>=0 && (h>0 || (h==0 && l>=0)); // monoclinic-B
      case 13: return l>=0 && ((k>=0 && h>0) || (h==0 && k==0)); // tetragonal-C, hexagonal-C
      case 15: return (k>=0 && h>0) || (h==0 && k==0 && l>=0); // trigonal-C
      case 18: return k>=0 && l>=0 && ((h>k && h>l) || (h==k && h>=l)); // cubic-T
      case 19: return h>=k && k>=l && l>=0; // cubic-O
    }
    unreachable();
  }

  const char* condition_str() const {
    switch (idx) {
      case 0: return "l>0 or (l=0 and (h>0 or (h=0 and k>=0)))";
      case 1: return "k>=0 and (l>0 or (l=0 and h>=0))";
      case 12:
      case 2: return "h>=0 and k>=0 and l>=0";
      case 3: return "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))";
      case 14:
      case 4: return "h>=k and k>=0 and l>=0";
      case 5: return "(h>=0 and k>0) or (h=0 and k=0 and l>=0)";
      case 16:
      case 6: return "h>=k and k>=0 and (k>0 or l>=0)";
      case 17:
      case 7: return "h>=k and k>=0 and (h>k or l>=0)";
      case 8: return "h>=0 and ((l>=h and k>h) or (l=h and k=h))";
      case 9: return "k>=l and l>=h and h>=0";
      case 10: return "k>0 or (k==0 and (h>0 or (h=0 and l>=0)))";
      case 11: return "k>=0 and (h>0 or (h=0 and l>=0))";
      case 13: return "l>=0 and ((k>=0 and h>0) or (h=0 and k==0))";
      case 15: return "(k>=0 and h>0) or (h=0 and k==0 and l>=0)";
      case 18: return "k>=0 and l>=0 and ((h>k and h>l) or (h=k and h>=l))";
      case 19: return "h>=k and k>=l and l>=0";
    }
    unreachable();
  }

  /// Returns hkl in asu and MTZ ISYM - 2*n-1 for reflections in the positive
  /// asu (I+ of a Friedel pair), 2*n for reflections in the negative asu (I-).
  std::pair<Op::Miller,int> to_asu(const Op::Miller& hkl, const std::vector<Op>& sym_ops) const {
    int isym = 0;
    for (const Op& op : sym_ops) {
      ++isym;
      Op::Miller new_hkl = op.apply_to_hkl_without_division(hkl);
      if (is_in(new_hkl))
        return {Op::divide_hkl_by_DEN(new_hkl), isym};
      ++isym;
      Op::Miller negated_new_hkl{{-new_hkl[0], -new_hkl[1], -new_hkl[2]}};
      if (is_in(negated_new_hkl))
        return {Op::divide_hkl_by_DEN(negated_new_hkl), isym};
    }
    fail("Oops, maybe inconsistent GroupOps?");
  }

  std::pair<Op::Miller,int> to_asu(const Op::Miller& hkl, const GroupOps& gops) const {
    return to_asu(hkl, gops.sym_ops);
  }

  /// Similar to to_asu(), but the second returned value is sign: true for + or centric
  std::pair<Op::Miller,bool> to_asu_sign(const Op::Miller& hkl, const GroupOps& gops) const {
    std::pair<Op::Miller,bool> neg = {{0,0,0}, true};
    for (const Op& op : gops.sym_ops) {
      Op::Miller new_hkl = op.apply_to_hkl_without_division(hkl);
      if (is_in(new_hkl))
        return {Op::divide_hkl_by_DEN(new_hkl), true};
      Op::Miller negated_new_hkl{{-new_hkl[0], -new_hkl[1], -new_hkl[2]}};
      if (is_in(negated_new_hkl))
        // don't return it yet, because for centric reflection we prefer (+)
        neg = {Op::divide_hkl_by_DEN(negated_new_hkl), false};
    }
    if (neg.second)
      fail("Oops, maybe inconsistent GroupOps?");
    return neg;
  }
};

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
