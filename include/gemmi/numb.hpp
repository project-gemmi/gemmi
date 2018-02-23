// Copyright 2017 Global Phasing Ltd.
//
// Utilities for parsing numbers: integer, real and so-called numb.
// Numb is the numeric type in CIF. It's a number with optional
// standard uncertainty (s.u.) in brackets: 1.23(8).
// mmcif file do not use s.u. though - they define own numeric categories.

#ifndef GEMMI_NUMB_HPP_
#define GEMMI_NUMB_HPP_

#include <cmath>   // for NAN
#include <cstddef> // for size_t
#include <tao/pegtl.hpp>
#include "cifdoc.hpp" // for is_null

namespace gemmi {
namespace cif {
using std::size_t;
namespace pegtl = tao::pegtl;

namespace numb_rules {
  using namespace pegtl;

  struct sign : opt<one<'+', '-'>> {};
  struct e : one<'e', 'E'> {};
  struct exponent : seq<sign, plus<digit>> {};
  struct uint_digit : digit {};
  struct fraction : plus<digit> {};
  struct base : if_then_else<one<'.'>, fraction,
                                       seq<plus<uint_digit>,
                                           opt<one<'.'>, opt<fraction>>>> {};
  // Error in brackets ,as per CIF spec. We ignore the value for now.
  struct err : seq<one<'('>, plus<digit>, one<')'>> {};
  struct numb : seq<sign, base, opt<e, exponent>, opt<err>, pegtl::eof> {};
}

// Actions for getting the number. For now we ignore s.u., so the actions
// are equivalent to locale-independent std::stod().
template<typename Rule> struct ActionNumb : pegtl::nothing<Rule> {};
template<> struct ActionNumb<numb_rules::uint_digit> {
  template<typename Input> static void apply(const Input& in, double& d) {
      d = d * 10 + (*in.begin() - '0');
  }
};
template<> struct ActionNumb<numb_rules::fraction> {
  template<typename Input> static void apply(const Input& in, double& d) {
    double mult = 0.1;
    for (const auto* p = in.begin(); p != in.end(); ++p, mult *= 0.1)
      d += mult * (*p - '0');
  }
};
template<> struct ActionNumb<numb_rules::exponent> {
  template<typename Input> static void apply(const Input& in, double& d) {
    int n = 0;
    bool neg = false;
    const auto* p = in.begin();
    if (*p == '-')
      neg = true;
    else if (*p != '+')
      n = *p - '0';
    for (++p; p != in.end(); ++p)
      n = n * 10 + (*p - '0');
    // We don't expect too many exponents in CIF files, let's have
    // only this small LUT.
    static const double e[] = { 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };
    double uexp = n <= 8 ? e[n] : std::pow(10, n);
    if (neg)
      d /= uexp;
    else
      d *= uexp;
  }
};
template<> struct ActionNumb<numb_rules::numb> {
  template<typename Input> static void apply(const Input& in, double& d) {
    if (*in.begin() == '-')
      d = -d;
  }
};

inline bool is_numb(const std::string& s) {
  pegtl::memory_input<> in(s, "");
  return pegtl::parse<numb_rules::numb, pegtl::nothing>(in);
}

inline double as_number(const std::string& s, double nan=NAN) {
  double d = 0;
  pegtl::memory_input<> in(s, "");
  if (pegtl::parse<numb_rules::numb, ActionNumb>(in, d))
    return d;
  return nan;
}


// Integers.

namespace int_rules {
  using namespace pegtl;
  struct sign : opt<one<'+', '-'>> {};
  struct int_ : seq<sign, plus<digit>, pegtl::eof> {};
}
template<typename Rule> struct ActionInt : pegtl::nothing<Rule> {};
template<> struct ActionInt<pegtl::digit> {
  template<typename Input> static void apply(const Input& in, int& n) {
      n = n * 10 + (*in.begin() - '0');
  }
};
template<> struct ActionInt<int_rules::int_> {
  template<typename Input> static void apply(const Input& in, int& n) {
    if (*in.begin() == '-')
      n = -n;
  }
};


// utility functions

inline int as_int(const std::string& s) {
  int n = 0;
  pegtl::memory_input<> in(s, "");
  if (pegtl::parse<int_rules::int_, ActionInt>(in, n))
    return n;
  throw std::runtime_error("not an integer number: " + s);
}

inline int as_int(const std::string& s, int default_) {
  return is_null(s) ? default_ : as_int(s);
}


} // namespace cif
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
