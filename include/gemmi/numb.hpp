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
#include "third_party/tao/pegtl.hpp"

namespace gemmi {
namespace cif {
using std::size_t;
namespace pegtl = tao::pegtl;

namespace numb_rules {
  struct sign : pegtl::opt<pegtl::one<'+', '-'>> {};
  struct e : pegtl::one<'e', 'E'> {};
  struct exponent : pegtl::seq<sign, pegtl::plus<pegtl::digit>> {};
  struct uint_digit : pegtl::digit {};
  struct fraction : pegtl::plus<pegtl::digit> {};
  struct full_base : pegtl::seq<pegtl::plus<uint_digit>,
                                pegtl::opt<pegtl::one<'.'>,
                                           pegtl::opt<fraction>>> {};
  struct base : pegtl::if_then_else<pegtl::one<'.'>, fraction, full_base> {};
  // Error in brackets ,as per CIF spec. We ignore the value for now.
  struct err : pegtl::seq<pegtl::one<'('>, pegtl::plus<pegtl::digit>,
                          pegtl::one<')'>> {};
  struct numb : pegtl::seq<sign, base, pegtl::opt<e, exponent>,
                           pegtl::opt<err>, pegtl::eof> {};
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

// for use in templates (see also as_any() functions in cifdoc.hpp)
inline float as_any(const std::string& s, float null) {
  return (float) as_number(s, null);
}
inline double as_any(const std::string& s, double null) {
  return as_number(s, null);
}

} // namespace cif
} // namespace gemmi
#endif
