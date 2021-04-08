// Copyright 2020 Global Phasing Ltd.
//
// Functions that convert string to floating-point number ignoring locale.
// Wrappers around https://github.com/fastfloat/fast_float/

#ifndef GEMMI_ATOF_HPP_
#define GEMMI_ATOF_HPP_

#include "atox.hpp"   // for is_space
#include "third_party/fast_float/fast_float.h"

namespace gemmi {

using fast_float::from_chars_result;

inline from_chars_result fast_from_chars(const char* start, const char* end, double& d) {
  while (is_space(*start))
    ++start;
  if (*start == '+')
    ++start;
  return fast_float::from_chars(start, end, d);
}

inline from_chars_result fast_from_chars(const char* start, double& d) {
  while (is_space(*start))
    ++start;
  if (*start == '+')
    ++start;
  return fast_float::from_chars(start, start + std::strlen(start), d);
}

} // namespace gemmi
#endif
