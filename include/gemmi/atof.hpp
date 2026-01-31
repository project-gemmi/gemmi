//! @file
//! @brief Locale-independent string-to-floating-point conversion.
//!
//! Functions that convert strings to floating-point numbers ignoring locale.
//! Simple wrappers around fastfloat::from_chars() with leading whitespace
//! and '+' sign handling.

// Copyright 2020 Global Phasing Ltd.

#ifndef GEMMI_ATOF_HPP_
#define GEMMI_ATOF_HPP_

#include "atox.hpp"   // for is_space
#include "third_party/fast_float.h"

namespace gemmi {

using fast_float::from_chars_result;

//! @brief Convert string to double (with range, skipping whitespace and '+').
//! @param start Start of string
//! @param end End of string
//! @param d Output double value
//! @return Conversion result with pointer to first unconverted character
inline from_chars_result fast_from_chars(const char* start, const char* end, double& d) {
  while (start < end && is_space(*start))
    ++start;
  if (start < end && *start == '+')
    ++start;
  return fast_float::from_chars(start, end, d);
}

//! @brief Convert C-string to double (skipping whitespace and '+').
//! @param start Start of null-terminated string
//! @param d Output double value
//! @return Conversion result with pointer to first unconverted character
inline from_chars_result fast_from_chars(const char* start, double& d) {
  while (is_space(*start))
    ++start;
  if (*start == '+')
    ++start;
  return fast_float::from_chars(start, start + std::strlen(start), d);
}

//! @brief Convert C-string to double (atof replacement).
//! @param p C-string to convert
//! @param endptr Optional pointer to store end of conversion
//! @return Converted double value
inline double fast_atof(const char* p, const char** endptr=nullptr) {
  double d = 0;
  auto result = fast_from_chars(p, d);
  if (endptr)
    *endptr = result.ptr;
  return d;
}

} // namespace gemmi
#endif
