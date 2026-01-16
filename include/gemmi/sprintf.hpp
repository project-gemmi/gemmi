//! @file
//! @brief Interface to stb_sprintf: snprintf_z, to_str(float|double).
//!
//! Provides locale-independent string formatting functions for floating-point numbers
//! and integers. Uses stb_sprintf internally for consistent cross-platform behavior.

// Copyright 2017 Global Phasing Ltd.

#ifndef GEMMI_SPRINTF_HPP_
#define GEMMI_SPRINTF_HPP_

#include <string>
#ifdef __has_include
# if __has_include(<charconv>) && !(defined(_MSVC_LANG) && _MSVC_LANG < 201703L)
#  include <charconv>
# endif
#endif

#if __cpp_lib_to_chars < 201611L
# include <algorithm> // for min
#endif

#include "fail.hpp"  // for GEMMI_DLL

namespace gemmi {

// On MinGW format(printf) doesn't support %zu.
#if (defined(__GNUC__) && !defined(__MINGW32__)) || defined(__clang__)
# define GEMMI_ATTRIBUTE_FORMAT(fmt,va) __attribute__((format(printf,fmt,va)))
#else
# define GEMMI_ATTRIBUTE_FORMAT(fmt,va)
#endif

//! @brief Locale-independent snprintf that is always zero-terminated.
//! @param buf Output buffer
//! @param count Buffer size
//! @param fmt Format string (printf-style)
//! @param ... Format arguments
//! @return Number of characters written (excluding null terminator)
//!
//! stb_snprintf in gemmi namespace - like snprintf, but ignores locale
//! and is always zero-terminated (hence _z).
GEMMI_DLL int snprintf_z(char *buf, int count, char const *fmt, ...)
                                                         GEMMI_ATTRIBUTE_FORMAT(3,4);

//! @brief Locale-independent sprintf.
//! @param buf Output buffer (must be large enough)
//! @param fmt Format string (printf-style)
//! @param ... Format arguments
//! @return Number of characters written
//!
//! stb_sprintf in gemmi namespace - ignores locale for consistent output.
GEMMI_DLL int sprintf_z(char *buf, char const *fmt, ...) GEMMI_ATTRIBUTE_FORMAT(2,3);

//! @brief Convert double to string with 9 significant digits.
//! @param d Double value
//! @return String representation
//!
//! Uses "%.9g" format for compact representation with full precision.
inline std::string to_str(double d) {
  char buf[24];
  int len = sprintf_z(buf, "%.9g", d);
  return std::string(buf, len > 0 ? len : 0);
}

//! @brief Convert float to string with 6 significant digits.
//! @param d Float value
//! @return String representation
//!
//! Uses "%.6g" format for compact representation appropriate for float precision.
inline std::string to_str(float d) {
  char buf[16];
  int len = sprintf_z(buf, "%.6g", d);
  return std::string(buf, len > 0 ? len : 0);
}

//! @brief Convert double to string with fixed precision.
//! @tparam Prec Number of decimal places (0-6)
//! @param d Double value
//! @return String representation
//!
//! Uses fixed-point format for values in range [-1e8, 1e8], otherwise uses %g.
template<int Prec>
std::string to_str_prec(double d) {
  static_assert(Prec >= 0 && Prec < 7, "unsupported precision");
  char buf[16];
  int len = d > -1e8 && d < 1e8 ? sprintf_z(buf, "%.*f", Prec, d)
                                : sprintf_z(buf, "%g", d);
  return std::string(buf, len > 0 ? len : 0);
}

//! @brief Convert integer to zero-terminated string in buffer.
//! @param first Start of buffer
//! @param last End of buffer
//! @param value Integer value
//! @return Pointer to null terminator
//!
//! Zero-terminated to_chars() - uses std::to_chars if available, otherwise snprintf_z.
inline char* to_chars_z(char* first, char* last, int value) {
#if __cpp_lib_to_chars >= 201611L
  auto result = std::to_chars(first, last-1, value);
  *result.ptr = '\0';
  return result.ptr;
#else
  int n = snprintf_z(first, int(last - first), "%d", value);
  return std::min(first + n, last - 1);
#endif
}

//! @brief Convert size_t to zero-terminated string in buffer.
//! @param first Start of buffer
//! @param last End of buffer
//! @param value Size value
//! @return Pointer to null terminator
//!
//! Zero-terminated to_chars() for size_t - uses std::to_chars if available.
inline char* to_chars_z(char* first, char* last, size_t value) {
#if __cpp_lib_to_chars >= 201611L
  auto result = std::to_chars(first, last-1, value);
  *result.ptr = '\0';
  return result.ptr;
#else
  int n = snprintf_z(first, int(last - first), "%zu", value);
  return std::min(first + n, last - 1);
#endif
}

} // namespace gemmi
#endif
