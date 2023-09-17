// Copyright 2017 Global Phasing Ltd.
//
// interface to stb_sprintf: snprintf_z, to_str(float|double)

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
/// stb_snprintf in gemmi namespace - like snprintf, but ignores locale
/// and is always zero-terminated (hence _z).
GEMMI_DLL int snprintf_z(char *buf, int count, char const *fmt, ...)
                                                         GEMMI_ATTRIBUTE_FORMAT(3,4);
/// stb_sprintf in gemmi namespace
GEMMI_DLL int sprintf_z(char *buf, char const *fmt, ...) GEMMI_ATTRIBUTE_FORMAT(2,3);

inline std::string to_str(double d) {
  char buf[24];
  int len = sprintf_z(buf, "%.9g", d);
  return std::string(buf, len > 0 ? len : 0);
}

inline std::string to_str(float d) {
  char buf[16];
  int len = sprintf_z(buf, "%.6g", d);
  return std::string(buf, len > 0 ? len : 0);
}

template<int Prec>
std::string to_str_prec(double d) {
  static_assert(Prec >= 0 && Prec < 7, "unsupported precision");
  char buf[16];
  int len = d > -1e8 && d < 1e8 ? sprintf_z(buf, "%.*f", Prec, d)
                                : sprintf_z(buf, "%g", d);
  return std::string(buf, len > 0 ? len : 0);
}

/// zero-terminated to_chars()
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
