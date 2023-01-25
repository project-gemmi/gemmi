// Copyright 2017 Global Phasing Ltd.
//
// interface to stb_sprintf: gstb_snprintf, to_str(float|double)

#ifndef GEMMI_SPRINTF_HPP_
#define GEMMI_SPRINTF_HPP_

#ifdef USE_STD_SNPRINTF  // for benchmarking and testing only
# include <cstdio>
# define gstb_snprintf std::snprintf
# define gstb_sprintf std::sprintf
#else
# define STB_SPRINTF_DECORATE(name) gstb_##name
// To use system stb_sprintf.h (not recommended, but some Linux distros
// don't like bundled libraries) just remove third_party/stb_sprintf.h.
# if defined(__has_include)
#  if !__has_include("third_party/stb_sprintf.h")
#   define GEMMI_USE_SYSTEM_STB 1
#  endif
# endif
# ifdef GEMMI_USE_SYSTEM_STB
#  pragma message("Using system stb_sprintf.h, not the bundled one. It may not work.")
#  include <stb/stb_sprintf.h>
# else
#  include "third_party/stb_sprintf.h"
# endif
#endif
#include <string>

namespace gemmi {

inline std::string to_str(double d) {
  char buf[24];
  int len = gstb_sprintf(buf, "%.9g", d);
  return std::string(buf, len > 0 ? len : 0);
}

inline std::string to_str(float d) {
  char buf[16];
  int len = gstb_sprintf(buf, "%.6g", d);
  return std::string(buf, len > 0 ? len : 0);
}

template<int Prec>
std::string to_str_prec(double d) {
  static_assert(Prec >= 0 && Prec < 7, "unsupported precision");
  char buf[16];
  int len = d > -1e8 && d < 1e8 ? gstb_sprintf(buf, "%.*f", Prec, d)
                                : gstb_sprintf(buf, "%g", d);
  return std::string(buf, len > 0 ? len : 0);
}

} // namespace gemmi
#endif
