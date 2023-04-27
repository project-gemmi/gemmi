
#include <gemmi/sprintf.hpp>
#ifndef USE_STD_SNPRINTF

#define STB_SPRINTF_IMPLEMENTATION
#define STB_SPRINTF_STATIC
#define STB_SPRINTF_NOUNALIGNED 1

// Making functions from stb_sprintf static may trigger warnings.
#if defined(__GNUC__)
# pragma GCC diagnostic ignored "-Wunused-function"
#endif

// To use system stb_sprintf.h (not recommended, but some Linux distros
// don't like bundled libraries) define GEMMI_USE_SYSTEM_STB or remove
// third_party/stb_sprintf.h.
# if defined(__has_include)
#  if !__has_include("gemmi/third_party/stb_sprintf.h")
#   define GEMMI_USE_SYSTEM_STB 1
#  endif
# endif
# ifdef GEMMI_USE_SYSTEM_STB
#  pragma message("Using system stb_sprintf.h, not the bundled one. It may not work.")
#  include <stb/stb_sprintf.h>
# else
#  include "gemmi/third_party/stb_sprintf.h"
# endif

namespace gemmi {

// We copy functions from sprintf.h only to have them declared with GEMMI_DLL.
int gstb_sprintf(char *buf, char const *fmt, ...) {
  int result;
  va_list va;
  va_start(va, fmt);
  result = STB_SPRINTF_DECORATE(vsprintfcb)(0, 0, buf, fmt, va);
  va_end(va);
  return result;
}
int gstb_snprintf(char *buf, int count, char const *fmt, ...) {
  int result;
  va_list va;
  va_start(va, fmt);
  result = STB_SPRINTF_DECORATE(vsnprintf)(buf, count, fmt, va);
  va_end(va);
  return result;
}

}  // namespace gemmi
#endif  // !USE_STD_SNPRINTF
