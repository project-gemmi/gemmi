// Copyright 2017 Global Phasing Ltd.
//
// fail(), unreachable() and __declspec/__attribute__ macros

#ifndef GEMMI_FAIL_HPP_
#define GEMMI_FAIL_HPP_

#include <cerrno>     // for errno
#include <stdexcept>  // for runtime_error
#include <system_error> // for system_error
#include <string>
#include <utility>    // for forward

#ifdef  __INTEL_COMPILER
// warning #2196: routine is both "inline" and "noinline"
# pragma warning disable 2196
#endif
#if defined(__GNUG__) && !defined(__clang__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wattributes"
#endif

#if defined(__GNUC__) || defined(__clang__)
# define GEMMI_COLD __attribute__((cold))
#elif defined(_MSC_VER)
# define GEMMI_COLD __declspec(noinline)
#else
# define GEMMI_COLD __attribute__((noinline))
#endif

#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
#  define GEMMI_LIKELY(x) (x) [[likely]]
#  define GEMMI_UNLIKELY(x) (x) [[unlikely]]
#elif defined(__GNUC__) || defined(__clang__)
#  define GEMMI_LIKELY(x) (__builtin_expect(!!(x), 1))
#  define GEMMI_UNLIKELY(x) (__builtin_expect(!!(x), 0))
#else
#  define GEMMI_LIKELY(x) (x)
#  define GEMMI_UNLIKELY(x) (x)
#endif

#if defined(_WIN32)
# if defined(GEMMI_SHARED)
#  if defined(GEMMI_BUILD)
#   define GEMMI_DLL __declspec(dllexport)
#  else
#   define GEMMI_DLL __declspec(dllimport)
#  endif  // GEMMI_BUILD
# else
#  define GEMMI_DLL
# endif  // GEMMI_SHARED
#else
# define GEMMI_DLL __attribute__((visibility("default")))
#endif

namespace gemmi {

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

template<typename T, typename... Args> [[noreturn]]
void fail(std::string&& str, T&& arg1, Args&&... args) {
  str += arg1;
  fail(std::move(str), std::forward<Args>(args)...);
}

[[noreturn]]
inline GEMMI_COLD void fail(const char* msg) { throw std::runtime_error(msg); }

[[noreturn]]
inline GEMMI_COLD void sys_fail(const std::string& msg) {
  throw std::system_error(errno, std::system_category(), msg);
}
[[noreturn]]
inline GEMMI_COLD void sys_fail(const char* msg) {
  throw std::system_error(errno, std::system_category(), msg);
}

// unreachable() is used to silence GCC -Wreturn-type and hint the compiler
[[noreturn]] inline void unreachable() {
#if defined(__GNUC__) || defined(__clang__)
  __builtin_unreachable();
#elif defined(_MSC_VER)
  __assume(0);
#endif
}

#if defined(__GNUG__) && !defined(__clang__)
# pragma GCC diagnostic pop
#endif

} // namespace gemmi
#endif
