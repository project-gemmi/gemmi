// Copyright 2019 Global Phasing Ltd.
//
// gemmi::tostr() - converts a list of arguments to string (uses ostringstream).

#ifndef GEMMI_TOSTR_HPP_
#define GEMMI_TOSTR_HPP_

#include <sstream>
#include <string>

namespace gemmi {

namespace impl {
inline void add_to_stream(std::ostringstream&) {}

template<typename T, typename... Args>
void add_to_stream(std::ostringstream& os, T&& value, Args&&... args) {
  os << std::forward<T>(value);
  add_to_stream(os, std::forward<Args>(args)...);
}

inline void add_to_string(std::string&) {}

template <typename T, typename... Args>
void add_to_string(std::string& out, const T& value, Args const&... args) {
  out += value;
  add_to_string(out, args...);
}
} // namespace impl

template<typename T, typename... Args>
std::string tostr(T&& value, Args&&... args) {
  std::ostringstream os;
  impl::add_to_stream(os, std::forward<T>(value), std::forward<Args>(args)...);
  return os.str();
}

template <class... Args>
std::string cat(Args const&... args) {
  std::string out;
  impl::add_to_string(out, args...);
  return out;
}

} // namespace gemmi
#endif
