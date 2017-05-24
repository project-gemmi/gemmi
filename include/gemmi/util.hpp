// Copyright 2017 Global Phasing Ltd.
//
// Utilities.

#ifndef GEMMI_UTIL_HPP_
#define GEMMI_UTIL_HPP_

#include <algorithm>
#include <string>
#include <vector>
#include <cctype>

namespace gemmi {

inline bool ends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl && str.compare(str.length() - sl, sl, suffix) == 0;
}

inline bool starts_with(const std::string& str, const std::string& prefix) {
  size_t sl = prefix.length();
  return str.length() >= sl && str.compare(0, sl, prefix) == 0;
}

// Case-insensitive version. Assumes the suffix is lowercase and ascii.
inline bool iends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl &&
         std::equal(std::begin(suffix), std::end(suffix), str.end() - sl,
                    [](char c1, char c2) { return c1 == std::tolower(c2); });
}

inline std::string path_basename(const std::string& path) {
  size_t pos = path.find_last_of("\\/");
  return pos == std::string::npos ? path : path.substr(pos + 1);
}

template <class T>
bool in_vector(const T& x, const std::vector<T>& v) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
