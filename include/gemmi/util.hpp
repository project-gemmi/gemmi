// Copyright 2017 Global Phasing Ltd.
//
// Utilities.

#ifndef GEMMI_UTIL_HPP_
#define GEMMI_UTIL_HPP_

#include <algorithm>  // for equal, find
#include <cctype>     // for tolower
#include <iterator>   // for begin, end, make_move_iterator
#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>

namespace gemmi {

inline bool starts_with(const std::string& str, const std::string& prefix) {
  size_t sl = prefix.length();
  return str.length() >= sl && str.compare(0, sl, prefix) == 0;
}

inline bool ends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl && str.compare(str.length() - sl, sl, suffix) == 0;
}

// Case-insensitive version. Assumes the prefix/suffix is lowercase and ascii.
inline bool istarts_with(const std::string& str, const std::string& prefix) {
  return str.length() >= prefix.length() &&
         std::equal(std::begin(prefix), std::end(prefix), str.begin(),
                    [](char c1, char c2) { return c1 == std::tolower(c2); });
}
inline bool iends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl &&
         std::equal(std::begin(suffix), std::end(suffix), str.end() - sl,
                    [](char c1, char c2) { return c1 == std::tolower(c2); });
}

inline bool giends_with(const std::string& str, const std::string& suffix) {
  return iends_with(str, suffix) || iends_with(str, suffix + ".gz");
}

inline std::string to_lower(std::string str) {
  for (char& c : str)
    if (c >= 'A' && c <= 'Z')
      c |= 0x20;
  return str;
}

inline std::string trim_str(const std::string& str)
{
  std::string ws = " \r\n\t";
  std::string::size_type first = str.find_first_not_of(ws);
  if (first == std::string::npos)
    return std::string{};
  std::string::size_type last = str.find_last_not_of(ws);
  return str.substr(first, last - first + 1);
}

namespace impl {
inline size_t length(char) { return 1; }
inline size_t length(const std::string& s) { return s.length(); }
}

template<typename S>
inline std::vector<std::string> split_str(const std::string& str, S sep) {
  std::vector<std::string> result;
  std::size_t start = 0, end;
  while ((end = str.find(sep, start)) != std::string::npos) {
    result.emplace_back(str, start, end - start);
    start = end + impl::length(sep);
  }
  result.emplace_back(str, start);
  return result;
}

template<typename T, typename S, typename F>
std::string join_str(const T& iterable, const S& sep, const F& getter) {
  std::string r;
  bool first = true;
  for (const auto& item : iterable) {
    if (!first)
      r += sep;
    r += getter(item);
    first = false;
  }
  return r;
}

template<typename T, typename S>
std::string join_str(const T& iterable, const S& sep) {
  return join_str(iterable, sep, [](const std::string& t) { return t; });
}

template <class T>
bool in_vector(const T& x, const std::vector<T>& v) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

template <class T>
void vector_move_extend(std::vector<T>& dst, std::vector<T>&& src) {
  if (dst.empty())
    dst = std::move(src);
  else
    dst.insert(dst.end(), std::make_move_iterator(src.begin()),
                          std::make_move_iterator(src.end()));
}

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
