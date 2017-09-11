// Copyright 2017 Global Phasing Ltd.
//
// Utilities.

#ifndef GEMMI_UTIL_HPP_
#define GEMMI_UTIL_HPP_

#include <algorithm>  // for equal, find
#include <cctype>     // for tolower
#include <cmath>      // for floor
#include <cstdio>     // for FILE, fopen, fclose
#include <iterator>   // for begin, end, make_move_iterator
#include <memory>     // for unique_ptr
#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>

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

inline std::vector<std::string> split_str(const std::string &str, char sep) {
  std::vector<std::string> result;
  std::size_t start = 0, end;
  while ((end = str.find(sep, start)) != std::string::npos) {
    result.emplace_back(str, start, end - start);
    start = end + 1;
  }
  result.emplace_back(str, start);
  return result;
}

template<typename T>
std::string join_str(const T &iterable, const std::string& sep) {
  std::string r;
  bool first = true;
  for (const std::string& i : iterable) {
    if (!first)
      r += sep;
    r += i;
    first = false;
  }
  return r;
}


inline std::string path_basename(const std::string& path) {
  size_t pos = path.find_last_of("\\/");
  return pos == std::string::npos ? path : path.substr(pos + 1);
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

inline int iround(double d) { return static_cast<int>(std::floor(d+0.5)); }

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }


// file operations
typedef std::unique_ptr<FILE, decltype(&std::fclose)> fileptr_t;

inline fileptr_t file_open(const char *path, const char *mode) {
  std::FILE* file;
#ifdef _WIN32
  if (::fopen_s(&file, path, mode) != 0)
#else
  if ((file = std::fopen(path, mode)) == nullptr)
#endif
    fail("Failed to open file: " + std::string(path));
  return fileptr_t(file, &std::fclose);
}

inline std::size_t file_size(FILE* f, const std::string& path) {
  if (std::fseek(f, 0, SEEK_END) != 0)
    fail(path + ": fseek failed");
  long length = std::ftell(f);
  if (length < 0)
    fail(path + ": ftell failed");
  if (std::fseek(f, 0, SEEK_SET) != 0)
    fail(path + ": fseek failed");
  return length;
}

// for transparent handling of stdin along filenames
class MaybeStdin {
public:
  explicit MaybeStdin(const std::string& path) : path_(path) {}
  bool is_stdin() const { return path_ == "-"; };
  const std::string& path() const { return path_; };
  size_t mem_size() const { return 0; };
  std::unique_ptr<char[]> memory() { return nullptr; }
private:
  std::string path_;
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
