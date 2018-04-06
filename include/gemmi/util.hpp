// Copyright 2017 Global Phasing Ltd.
//
// Utilities.

#ifndef GEMMI_UTIL_HPP_
#define GEMMI_UTIL_HPP_

#include <algorithm>  // for equal, find
#include <cctype>     // for tolower
#include <cstdio>     // for FILE, fopen, fclose, fgetc, fgets
#include <iterator>   // for begin, end, make_move_iterator
#include <memory>     // for unique_ptr
#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>

#if defined(_WIN32) && defined(GEMMI_WINDOWS_PATHS_IN_UTF8)
#include <locale>
#include <codecvt>
#endif

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
inline std::vector<std::string> split_str(const std::string &str, S sep) {
  std::vector<std::string> result;
  std::size_t start = 0, end;
  while ((end = str.find(sep, start)) != std::string::npos) {
    result.emplace_back(str, start, end - start);
    start = end + impl::length(sep);
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

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }


// file operations
typedef std::unique_ptr<FILE, decltype(&std::fclose)> fileptr_t;

inline fileptr_t file_open(const char *path, const char *mode) {
  std::FILE* file;
#if defined(_WIN32) && defined(GEMMI_WINDOWS_PATHS_IN_UTF8)
  std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> convert;
  std::wstring wpath = convert.from_bytes(path);
  std::wstring wmode = convert.from_bytes(mode);
  if ((file = ::_wfopen(wpath.c_str(), wmode.c_str())) == nullptr)
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
  bool get_line_stream() const { return false; }
private:
  std::string path_;
};

inline bool is_pdb_code(const std::string& str) {
  return str.length() == 4 && std::isdigit(str[0]) && std::isalnum(str[1]) &&
                              std::isalnum(str[2]) && std::isalnum(str[3]);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
