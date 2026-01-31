//! @file
//! @brief Utilities for working with strings and vectors.
//!
//! Provides helper functions for string manipulation (concatenation, case conversion,
//! trimming, splitting, joining), case-insensitive comparisons, and vector operations.
//! These utilities are used throughout GEMMI for text processing and collection management.

// Copyright 2017 Global Phasing Ltd.

#ifndef GEMMI_UTIL_HPP_
#define GEMMI_UTIL_HPP_

#include <cassert>
#include <cctype>     // for isspace
#include <cstring>    // for strncmp
#include <algorithm>  // for equal, find, remove_if
#include <iterator>   // for begin, end, make_move_iterator
#include <string>
#include <vector>

namespace gemmi {

//   #####   string helpers   #####

//! @brief Append integer to string.
//! @param out String to append to
//! @param v Integer value
inline void append_to_str(std::string& out, int v) { out += std::to_string(v); }

//! @brief Append size_t to string.
//! @param out String to append to
//! @param v Size value
inline void append_to_str(std::string& out, size_t v) { out += std::to_string(v); }

//! Double overload is deleted to prevent accidental use without proper formatting.
void append_to_str(std::string& out, double) = delete;

//! @brief Generic append to string (for string-like types).
//! @tparam T Type with += operator
//! @param out String to append to
//! @param v Value to append
template<typename T>
void append_to_str(std::string& out, const T& v) { out += v; }

//! @brief Base case for cat_to recursion.
//! @param out Unused output string
inline void cat_to(std::string&) {}

//! @brief Concatenate multiple values to a string (variadic).
//! @tparam T Type of first value
//! @tparam Args Types of remaining values
//! @param out String to append to
//! @param value First value
//! @param args Remaining values
template <typename T, typename... Args>
void cat_to(std::string& out, const T& value, Args const&... args) {
  append_to_str(out, value);
  cat_to(out, args...);
}

//! @brief Concatenate multiple values into a new string.
//! @tparam Args Types of values to concatenate
//! @param args Values to concatenate
//! @return Concatenated string
template <class... Args>
std::string cat(Args const&... args) {
  std::string out;
  cat_to(out, args...);
  return out;
}

//! @brief Check if string starts with prefix.
//! @param str String to check
//! @param prefix Prefix to search for
//! @return True if str starts with prefix
inline bool starts_with(const std::string& str, const std::string& prefix) {
  size_t sl = prefix.length();
  return str.length() >= sl && str.compare(0, sl, prefix) == 0;
}

//! @brief Check if C-string starts with literal.
//! @tparam N Size of literal array
//! @param a C-string to check
//! @param b String literal
//! @return True if a starts with b
template<size_t N> bool starts_with(const char* a, const char (&b)[N]) {
  return std::strncmp(a, b, N-1) == 0;
}

//! @brief Check if string ends with suffix.
//! @param str String to check
//! @param suffix Suffix to search for
//! @return True if str ends with suffix
inline bool ends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl && str.compare(str.length() - sl, sl, suffix) == 0;
}

//! @brief Convert character to lowercase (ASCII only).
//! @param c Character to convert
//! @return Lowercase character
//!
//! Can be faster than std::tolower() because it takes char not int.
inline char lower(char c) {
  if (c >= 'A' && c <= 'Z')
    return c | 0x20;
  return c;
}

//! @brief Convert character to uppercase (ASCII only).
//! @param c Character to convert
//! @return Uppercase character
//!
//! Works as expected only for a-zA-Z.
inline char alpha_up(char c) { return c & ~0x20; }

//! @brief Convert string to lowercase (ASCII only).
//! @param str String to convert
//! @return Lowercase string
inline std::string to_lower(std::string str) {
  for (char& c : str)
    if (c >= 'A' && c <= 'Z')
      c |= 0x20;
  return str;
}

//! @brief Convert string to uppercase (ASCII only).
//! @param str String to convert
//! @return Uppercase string
inline std::string to_upper(std::string str) {
  for (char& c : str)
    if (c >= 'a' && c <= 'z')
      c &= ~0x20;
  return str;
}

//! @brief Case-insensitive character comparison (ASCII only).
//! @param a First character
//! @param b Second character
//! @return True if characters are equal (case-insensitive)
inline bool isame(char a, char b) {
  return a == b || ((a^b) == 0x20 && (a|0x20) >= 'a' && (a|0x20) <= 'z');
}

// Case-insensitive comparisons. The second arg must be lowercase.

//! @brief Case-insensitive string comparison from offset.
//! @param str String to check
//! @param offset Offset into str
//! @param low Lowercase string to compare against
//! @return True if strings match (case-insensitive)
inline bool iequal_from(const std::string& str, size_t offset, const std::string& low) {
  return str.length() == low.length() + offset &&
         std::equal(std::begin(low), std::end(low), str.begin() + offset,
                    [](char c1, char c2) { return c1 == lower(c2); });
}

//! @brief Case-insensitive string comparison.
//! @param str String to check
//! @param low Lowercase string to compare against
//! @return True if strings match (case-insensitive)
inline bool iequal(const std::string& str, const std::string& low) {
  return iequal_from(str, 0, low);
}

//! @brief Case-insensitive prefix check.
//! @param str String to check
//! @param prefix Lowercase prefix to search for
//! @return True if str starts with prefix (case-insensitive)
inline bool istarts_with(const std::string& str, const std::string& prefix) {
  return str.length() >= prefix.length() &&
         std::equal(std::begin(prefix), std::end(prefix), str.begin(),
                    [](char c1, char c2) { return c1 == lower(c2); });
}

//! @brief Case-insensitive suffix check.
//! @param str String to check
//! @param suffix Lowercase suffix to search for
//! @return True if str ends with suffix (case-insensitive)
inline bool iends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl &&
         std::equal(std::begin(suffix), std::end(suffix), str.end() - sl,
                    [](char c1, char c2) { return c1 == lower(c2); });
}

//! @brief Check if string ends with suffix or suffix.gz.
//! @param str String to check
//! @param suffix Lowercase suffix to search for
//! @return True if str ends with suffix or suffix.gz (case-insensitive)
inline bool giends_with(const std::string& str, const std::string& suffix) {
  return iends_with(str, suffix) || iends_with(str, suffix + ".gz");
}

//! @brief Trim whitespace from both ends of string.
//! @param str String to trim
//! @return Trimmed string
inline std::string trim_str(const std::string& str) {
  const std::string ws = " \r\n\t";
  std::string::size_type first = str.find_first_not_of(ws);
  if (first == std::string::npos)
    return std::string{};
  std::string::size_type last = str.find_last_not_of(ws);
  return str.substr(first, last - first + 1);
}

//! @brief Trim whitespace from end of string.
//! @param str String to trim
//! @return Right-trimmed string
inline std::string rtrim_str(const std::string& str) {
  std::string::size_type last = str.find_last_not_of(" \r\n\t");
  return str.substr(0, last == std::string::npos ? 0 : last + 1);
}

//! @brief Find end of trimmed C-string (excluding trailing whitespace).
//! @param start Start of string
//! @param end End of string (or nullptr to search for null terminator)
//! @return Pointer after last non-whitespace character
//!
//! end is after the last character of the string (typically \0).
inline const char* rtrim_cstr(const char* start, const char* end=nullptr) {
  if (!start)
    return nullptr;
  if (!end) {
    end = start;
    while (*end != '\0')
      ++end;
  }
  while (end > start && std::isspace(end[-1]))
    --end;
  return end;
}

namespace impl {
inline size_t length(char) { return 1; }
inline size_t length(const std::string& s) { return s.length(); }
}

//! @brief Split string by separator into vector.
//! @tparam S Separator type (char or string)
//! @param str String to split
//! @param sep Separator
//! @param result Vector to append results to
//!
//! Takes a single separator (usually char or string). May return empty fields.
template<typename S>
void split_str_into(const std::string& str, S sep,
                    std::vector<std::string>& result) {
  std::size_t start = 0, end;
  while ((end = str.find(sep, start)) != std::string::npos) {
    result.emplace_back(str, start, end - start);
    start = end + impl::length(sep);
  }
  result.emplace_back(str, start);
}

//! @brief Split string by separator.
//! @tparam S Separator type (char or string)
//! @param str String to split
//! @param sep Separator
//! @return Vector of substrings
template<typename S>
std::vector<std::string> split_str(const std::string& str, S sep) {
  std::vector<std::string> result;
  split_str_into(str, sep, result);
  return result;
}

//! @brief Split string by multiple separators into vector.
//! @param str String to split
//! @param seps String of separator characters
//! @param result Vector to append results to
//!
//! Takes multiple 1-char separators as a string. Discards empty fields.
inline void split_str_into_multi(const std::string& str, const char* seps,
                                 std::vector<std::string>& result) {
  std::size_t start = str.find_first_not_of(seps);
  while (start != std::string::npos) {
    std::size_t end = str.find_first_of(seps, start);
    result.emplace_back(str, start, end - start);
    start = str.find_first_not_of(seps, end);
  }
}

//! @brief Split string by multiple separators.
//! @param str String to split
//! @param seps String of separator characters
//! @return Vector of substrings (empty fields discarded)
inline std::vector<std::string> split_str_multi(const std::string& str,
                                                const char* seps=" \t") {
  std::vector<std::string> result;
  split_str_into_multi(str, seps, result);
  return result;
}

//! @brief Join iterator range to string with separator and getter.
//! @tparam T Iterator type
//! @tparam S Separator type
//! @tparam F Getter function type
//! @param begin Start iterator
//! @param end End iterator
//! @param sep Separator
//! @param getter Function to extract string from element
//! @return Joined string
template<typename T, typename S, typename F>
std::string join_str(T begin, T end, const S& sep, const F& getter) {
  std::string r;
  bool first = true;
  for (T i = begin; i != end; ++i) {
    if (!first)
      r += sep;
    r += getter(*i);
    first = false;
  }
  return r;
}

//! @brief Join iterator range to string with separator.
//! @tparam T Iterator type
//! @tparam S Separator type
//! @param begin Start iterator
//! @param end End iterator
//! @param sep Separator
//! @return Joined string
template<typename T, typename S>
std::string join_str(T begin, T end, const S& sep) {
  return join_str(begin, end, sep, [](const std::string& t) { return t; });
}

//! @brief Join iterable container to string with separator and getter.
//! @tparam T Container type
//! @tparam S Separator type
//! @tparam F Getter function type
//! @param iterable Container
//! @param sep Separator
//! @param getter Function to extract string from element
//! @return Joined string
template<typename T, typename S, typename F>
std::string join_str(const T& iterable, const S& sep, const F& getter) {
  return join_str(iterable.begin(), iterable.end(), sep, getter);
}

//! @brief Join iterable container to string with separator.
//! @tparam T Container type
//! @tparam S Separator type
//! @param iterable Container
//! @param sep Separator
//! @return Joined string
template<typename T, typename S>
std::string join_str(const T& iterable, const S& sep) {
  return join_str(iterable.begin(), iterable.end(), sep);
}

//! @brief Append item to string with separator (only if string is non-empty).
//! @tparam T Item type
//! @tparam S Separator type
//! @param str String to append to
//! @param sep Separator
//! @param item Item to append
template<typename T, typename S>
void string_append_sep(std::string& str, S sep, const T& item) {
  if (!str.empty())
    str += sep;
  str += item;
}

//! @brief Replace all occurrences of substring.
//! @param s String to modify
//! @param old Substring to replace
//! @param new_ Replacement substring
inline void replace_all(std::string &s,
                        const std::string &old, const std::string &new_) {
  std::string::size_type pos = 0;
  while ((pos = s.find(old, pos)) != std::string::npos) {
    s.replace(pos, old.size(), new_);
    pos += new_.size();
  }
}

//! @brief Check if name is in comma-separated list.
//! @param name Name to search for
//! @param list Comma-separated string
//! @param sep Separator character
//! @return True if name is in list
//!
//! list is a comma separated string.
inline bool is_in_list(const std::string& name, const std::string& list,
                       char sep=',') {
  if (name.length() >= list.length())
    return name == list;
  for (size_t start=0, end=0; end != std::string::npos; start=end+1) {
    end = list.find(sep, start);
    if (list.compare(start, end - start, name) == 0)
      return true;
  }
  return false;
}

//   #####   vector helpers   #####

//! @brief Check if value is in vector.
//! @tparam T Element type
//! @param x Value to search for
//! @param v Vector to search
//! @return True if x is in v
template <class T>
bool in_vector(const T& x, const std::vector<T>& v) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

//! @brief Check if predicate matches any vector element.
//! @tparam F Predicate function type
//! @tparam T Element type
//! @param f Predicate function
//! @param v Vector to search
//! @return True if f matches any element
template <typename F, typename T>
bool in_vector_f(F f, const std::vector<T>& v) {
  return std::find_if(v.begin(), v.end(), f) != v.end();
}

//! @brief Get pointer to one-past-end of vector.
//! @tparam T Element type
//! @param v Vector
//! @return Pointer to end
template <class T>
T* vector_end_ptr(std::vector<T>& v) { return v.data() + v.size(); }

//! @brief Get const pointer to one-past-end of vector.
//! @tparam T Element type
//! @param v Vector
//! @return Const pointer to end
template <class T>
const T* vector_end_ptr(const std::vector<T>& v) { return v.data() + v.size(); }

//! @brief Move-extend destination vector with source vector.
//! @tparam T Element type
//! @param dst Destination vector
//! @param src Source vector (will be moved from)
template <class T>
void vector_move_extend(std::vector<T>& dst, std::vector<T>&& src) {
  if (dst.empty())
    dst = std::move(src);
  else
    dst.insert(dst.end(), std::make_move_iterator(src.begin()),
                          std::make_move_iterator(src.end()));
}

//! @brief Remove elements matching condition (erase-remove idiom).
//! @tparam T Element type
//! @tparam F Predicate function type
//! @param v Vector to modify
//! @param condition Predicate function
//!
//! Wrapper around the erase-remove idiom.
template <class T, typename F>
void vector_remove_if(std::vector<T>& v, F&& condition) {
  v.erase(std::remove_if(v.begin(), v.end(), condition), v.end());
}

//! @brief Insert columns into 2D array stored in vector.
//! @tparam T Element type
//! @param data 2D array (old_width x length) in a vector
//! @param old_width Number of columns before insertion
//! @param length Number of rows
//! @param n Number of columns to insert
//! @param pos Column position for insertion
//! @param new_value Value for new columns
//!
//! Insert n new columns at position pos.
template <class T>
void vector_insert_columns(std::vector<T>& data, size_t old_width,
                           size_t length, size_t n, size_t pos, const T& new_value) {
  assert(data.size() == old_width * length);
  assert(pos <= old_width);
  data.resize(data.size() + n * length);
  typename std::vector<T>::iterator dst = data.end();
  for (size_t i = length; i-- != 0; ) {
    for (size_t j = old_width; j-- != pos; )
      *--dst = data[i * old_width + j];
    for (size_t j = n; j-- != 0; )
      *--dst = new_value;
    for (size_t j = pos; j-- != 0; )
      *--dst = data[i * old_width + j];
  }
  assert(dst == data.begin());
}

//! @brief Remove column from 2D array stored in vector.
//! @tparam T Element type
//! @param data 2D array with new_width+1 columns, in a vector
//! @param new_width Number of columns after removal
//! @param pos Column position to remove
//!
//! Remove column at position pos.
template <class T>
void vector_remove_column(std::vector<T>& data, size_t new_width, size_t pos) {
  assert(pos <= new_width);
  for (size_t source = pos + 1; source < data.size(); ++source)
    for (size_t i = 0; i < new_width && source < data.size(); ++i)
      data[pos++] = data[source++];
  data.resize(pos);
}


//   #####   other helpers   #####

//! @brief Numeric ID for case-insensitive comparison of 4 letters.
//! @param s String with 4 chars or 3 chars + NUL
//! @return Numeric ID
//!
//! Numeric ID used for case-insensitive comparison of 4 letters.
//! s must have 4 chars or 3 chars + NUL, ' ' and NUL are equivalent in s.
constexpr int ialpha4_id(const char* s) {
  return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
}

//! @brief Numeric ID for case-insensitive comparison of 3 letters.
//! @param s String with 3 chars
//! @return Numeric ID
//!
//! Numeric ID used for case-insensitive comparison of 3 letters.
constexpr int ialpha3_id(const char* s) {
  return (s[0] << 16 | s[1] << 8 | s[2]) & ~0x20202020;
}

} // namespace gemmi
#endif
