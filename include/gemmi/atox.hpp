//! @file
//! @brief Locale-independent string to integer conversion.
//!
//! Locale-independent functions that convert strings to integers,
//! equivalents of standard isspace and isdigit, and a few helper functions.
//!
//! This file is named similarly to the standard functions atoi() and atof().
//! But the functions here are not meant to be equivalent to the standard
//! library functions. They are locale-independent (a good thing when reading
//! numbers from files). They don't set errno, don't signal overflow and
//! underflow. Due to the limited scope these functions tend to be faster
//! than the standard-library ones.

// Copyright 2018 Global Phasing Ltd.
//
// Locale-independent functions that convert strings to integers,
// equivalents of standard isspace and isdigit, and a few helper functions.
//
// This file is named similarly to the standard functions atoi() and atof().
// But the functions here are not meant to be equivalent to the standard
// library functions. They are locale-independent (a good thing when reading
// numbers from files). They don't set errno, don't signal overflow and
// underflow. Due to the limited scope these functions tend to be faster
// than the standard-library ones.

#ifndef GEMMI_ATOX_HPP_
#define GEMMI_ATOX_HPP_

#include <cstdint>
#include <stdexcept>  // for invalid_argument
#include <string>

namespace gemmi {

//! @brief Check if character is whitespace (C locale).
//! @param c Character to test
//! @return True if c is whitespace (9-13 or 32)
//!
//! Equivalent of std::isspace for C locale (no handling of EOF).
inline bool is_space(char c) {
  static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
  };
  return table[(std::uint8_t)c] != 0;
}

//! @brief Check if character is blank (C locale).
//! @param c Character to test
//! @return True if c is space or tab
//!
//! Equivalent of std::isblank for C locale (no handling of EOF).
inline bool is_blank(char c) {
  return c == ' ' || c == '\t';
}

//! @brief Check if character is digit (C locale).
//! @param c Character to test
//! @return True if c is '0'-'9'
//!
//! Equivalent of std::isdigit for C locale (no handling of EOF).
inline bool is_digit(char c) {
  return c >= '0' && c <= '9';
}

//! @brief Skip blank characters.
//! @param p Pointer to string
//! @return Pointer past blank characters
inline const char* skip_blank(const char* p) {
  if (p)
    while (is_blank(*p))
      ++p;
  return p;
}

//! @brief Skip non-whitespace characters.
//! @param p Pointer to string
//! @return Pointer to next whitespace or null terminator
inline const char* skip_word(const char* p) {
  if (p)
    while (*p != '\0' && !is_space(*p))
      ++p;
  return p;
}

//! @brief Read word from line.
//! @param line Input line
//! @return Word (whitespace-delimited)
inline std::string read_word(const char* line) {
  line = skip_blank(line);
  return std::string(line, skip_word(line));
}

//! @brief Read word from line with end pointer.
//! @param line Input line
//! @param endptr Output: pointer past the word
//! @return Word (whitespace-delimited)
inline std::string read_word(const char* line, const char** endptr) {
  line = skip_blank(line);
  *endptr = skip_word(line);
  return std::string(line, *endptr);
}

//! @brief Convert string to integer.
//! @param p String pointer
//! @param checked If true, validate full string is consumed
//! @param length Maximum length to parse (0=unlimited)
//! @return Parsed integer
//! @throws std::invalid_argument if checked=true and parsing fails
//!
//! No checking for overflow.
inline int string_to_int(const char* p, bool checked, size_t length=0) {
  int mult = -1;
  int n = 0;
  size_t i = 0;
  while ((length == 0 || i < length) && is_space(p[i]))
    ++i;
  if (p[i] == '-') {
    mult = 1;
    ++i;
  } else if (p[i] == '+') {
    ++i;
  }
  bool has_digits = false;
  // use negative numbers because INT_MIN < -INT_MAX
  for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    n = n * 10 - (p[i] - '0');
    has_digits = true;
  }
  if (checked) {
    while ((length == 0 || i < length) && is_space(p[i]))
      ++i;
    if (!has_digits || p[i] != '\0')
      throw std::invalid_argument("not an integer: " +
                                  std::string(p, length ? length : i+1));
  }
  return mult * n;
}

//! @brief Convert string to integer.
//! @param str String to parse
//! @param checked If true, validate full string is consumed
//! @return Parsed integer
//! @throws std::invalid_argument if checked=true and parsing fails
inline int string_to_int(const std::string& str, bool checked) {
  return string_to_int(str.c_str(), checked);
}

//! @brief Simple atoi implementation.
//! @param p String pointer
//! @param endptr Output: pointer past parsed integer
//! @return Parsed integer
inline int simple_atoi(const char* p, const char** endptr=nullptr) {
  int mult = -1;
  int n = 0;
  while (is_space(*p))
    ++p;
  if (*p == '-') {
    mult = 1;
    ++p;
  } else if (*p == '+') {
    ++p;
  }
  for (; is_digit(*p); ++p)
    n = n * 10 - (*p - '0'); // use negative numbers because INT_MIN < -INT_MAX
  if (endptr)
    *endptr = p;
  return mult * n;
}

//! @brief Parse unsigned integer (no sign allowed).
//! @param p String pointer
//! @param endptr Output: pointer past parsed integer
//! @return Parsed integer
inline int no_sign_atoi(const char* p, const char** endptr=nullptr) {
  int n = 0;
  while (is_space(*p))
    ++p;
  for (; is_digit(*p); ++p)
    n = n * 10 + (*p - '0');
  if (endptr)
    *endptr = p;
  return n;
}

} // namespace gemmi
#endif
