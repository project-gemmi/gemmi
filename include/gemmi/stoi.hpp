// Copyright 2018 Global Phasing Ltd.
//
// string to integer

#ifndef GEMMI_STOI_HPP_
#define GEMMI_STOI_HPP_

#include <stdexcept>  // for invalid_argument
#include <string>

namespace gemmi {

// equivalent of std::isspace for C locale
inline bool isspace_c(char c) {
  return c == ' ' || c == '\t' || c == '\n'|| c == '\f' ||
         c == '\v' || c == '\r' ;
}

// no checking for overflow
inline int string_to_int(const char* p, bool checked, size_t length=0) {
  int mult = -1;
  int n = 0;
  size_t i = 0;
  while ((length == 0 || i < length) && isspace_c(p[i]))
    ++i;
  if (p[i] == '-') {
    mult = 1;
    ++i;
  } else if (p[i] == '+') {
    ++i;
  }
  bool has_digits = false;
  // use negative numbers because INT_MIN < -INT_MAX
  for (; (length == 0 || i < length) && p[i] >= '0' && p[i] <= '9'; ++i) {
    n = n * 10 - (p[i] - '0');
    has_digits = true;
  }
  if (checked) {
    while ((length == 0 || i < length) && isspace_c(p[i]))
      ++i;
    if (!has_digits || p[i] != '\0')
      throw std::invalid_argument("not an integer: " +
                                   std::string(p, length ? length : i+1));
  }
  return mult * n;
}

inline int string_to_int(const std::string& str, bool checked) {
  return string_to_int(str.c_str(), checked);
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
