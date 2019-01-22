// Copyright 2018 Global Phasing Ltd.
//
// Locale-independent C-string to double conversion.

#ifndef GEMMI_ATOF_HPP_
#define GEMMI_ATOF_HPP_

namespace gemmi {

// no checking for overflow
// no support for scientific notation
inline double simple_atof(const char* p, const char** endptr=nullptr) {
  // equivalent of std::isspace for C locale
  while (*p  == ' ' || *p == '\t' || *p == '\n'||
         *p == '\f' || *p == '\v' || *p == '\r')
    ++p;
  int sign = 1;
  if (*p == '-') {
    ++p;
    sign = -1;
  } else if (*p == '+') {
    ++p;
  }
  double d = 0;
  for (; *p >= '0' && *p <= '9'; ++p)
    d = d * 10 + (*p - '0');
  if (*p == '.')
    for (double mult = 0.1; *++p >= '0' && *p <= '9'; mult *= 0.1)
      d += mult * (*p - '0');
  if (endptr)
    *endptr = p;
  return sign * d;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
