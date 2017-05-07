// Copyright 2017 Global Phasing Ltd.
//
// Utilities.

#ifndef GEMMI_UTIL_HH_
#define GEMMI_UTIL_HH_

#include <string>

namespace gemmi {

bool ends_with(const std::string& str, const std::string& suffix) {
  size_t sl = suffix.length();
  return str.length() >= sl && str.compare(str.length() - sl, sl, suffix) == 0;
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
