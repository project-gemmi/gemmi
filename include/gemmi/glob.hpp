//! @file
//! @brief Glob pattern matching.

// Copyright Global Phasing Ltd.
//
// Glob pattern matching

#ifndef GEMMI_GLOB_HPP_
#define GEMMI_GLOB_HPP_

#include <string>

namespace gemmi {

//! @brief Match string against glob pattern.
//! @param pattern Glob pattern ('*' and '?' wildcards)
//! @param str String to match
//! @return True if string matches pattern
//!
//! Linear-time glob matching: https://research.swtch.com/glob
inline bool glob_match(const std::string& pattern, const std::string& str) {
  size_t pat_next = 0;
  size_t str_next = std::string::npos;
  size_t pat_pos = 0;
  size_t str_pos = 0;
  while (pat_pos < pattern.size() || str_pos < str.size()) {
    if (pat_pos < pattern.size()) {
      char c = pattern[pat_pos];
      if (c == '*') {
        pat_next = pat_pos;
        str_next = str_pos + 1;
        pat_pos++;
        continue;
      }
      if (str_pos < str.size() && (c == '?' || c == str[str_pos])) {
        pat_pos++;
        str_pos++;
        continue;
      }
    }
    if (str_next > str.size())
      return false;
    pat_pos = pat_next;
    str_pos = str_next;
  }
  return true;
}

} // namespace gemmi

#endif
