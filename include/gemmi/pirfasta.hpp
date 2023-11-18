// Copyright 2020 Global Phasing Ltd.
//
// Read sequence from PIR or (multi-)FASTA format.

#ifndef GEMMI_PIRFASTA_HPP_
#define GEMMI_PIRFASTA_HPP_

#include <cctype>   // for isspace
#include <vector>
#include <algorithm>  // for min
#include "fail.hpp"

namespace gemmi {

struct FastaSeq {
  std::string header;
  std::string seq;
};

// PIR format starts with one of: >P1; >F1; >DL; >DC; >RL; >RC; >XX;
inline bool is_pir_format(const std::string& s) {
  return s.length() > 4 && s[0] == '>' && s[3] == ';' && (
        ((s[1] == 'P' || s[1] == 'F') && s[2] == '1') ||
        ((s[1] == 'D' || s[1] == 'R') && (s[2] == 'L' || s[2] == 'C')) ||
        (s[1] == 'X' && s[2] == 'X'));
}

inline std::vector<FastaSeq> read_pir_or_fasta(const std::string& str) {
  if (str[0] != '>')
    fail("PIR/FASTA files start with '>'");
  bool pir = is_pir_format(str);
  std::vector<FastaSeq> r;
  int blank_lines = 0;
  int paren_level = 0;
  bool ended = false;
  for (size_t pos=0, end=0; end != std::string::npos; pos = end + 1) {
    end = str.find('\n', pos);
    if (str[pos] == '>') {
      ended = false;
      if (paren_level != 0)
        break;
      r.emplace_back();
      if (pir && end != std::string::npos)
        end = str.find('\n', end+1);
      r.back().header = str.substr(pos+1, end-(pos+1));
    } else {
      std::string& seq = r.back().seq;
      for (size_t i = pos; i < std::min(end, str.size()); ++i) {
        char c = str[i];
        if (std::isspace(c)) {
          if (c == '\n')
            ++blank_lines;
          continue;
        }
        // handle non-blank characters
        if (ended)
          fail("'*' is interpreted as sequence terminator, here it is followed by: ", c);
        if (blank_lines >= 2)
          fail("blank lines can be followed only by a line starting with '>'");
        blank_lines = 0;
        if (('a' <= (c | 0x20) && (c | 0x20) <= 'z') || c == '-' ||
            (paren_level != 0 && '0' <= c && c <= '9')) {
          // good character, nothing to be done here
        } else if (c == '*') {
          ended = true;
          continue;
        } else if (c == '(') {
          if (++paren_level > 1)
            fail("nested parentheses are not allowed");
        } else if (c == ')') {
          if (--paren_level < 0)
            fail("')' without matching '('");
        } else {
          fail("unexpected character in sequence: ", c);
        }
        seq += c;
      }
    }
  }
  if (paren_level != 0)
    fail("unmatched '('");
  return r;
}

} // namespace gemmi
#endif
