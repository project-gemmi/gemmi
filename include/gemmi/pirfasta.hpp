// Copyright 2020 Global Phasing Ltd.
//
// Read sequences from PIR or (multi-)FASTA formats.

/// @file
/// @brief Parsing of sequence data in PIR and FASTA formats.

#ifndef GEMMI_PIRFASTA_HPP_
#define GEMMI_PIRFASTA_HPP_

#include <cctype>   // for isspace
#include <vector>
#include <algorithm>  // for min
#include "fail.hpp"

namespace gemmi {

/// @brief A sequence record with header and body.
///
/// Represents a single sequence entry from a FASTA or PIR file.
struct FastaSeq {
  /// @brief Sequence header/description line (without the '>' prefix).
  std::string header;
  /// @brief The sequence itself (letters and allowed special characters).
  std::string seq;
};

/// @brief Determine whether a file starts in PIR format.
///
/// PIR format is identified by a sequence code at position 1-3 of the first line,
/// followed by a semicolon at position 3:
/// - >P1; — Protein
/// - >F1; — Protein (structure known)
/// - >DL; — DNA
/// - >DC; — DNA (structure known)
/// - >RL; — RNA
/// - >RC; — RNA (structure known)
/// - >XX; — Generic/Unknown
///
/// FASTA format uses '>something' (arbitrary text after '>').
///
/// @param s The first line of a sequence file (expected to start with '>').
/// @return true if the line matches PIR format syntax, false otherwise.
inline bool is_pir_format(const std::string& s) {
  return s.length() > 4 && s[0] == '>' && s[3] == ';' && (
        ((s[1] == 'P' || s[1] == 'F') && s[2] == '1') ||
        ((s[1] == 'D' || s[1] == 'R') && (s[2] == 'L' || s[2] == 'C')) ||
        (s[1] == 'X' && s[2] == 'X'));
}

/// @brief Parse a string containing PIR or FASTA format sequences.
///
/// Supports both FASTA and PIR formats:
/// - **FASTA format:** Lines start with '>' followed by a description, then sequence lines.
/// - **PIR format:** Like FASTA, but has a second header line immediately after the first.
///   The sequence starts on the second line after the primary header.
///
/// Sequence parsing rules:
/// - Alphabetic characters (A-Z, case-insensitive) are included in the sequence.
/// - Hyphens ('-') represent gaps/insertions.
/// - Asterisks ('*') mark sequence terminators; content after '*' is not allowed.
/// - Parentheses '(...)' are allowed (for numbering in PIR); nesting is not allowed.
/// - Digits are allowed only within parentheses.
/// - Whitespace is ignored.
/// - Two or more blank lines separate records.
///
/// @param str The complete file content as a string.
/// @return A vector of FastaSeq records, each with header and sequence.
/// @throws gemmi::fail() If the input is malformed (no leading '>', unexpected
///   characters, unmatched parentheses, '*' not at end, etc.).
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
