// Copyright 2017 Global Phasing Ltd.
//
// Read PDB format into a Structure from model.hh.
// https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

#ifndef GEMMI_READ_PDB_HH_
#define GEMMI_READ_PDB_HH_

#include <string>
#include <cstdio>
#include <memory>
#include "model.hh"

namespace gemmi {
namespace mol {

inline int read_pdb_int(const char* p, int field_length) {
  (void) p;
  (void) field_length;
  //TODO
  return 0;
}

inline double read_pdb_number(const char* p, int field_length) {
  (void) p;
  (void) field_length;
  //TODO
  return 0.;
}

inline std::string read_pdb_string(const char* p, int field_length) {
  (void) p;
  (void) field_length;
  //TODO trim
  return std::string(p, field_length);
}

// Compare the first 4 letters of s, ignoring case, with uppercase record.
// Both args must have at least 3+1 chars. ' ' and NUL are equivalent in s.
inline bool is_record_type(const char* s, const char* record) {
  return ((s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020) ==
          (record[0] << 24 | record[1] << 16 | record[2] << 8 | record[3]);
}

inline void read_pdb_line(const char* line, Structure* st) {
  if (is_record_type(line, "ATOM") || is_record_type(line, "HETATM")) {
  } else if (is_record_type(line, "ANISOU")) {
  } else if (is_record_type(line, "REMARK")) {
  } else if (is_record_type(line, "CONECT")) {

  } else if (is_record_type(line, "HEADER")) {
  } else if (is_record_type(line, "TITLE")) {
    printf("title %s", line);
  } else if (is_record_type(line, "CRYST1")) {
  } else if (is_record_type(line, "MTRIXn")) {
  } else if (is_record_type(line, "MODEL")) {
  //} else if (is_record_type(line, "ENDMDL")) {
  //} else if (is_record_type(line, "TER")) {
  //} else if (is_record_type(line, "SCALEn")) {
  } else if (is_record_type(line, "END")) {  // NUL == ' ' & ~0x20
  }
}

inline Structure read_pdb_from_cstream(FILE* f) {
  Structure st;
  char buf[88] = {0};
  while(fgets(buf, 82, f))
    read_pdb_line(buf, &st);
  return st;
}

inline Structure read_pdb(const std::string& path) {
  std::unique_ptr<FILE, decltype(&std::fclose)> f(std::fopen(path.c_str(), "r"),
                                                  &std::fclose);
  if (!f)
    throw std::runtime_error("Failed to open file: " + path);
  return read_pdb_from_cstream(f.get());
}

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
