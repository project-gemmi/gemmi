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
#include <iostream> // temporary

namespace gemmi {
namespace mol {

class CstreamLineInput {
public:
  CstreamLineInput(FILE* f) : f_(f) {}

  size_t copy_line(char* line) {
    if (!fgets(line, 82, f_))
      return 0;
    size_t len = std::strlen(line);
    // We don't expect lines longer than 80 characters, but if one is found,
    // just discard the rest of the line.
    if (len > 0 && line[len-1] != '\n')
      for (int c = fgetc(f_); c != 0 && c != '\n'; c = fgetc(f_))
        continue;
    return len;
  }
private:
  FILE* f_;
  size_t line_number_ = 0;
};

inline std::string rtrimmed(std::string s) {
  auto p = std::find_if_not(s.rbegin(), s.rend(),
                            [](int c) { return std::isspace(c); });
  s.erase(p.base(), s.end());
  return s;
}

inline int read_pdb_int(const char* p, int field_length) {
  int sign = 1;
  int n = 0;
  int i = 0;
  while (std::isspace(p[i]) && i < field_length)
    ++i;
  if (p[i] == '-') {
    ++i;
    sign = -1;
  } else if (p[i] == '+') {
    ++i;
  }
  for (; p[i] >= '0' && p[i] <= '9' && i < field_length; ++i) {
    n = n * 10 + (p[i] - '0');
  }
  return sign * n;
}

inline double read_pdb_number(const char* p, int field_length) {
  int sign = 1;
  double d = 0;
  int i = 0;
  while (std::isspace(p[i]) && i < field_length)
    ++i;
  if (p[i] == '-') {
    ++i;
    sign = -1;
  } else if (p[i] == '+') {
    ++i;
  }
  for (; i < field_length && p[i] >= '0' && p[i] <= '9'; ++i)
    d = d * 10 + (p[i] - '0');
  if (i < field_length && p[i] == '.') {
    double mult = 0.1;
    for (++i; i < field_length && p[i] >= '0' && p[i] <= '9'; ++i, mult *= 0.1)
      d += mult * (p[i] - '0');
  }
  return sign * d;
}

inline std::string read_pdb_string(const char* p, int field_length) {
  // left trim
  while (std::isspace(*p)) {
    ++p;
    --field_length;
    if (field_length == 0)
      break;
  }
  // right trim
  while (field_length != 0 && std::isspace(p[field_length-1]))
    --field_length;

  return std::string(p, field_length);
}

// Compare the first 4 letters of s, ignoring case, with uppercase record.
// Both args must have at least 3+1 chars. ' ' and NUL are equivalent in s.
inline bool is_record_type(const char* s, const char* record) {
  return ((s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020) ==
          (record[0] << 24 | record[1] << 16 | record[2] << 8 | record[3]);
}

template<typename InputType> Structure read_pdb_from_input(InputType&& in) {
  Structure st;
  char line[88] = {0};
  size_t len;
  while ((len = in.copy_line(line))) {
    if (is_record_type(line, "ATOM") || is_record_type(line, "HETATM")) {
    } else if (is_record_type(line, "ANISOU")) {
    } else if (is_record_type(line, "REMARK")) {
    } else if (is_record_type(line, "CONECT")) {
    } else if (is_record_type(line, "HEADER")) {
      if (len > 50)
        st.info["_struct_keywords.pdbx_keywords"] =
                                    rtrimmed(std::string(line+10, 40));
      if (len > 59) { // date in PDB has format 28-MAR-07
        std::string date(line+50, 9);
        const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
                              "JUL07AUG08SEP09OCT10NOV11DEC122222";
        const char* m = strstr(months, date.substr(3, 3).c_str());
        st.info["_database_PDB_rev.date_original"] =
          (date[7] > '6' ? "19" : "20") + date.substr(7, 2) + "-" +
          (m ? std::string(m+3, 2) : "??") + "-" + date.substr(0, 2);
      }
      if (len > 66)
        st.info["_entry.id"] = std::string(line+62, 4);

    } else if (is_record_type(line, "TITLE")) {
      if (len > 10)
        st.info["_struct.title"] += rtrimmed(std::string(line+10, len-10-1));
    } else if (is_record_type(line, "KEYWDS")) {
      if (len > 10)
        st.info["_struct_keywords.text"] +=
                                    rtrimmed(std::string(line+10, len-10-1));
    } else if (is_record_type(line, "CRYST1")) {
      if (len > 54) {
        st.cell.a = read_pdb_number(line+6, 9);
        st.cell.b = read_pdb_number(line+15, 9);
        st.cell.c = read_pdb_number(line+24, 9);
        st.cell.alpha = read_pdb_number(line+33, 7);
        st.cell.beta = read_pdb_number(line+40, 7);
        st.cell.gamma = read_pdb_number(line+47, 7);
        st.cell.calculate_matrices();
      }
      if (len > 56)
        st.sg_hm = read_pdb_string(line+55, 11);
      if (len > 67)
        st.info["_cell.Z_PDB"] = read_pdb_string(line+66, 4);
    } else if (is_record_type(line, "MTRIXn")) {
    } else if (is_record_type(line, "MODEL")) {
    //} else if (is_record_type(line, "ENDMDL")) {
    //} else if (is_record_type(line, "TER")) {
    //} else if (is_record_type(line, "SCALEn")) {
    } else if (is_record_type(line, "END")) {  // NUL == ' ' & ~0x20
      break;
    }
  }
  return st;
}

inline Structure read_pdb_from_cstream(FILE* f) {
  CstreamLineInput input(f);
  return read_pdb_from_input(input);
}

/*
inline Structure read_pdb_from_istream(std::istream &is) {
  IstreamLineInput input(is);
  return read_pdb_from_input(input);
}

inline Structure read_pdb_from_memory(const char* data, size_t size) {
  MemoryLineInput input(data, size);
  return read_pdb_from_input(input);
}
*/


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
