// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format.

#ifndef GEMMI_TO_PDB_HH_
#define GEMMI_TO_PDB_HH_

#include <cstring>
#include <algorithm>
#include <ostream>
#include <stb_sprintf.h>
#include "model.hh"

namespace gemmi {
namespace mol {

#define WRITE(...) do { \
    stbsp_snprintf(buf, 82, __VA_ARGS__); \
    os.write(buf, 81); \
    } while(0)

#define WRITEU(...) do { \
    stbsp_snprintf(buf, 82, __VA_ARGS__); \
    for (int i_ = 0; i_ != 80; i_++) \
      if (buf[i_] >= 'a' && buf[i_] <= 'z') buf[i_] -= 0x20; \
    os.write(buf, 81); \
    } while(0)

const char* find_last_break(const char *str, int max_len) {
  int last_break = 0;
  for (int i = 0; i < max_len; i++) {
    if (str[i] == '\0')
      return str + i;
    if (str[i] == ' ' || str[i] == '-')
      last_break = i + 1;
  }
  return str + (last_break != 0 ? last_break : max_len);
}

// Write record with possible continuation lines, with the format:
// 1-6 record name, 8-10 continuation, 11-79 string.
inline void write_multiline(std::ostream& os, const char* record_name,
                            const char* text) {
  if (text == nullptr)
    return;
  char buf[83];
  const char *end = find_last_break(text, 70);
  WRITEU("%-6s    %-70.*s\n", record_name, end-text, text);
  for (int n = 2; n < 1000 && *end != '\0'; ++n) {
    const char *start = end;
    end = find_last_break(start, 69);
    WRITEU("%-6s %3d %-69.*s\n", record_name, n, end-start, start);
  }
}

inline double avoid_neg_zero(double x, double prec) {
  return x > 0 || x < -0.5 * prec ? x : 0.0;
}

inline void write_pdb(const Structure& st, std::ostream& os) {
  const char* months = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC???";
  char buf[83] = {0};

  const char* date = st.get_info("_database_PDB_rev.date_original");
  std::string pdb_date;
  if (date && std::strlen(date) == 10) {
    unsigned month_idx = 10 * (date[5] - '0') + date[6] - '0' - 1;
    pdb_date = std::string(date + 8, 2) + "-" +
               std::string(months + 3 * std::min(month_idx, 13u), 3) +
               "-" + std::string(date + 2, 2);
  }
  WRITEU("HEADER    %-40s%-9s   %-18s\n",
         // "classification" in PDB == _struct_keywords.pdbx_keywords in mmCIF
         st.get_info("_struct_keywords.pdbx_keywords", ""),
         pdb_date.c_str(), st.get_info("_entry.id"));
  write_multiline(os, "TITLE", st.get_info("_struct.title"));
  if (st.models.size() > 1)
    WRITE("NUMMDL    %-6jd %63s\n", st.models.size(), "");
  // TODO: SEQRES
  WRITE("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s          \n",
        st.cell.a, st.cell.b, st.cell.c,
        st.cell.alpha, st.cell.beta, st.cell.gamma,
        st.sg_hm.empty() ? "P 1" : st.sg_hm.c_str(),
        st.get_info("_cell.Z_PDB", "1"));
  // TODO: SCALE
  //TODO: special handling of large structures (>62 chains or >=1M atoms)
  for (const mol::Model& model : st.models) {
    int serial = 0;
    if (st.models.size() > 1)
      WRITE("MODEL %8s %65s\n", model.name.c_str(), "");
    for (const mol::Chain& chain : model.chains) {
      if (chain.auth_name.empty())
        throw std::runtime_error("empty chain name");
      if (chain.auth_name.length() > 1)
        throw std::runtime_error("long chain name: " + chain.auth_name);
      for (const mol::Residue& res : chain.residues) {
        bool standard = res.has_standard_pdb_name();
        for (const mol::Atom& a : res.atoms) {
          //  1- 6  6s  record name
          //  7-11  5d  integer serial
          // 12     1   -
          // 13-16  4s  atom name (from 13 only if 4-char or 2-char symbol)
          // 17     1c  altloc
          // 18-20  3s  residue name
          // 21     1   -
          // 22     1s  chain
          // 23-26  4d  integer residue sequence number
          // 27     1c  insertion code
          // 28-30  3   - 
          // 31-38  8f  x (8.3)
          // 39-46  8f  y
          // 47-54  8f  z
          // 55-60  6f  occupancy (6.2)
          // 61-66  6f  temperature factor (6.2)
          // 67-76 10   -
          // 77-78  2s  element symbol, right-justified
          // 79-80  2s  charge
          bool empty13 = (a.element.uname()[1] == '\0' && a.name.size() < 4);
          WRITE("%-6s%5d %c%-3s%c%3s"
                " %1s%4d%c"
                "   %8.3f%8.3f%8.3f"
                "%6.2f%6.2f          %2s%c%c\n",
                standard ? "ATOM" : "HETATM",
                ++serial,
                empty13 ? ' ' : a.name[0],
                a.name.c_str() + (empty13 || a.name.empty() ? 0 : 1),
                a.altloc ? a.altloc : ' ',
                res.name.c_str(),
                chain.auth_name.c_str(),
                res.seq_id_for_pdb(),
                res.ins_code ? res.ins_code : ' ',
                avoid_neg_zero(a.x, 1e-3),
                avoid_neg_zero(a.y, 1e-3),
                avoid_neg_zero(a.z, 1e-3),
                a.occ, a.b_iso,
                a.element.uname(),
                // charge is written as 1+ or 2-, etc, or just empty space
                a.charge ? a.charge > 0 ? '0'+a.charge : '0'-a.charge : ' ',
                a.charge ? a.charge > 0 ? '+' : '-' : ' ');
          if (a.u11 != 0.0f) {
            // re-using part of the buffer
            memcpy(buf, "ANISOU", 6);
            stbsp_snprintf(buf+28, 43, "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
                           avoid_neg_zero(a.u11*1e4, 0),
                           avoid_neg_zero(a.u22*1e4, 0),
                           avoid_neg_zero(a.u33*1e4, 0),
                           avoid_neg_zero(a.u12*1e4, 0),
                           avoid_neg_zero(a.u13*1e4, 0),
                           avoid_neg_zero(a.u23*1e4, 0));
            buf[28+42] = ' ';
            os.write(buf, 81);
          }
        }
      }
      // TODO: proper detection when polymer ends
      if (chain.residues.back().seq_id != mol::Residue::UnknownId) {
        // re-using part of the buffer in the middle, e.g.:
        // TER    4153      LYS B 286
        stbsp_sprintf(buf, "TER   %5d", ++serial);
        std::memset(buf+11, ' ', 6);
        std::memset(buf+30, ' ', 50);
        os.write(buf, 81);
      }
    }
    if (st.models.size() > 1)
      WRITE("%-80s\n", "ENDMDL");
  }
  WRITE("%-80s\n", "END");
}

#undef WRITE

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
