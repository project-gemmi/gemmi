// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format.

#ifndef GEMMI_TO_PDB_HH_
#define GEMMI_TO_PDB_HH_

#include <cctype>
#include <cstring>
#include <algorithm>
#include <ostream>
#ifdef USE_STD_SNPRINTF
# include <cstdio>
#else
# include <stb_sprintf.h>
#endif
#include "model.hh"

namespace gemmi {
namespace mol {

#ifdef USE_STD_SNPRINTF  // for benchmarking and testing only
#define stbsp_snprintf  std::snprintf
#endif

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
// 1-6 record name, 8-10 continuation, 11-lastcol string.
inline void write_multiline(std::ostream& os, const char* record_name,
                            const char* text, int lastcol) {
  if (text == nullptr)
    return;
  char buf[88]; // a few bytes extra, just in case
  const char *end = find_last_break(text, lastcol-10);
  WRITEU("%-6s    %-70.*s\n", record_name, static_cast<int>(end-text), text);
  for (int n = 2; n < 1000 && *end != '\0'; ++n) {
    const char *start = end;
    end = find_last_break(start, lastcol-11);
    int len = end - start;
    WRITEU("%-6s %3d %-69.*s\n", record_name, n, len, start);
  }
}

inline void write_pdb(const Structure& st, std::ostream& os) {
  const char* months = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC???";
  char buf[88];

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
  write_multiline(os, "TITLE", st.get_info("_struct.title"), 80);
  write_multiline(os, "KEYWDS", st.get_info("_struct_keywords.text"), 79);
  if (st.models.size() > 1)
    WRITE("NUMMDL    %-6jd %63s\n", st.models.size(), "");
  // TODO: SEQRES
  const UnitCell& cell = st.cell;
  WRITE("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s          \n",
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma,
        st.sg_hm.empty() ? "P 1" : st.sg_hm.c_str(),
        st.get_info("_cell.Z_PDB", "1"));
  // We add a small number to avoid negative 0.
  WRITE("SCALE1 %13.6f%10.6f%10.6f %14.5f %24s\n",
        cell.frac.a11+1e-15, cell.frac.a12+1e-15, cell.frac.a13+1e-15, 0.0, "");
  WRITE("SCALE2 %13.6f%10.6f%10.6f %14.5f %24s\n",
        cell.frac.a21+1e-15, cell.frac.a22+1e-15, cell.frac.a23+1e-15, 0.0, "");
  WRITE("SCALE3 %13.6f%10.6f%10.6f %14.5f %24s\n",
        cell.frac.a31+1e-15, cell.frac.a32+1e-15, cell.frac.a33+1e-15, 0.0, "");

  for (size_t i = 0; i != st.ncs.size(); i++) {
    const NcsOp& op = st.ncs[i];
    char g = op.given ? '1' : ' ';
    WRITE("MTRIX%d %3jd%10.6f%10.6f%10.6f %14.5f    %-21c\n",
          1, i+1, op.rot.a11, op.rot.a12, op.rot.a13, op.tran.x, g);
    WRITE("MTRIX%d %3jd%10.6f%10.6f%10.6f %14.5f    %-21c\n",
          2, i+1, op.rot.a21, op.rot.a22, op.rot.a23, op.tran.y, g);
    WRITE("MTRIX%d %3jd%10.6f%10.6f%10.6f %14.5f    %-21c\n",
          3, i+1, op.rot.a31, op.rot.a32, op.rot.a33, op.tran.z, g);
  }
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
        bool standard = chain.entity_type != EntityType::NonPolymer &&
                        res.has_standard_pdb_name();
        for (const mol::Atom& a : res.atoms) {
          if (serial == 1000000)
            throw std::runtime_error("Too many atoms for PDB file.");
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
                a.altloc ? std::toupper(a.altloc) : ' ',
                res.name.c_str(),
                chain.auth_name.c_str(),
                res.seq_id_for_pdb(),
                res.ins_code ? res.ins_code : ' ',
                // We want to avoid negative zero and round them numbers up
                // if they originally had one digit more and that digit was 5.
                a.pos.x > -5e-4 && a.pos.x < 0 ? 0 : a.pos.x + 1e-10,
                a.pos.y > -5e-4 && a.pos.y < 0 ? 0 : a.pos.y + 1e-10,
                a.pos.z > -5e-4 && a.pos.z < 0 ? 0 : a.pos.z + 1e-10,
                // Occupancy is stored as single prec, but we know it's <= 1,
                // so no precision is lost even if it had 6 digits after dot.
                a.occ + 1e-6,
                // B is harder to get rounded right. It is stored as float,
                // and may be given with more than single precision in mmCIF
                // If it was originally %.5f (5TIS) we need to add 0.5 * 10^-5.
                a.b_iso + 0.5e-5,
                a.element.uname(),
                // Charge is written as 1+ or 2-, etc, or just empty space.
                // Sometimes PDB files have explicit 0s (5M05); we ignore them.
                a.charge ? a.charge > 0 ? '0'+a.charge : '0'-a.charge : ' ',
                a.charge ? a.charge > 0 ? '+' : '-' : ' ');
          if (a.u11 != 0.0f) {
            // re-using part of the buffer
            memcpy(buf, "ANISOU", 6);
            stbsp_snprintf(buf+28, 43, "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
                           a.u11*1e4 + 1e-6, a.u22*1e4 + 1e-6, a.u33*1e4 + 1e-6,
                           a.u12*1e4 + 1e-6, a.u13*1e4 + 1e-6, a.u23*1e4 + 1e-6);
            buf[28+42] = ' ';
            os.write(buf, 81);
          }
        }
      }
      // TODO: proper detection when polymer ends
      if (chain.residues.back().seq_id != mol::Residue::UnknownId) {
        // re-using part of the buffer in the middle, e.g.:
        // TER    4153      LYS B 286
        stbsp_snprintf(buf, 82, "TER   %5d", ++serial);
        std::memset(buf+11, ' ', 6);
        std::memset(buf+28, ' ', 52);
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
