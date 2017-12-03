// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format (Structure -> pdb file).

#ifndef GEMMI_TO_PDB_HPP_
#define GEMMI_TO_PDB_HPP_

#include <cassert>
#include <cctype>
#include <cstring>
#include <algorithm>
#include <ostream>
#include "sprintf.hpp"
#include "model.hpp"
#include "util.hpp"

namespace gemmi {

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

namespace impl {

// works for non-negative values only
inline char *base36_encode(char* buffer, int width, int value) {
  const char base36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  buffer[width] = '\0';
  do {
    buffer[--width] = base36[value % 36];
    value /= 36;
  } while (value != 0 && width != 0);
  while (width != 0)
    buffer[--width] = ' ';
  return buffer;
}

// based on http://cci.lbl.gov/hybrid_36/
inline char* encode_serial_in_hybrid36(char* str, int serial) {
  assert(serial >= 0);
  if (serial < 100000) {
    stbsp_sprintf(str, "%5d", serial);
    return str;
  }
  return base36_encode(str, 5, serial - 100000 + 10 * 36 * 36 * 36 * 36);
}

// based on http://cci.lbl.gov/hybrid_36/
inline char* encode_seq_id_in_hybrid36(char* str, int seq_id) {
  if (seq_id > -1000 && seq_id < 10000) {
    stbsp_sprintf(str, "%4d", seq_id);
    return str;
  }
  return base36_encode(str, 4, seq_id - 10000 + 10 * 36 * 36 * 36);
}

inline char* write_seq_id(char* str, const Residue& res) {
  encode_seq_id_in_hybrid36(str, res.seq_id_for_pdb());
  str[4] = res.snic.ins_code ? res.snic.ins_code : ' ';
  str[5] = '\0';
  return str;
}

inline const char* find_last_break(const char *str, int max_len) {
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

} // namespace impl

inline void write_pdb(const Structure& st, std::ostream& os,
                      bool iotbx_compat=false) {
  const char* months = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC???";
  char buf[88];
  char buf8[8];
  char buf8a[8];

  const char* date =
    st.get_info("_pdbx_database_status.recvd_initial_deposition_date");
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
         pdb_date.c_str(), st.get_info("_entry.id", ""));
  impl::write_multiline(os, "TITLE", st.get_info("_struct.title"), 80);
  impl::write_multiline(os, "KEYWDS", st.get_info("_struct_keywords.text"), 79);
  impl::write_multiline(os, "EXPDTA", st.get_info("_exptl.method"), 79);
  if (st.models.size() > 1)
    WRITE("NUMMDL    %-6zu %63s\n", st.models.size(), "");

  if (!st.models.empty() && !iotbx_compat) {
    // SEQRES
    for (const Chain& ch : st.models[0].chains)
      if (ch.entity && ch.entity->type == EntityType::Polymer) {
        const std::string& chain_name = ch.name_for_pdb();
        int seq_len = 0;
        int prev_seq_num = -1;
        for (const SequenceItem& si : ch.entity->sequence)
          if (si.num < 0 || si.num != prev_seq_num) {
            ++seq_len;
            prev_seq_num = si.num;
          }
        prev_seq_num = -1;
        int row = 0;
        int col = 0;
        for (const SequenceItem& si : ch.entity->sequence) {
          if (si.num >= 0 && si.num == prev_seq_num)
            continue;
          prev_seq_num = si.num;
          if (col == 0)
            stbsp_snprintf(buf, 82, "SEQRES%4d%2s%5d %62s\n",
                           ++row, chain_name.c_str(), seq_len, "");
          const std::string& res = si.mon;
          memcpy(buf + 18 + 4*col + 4-res.length(), res.c_str(), res.length());
          if (++col == 13) {
            os.write(buf, 81);
            col = 0;
          }
        }
        if (col != 0)
          os.write(buf, 81);
      }

    // CISPEP (note: we use only the first conformation)
    int counter = 0;
    for (const Model& model : st.models)
      for (const Chain& chain : model.chains) {
        const char* cname = chain.name_for_pdb().c_str();
        for (const Residue& res : chain.residues)
          if (res.is_cis)
            if (const Residue* next = res.next_bonded_aa()) {
              WRITE("CISPEP%4d %3s%2s %5s   %3s%2s %5s %9s %12.2f %20s\n",
                  ++counter,
                  res.name.c_str(), cname, impl::write_seq_id(buf8, res),
                  next->name.c_str(), cname, impl::write_seq_id(buf8a, *next),
                  st.models.size() > 1 ? model.name.c_str() : "0",
                  res.calculate_omega(*next) * (180. / 3.14159265358979323846),
                  "");
            }
      }
  }

  const UnitCell& cell = st.cell;
  WRITE("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s          \n",
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma,
        st.sg_hm.empty() ? "P 1" : st.sg_hm.c_str(),
        st.get_info("_cell.Z_PDB", ""));
  for (int i = 0; i < 3; ++i)
    WRITE("ORIGX%d %13.6f%10.6f%10.6f %14.5f %24s\n",
          i+1, st.origx.x[i], st.origx.y[i], st.origx.z[i], st.origx.w[i], "");
  const Matrix33& frac = cell.frac;
  // We add a small number to avoid negative 0.
  WRITE("SCALE1 %13.6f%10.6f%10.6f %14.5f %24s\n",
        frac.a11+1e-15, frac.a12+1e-15, frac.a13+1e-15, cell.shift.x+1e-15, "");
  WRITE("SCALE2 %13.6f%10.6f%10.6f %14.5f %24s\n",
        frac.a21+1e-15, frac.a22+1e-15, frac.a23+1e-15, cell.shift.y+1e-15, "");
  WRITE("SCALE3 %13.6f%10.6f%10.6f %14.5f %24s\n",
        frac.a31+1e-15, frac.a32+1e-15, frac.a33+1e-15, cell.shift.z+1e-15, "");

  for (const NcsOp& op : st.ncs)
    for (int j = 0; j < 3; ++j) {
      auto r = op.transform.row(j);
      WRITE("MTRIX%d %3.3s%10.6f%10.6f%10.6f %14.5f    %-21c\n",
            j + 1, op.id.c_str(), r.x, r.y, r.z, r.w, op.given ? '1' : ' ');
    }

  for (const Model& model : st.models) {
    int serial = 0;
    if (st.models.size() > 1)
      WRITE("MODEL %8s %65s\n", model.name.c_str(), "");
    for (const Chain& chain : model.chains) {
      if (chain.force_pdb_serial)
        serial = chain.force_pdb_serial - 1;
      const std::string& chain_name = chain.name_for_pdb();
      if (chain_name.empty())
        gemmi::fail("empty chain name");
      if (chain_name.length() > 2)
        gemmi::fail("long chain name: " + chain_name);
      for (const Residue& res : chain.residues) {
        bool standard = res.get_info().pdb_standard && !(chain.entity &&
                                 chain.entity->type == EntityType::NonPolymer);
        for (const Atom& a : res.atoms) {
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
          // 67-76  6   -
          // 73-76      segment identifier, left-justified (non-standard)
          // 77-78  2s  element symbol, right-justified
          // 79-80  2s  charge
          bool empty13 = (a.element.uname()[1] == '\0' && a.name.size() < 4);
          WRITE("%-6s%5s %c%-3s%c%3s"
                "%2s%4s%c"
                "   %8.3f%8.3f%8.3f"
                "%6.2f%6.2f      %-4.4s%2s%c%c\n",
                standard ? "ATOM" : "HETATM",
                impl::encode_serial_in_hybrid36(buf8, ++serial),
                empty13 ? ' ' : a.name[0],
                a.name.c_str() + (empty13 || a.name.empty() ? 0 : 1),
                a.altloc ? std::toupper(a.altloc) : ' ',
                res.name.c_str(),
                chain_name.c_str(),
                impl::encode_seq_id_in_hybrid36(buf8a, res.seq_id_for_pdb()),
                res.snic.ins_code ? res.snic.ins_code : ' ',
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
                res.segment.c_str(),
                a.element.uname(),
                // Charge is written as 1+ or 2-, etc, or just empty space.
                // Sometimes PDB files have explicit 0s (5M05); we ignore them.
                a.charge ? a.charge > 0 ? '0'+a.charge : '0'-a.charge : ' ',
                a.charge ? a.charge > 0 ? '+' : '-' : ' ');
          if (a.u11 != 0.0f) {
            // re-using part of the buffer
            memcpy(buf, "ANISOU", 6);
            const double eps = 1e-6;
            stbsp_snprintf(buf+28, 43, "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
                           a.u11*1e4 + eps, a.u22*1e4 + eps, a.u33*1e4 + eps,
                           a.u12*1e4 + eps, a.u13*1e4 + eps, a.u23*1e4 + eps);
            buf[28+42] = ' ';
            os.write(buf, 81);
          }
        }
      }
      if (chain.entity && chain.entity->type == EntityType::Polymer) {
        if (iotbx_compat) {
          WRITE("%-80s\n", "TER");
        } else {
          // re-using part of the buffer in the middle, e.g.:
          // TER    4153      LYS B 286
          stbsp_snprintf(buf, 82, "TER   %5s",
                         impl::encode_serial_in_hybrid36(buf8, ++serial));
          std::memset(buf+11, ' ', 6);
          std::memset(buf+28, ' ', 52);
          os.write(buf, 81);
        }
      }
    }
    if (st.models.size() > 1)
      WRITE("%-80s\n", "ENDMDL");
  }
  WRITE("%-80s\n", "END");
}

#undef WRITE

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
