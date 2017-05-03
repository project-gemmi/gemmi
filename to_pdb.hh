// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format.

#ifndef GEMMI_TO_PDB_HH_
#define GEMMI_TO_PDB_HH_

#include <ostream>
#include <cstring>
#include <algorithm>
#include <stb_sprintf.h>
#include "model.hh"

namespace gemmi {
namespace mol {

inline void write_pdb(const Structure& st, std::ostream& os) {
  const char* months = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC???";
  char buf[83] = {0};

#define WRITE(...) do { \
    stbsp_snprintf(buf, 82, __VA_ARGS__); \
    os.write(buf, 81); \
    } while(0)

  const char* date = st.get_info("_database_PDB_rev.date_original");
  std::string pdb_date;
  if (date && std::strlen(date) == 10) {
    unsigned month_idx = 10 * (date[5] - '0') + date[6] - '0' - 1;
    pdb_date = std::string(date + 8, 2) + "-" +
               std::string(months + 3 * std::min(month_idx, 13u), 3) +
               "-" + std::string(date + 2, 2);
  }
  WRITE("HEADER    %-40s%-9s   %-18s\n",
        // 11-50 is called classification in the PDB spec, but in mmCIF it is:
        st.get_info("_struct_keywords.pdbx_keywords", ""),
        pdb_date.c_str(), st.get_info("_entry.id"));
  const char* title = st.get_info("_struct.title");
  if (title) // TODO: line breaking for long strings
    WRITE("TITLE     %-70.60s\n", title);
  // TODO: SEQRES
  WRITE("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s          \n",
        st.cell.a, st.cell.b, st.cell.c,
        st.cell.alpha, st.cell.beta, st.cell.gamma,
        st.sg_hm.c_str(), st.get_info("_cell.Z_PDB", ""));
  // TODO: SCALE
  //TODO: special handling of large structures (>62 chains or >=1M atoms)
  for (const mol::Model& model : st.models) {
    int serial = 0;
    if (st.models.size() > 1)
      WRITE("%-10s%-70s\n", "MODEL", model.name.c_str());
    char prev_chain = 0;
    for (const mol::Chain& chain : model.chains) {
      if (chain.auth_name.empty())
        throw std::runtime_error("empty chain name");
      if (chain.auth_name.length() > 1)
        throw std::runtime_error("long chain name: " + chain.auth_name);
      // TODO: proper detection when polymer ends
      if (prev_chain && (prev_chain != chain.auth_name[0] ||
                         chain.residues[0].seq_id == mol::Residue::UnknownId)) {
        // re-using part of the buffer in the middle, e.g.:
        // TER    4153      LYS B 286
        stbsp_sprintf(buf, "TER   %5d", ++serial);
        std::memset(buf+11, ' ', 6);
        std::memset(buf+30, ' ', 50);
        os.write(buf, 81);
        prev_chain = 0;
      }
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
          WRITE("%-6s%5d %c%-3s%c%-3s"
                " %1s%4d%c"
                "   %8.3f%8.3f%8.3f"
                "%6.2f%6.2f          %2s%c%c\n",
                standard ? "ATOM" : "HETATM",
                ++serial,
                empty13 ? ' ' : a.name[0],
                a.name.c_str() + (empty13 || a.name.empty() ? 0 : 1),
                a.altloc ? a.altloc : ' ',
                res.name.c_str(),
                chain.auth_name.c_str(), res.seq_id_for_pdb(),
                res.ins_code ? res.ins_code : ' ',
                a.x, a.y, a.z,
                a.occ, a.b_iso, a.element.uname(),
                a.charge >= 0 ? ' ' : '-',
                a.charge != 0 ? a.charge + '0' : ' ');
        }
        if (res.seq_id != mol::Residue::UnknownId)
          prev_chain = chain.auth_name[0];
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
