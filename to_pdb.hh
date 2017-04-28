// Copyright 2017 Global Phasing Ltd.
//
// Writing PDB file format.

#ifndef GEMMI_TO_PDB_HH_
#define GEMMI_TO_PDB_HH_

#include <ostream>
#include <stb_sprintf.h>
#include "model.hh"

namespace gemmi {
namespace mol {

void write_pdb(const Structure& st, std::ostream& os) {
  char buf[83] = {0};

#define WRITE(...) do { \
    stbsp_snprintf(buf, 82, __VA_ARGS__); \
    os.write(buf, 81); \
    } while(0)

  int serial = 0;
  for (const mol::Model& model : st.models) {
    if (st.models.size() > 1)
      WRITE("%-10s%-70s\n", "MODEL", model.name.c_str());
    for (const mol::Chain& chain : model.chains) {
      for (const mol::Residue& res : chain.residues) {
        for (const mol::Atom& a : res.atoms) {
          //  1- 6  6s  record name
          //  7-11  5d  integer serial
          // 12     1   -
          // 13-16  4s  atom name
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
          WRITE("%-6s%5d %-4s%c%-3s"
                " %1s%4d%c"
                "   %8.3f%8.3f%8.3f"
                "%6.2f%6.2f          %2s%c%c\n",
                "ATOM", ++serial, a.name.c_str(),
                a.altloc ? a.altloc : ' ',
                res.name.c_str(),
                chain.name.c_str(), res.seq_id,
                res.ins_code ? res.ins_code : ' ',
                a.x, a.y, a.z,
                a.occ, a.b_iso, a.element.uname(),
                a.charge >= 0 ? ' ' : '-',
                a.charge != 0 ? a.charge + '0' : ' ');
        }
      }
      WRITE("%-80s\n", "TER");
    }
    if (st.models.size() > 1)
      WRITE("%-80s\n", "ENDMDL");
  }
}

#undef WRITE

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
