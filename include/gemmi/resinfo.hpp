// Copyright 2018 Global Phasing Ltd.
//
// List of common residues with basic data.

#ifndef GEMMI_RESINFO_HPP_
#define GEMMI_RESINFO_HPP_

#include <cstdint>  // for uint8_t
#include <string>  // for uint8_t

namespace gemmi {

struct ResidueInfo {
  // simple classification, ELS is something else
  enum Kind : char {
    UNKNOWN=0, AA=1, NA=2, RNA=(2|4), DNA=(2|8), SOL=16, HOH=(16|32), ELS=64
  };
  Kind kind;
  // PDB format has non-standard residues (modified AA) marked as HETATM.
  bool pdb_standard;
  // rough count of hydrogens used to estimate mass with implicit hydrogens
  std::uint8_t hydrogen_count;

  bool is_water() const { return kind == HOH; }
  bool is_nucleic() const { return kind & NA; }
  bool is_amino() const { return kind == AA; }
};

// hyderogen_count needs to be verified
inline const ResidueInfo find_tabulated_residue(const std::string& name) {
  if (name.size() == 3) {
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(name.c_str())) {
      case ID("ALA"): return { ResidueInfo::AA,  true,  5 };
      case ID("ARG"): return { ResidueInfo::AA,  true, 13 };
      case ID("ASN"): return { ResidueInfo::AA,  true,  6 };
      case ID("ABA"): return { ResidueInfo::AA,  false, 7 };
      case ID("ASP"): return { ResidueInfo::AA,  true,  4 };
      case ID("ASX"): return { ResidueInfo::AA,  true,  4 };
      case ID("CYS"): return { ResidueInfo::AA,  true,  4 };
      case ID("CSH"): return { ResidueInfo::AA,  false, 5 };
      case ID("GLN"): return { ResidueInfo::AA,  true,  8 };
      case ID("GLU"): return { ResidueInfo::AA,  true,  6 };
      case ID("GLX"): return { ResidueInfo::AA,  true,  8 };
      case ID("GLY"): return { ResidueInfo::AA,  true,  3 };
      case ID("HIS"): return { ResidueInfo::AA,  true,  8 };
      case ID("ILE"): return { ResidueInfo::AA,  true, 11 };
      case ID("LEU"): return { ResidueInfo::AA,  true, 11 };
      case ID("LYS"): return { ResidueInfo::AA,  true, 13 };
      case ID("MET"): return { ResidueInfo::AA,  true,  9 };
      case ID("MSE"): return { ResidueInfo::AA,  false, 9 };
      case ID("ORN"): return { ResidueInfo::AA,  false,10 };
      case ID("PHE"): return { ResidueInfo::AA,  true,  9 };
      case ID("PRO"): return { ResidueInfo::AA,  true,  7 };
      case ID("SER"): return { ResidueInfo::AA,  true,  5 };
      case ID("THR"): return { ResidueInfo::AA,  true,  7 };
      case ID("TRP"): return { ResidueInfo::AA,  true, 10 };
      case ID("TYR"): return { ResidueInfo::AA,  true,  9 };
      case ID("UNK"): return { ResidueInfo::AA,  true,  0 };
      case ID("VAL"): return { ResidueInfo::AA,  true,  9 };
      case ID("SEC"): return { ResidueInfo::AA,  true,  6 };
      case ID("PYL"): return { ResidueInfo::AA,  true, 19 };
      case ID("HOH"): return { ResidueInfo::HOH, false, 2 };
      case ID("WAT"): return { ResidueInfo::HOH, false, 2 };
      case ID("H20"): return { ResidueInfo::HOH, false, 2 };
      case ID("DOD"): return { ResidueInfo::HOH, false, 2 };
      case ID("HEM"): return { ResidueInfo::ELS, false,30 };
      case ID("SO4"): return { ResidueInfo::ELS, false, 0 };
      case ID("SUL"): return { ResidueInfo::ELS, false, 0 };
    }
#undef ID
  } else if (name.size() == 1) {
    switch (name[0]) {
      case 'A': return { ResidueInfo::RNA, true, 13 };
      case 'C': return { ResidueInfo::RNA, true, 13 };
      case 'G': return { ResidueInfo::RNA, true, 13 };
      case 'I': return { ResidueInfo::RNA, true, 12 };
      case 'T': return { ResidueInfo::RNA, true, 14 };
      case 'U': return { ResidueInfo::RNA, true, 12 };
    }
  } else if (name.size() == 2) {
    if (name[0] == 'D' || name[0] == '+')
      switch (name[1]) {
        case 'A': return { ResidueInfo::DNA, true, 13 };
        case 'C': return { ResidueInfo::DNA, true, 13 };
        case 'G': return { ResidueInfo::DNA, true, 13 };
        case 'I': return { ResidueInfo::DNA, true, 12 };
        case 'T': return { ResidueInfo::DNA, true, 14 };
        case 'U': return { ResidueInfo::DNA, true, 12 };
      }
  }
  return { ResidueInfo::UNKNOWN, false, 0 };
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
