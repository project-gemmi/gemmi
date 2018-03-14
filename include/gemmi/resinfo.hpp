// Copyright 2018 Global Phasing Ltd.
//
// List of common residues with basic data.

#ifndef GEMMI_RESINFO_HPP_
#define GEMMI_RESINFO_HPP_

#include <cstdint>  // for uint8_t
#include <string>  // for uint8_t

namespace gemmi {

struct ResidueInfo {
  // Simple approximate classification. ELS is something else (ligand).
  enum Kind : unsigned char {
    UNKNOWN=0, AA=1, NA=2, RNA=(2|4), DNA=(2|8), SOL=16, HOH=(16|32),
    PYR=64, ELS=128
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
inline ResidueInfo find_tabulated_residue(const std::string& name) {
  if (name.size() == 3) {
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(name.c_str())) {
      case ID("ALA"): return { ResidueInfo::AA,  true,   5 };
      case ID("ARG"): return { ResidueInfo::AA,  true,  13 };
      case ID("ASN"): return { ResidueInfo::AA,  true,   6 };
      case ID("ABA"): return { ResidueInfo::AA,  false,  7 };
      case ID("ASP"): return { ResidueInfo::AA,  true,   5 };
      case ID("ASX"): return { ResidueInfo::AA,  true,   4 };
      case ID("CYS"): return { ResidueInfo::AA,  true,   5 };
      case ID("CSH"): return { ResidueInfo::AA,  false,  15 };
      case ID("GLN"): return { ResidueInfo::AA,  true,   8 };
      case ID("GLU"): return { ResidueInfo::AA,  true,   7 };
      case ID("GLX"): return { ResidueInfo::AA,  true,   6 };
      case ID("GLY"): return { ResidueInfo::AA,  true,   3 };
      case ID("HIS"): return { ResidueInfo::AA,  true,   8 };
      case ID("ILE"): return { ResidueInfo::AA,  true,  11 };
      case ID("LEU"): return { ResidueInfo::AA,  true,  11 };
      case ID("LYS"): return { ResidueInfo::AA,  true,  13 };
      case ID("MET"): return { ResidueInfo::AA,  true,   9 };
      case ID("MSE"): return { ResidueInfo::AA,  false,  9 };
      case ID("ORN"): return { ResidueInfo::AA,  false, 10 };
      case ID("PHE"): return { ResidueInfo::AA,  true,   9 };
      case ID("PRO"): return { ResidueInfo::AA,  true,   7 };
      case ID("SER"): return { ResidueInfo::AA,  true,   5 };
      case ID("THR"): return { ResidueInfo::AA,  true,   7 };
      case ID("TRY"): // fall-through - synonym for TRP
      case ID("TRP"): return { ResidueInfo::AA,  true,  10 };
      case ID("TYR"): return { ResidueInfo::AA,  true,   9 };
      case ID("UNK"): return { ResidueInfo::AA,  true,   7 };
      case ID("VAL"): return { ResidueInfo::AA,  true,   9 };
      case ID("SEC"): return { ResidueInfo::AA,  true,   5 };
      case ID("PYL"): return { ResidueInfo::AA,  true,  19 };

      case ID("SEP"): return { ResidueInfo::AA,  false,  6 };
      case ID("TPO"): return { ResidueInfo::AA,  false,  8 };
      case ID("PCA"): return { ResidueInfo::AA,  false,  5 };
      case ID("CSO"): return { ResidueInfo::AA,  false,  5 };
      case ID("PTR"): return { ResidueInfo::AA,  false, 10 };
      case ID("KCX"): return { ResidueInfo::AA,  false, 12 };
      case ID("CSD"): return { ResidueInfo::AA,  false,  5 };
      case ID("LLP"): return { ResidueInfo::AA,  false, 20 };
      case ID("CME"): return { ResidueInfo::AA,  false,  9 };
      case ID("MLY"): return { ResidueInfo::AA,  false, 16 };
      case ID("DAL"): return { ResidueInfo::AA,  false,  5 };
      case ID("TYS"): return { ResidueInfo::AA,  false,  9 };
      case ID("OCS"): return { ResidueInfo::AA,  false,  5 };
      case ID("M3L"): return { ResidueInfo::AA,  false, 19 };
      case ID("FME"): return { ResidueInfo::AA,  false,  9 };
      case ID("ALY"): return { ResidueInfo::AA,  false, 14 };
      case ID("HYP"): return { ResidueInfo::AA,  false,  7 };
      case ID("CAS"): return { ResidueInfo::AA,  false, 10 };
      case ID("CRO"): return { ResidueInfo::AA,  false, 15 };
      case ID("CSX"): return { ResidueInfo::AA,  false,  5 };
      case ID("DPR"): return { ResidueInfo::AA,  false,  7 };
      case ID("DGL"): return { ResidueInfo::AA,  false,  7 };
      case ID("DVA"): return { ResidueInfo::AA,  false,  9 };
      case ID("CSS"): return { ResidueInfo::AA,  false,  5 };
      case ID("DPN"): return { ResidueInfo::AA,  false,  9 };
      case ID("DSN"): return { ResidueInfo::AA,  false,  5 };
      case ID("DLE"): return { ResidueInfo::AA,  false, 11 };
      case ID("HIC"): return { ResidueInfo::AA,  false,  9 };
      case ID("NLE"): return { ResidueInfo::AA,  false, 11 };
      case ID("MVA"): return { ResidueInfo::AA,  false, 11 };
      case ID("MLZ"): return { ResidueInfo::AA,  false, 14 };
      case ID("CR2"): return { ResidueInfo::AA,  false, 11 };
      case ID("SAR"): return { ResidueInfo::AA,  false,  5 };
      case ID("DAR"): return { ResidueInfo::AA,  false, 13 };
      case ID("DLY"): return { ResidueInfo::AA,  false, 12 };
      case ID("YCM"): return { ResidueInfo::AA,  false,  8 };
      case ID("NRQ"): return { ResidueInfo::AA,  false, 15 };
      case ID("CGU"): return { ResidueInfo::AA,  false,  7 };
      case ID("0TD"): return { ResidueInfo::AA,  false,  7 };
      case ID("MLE"): return { ResidueInfo::AA,  false, 13 };
      case ID("DAS"): return { ResidueInfo::AA,  false,  5 };
      case ID("DTR"): return { ResidueInfo::AA,  false, 10 };
      case ID("CXM"): return { ResidueInfo::AA,  false,  9 };
      case ID("TPQ"): return { ResidueInfo::AA,  false,  7 };
      case ID("DCY"): return { ResidueInfo::AA,  false,  5 };
      case ID("DSG"): return { ResidueInfo::AA,  false,  6 };
      case ID("DTY"): return { ResidueInfo::AA,  false,  9 };
      case ID("DHI"): return { ResidueInfo::AA,  false,  8 };
      case ID("MEN"): return { ResidueInfo::AA,  false,  8 };
      case ID("DTH"): return { ResidueInfo::AA,  false,  7 };
      case ID("SAC"): return { ResidueInfo::AA,  false,  7 };
      case ID("DGN"): return { ResidueInfo::AA,  false,  8 };
      case ID("AIB"): return { ResidueInfo::AA,  false,  7 };
      case ID("SMC"): return { ResidueInfo::AA,  false,  7 };
      case ID("BMT"): return { ResidueInfo::AA,  false, 17 };
      case ID("DIL"): return { ResidueInfo::AA,  false, 11 };
      case ID("PSU"): return { ResidueInfo::RNA, false, 13 };
      case ID("5MU"): return { ResidueInfo::RNA, false, 15 };
      case ID("7MG"): return { ResidueInfo::RNA, false, 18 };
      case ID("OMG"): return { ResidueInfo::RNA, false, 16 };
      case ID("UR3"): return { ResidueInfo::RNA, false, 15 };
      case ID("OMC"): return { ResidueInfo::RNA, false, 16 };
      case ID("2MG"): return { ResidueInfo::RNA, false, 16 };
      case ID("H2U"): return { ResidueInfo::RNA, false, 15 };
      case ID("4SU"): return { ResidueInfo::RNA, false, 13 };
      case ID("OMU"): return { ResidueInfo::RNA, false, 15 };
      case ID("4OC"): return { ResidueInfo::RNA, false, 18 };
      case ID("MA6"): return { ResidueInfo::RNA, false, 18 };
      case ID("M2G"): return { ResidueInfo::RNA, false, 18 };
      case ID("1MA"): return { ResidueInfo::RNA, false, 16 };
      case ID("6MZ"): return { ResidueInfo::RNA, false, 16 };
      case ID("CCC"): return { ResidueInfo::RNA, false, 13 };
      case ID("2MA"): return { ResidueInfo::RNA, false, 16 };
      case ID("1MG"): return { ResidueInfo::RNA, false, 16 };
      case ID("5BU"): return { ResidueInfo::RNA, false, 12 };
      case ID("MIA"): return { ResidueInfo::RNA, false, 26 };
      case ID("DOC"): return { ResidueInfo::DNA, false, 14 };
      case ID("8OG"): return { ResidueInfo::DNA, false, 14 };
      case ID("5CM"): return { ResidueInfo::DNA, false, 16 };
      case ID("3DR"): return { ResidueInfo::DNA, false, 11 };
      case ID("BRU"): return { ResidueInfo::DNA, false, 12 };
      case ID("CBR"): return { ResidueInfo::DNA, false, 13 };

      case ID("H2O"): // fall-through - synonym for HOH
      case ID("WAT"): // fall-through - synonym for HOH
      case ID("HOH"): return { ResidueInfo::HOH, false,  2 };
      case ID("DOD"): return { ResidueInfo::HOH, false,  2 };
      case ID("HEM"): return { ResidueInfo::ELS, false, 30 };
      case ID("SUL"): // fall-through - synonym for SO4
      case ID("SO4"): return { ResidueInfo::ELS, false,  0 };
      case ID("GOL"): return { ResidueInfo::ELS, false,  8 };
      case ID("EDO"): return { ResidueInfo::ELS, false,  6 };
      case ID("NAG"): return { ResidueInfo::PYR, false, 15 };
      case ID("PO4"): return { ResidueInfo::ELS, false,  0 };
      case ID("ACT"): return { ResidueInfo::ELS, false,  3 };
      case ID("PEG"): return { ResidueInfo::ELS, false, 10 };
      case ID("MAN"): return { ResidueInfo::PYR, false, 12 };
      case ID("FAD"): return { ResidueInfo::ELS, false, 33 };
      case ID("BMA"): return { ResidueInfo::PYR, false, 12 };
      case ID("ADP"): return { ResidueInfo::ELS, false, 15 };
      case ID("DMS"): return { ResidueInfo::ELS, false,  6 };
      // ACE is a non-polymer that occurs primarily in polymers
      case ID("ACE"): return { ResidueInfo::ELS, false,  4 };
      case ID("MPD"): return { ResidueInfo::ELS, false, 14 };
      case ID("MES"): return { ResidueInfo::ELS, false, 13 };
      case ID("NAD"): return { ResidueInfo::ELS, false, 27 };
      case ID("NAP"): return { ResidueInfo::ELS, false, 28 };
      case ID("TRS"): return { ResidueInfo::ELS, false, 12 };
      case ID("ATP"): return { ResidueInfo::ELS, false, 16 };
      case ID("PG4"): return { ResidueInfo::ELS, false, 18 };
      case ID("GDP"): return { ResidueInfo::ELS, false, 15 }; // RNA in CCD
      case ID("FUC"): return { ResidueInfo::PYR, false, 12 };
      case ID("FMT"): return { ResidueInfo::ELS, false,  2 };
      case ID("NH2"): return { ResidueInfo::ELS, false,  2 }; // ?
      case ID("GAL"): return { ResidueInfo::PYR, false, 12 };
      case ID("PGE"): return { ResidueInfo::ELS, false, 14 };
      case ID("FMN"): return { ResidueInfo::ELS, false, 21 };
      case ID("PLP"): return { ResidueInfo::ELS, false, 10 };
      case ID("EPE"): return { ResidueInfo::ELS, false, 18 };
      case ID("SF4"): return { ResidueInfo::ELS, false,  0 };
      case ID("BME"): return { ResidueInfo::ELS, false,  6 };
      case ID("CIT"): return { ResidueInfo::ELS, false,  8 };
    }
#undef ID
  } else if (name.size() == 1) {
    switch (name[0]) {
      case 'A': return { ResidueInfo::RNA, true, 14 };
      case 'C': return { ResidueInfo::RNA, true, 14 };
      case 'G': return { ResidueInfo::RNA, true, 14 };
      case 'I': return { ResidueInfo::RNA, true, 13 };
      case 'U': return { ResidueInfo::RNA, true, 13 };

      case 'K': return { ResidueInfo::ELS, false, 0 };
    }
  } else if (name.size() == 2) {
    if (name[0] == 'D' || name[0] == '+')
      switch (name[1]) {
        case 'A': return { ResidueInfo::DNA, true, 14 };
        case 'C': return { ResidueInfo::DNA, true, 14 };
        case 'G': return { ResidueInfo::DNA, true, 14 };
        case 'I': return { ResidueInfo::DNA, true, 13 };
        case 'T': return { ResidueInfo::DNA, true, 15 };
        case 'U': return { ResidueInfo::DNA, true, 13 };
      }
    else
#define ID(s) (s[0] << 8 | s[1])
      switch (ID(name.c_str())) {
        case ID("ZN"): return { ResidueInfo::ELS, false, 0 };
        case ID("MG"): return { ResidueInfo::ELS, false, 0 };
        case ID("CL"): return { ResidueInfo::ELS, false, 0 };
        case ID("CA"): return { ResidueInfo::ELS, false, 0 };
        case ID("NA"): return { ResidueInfo::ELS, false, 0 };
        case ID("MN"): return { ResidueInfo::ELS, false, 0 };
        case ID("FE"): return { ResidueInfo::ELS, false, 0 };
        case ID("NI"): return { ResidueInfo::ELS, false, 0 };
        case ID("CU"): return { ResidueInfo::ELS, false, 0 };
        case ID("CD"): return { ResidueInfo::ELS, false, 0 };
        case ID("CO"): return { ResidueInfo::ELS, false, 0 };
        case ID("HG"): return { ResidueInfo::ELS, false, 0 };
        case ID("BR"): return { ResidueInfo::ELS, false, 0 };
        case ID("BA"): return { ResidueInfo::ELS, false, 0 };
        case ID("SR"): return { ResidueInfo::ELS, false, 0 };
        case ID("NO"): return { ResidueInfo::ELS, false, 0 };
      }
#undef ID
  }
  return { ResidueInfo::UNKNOWN, false, 0 };
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
