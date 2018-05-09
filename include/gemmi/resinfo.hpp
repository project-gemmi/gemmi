// Copyright 2018 Global Phasing Ltd.
//
// List of common residues with basic data.

#ifndef GEMMI_RESINFO_HPP_
#define GEMMI_RESINFO_HPP_

#include <cstdint>  // for uint8_t
#include <string>

namespace gemmi {

struct ResidueInfo {
  // Simple approximate classification. ELS is something else (ligand).
  enum Kind : unsigned char {
    UNKNOWN=0, AA=1, NA=2, RNA=(2|4), DNA=(2|8), SOL=16, HOH=(16|32),
    PYR=64, ELS=128
  };
  Kind kind;
  // one-letter code or space (uppercase iff it is a standard residues)
  char one_letter_code;
  // rough count of hydrogens used to estimate mass with implicit hydrogens
  std::uint8_t hydrogen_count;

  bool found() const { return kind != UNKNOWN; }
  bool is_water() const { return kind == HOH; }
  bool is_nucleic() const { return kind & NA; }
  bool is_amino() const { return kind == AA; }
  // PDB format has non-standard residues (modified AA) marked as HETATM.
  bool is_standard() const { return (one_letter_code & 0x20) == 0; }
};

// hyderogen_count needs to be verified
inline ResidueInfo find_tabulated_residue(const std::string& name) {
  if (name.size() == 3) {
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(name.c_str())) {
      case ID("ALA"): return { ResidueInfo::AA,  'A',   7 };
      case ID("ARG"): return { ResidueInfo::AA,  'R',  15 };
      case ID("ASN"): return { ResidueInfo::AA,  'N',   8 };
      case ID("ABA"): return { ResidueInfo::AA,  'a',   9 };
      case ID("ASP"): return { ResidueInfo::AA,  'D',   7 };
      case ID("ASX"): return { ResidueInfo::AA,  'B',   6 };
      case ID("CYS"): return { ResidueInfo::AA,  'C',   7 };
      case ID("CSH"): return { ResidueInfo::AA,  's',  17 };
      case ID("GLN"): return { ResidueInfo::AA,  'Q',  10 };
      case ID("GLU"): return { ResidueInfo::AA,  'E',   9 };
      case ID("GLX"): return { ResidueInfo::AA,  'Z',   8 };
      case ID("GLY"): return { ResidueInfo::AA,  'G',   5 };
      case ID("HIS"): return { ResidueInfo::AA,  'H',  10 };
      case ID("ILE"): return { ResidueInfo::AA,  'I',  13 };
      case ID("LEU"): return { ResidueInfo::AA,  'L',  13 };
      case ID("LYS"): return { ResidueInfo::AA,  'K',  15 };
      case ID("MET"): return { ResidueInfo::AA,  'M',  11 };
      case ID("MSE"): return { ResidueInfo::AA,  'm',  11 };
      case ID("ORN"): return { ResidueInfo::AA,  'a',  12 };
      case ID("PHE"): return { ResidueInfo::AA,  'F',  11 };
      case ID("PRO"): return { ResidueInfo::AA,  'P',   9 };
      case ID("SER"): return { ResidueInfo::AA,  'S',   7 };
      case ID("THR"): return { ResidueInfo::AA,  'T',   9 };
      case ID("TRY"): // fall-through - synonym for TRP
      case ID("TRP"): return { ResidueInfo::AA,  'W',  12 };
      case ID("TYR"): return { ResidueInfo::AA,  'Y',  11 };
      case ID("UNK"): return { ResidueInfo::AA,  'X',   9 };
      case ID("VAL"): return { ResidueInfo::AA,  'V',  11 };
      case ID("SEC"): return { ResidueInfo::AA,  'U',   7 };
      case ID("PYL"): return { ResidueInfo::AA,  'O',  21 };

      case ID("SEP"): return { ResidueInfo::AA,  's',   8 };
      case ID("TPO"): return { ResidueInfo::AA,  't',  10 };
      case ID("PCA"): return { ResidueInfo::AA,  'e',   7 };
      case ID("CSO"): return { ResidueInfo::AA,  'c',   7 };
      case ID("PTR"): return { ResidueInfo::AA,  'y',  12 };
      case ID("KCX"): return { ResidueInfo::AA,  'k',  14 };
      case ID("CSD"): return { ResidueInfo::AA,  'c',   7 };
      case ID("LLP"): return { ResidueInfo::AA,  'k',  22 };
      case ID("CME"): return { ResidueInfo::AA,  'c',  11 };
      case ID("MLY"): return { ResidueInfo::AA,  'k',  18 };
      case ID("DAL"): return { ResidueInfo::AA,  'a',   7 };
      case ID("TYS"): return { ResidueInfo::AA,  'y',  11 };
      case ID("OCS"): return { ResidueInfo::AA,  'c',   7 };
      case ID("M3L"): return { ResidueInfo::AA,  'k',  21 };
      case ID("FME"): return { ResidueInfo::AA,  'm',  11 };
      case ID("ALY"): return { ResidueInfo::AA,  'k',  16 };
      case ID("HYP"): return { ResidueInfo::AA,  'p',   9 };
      case ID("CAS"): return { ResidueInfo::AA,  'c',  12 };
      case ID("CRO"): return { ResidueInfo::AA,  't',  17 };
      case ID("CSX"): return { ResidueInfo::AA,  'c',   7 };
      case ID("DPR"): return { ResidueInfo::AA,  'p',   9 };
      case ID("DGL"): return { ResidueInfo::AA,  'e',   9 };
      case ID("DVA"): return { ResidueInfo::AA,  'v',  11 };
      case ID("CSS"): return { ResidueInfo::AA,  'c',   7 };
      case ID("DPN"): return { ResidueInfo::AA,  'f',  11 };
      case ID("DSN"): return { ResidueInfo::AA,  's',   7 };
      case ID("DLE"): return { ResidueInfo::AA,  'l',  13 };
      case ID("HIC"): return { ResidueInfo::AA,  'h',  11 };
      case ID("NLE"): return { ResidueInfo::AA,  'l',  13 };
      case ID("MVA"): return { ResidueInfo::AA,  'v',  13 };
      case ID("MLZ"): return { ResidueInfo::AA,  'k',  16 };
      case ID("CR2"): return { ResidueInfo::AA,  'g',  13 };
      case ID("SAR"): return { ResidueInfo::AA,  'g',   7 };
      case ID("DAR"): return { ResidueInfo::AA,  'r',  15 };
      case ID("DLY"): return { ResidueInfo::AA,  'k',  14 };
      case ID("YCM"): return { ResidueInfo::AA,  'c',  10 };
      case ID("NRQ"): return { ResidueInfo::AA,  'm',  17 };
      case ID("CGU"): return { ResidueInfo::AA,  'e',   9 };
      case ID("0TD"): return { ResidueInfo::AA,  'd',   9 };
      case ID("MLE"): return { ResidueInfo::AA,  'l',  15 };
      case ID("DAS"): return { ResidueInfo::AA,  'd',   7 };
      case ID("DTR"): return { ResidueInfo::AA,  'w',  12 };
      case ID("CXM"): return { ResidueInfo::AA,  'm',  11 };
      case ID("TPQ"): return { ResidueInfo::AA,  'y',   9 };
      case ID("DCY"): return { ResidueInfo::AA,  'c',   7 };
      case ID("DSG"): return { ResidueInfo::AA,  'n',   8 };
      case ID("DTY"): return { ResidueInfo::AA,  'y',  11 };
      case ID("DHI"): return { ResidueInfo::AA,  'h',  10 };
      case ID("MEN"): return { ResidueInfo::AA,  'n',  10 };
      case ID("DTH"): return { ResidueInfo::AA,  't',   9 };
      case ID("SAC"): return { ResidueInfo::AA,  's',   9 };
      case ID("DGN"): return { ResidueInfo::AA,  'q',  10 };
      case ID("AIB"): return { ResidueInfo::AA,  'a',   9 };
      case ID("SMC"): return { ResidueInfo::AA,  'c',   9 };
      case ID("IAS"): return { ResidueInfo::AA,  'd',   7 };
      case ID("CIR"): return { ResidueInfo::AA,  'r',  13 };
      case ID("BMT"): return { ResidueInfo::AA,  't',  19 };
      case ID("DIL"): return { ResidueInfo::AA,  'i',  13 };
      case ID("FGA"): return { ResidueInfo::AA,  'e',   9 };
      case ID("PHI"): return { ResidueInfo::AA,  'f',  10 };
      case ID("CRQ"): return { ResidueInfo::AA,  'q',  16 };
      case ID("SME"): return { ResidueInfo::AA,  'm',  11 };
      case ID("GHP"): return { ResidueInfo::AA,  'g',   9 };
      case ID("MHO"): return { ResidueInfo::AA,  'm',  11 };
      case ID("NEP"): return { ResidueInfo::AA,  'h',  10 };
      case ID("TRQ"): return { ResidueInfo::AA,  'w',  10 };
      case ID("TOX"): return { ResidueInfo::AA,  'w',  12 };
      case ID("ALC"): return { ResidueInfo::AA,  'a',  17 };
      case ID("3FG"): return { ResidueInfo::AA,  ' ',   9 };
      case ID("SCH"): return { ResidueInfo::AA,  'c',   9 };
      case ID("MDO"): return { ResidueInfo::AA,  'a',  11 };
      case ID("MAA"): return { ResidueInfo::AA,  'a',   9 };
      case ID("GYS"): return { ResidueInfo::AA,  's',  15 };
      case ID("MK8"): return { ResidueInfo::AA,  'l',  15 };
      case ID("CR8"): return { ResidueInfo::AA,  'h',  16 };
      case ID("KPI"): return { ResidueInfo::AA,  'k',  16 };
      case ID("SCY"): return { ResidueInfo::AA,  'c',   9 };
      case ID("DHA"): return { ResidueInfo::AA,  's',   5 };
      case ID("OMY"): return { ResidueInfo::AA,  'y',  10 };
      case ID("CAF"): return { ResidueInfo::AA,  'c',  12 };
      case ID("0AF"): return { ResidueInfo::AA,  'w',  12 };
      case ID("SNN"): return { ResidueInfo::AA,  'n',   6 };
      case ID("MHS"): return { ResidueInfo::AA,  'h',  11 };
      case ID("MLU"): return { ResidueInfo::AA,  ' ',  15 };
      case ID("SNC"): return { ResidueInfo::AA,  'c',   6 };
      case ID("PHD"): return { ResidueInfo::AA,  'd',   8 };
      case ID("B3E"): return { ResidueInfo::AA,  'e',  11 };
      case ID("MEA"): return { ResidueInfo::AA,  'f',  13 };
      case ID("MED"): return { ResidueInfo::AA,  'm',  11 };
      case ID("OAS"): return { ResidueInfo::AA,  's',   9 };
      case ID("GL3"): return { ResidueInfo::AA,  'g',   5 };
      case ID("FVA"): return { ResidueInfo::AA,  'v',  11 };
      case ID("PHL"): return { ResidueInfo::AA,  'f',  13 };
      case ID("CRF"): return { ResidueInfo::AA,  't',  18 };
      case ID("OMZ"): return { ResidueInfo::AA,  ' ',  10 };
      case ID("BFD"): return { ResidueInfo::AA,  'd',   6 };
      case ID("MEQ"): return { ResidueInfo::AA,  'q',  12 };
      case ID("DAB"): return { ResidueInfo::AA,  'a',  10 };
      case ID("AGM"): return { ResidueInfo::AA,  'r',  17 };

      case ID("PSU"): return { ResidueInfo::RNA, 'u',  13 };
      case ID("5MU"): return { ResidueInfo::RNA, 'u',  15 };
      case ID("7MG"): return { ResidueInfo::RNA, 'g',  18 };
      case ID("OMG"): return { ResidueInfo::RNA, 'g',  16 };
      case ID("UR3"): return { ResidueInfo::RNA, 'u',  15 };
      case ID("OMC"): return { ResidueInfo::RNA, 'c',  16 };
      case ID("2MG"): return { ResidueInfo::RNA, 'g',  16 };
      case ID("H2U"): return { ResidueInfo::RNA, 'u',  15 };
      case ID("4SU"): return { ResidueInfo::RNA, 'u',  13 };
      case ID("OMU"): return { ResidueInfo::RNA, 'u',  15 };
      case ID("4OC"): return { ResidueInfo::RNA, 'c',  18 };
      case ID("MA6"): return { ResidueInfo::RNA, 'a',  18 };
      case ID("M2G"): return { ResidueInfo::RNA, 'g',  18 };
      case ID("1MA"): return { ResidueInfo::RNA, 'a',  16 };
      case ID("6MZ"): return { ResidueInfo::RNA, 'a',  16 };
      case ID("CCC"): return { ResidueInfo::RNA, 'c',  13 };
      case ID("2MA"): return { ResidueInfo::RNA, 'a',  16 };
      case ID("1MG"): return { ResidueInfo::RNA, 'g',  16 };
      case ID("5BU"): return { ResidueInfo::RNA, 'u',  12 };
      case ID("MIA"): return { ResidueInfo::RNA, 'a',  26 };
      case ID("DOC"): return { ResidueInfo::DNA, 'c',  14 };
      case ID("8OG"): return { ResidueInfo::DNA, 'g',  14 };
      case ID("5CM"): return { ResidueInfo::DNA, 'c',  16 };
      case ID("3DR"): return { ResidueInfo::DNA, ' ',  11 };
      case ID("BRU"): return { ResidueInfo::DNA, 'u',  12 };
      case ID("CBR"): return { ResidueInfo::DNA, 'c',  13 };

      case ID("H2O"): // fall-through - synonym ' ',for HOH
      case ID("WAT"): // fall-through - synonym for HOH
      case ID("HOH"): return { ResidueInfo::HOH, ' ',   2 };
      case ID("DOD"): return { ResidueInfo::HOH, ' ',   2 };
      case ID("HEM"): return { ResidueInfo::ELS, ' ',  32 };
      case ID("SUL"): // fall-through - synonym for SO4
      case ID("SO4"): return { ResidueInfo::ELS, ' ',   0 };
      case ID("GOL"): return { ResidueInfo::ELS, ' ',   8 };
      case ID("EDO"): return { ResidueInfo::ELS, ' ',   6 };
      case ID("NAG"): return { ResidueInfo::PYR, ' ',  15 };
      case ID("PO4"): return { ResidueInfo::ELS, ' ',   0 };
      case ID("ACT"): return { ResidueInfo::ELS, ' ',   3 };
      case ID("PEG"): return { ResidueInfo::ELS, ' ',  10 };
      case ID("MAN"): return { ResidueInfo::PYR, ' ',  12 };
      case ID("FAD"): return { ResidueInfo::ELS, ' ',  33 };
      case ID("BMA"): return { ResidueInfo::PYR, ' ',  12 };
      case ID("ADP"): return { ResidueInfo::ELS, ' ',  15 };
      case ID("DMS"): return { ResidueInfo::ELS, ' ',   6 };
      // ACE is a non-polymer that occurs primarily in polymers
      case ID("ACE"): return { ResidueInfo::ELS, ' ',   4 };
      case ID("MPD"): return { ResidueInfo::ELS, ' ',  14 };
      case ID("MES"): return { ResidueInfo::ELS, ' ',  13 };
      case ID("NAD"): return { ResidueInfo::ELS, ' ',  27 };
      case ID("NAP"): return { ResidueInfo::ELS, ' ',  28 };
      case ID("TRS"): return { ResidueInfo::ELS, ' ',  12 };
      case ID("ATP"): return { ResidueInfo::ELS, ' ',  16 };
      case ID("PG4"): return { ResidueInfo::ELS, ' ',  18 };
      case ID("GDP"): return { ResidueInfo::ELS, 'g',  15 }; // RNA in CCD
      case ID("FUC"): return { ResidueInfo::PYR, ' ',  12 };
      case ID("FMT"): return { ResidueInfo::ELS, ' ',   2 };
      case ID("NH2"): return { ResidueInfo::ELS, ' ',   2 }; // ?
      case ID("GAL"): return { ResidueInfo::PYR, ' ',  12 };
      case ID("PGE"): return { ResidueInfo::ELS, ' ',  14 };
      case ID("FMN"): return { ResidueInfo::ELS, ' ',  21 };
      case ID("PLP"): return { ResidueInfo::ELS, ' ',  10 };
      case ID("EPE"): return { ResidueInfo::ELS, ' ',  18 };
      case ID("SF4"): return { ResidueInfo::ELS, ' ',   0 };
      case ID("BME"): return { ResidueInfo::ELS, ' ',   6 };
      case ID("CIT"): return { ResidueInfo::ELS, ' ',   8 };
    }
#undef ID
  } else if (name.size() == 1) {
    switch (name[0]) {
      case 'A': return { ResidueInfo::RNA, 'A',  14 };
      case 'C': return { ResidueInfo::RNA, 'C',  14 };
      case 'G': return { ResidueInfo::RNA, 'G',  14 };
      case 'I': return { ResidueInfo::RNA, 'I',  13 };
      case 'U': return { ResidueInfo::RNA, 'U',  13 };

      case 'K': return { ResidueInfo::ELS, ' ',   0 };
    }
  } else if (name.size() == 2) {
    if (name[0] == 'D' || name[0] == '+')
      switch (name[1]) {
        case 'A': return { ResidueInfo::DNA, 'A',  14 };
        case 'C': return { ResidueInfo::DNA, 'C',  14 };
        case 'G': return { ResidueInfo::DNA, 'G',  14 };
        case 'I': return { ResidueInfo::DNA, 'I',  13 };
        case 'T': return { ResidueInfo::DNA, 'T',  15 };
        case 'U': return { ResidueInfo::DNA, 'U',  13 };
      }
    else
#define ID(s) (s[0] << 8 | s[1])
      switch (ID(name.c_str())) {
        case ID("ZN"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("MG"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("CL"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("CA"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("NA"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("MN"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("FE"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("NI"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("CU"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("CD"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("CO"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("HG"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("BR"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("BA"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("SR"): return { ResidueInfo::ELS, ' ', 0 };
        case ID("NO"): return { ResidueInfo::ELS, ' ', 0 };
      }
#undef ID
  }
  return { ResidueInfo::UNKNOWN, ' ', 0 };
}

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
