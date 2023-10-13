// Copyright 2018-2022 Global Phasing Ltd.

#include <gemmi/resinfo.hpp>
#include <gemmi/atox.hpp>  // for is_space

namespace gemmi {

// hydrogen_count needs to be verified
ResidueInfo find_tabulated_residue(const std::string& name) {
  using RI = ResidueKind;
  if (name.size() == 3) {
#define ID(s) (s[0] << 16 | s[1] << 8 | s[2])
    switch (ID(name.c_str())) {
      case ID("ALA"): return { RI::AA,  1, 'A',   7, 89.0932f };
      case ID("ARG"): return { RI::AA,  1, 'R',  15, 175.209f };
      case ID("ASN"): return { RI::AA,  1, 'N',   8, 132.118f };
      case ID("ABA"): return { RI::AA,  1, 'a',   9, 103.120f };
      case ID("ASP"): return { RI::AA,  1, 'D',   7, 133.103f };
      case ID("ASX"): return { RI::AA,  1, 'B',   6, 100.096f };
      case ID("CYS"): return { RI::AA,  1, 'C',   7, 121.158f };  // also BUF
      case ID("CSH"): return { RI::AA,  1, 's',  17, 283.284f };
      case ID("GLN"): return { RI::AA,  1, 'Q',  10, 146.144f };
      case ID("GLU"): return { RI::AA,  1, 'E',   9, 147.129f };
      case ID("GLX"): return { RI::AA,  1, 'Z',   8, 114.123f };
      case ID("GLY"): return { RI::AA,  1, 'G',   5, 75.0666f };  // also BUF
      case ID("HIS"): return { RI::AA,  1, 'H',  10, 156.162f };
      case ID("ILE"): return { RI::AA,  1, 'I',  13, 131.173f };
      case ID("LEU"): return { RI::AA,  1, 'L',  13, 131.173f };
      case ID("LYS"): return { RI::AA,  1, 'K',  15, 147.196f };
      case ID("MET"): return { RI::AA,  1, 'M',  11, 149.211f };
      case ID("MSE"): return { RI::AA,  1, 'm',  11, 196.106f };
      case ID("ORN"): return { RI::AA,  1, 'a',  12, 132.161f };
      case ID("PHE"): return { RI::AA,  1, 'F',  11, 165.189f };
      case ID("PRO"): return { RI::AA,  1, 'P',   9, 115.130f };
      case ID("SER"): return { RI::AA,  1, 'S',   7, 105.093f };
      case ID("THR"): return { RI::AA,  1, 'T',   9, 119.119f };
      case ID("TRY"): // fall-through - synonym for TRP
      case ID("TRP"): return { RI::AA,  1, 'W',  12, 204.225f };
      case ID("TYR"): return { RI::AA,  1, 'Y',  11, 181.189f };
      case ID("UNK"): return { RI::AA,  1, 'X',   9, 103.120f };
      case ID("VAL"): return { RI::AA,  1, 'V',  11, 117.146f };
      case ID("SEC"): return { RI::AA,  1, 'U',   7, 168.053f };
      case ID("PYL"): return { RI::AA,  1, 'O',  21, 255.313f };

      case ID("SEP"): return { RI::AA,  1, 's',   8, 185.072f };
      case ID("TPO"): return { RI::AA,  1, 't',  10, 199.099f };
      case ID("PCA"): return { RI::AA,  1, 'e',   7, 129.114f };
      case ID("CSO"): return { RI::AA,  1, 'c',   7, 137.158f };
      case ID("PTR"): return { RI::AA,  1, 'y',  12, 261.168f };
      case ID("KCX"): return { RI::AA,  1, 'k',  14, 190.197f };
      case ID("CSD"): return { RI::AA,  1, 'c',   7, 153.157f };
      case ID("LLP"): return { RI::AA,  1, 'k',  22, 375.314f };
      case ID("CME"): return { RI::AA,  1, 'c',  11, 197.276f };
      case ID("MLY"): return { RI::AA,  1, 'k',  18, 174.241f };
      case ID("DAL"): return { RI::AAD, 1, 'a',   7, 89.0932f };
      case ID("TYS"): return { RI::AA,  1, 'y',  11, 261.252f };
      case ID("OCS"): return { RI::AA,  1, 'c',   7, 169.156f };
      case ID("M3L"): return { RI::AA,  1, 'k',  21, 189.275f };
      case ID("FME"): return { RI::AA,  1, 'm',  11, 177.221f };
      case ID("ALY"): return { RI::AA,  1, 'k',  16, 188.224f };
      case ID("HYP"): return { RI::AA,  1, 'p',   9, 131.130f };
      case ID("CAS"): return { RI::AA,  1, 'c',  12, 225.141f };
      case ID("CRO"): return { RI::AA,  1, 't',  17, 319.313f };
      case ID("CSX"): return { RI::AA,  1, 'c',   7, 137.158f };
      case ID("DPR"): return { RI::AAD, 1, 'p',   9, 115.130f };  // also BUF
      case ID("DGL"): return { RI::AAD, 1, 'e',   9, 147.129f };
      case ID("DVA"): return { RI::AAD, 1, 'v',  11, 117.146f };
      case ID("CSS"): return { RI::AA,  1, 'c',   7, 153.223f };
      case ID("DPN"): return { RI::AAD, 1, 'f',  11, 165.189f };
      case ID("DSN"): return { RI::AAD, 1, 's',   7, 105.093f };
      case ID("DLE"): return { RI::AAD, 1, 'l',  13, 131.173f };
      case ID("HIC"): return { RI::AA,  1, 'h',  11, 169.181f };
      case ID("NLE"): return { RI::AA,  1, 'l',  13, 131.173f };
      case ID("MVA"): return { RI::AA,  1, 'v',  13, 131.173f };
      case ID("MLZ"): return { RI::AA,  1, 'k',  16, 160.214f };
      case ID("CR2"): return { RI::AA,  1, 'g',  13, 275.260f };
      case ID("SAR"): return { RI::AA,  1, 'g',   7, 89.0932f };
      case ID("DAR"): return { RI::AAD, 1, 'r',  15, 175.209f };
      case ID("DLY"): return { RI::AAD, 1, 'k',  14, 146.188f };
      case ID("YCM"): return { RI::AA,  1, 'c',  10, 178.209f };
      case ID("NRQ"): return { RI::AA,  1, 'm',  17, 347.389f };
      case ID("CGU"): return { RI::AA,  1, 'e',   9, 191.139f };
      case ID("0TD"): return { RI::AA,  1, 'd',   9, 179.194f };
      case ID("MLE"): return { RI::AA,  1, 'l',  15, 145.200f };
      case ID("DAS"): return { RI::AAD, 1, 'd',   7, 133.103f };
      case ID("DTR"): return { RI::AAD, 1, 'w',  12, 204.225f };
      case ID("CXM"): return { RI::AA,  1, 'm',  11, 193.221f };
      case ID("TPQ"): return { RI::AA,  1, 'y',   9, 211.171f };
      case ID("DCY"): return { RI::AAD, 1, 'c',   7, 121.158f };
      case ID("DSG"): return { RI::AAD, 1, 'n',   8, 132.118f };
      case ID("DTY"): return { RI::AAD, 1, 'y',  11, 181.189f };
      case ID("DHI"): return { RI::AAD, 1, 'h',  10, 156.162f };
      case ID("MEN"): return { RI::AA,  1, 'n',  10, 146.144f };
      case ID("DTH"): return { RI::AAD, 1, 't',   9, 119.119f };
      case ID("SAC"): return { RI::AA,  1, 's',   9, 147.129f };
      case ID("DGN"): return { RI::AAD, 1, 'q',  10, 146.144f };
      case ID("AIB"): return { RI::AA,  1, 'a',   9, 103.120f };
      case ID("SMC"): return { RI::AA,  1, 'c',   9, 135.185f };
      case ID("IAS"): return { RI::AA,  1, 'd',   7, 133.103f };
      case ID("CIR"): return { RI::AA,  1, 'r',  13, 175.186f };
      case ID("BMT"): return { RI::AA,  1, 't',  19, 201.263f };
      case ID("DIL"): return { RI::AAD, 1, 'i',  13, 131.173f };
      case ID("FGA"): return { RI::AA,  1, 'e',   9, 147.129f };
      case ID("PHI"): return { RI::AA,  1, 'f',  10, 291.086f };
      case ID("CRQ"): return { RI::AA,  1, 'q',  16, 344.322f };
      case ID("SME"): return { RI::AA,  1, 'm',  11, 165.211f };
      case ID("GHP"): return { RI::AA,  1, 'g',   9, 167.162f };  // d-peptide in CCD
      case ID("MHO"): return { RI::AA,  1, 'm',  11, 165.211f };
      case ID("NEP"): return { RI::AA,  1, 'h',  10, 235.134f };
      case ID("TRQ"): return { RI::AA,  1, 'w',  10, 234.208f };
      case ID("TOX"): return { RI::AA,  1, 'w',  12, 236.224f };
      case ID("ALC"): return { RI::AA,  1, 'a',  17, 171.237f };
      case ID("3FG"): return { RI::AA,  1, ' ',   9, 183.161f };
      case ID("SCH"): return { RI::AA,  1, 'c',   9, 167.250f };
      case ID("MDO"): return { RI::AA,  1, 'a',  11, 197.191f };
      case ID("MAA"): return { RI::AA,  1, 'a',   9, 103.120f };
      case ID("GYS"): return { RI::AA,  1, 's',  15, 305.286f };
      case ID("MK8"): return { RI::AA,  1, 'l',  15, 145.200f };
      case ID("CR8"): return { RI::AA,  1, 'h',  16, 354.340f };
      case ID("KPI"): return { RI::AA,  1, 'k',  16, 216.234f };
      case ID("SCY"): return { RI::AA,  1, 'c',   9, 163.195f };
      case ID("DHA"): return { RI::AA,  1, 's',   5, 87.0773f };
      case ID("OMY"): return { RI::AA,  1, 'y',  10, 231.633f };
      case ID("CAF"): return { RI::AA,  1, 'c',  12, 241.140f };
      case ID("0AF"): return { RI::AA,  1, 'w',  12, 220.225f };
      case ID("SNN"): return { RI::AA,  1, 'n',   6, 114.103f };
      case ID("MHS"): return { RI::AA,  1, 'h',  11, 169.181f };
      case ID("MLU"): return { RI::AAD, 1, ' ',  15, 145.200f };
      case ID("SNC"): return { RI::AA,  1, 'c',   6, 150.156f };
      case ID("PHD"): return { RI::AA,  1, 'd',   8, 213.083f };
      case ID("B3E"): return { RI::AA,  1, 'e',  11, 161.156f };
      case ID("MEA"): return { RI::AA,  1, 'f',  13, 179.216f };
      case ID("MED"): return { RI::AAD, 1, 'm',  11, 149.211f };
      case ID("OAS"): return { RI::AA,  1, 's',   9, 147.129f };
      case ID("GL3"): return { RI::AA,  1, 'g',   5, 91.1322f };
      case ID("FVA"): return { RI::AA,  1, 'v',  11, 145.156f };
      case ID("PHL"): return { RI::AA,  1, 'f',  13, 151.206f };
      case ID("CRF"): return { RI::AA,  1, 't',  18, 342.349f };
      case ID("OMZ"): return { RI::AA,  1, ' ',  10, 231.633f };  // d-peptide in CCD
      case ID("BFD"): return { RI::AA,  1, 'd',   6, 198.102f };
      case ID("MEQ"): return { RI::AA,  1, 'q',  12, 160.171f };
      case ID("DAB"): return { RI::AA,  1, 'a',  10, 118.134f };
      case ID("AGM"): return { RI::AA,  1, 'r',  17, 189.235f };

      case ID("PSU"): return { RI::RNA, 2, 'u',  13, 324.181f };
      case ID("5MU"): return { RI::RNA, 2, 'u',  15, 338.208f };
      case ID("7MG"): return { RI::RNA, 2, 'g',  18, 379.263f };
      case ID("OMG"): return { RI::RNA, 2, 'g',  16, 377.247f };
      case ID("UR3"): return { RI::RNA, 2, 'u',  15, 338.208f };
      case ID("OMC"): return { RI::RNA, 2, 'c',  16, 337.223f };
      case ID("2MG"): return { RI::RNA, 2, 'g',  16, 377.247f };
      case ID("H2U"): return { RI::RNA, 2, 'u',  15, 326.197f };
      case ID("4SU"): return { RI::RNA, 2, 'u',  13, 340.247f };
      case ID("OMU"): return { RI::RNA, 2, 'u',  15, 338.208f };
      case ID("4OC"): return { RI::RNA, 2, 'c',  18, 351.250f };
      case ID("MA6"): return { RI::RNA, 2, 'a',  18, 375.274f };
      case ID("M2G"): return { RI::RNA, 2, 'g',  18, 391.274f };
      case ID("1MA"): return { RI::RNA, 2, 'a',  16, 361.248f };
      case ID("6MZ"): return { RI::RNA, 2, 'a',  16, 361.248f };
      case ID("CCC"): return { RI::RNA, 2, 'c',  13, 385.161f };
      case ID("2MA"): return { RI::RNA, 2, 'a',  16, 361.248f };
      case ID("1MG"): return { RI::RNA, 2, 'g',  16, 377.247f };
      case ID("5BU"): return { RI::RNA, 2, 'u',  12, 403.077f };
      case ID("MIA"): return { RI::RNA, 2, 'a',  24, 461.430f };
      case ID("DOC"): return { RI::DNA, 2, 'c',  14, 291.198f };
      case ID("8OG"): return { RI::DNA, 2, 'g',  14, 363.221f };
      case ID("5CM"): return { RI::DNA, 2, 'c',  16, 321.224f };
      case ID("3DR"): return { RI::DNA, 2, ' ',  11, 198.111f };
      case ID("BRU"): return { RI::DNA, 2, 'u',  12, 387.078f };
      case ID("CBR"): return { RI::DNA, 2, 'c',  13, 386.093f };

      case ID("H2O"): // fall-through - synonym ' ',for HOH
      case ID("WAT"): // fall-through - synonym for HOH
      case ID("HOH"): return { RI::HOH, 0, ' ',   2, 18.0153f };
      case ID("DOD"): return { RI::HOH, 0, ' ',   2, 20.0276f };
      case ID("HEM"): return { RI::ELS, 0, ' ',  32, 616.487f };
      case ID("SUL"): // fall-through - synonym for SO4
      case ID("SO4"): return { RI::BUF, 0, ' ',   0, 96.0626f };
      case ID("GOL"): return { RI::BUF, 0, ' ',   8, 92.0938f };
      case ID("EDO"): return { RI::BUF, 0, ' ',   6, 62.0678f };
      case ID("NAG"): return { RI::PYR, 0, ' ',  15, 221.208f };
      case ID("PO4"): return { RI::ELS, 0, ' ',   0, 94.9714f };
      case ID("ACT"): return { RI::BUF, 0, ' ',   3, 59.0440f };
      case ID("PEG"): return { RI::ELS, 0, ' ',  10, 106.120f };
      case ID("MAN"): return { RI::PYR, 0, ' ',  12, 180.156f };  // also BUF
      case ID("FAD"): return { RI::ELS, 0, ' ',  33, 785.550f };
      case ID("BMA"): return { RI::PYR, 0, ' ',  12, 180.156f };  // also BUF
      case ID("ADP"): return { RI::ELS, 0, ' ',  15, 427.201f };
      case ID("DMS"): return { RI::BUF, 0, ' ',   6, 78.1334f };
      // ACE and NH2 are non-polymers present primarily in peptide chains
      case ID("ACE"): return { RI::ELS, 1, ' ',   4, 44.0526f };
      case ID("NH2"): return { RI::ELS, 1, ' ',   2, 16.0226f };  // ?
      case ID("MPD"): return { RI::BUF, 0, ' ',  14, 118.174f };
      case ID("MES"): return { RI::ELS, 0, ' ',  13, 195.237f };
      case ID("NAD"): return { RI::ELS, 0, ' ',  27, 663.425f };
      case ID("NAP"): return { RI::ELS, 0, ' ',  28, 743.405f };
      case ID("TRS"): return { RI::BUF, 0, ' ',  12, 122.143f };
      case ID("ATP"): return { RI::ELS, 0, ' ',  16, 507.181f };
      case ID("PG4"): return { RI::ELS, 0, ' ',  18, 194.226f };
      case ID("GDP"): return { RI::ELS, 2, 'g',  15, 443.201f };  // RNA in CCD
      case ID("FUC"): return { RI::PYR, 0, ' ',  12, 164.156f };
      case ID("FMT"): return { RI::BUF, 0, ' ',   2, 46.0254f };
      case ID("GAL"): return { RI::PYR, 0, ' ',  12, 180.156f };
      case ID("PGE"): return { RI::BUF, 0, ' ',  14, 150.173f };
      case ID("FMN"): return { RI::ELS, 0, ' ',  21, 456.344f };
      case ID("PLP"): return { RI::ELS, 0, ' ',  10, 247.142f };
      case ID("EPE"): return { RI::ELS, 0, ' ',  18, 238.305f };
      case ID("SF4"): return { RI::ELS, 0, ' ',   0, 351.640f };
      case ID("BME"): return { RI::ELS, 0, ' ',   6, 78.1334f };
      case ID("CIT"): return { RI::BUF, 0, ' ',   8, 192.124f };

      case ID("BE7"): return { RI::BUF, 0, ' ',   5, 357.156f };
      case ID("MRD"): return { RI::BUF, 0, ' ',  14, 118.174f };
      case ID("MHA"): return { RI::BUF, 0, ' ',  10, 190.154f };
      case ID("BU3"): return { RI::BUF, 0, ' ',  10, 90.1210f };
      case ID("PGO"): return { RI::BUF, 0, ' ',   8, 76.0944f };
      case ID("BU2"): return { RI::BUF, 0, ' ',  10, 90.1210f };
      case ID("PDO"): return { RI::BUF, 0, ' ',   8, 76.0944f };
      case ID("BU1"): return { RI::BUF, 0, ' ',  10, 90.1210f };
      case ID("PG6"): return { RI::BUF, 0, ' ',  26, 266.331f };
      case ID("1BO"): return { RI::BUF, 0, ' ',  10, 74.1216f };
      case ID("PE7"): return { RI::BUF, 0, ' ',  30, 342.449f };
      case ID("PG5"): return { RI::BUF, 0, ' ',  18, 178.226f };
      case ID("TFP"): return { RI::BUF, 0, ' ',  24, 407.496f };
      case ID("DHD"): return { RI::BUF, 0, ' ',   4, 160.082f };
      case ID("PEU"): return { RI::BUF, 0, ' ', 112, 1221.46f };
      case ID("TAU"): return { RI::BUF, 0, ' ',   7, 125.147f };
      case ID("SBT"): return { RI::BUF, 0, ' ',  10, 74.1216f };
      case ID("SAL"): return { RI::BUF, 0, ' ',   6, 138.121f };
      case ID("IOH"): return { RI::BUF, 0, ' ',   8, 60.0950f };
      case ID("IPA"): return { RI::BUF, 0, ' ',   8, 60.0950f };
      case ID("PIG"): return { RI::BUF, 0, ' ',  14, 150.173f };
      case ID("B3P"): return { RI::BUF, 0, ' ',  26, 282.334f };
      case ID("BTB"): return { RI::BUF, 0, ' ',  19, 209.240f };
      case ID("NHE"): return { RI::BUF, 0, ' ',  17, 207.290f };
      case ID("C8E"): return { RI::BUF, 0, ' ',  34, 306.438f };
      case ID("OTE"): return { RI::BUF, 0, ' ',  34, 306.438f };
      case ID("PE4"): return { RI::BUF, 0, ' ',  34, 354.436f };
      case ID("XPE"): return { RI::BUF, 0, ' ',  42, 458.541f };
      case ID("PE8"): return { RI::BUF, 0, ' ',  34, 370.436f };
      case ID("P33"): return { RI::BUF, 0, ' ',  30, 326.383f };
      case ID("N8E"): return { RI::BUF, 0, ' ',  38, 350.491f };
      case ID("2OS"): return { RI::BUF, 0, ' ',  36, 468.493f };
      case ID("1PS"): return { RI::BUF, 0, ' ',  11, 201.243f };
      case ID("CPS"): return { RI::BUF, 0, ' ',  58, 614.877f };
      case ID("DMX"): return { RI::BUF, 0, ' ',  19, 257.349f };
      case ID("MPO"): return { RI::BUF, 0, ' ',  15, 209.263f };
      case ID("GCD"): return { RI::PYR, 0, ' ',   8, 176.124f };
      case ID("DXG"): return { RI::BUF, 0, ' ',   8, 192.124f };
      case ID("CM5"): return { RI::BUF, 0, ' ',  42, 494.573f };
      case ID("ACA"): return { RI::BUF, 1, ' ',  13, 131.173f }; // peptide linking
      case ID("ACN"): return { RI::BUF, 0, ' ',   6, 58.0791f };
      case ID("CCN"): return { RI::BUF, 0, ' ',   3, 41.0519f };
      case ID("GLC"): return { RI::PYR, 0, ' ',  12, 180.156f };
      case ID("DR6"): return { RI::BUF, 0, ' ', 142, 1527.90f };
      case ID("NH4"): return { RI::BUF, 0, ' ',   4, 18.0385f };
      case ID("AZI"): return { RI::BUF, 0, ' ',   0, 42.0201f };
      case ID("BNG"): return { RI::PYR, 0, ' ',  30, 306.395f };
      case ID("BOG"): return { RI::PYR, 0, ' ',  28, 292.369f };
      case ID("BGC"): return { RI::PYR, 0, ' ',  12, 180.156f };
      case ID("BCN"): return { RI::BUF, 0, ' ',  13, 163.172f };
      case ID("BRO"): return { RI::BUF, 0, ' ',   0, 79.9040f };
      case ID("CAC"): return { RI::BUF, 0, ' ',   6, 136.989f };
      case ID("CBX"): return { RI::BUF, 0, ' ',   2, 46.0254f };
      case ID("ACY"): return { RI::BUF, 0, ' ',   4, 60.0520f };
      case ID("CBM"): return { RI::BUF, 0, ' ',   4, 60.0520f };
      case ID("CLO"): return { RI::BUF, 0, ' ',   0, 35.4530f };
      case ID("3CO"): return { RI::BUF, 0, ' ',   0, 58.9332f };
      case ID("NCO"): return { RI::BUF, 0, ' ',  18, 161.116f };
      case ID("CU1"): return { RI::BUF, 0, ' ',   0, 63.5460f };
      case ID("CYN"): return { RI::BUF, 0, ' ',   0, 26.0174f };
      case ID("MA4"): return { RI::BUF, 0, ' ',  44, 508.600f };
      case ID("TAR"): return { RI::BUF, 0, ' ',   6, 150.087f };
      case ID("GLO"): return { RI::BUF, 0, ' ',  12, 180.156f };  // d-saccharide
      case ID("MTL"): return { RI::BUF, 0, ' ',  14, 182.172f };
      case ID("SOR"): return { RI::BUF, 0, ' ',  14, 182.172f };
      case ID("DMU"): return { RI::BUF, 0, ' ',  42, 482.562f };  // d-saccharide
      case ID("DDQ"): return { RI::BUF, 0, ' ',  27, 201.349f };
      case ID("DMF"): return { RI::BUF, 0, ' ',   7, 73.0938f };
      case ID("DIO"): return { RI::BUF, 0, ' ',   8, 88.1051f };
      case ID("DOX"): return { RI::BUF, 0, ' ',   8, 88.1051f };
      case ID("12P"): return { RI::BUF, 0, ' ',  50, 546.646f };
      case ID("SDS"): return { RI::BUF, 0, ' ',  26, 266.397f };
      case ID("LMT"): return { RI::BUF, 0, ' ',  46, 510.615f };  // d-saccharide
      case ID("EOH"): return { RI::BUF, 0, ' ',   6, 46.0684f };
      case ID("EEE"): return { RI::BUF, 0, ' ',   8, 88.1051f };
      case ID("EGL"): return { RI::BUF, 0, ' ',   6, 62.0678f };
      case ID("FLO"): return { RI::BUF, 0, ' ',   0, 18.9984f };
      case ID("TRT"): return { RI::BUF, 0, ' ',  36, 352.508f };
      case ID("FCY"): return { RI::BUF, 0, ' ',   7, 121.158f };
      case ID("FRU"): return { RI::BUF, 0, ' ',  12, 180.156f };  // saccharide
      case ID("GBL"): return { RI::BUF, 0, ' ',   6, 86.0892f };
      case ID("GPX"): return { RI::BUF, 0, ' ',  14, 505.165f };
      case ID("HTO"): return { RI::BUF, 0, ' ',  16, 148.200f };
      case ID("HTG"): return { RI::BUF, 0, ' ',  26, 294.408f };
      case ID("B7G"): return { RI::BUF, 0, ' ',  26, 278.342f };
      case ID("C10"): return { RI::BUF, 0, ' ',  46, 422.596f };
      case ID("16D"): return { RI::BUF, 0, ' ',  16, 116.205f };
      case ID("HEZ"): return { RI::BUF, 0, ' ',  14, 118.174f };
      case ID("IOD"): return { RI::BUF, 0, ' ',   0, 126.904f };
      case ID("IDO"): return { RI::BUF, 0, ' ',   0, 126.904f };
      case ID("ICI"): return { RI::BUF, 0, ' ',   8, 192.124f };
      case ID("ICT"): return { RI::BUF, 0, ' ',   8, 192.124f };
      case ID("TLA"): return { RI::BUF, 0, ' ',   6, 150.087f };
      case ID("LAT"): return { RI::BUF, 0, ' ',  22, 342.296f };  // saccharide
      case ID("LBT"): return { RI::BUF, 0, ' ',  22, 342.296f };  // saccharide
      case ID("LDA"): return { RI::BUF, 0, ' ',  31, 229.402f };
      case ID("MN3"): return { RI::BUF, 0, ' ',   0, 54.9380f };
      case ID("MRY"): return { RI::BUF, 0, ' ',  10, 122.120f };
      case ID("MOH"): return { RI::BUF, 0, ' ',   4, 32.0419f };
      case ID("BEQ"): return { RI::BUF, 0, ' ',  38, 342.517f };
      case ID("C15"): return { RI::BUF, 0, ' ',  38, 336.554f };
      case ID("MG8"): return { RI::BUF, 0, ' ',  31, 321.410f };
      case ID("POL"): return { RI::BUF, 0, ' ',   8, 60.0950f };
      case ID("NO3"): return { RI::BUF, 0, ' ',   0, 62.0049f };
      case ID("JEF"): return { RI::BUF, 0, ' ',  63, 597.822f };
      case ID("P4C"): return { RI::BUF, 0, ' ',  28, 324.367f };
      case ID("CE1"): return { RI::BUF, 0, ' ',  58, 538.755f };
      case ID("DIA"): return { RI::BUF, 0, ' ',  20, 144.258f };
      case ID("CXE"): return { RI::BUF, 0, ' ',  42, 378.544f };
      case ID("IPH"): return { RI::BUF, 0, ' ',   6, 94.1112f };
      case ID("PIN"): return { RI::BUF, 0, ' ',  18, 302.368f };
      case ID("15P"): return { RI::BUF, 0, ' ', 140, 1529.83f };
      case ID("CRY"): return { RI::BUF, 0, ' ',   8, 92.0938f };
      case ID("PGR"): return { RI::BUF, 0, ' ',   8, 76.0944f };
      case ID("PGQ"): return { RI::BUF, 0, ' ',   8, 76.0944f };
      case ID("SPD"): return { RI::BUF, 0, ' ',  19, 145.246f };
      case ID("SPK"): return { RI::BUF, 0, ' ',  30, 206.372f };
      case ID("SPM"): return { RI::BUF, 0, ' ',  26, 202.340f };
      case ID("SUC"): return { RI::PYR, 0, ' ',  22, 342.296f };
      case ID("TBU"): return { RI::BUF, 0, ' ',  10, 74.1216f };
      case ID("TMA"): return { RI::BUF, 0, ' ',  12, 74.1448f };
      case ID("TEP"): return { RI::BUF, 0, ' ',   8, 180.164f };
      case ID("SCN"): return { RI::BUF, 0, ' ',   0, 58.0824f };
      case ID("TRE"): return { RI::PYR, 0, ' ',  22, 342.296f };
      case ID("ETF"): return { RI::BUF, 0, ' ',   3, 100.040f };
      case ID("144"): return { RI::BUF, 0, ' ',  12, 122.143f };
      case ID("UMQ"): return { RI::BUF, 0, ' ',  44, 496.589f };
      case ID("URE"): return { RI::BUF, 0, ' ',   4, 60.0553f };
      case ID("YT3"): return { RI::BUF, 0, ' ',   0, 88.9059f };
      case ID("ZN2"): return { RI::BUF, 0, ' ',   0, 65.3800f };
      case ID("FE2"): return { RI::BUF, 0, ' ',   0, 55.8450f };
      case ID("3NI"): return { RI::BUF, 0, ' ',   0, 58.6934f };
    }
#undef ID
  } else if (name.size() == 1) {
    switch (name[0]) {
      case 'A': return { RI::RNA, 2, 'A',  14, 347.221f };
      case 'C': return { RI::RNA, 2, 'C',  14, 323.197f };
      case 'G': return { RI::RNA, 2, 'G',  14, 363.221f };
      case 'I': return { RI::RNA, 2, 'I',  13, 348.206f };
      case 'U': return { RI::RNA, 2, 'U',  13, 324.181f };
      case 'N': return { RI::RNA, 2, 'N',  11, 214.11f };

      case 'F': return { RI::BUF, 0, ' ',   0, 18.9984f };
      case 'K': return { RI::BUF, 0, ' ',   0, 39.0983f };
    }
  } else if (name.size() == 2) {
    if (name[0] == 'D' || name[0] == '+')
      switch (name[1]) {
        case 'A': return { RI::DNA, 2, 'A',  14, 331.222f };
        case 'C': return { RI::DNA, 2, 'C',  14, 307.197f };
        case 'G': return { RI::DNA, 2, 'G',  14, 347.221f };
        case 'I': return { RI::DNA, 2, 'I',  13, 332.207f };
        case 'T': return { RI::DNA, 2, 'T',  15, 322.208f };
        case 'U': return { RI::DNA, 2, 'U',  13, 308.182f };
        case 'N': return { RI::DNA, 2, 'N',  14, 198.111f };  // unknown DNA
      }
    else
#define ID(s) (s[0] << 8 | s[1])
      switch (ID(name.c_str())) {
        case ID("AG"): return { RI::BUF, 0, ' ',   0, 107.868f };
        case ID("AL"): return { RI::BUF, 0, ' ',   0, 26.9815f };
        case ID("BA"): return { RI::BUF, 0, ' ',   0, 137.327f };
        case ID("BR"): return { RI::BUF, 0, ' ',   0, 79.9040f };
        case ID("CA"): return { RI::BUF, 0, ' ',   0, 40.0780f };
        case ID("CD"): return { RI::BUF, 0, ' ',   0, 112.411f };
        case ID("CL"): return { RI::BUF, 0, ' ',   0, 35.4530f };
        case ID("CM"): return { RI::BUF, 0, ' ',   4, 60.0520f };
        case ID("CN"): return { RI::BUF, 0, ' ',   0, 27.0253f };
        case ID("CO"): return { RI::BUF, 0, ' ',   0, 58.9332f };
        case ID("CS"): return { RI::BUF, 0, ' ',   0, 132.905f };
        case ID("CU"): return { RI::BUF, 0, ' ',   0, 63.5460f };
        case ID("FE"): return { RI::BUF, 0, ' ',   0, 55.8450f };
        case ID("HG"): return { RI::BUF, 0, ' ',   0, 200.590f };
        case ID("LI"): return { RI::BUF, 0, ' ',   0, 6.94100f };
        case ID("MG"): return { RI::BUF, 0, ' ',   0, 24.3050f };
        case ID("MN"): return { RI::BUF, 0, ' ',   0, 54.9380f };
        case ID("NA"): return { RI::BUF, 0, ' ',   0, 22.9898f };
        case ID("NI"): return { RI::BUF, 0, ' ',   0, 58.6934f };
        case ID("NO"): return { RI::ELS, 0, ' ',   0, 30.0061f };
        case ID("PB"): return { RI::BUF, 0, ' ',   0, 207.200f };
        case ID("RB"): return { RI::BUF, 0, ' ',   0, 85.4678f };
        case ID("SR"): return { RI::BUF, 0, ' ',   0, 87.6200f };
        case ID("Y1"): return { RI::BUF, 0, ' ',   0, 88.9059f };
        case ID("ZN"): return { RI::BUF, 0, ' ',   0, 65.3800f };
      }
#undef ID
  }
  return { RI::UNKNOWN, 0, ' ', 0, 0.0f };
}

std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
                                                    ResidueKind kind) {
  std::vector<std::string> r;
  r.reserve(seq.size());
  auto kind_str = [&]() {
    switch (kind) {
      case ResidueKind::AA: return "peptide";
      case ResidueKind::RNA: return "RNA";
      case ResidueKind::DNA: return "DNA";
      default: return "unknown";
    }
  };
  for (size_t i = 0; i != seq.size(); ++i) {
    char c = seq[i];
    if (is_space(c))
      continue;
    if (c == '(') { // special case, e.g. (MSE)
      size_t start = i + 1;
      i = seq.find(')', start);
      if (i == std::string::npos)
        fail("unmatched '(' in sequence");
      r.emplace_back(seq, start, i - start);
    } else {
      const char* str = expand_one_letter(c, kind);
      if (str == nullptr)
        fail("unexpected letter in ", kind_str(), " sequence: ", c,
             " (", std::to_string(int(c)), ')');
      r.emplace_back(str);
    }
  }
  return r;
}

} // namespace gemmi
