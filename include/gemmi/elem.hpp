// Copyright 2017 Global Phasing Ltd.
//
// Elements from the periodic table.

#ifndef GEMMI_ELEM_HPP_
#define GEMMI_ELEM_HPP_

#include <string>

namespace gemmi {

// elements
enum class El : unsigned char {
  X=0,  // unknown element is marked as X in PDB entries
  H=1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,  // 1-3
  K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,  // 4
  Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,  // 5
  Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,  // 6..
  Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,  // ..6
  Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,  // 7..
  Rf, Db, Sg, Bh, Hs, Mt,  Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og, // ..7
  D, // heh, what should we do with Deuterium?
  END
};

inline bool is_hydrogen(El el) { return el == El::H || el == El::D; }

// somewhat arbitrary division into metals and non-metals
inline bool is_metal(El el) {
  static constexpr bool table[] = {
    // X     H     He
    false, false, false,
    // Li  Be     B      C      N      O      F     Ne
    true, true, false, false, false, false, false, false,
    // Na  Mg    Al     Si     P      S      Cl     Ar
    true, true, true, false, false, false, false, false,
    // K   Ca    Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
    true, true, true, true, true, true, true, true, true, true, true, true,
    // Ga  Ge    As     Se     Br     Kr
    true, true, false, false, false, false,
    // Rb  Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
    true, true, true, true, true, true, true, true, true, true, true, true,
    // In  Sn    Sb    Te      I     Xe
    true, true, true, false, false, false,
    // Cs  Ba    La    Ce    Pr    Nd    Pm    Sm    Eu    Gd    Tb    Dy
    true, true, true, true, true, true, true, true, true, true, true, true,
    // Ho  Er    Tm    Yb    Lu    Hf    Ta    W     Re    Os    Ir    Pt
    true, true, true, true, true, true, true, true, true, true, true, true,
    // Au  Hg    Tl    Pb    Bi    Po    At     Rn
    true, true, true, true, true, true, false, false,
    // Fr  Ra    Ac    Th    Pa    U     Np    Pu    Am    Cm    Bk    Cf
    true, true, true, true, true, true, true, true, true, true, true, true,
    // Es  Fm    Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
    true, true, true, true, true, true, true, true, true, true, true, true,
    // Rg  Cn    Nh    Fl    Mc    Lv    Ts     Og
    true, true, true, true, true, true, false, false,
    // D    END
    false, false
  };
  static_assert(table[static_cast<int>(El::D)] == false, "Hmm");
  static_assert(sizeof(table) / sizeof(table[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return table[static_cast<int>(el)];
}


inline double molecular_weight(El el) {
  static constexpr double weights[] = {
    /*X*/ 0.0,
    /*H*/ 1.00794, /*He*/ 4.0026,
    /*Li*/ 6.941, /*Be*/ 9.012182, /*B*/ 10.811, /*C*/ 12.0107,
    /*N*/ 14.0067, /*O*/ 15.9994, /*F*/ 18.998403, /*Ne*/ 20.1797,
    /*Na*/ 22.98977, /*Mg*/ 24.305, /*Al*/ 26.981539, /*Si*/ 28.0855,
    /*P*/ 30.973761, /*S*/ 32.065, /*Cl*/ 35.453, /*Ar*/ 39.948,
    /*K*/ 39.0983, /*Ca*/ 40.078, /*Sc*/ 44.95591, /*Ti*/ 47.867,
    /*V*/ 50.9415, /*Cr*/ 51.9961, /*Mn*/ 54.93805, /*Fe*/ 55.845,
    /*Co*/ 58.9332, /*Ni*/ 58.6934, /*Cu*/ 63.546, /*Zn*/ 65.38,
    /*Ga*/ 69.723, /*Ge*/ 72.64, /*As*/ 74.9216, /*Se*/ 78.96,
    /*Br*/ 79.904, /*Kr*/ 83.798, /*Rb*/ 85.4678, /*Sr*/ 87.62,
    /*Y*/ 88.90585, /*Zr*/ 91.224, /*Nb*/ 92.9064,
    /*Mo*/ 95.95, /*Tc*/ 98, /*Ru*/ 101.07, /*Rh*/ 102.9055, /*Pd*/ 106.42,
    /*Ag*/ 107.8682, /*Cd*/ 112.411, /*In*/ 114.818, /*Sn*/ 118.71,
    /*Sb*/ 121.76, /*Te*/ 127.6, /*I*/ 126.90447, /*Xe*/ 131.293,
    /*Cs*/ 132.905, /*Ba*/ 137.327, /*La*/ 138.905, /*Ce*/ 140.116,
    /*Pr*/ 140.908, /*Nd*/ 144.24, /*Pm*/ 145, /*Sm*/ 150.36,
    /*Eu*/ 151.964, /*Gd*/ 157.25, /*Tb*/ 158.925, /*Dy*/ 162.5,
    /*Ho*/ 164.93, /*Er*/ 167.259, /*Tm*/ 168.934, /*Yb*/ 173.05,
    /*Lu*/ 174.967, /*Hf*/ 178.49, /*Ta*/ 180.948, /*W*/ 183.84,
    /*Re*/ 186.207, /*Os*/ 190.23, /*Ir*/ 192.217, /*Pt*/ 195.084,
    /*Au*/ 196.967, /*Hg*/ 200.59, /*Tl*/ 204.383,
    /*Pb*/ 207.2, /*Bi*/ 208.98, /*Po*/ 209, /*At*/ 210, /*Rn*/ 222,
    /*Fr*/ 223, /*Ra*/ 226, /*Ac*/ 227, /*Th*/ 232.038, /*Pa*/ 231.036,
    /*U*/ 238.029, /*Np*/ 237, /*Pu*/ 244, /*Am*/ 243, /*Cm*/ 247,
    /*Bk*/ 247, /*Cf*/ 251, /*Es*/ 252, /*Fm*/ 257, /*Md*/ 258,
    /*No*/ 259, /*Lr*/ 262, /*Rf*/ 267, /*Db*/ 268, /*Sg*/ 271,
    /*Bh*/ 272, /*Hs*/ 270, /*Mt*/ 276, /*Ds*/ 281, /*Rg*/ 280, /*Cn*/ 285,
    /*Nh*/ 284, /*Fl*/ 289, /*Mc*/ 288, /*Lv*/ 293, /*Ts*/ 294, /*Og*/ 294,
    /*D*/ 2.0141, /*END*/ 0.0
  };
  static_assert(weights[static_cast<int>(El::D)] == 2.0141, "Hmm");
  static_assert(sizeof(weights) / sizeof(weights[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return weights[static_cast<int>(el)];
}

inline float covalent_radius(El el) {
  static constexpr float radii[] = {
    /*X*/ 0.50f,
    /*H*/ 0.23f, /*He*/ 0.93f,
    /*Li*/ 0.68f, /*Be*/ 0.35f, /*B*/ 0.83f, /*C*/ 0.68f, /*N*/ 0.68f,
    /*O*/ 0.68f, /*F*/ 0.64f, /*Ne*/ 1.12f,
    /*Na*/ 0.97f, /*Mg*/ 1.10f, /*Al*/ 1.35f, /*Si*/ 1.20f, /*P*/ 0.75f,
    /*S*/ 1.02f, /*Cl*/ 0.99f, /*Ar*/ 1.57f,
    /*K*/ 1.33f, /*Ca*/ 0.99f, /*Sc*/ 1.44f, /*Ti*/ 1.47f, /*V*/ 1.33f,
    /*Cr*/ 1.35f, /*Mn*/ 1.35f, /*Fe*/ 1.34f, /*Co*/ 1.33f, /*Ni*/ 1.50f,
    /*Cu*/ 1.52f, /*Zn*/ 1.45f, /*Ga*/ 1.22f, /*Ge*/ 1.17f, /*As*/ 1.21f,
    /*Se*/ 1.22f, /*Br*/ 1.21f, /*Kr*/ 1.91f,
    /*Rb*/ 1.47f, /*Sr*/ 1.12f, /*Y*/ 1.78f, /*Zr*/ 1.56f, /*Nb*/ 1.48f,
    /*Mo*/ 1.47f, /*Tc*/ 1.35f, /*Ru*/ 1.40f, /*Rh*/ 1.45f, /*Pd*/ 1.50f,
    /*Ag*/ 1.59f, /*Cd*/ 1.69f, /*In*/ 1.63f, /*Sn*/ 1.46f, /*Sb*/ 1.46f,
    /*Te*/ 1.47f, /*I*/ 1.40f, /*Xe*/ 1.98f,
    /*Cs*/ 1.67f, /*Ba*/ 1.34f, /*La*/ 1.87f, /*Ce*/ 1.83f, /*Pr*/ 1.82f,
    /*Nd*/ 1.81f, /*Pm*/ 1.80f, /*Sm*/ 1.80f, /*Eu*/ 1.99f, /*Gd*/ 1.79f,
    /*Tb*/ 1.76f, /*Dy*/ 1.75f, /*Ho*/ 1.74f, /*Er*/ 1.73f, /*Tm*/ 1.72f,
    /*Yb*/ 1.94f, /*Lu*/ 1.72f, /*Hf*/ 1.57f, /*Ta*/ 1.43f, /*W*/ 1.37f,
    /*Re*/ 1.35f, /*Os*/ 1.37f, /*Ir*/ 1.32f, /*Pt*/ 1.50f, /*Au*/ 1.50f,
    /*Hg*/ 1.70f, /*Tl*/ 1.55f, /*Pb*/ 1.54f, /*Bi*/ 1.54f, /*Po*/ 1.68f,
    /*At*/ 1.70f, /*Rn*/ 2.40f,
    /*Fr*/ 2.00f, /*Ra*/ 1.90f, /*Ac*/ 1.88f, /*Th*/ 1.79f, /*Pa*/ 1.61f,
    /*U*/ 1.58f, /*Np*/ 1.55f, /*Pu*/ 1.53f, /*Am*/ 1.51f, /*Cm*/ 1.50f,
    /*Bk*/ 1.50f, /*Cf*/ 1.50f, /*Es*/ 1.50f, /*Fm*/ 1.50f, /*Md*/ 1.50f,
    /*No*/ 1.50f, /*Lr*/ 1.50f, /*Rf*/ 1.57f, /*Db*/ 1.49f, /*Sg*/ 1.43f,
    /*Bh*/ 1.41f, /*Hs*/ 1.34f, /*Mt*/ 1.29f, /*Ds*/ 1.28f, /*Rg*/ 1.21f,
    /*Cn*/ 1.22f, /*Nh*/ 1.50f, /*Fl*/ 1.50f, /*Mc*/ 1.50f, /*Lv*/ 1.50f,
    /*Ts*/ 1.50f, /*Og*/ 1.50f,
    /*D*/ 0.23f, /*END*/ 0.0f
  };
  static_assert(radii[static_cast<int>(El::D)] == 0.23f, "Hmm");
  static_assert(sizeof(radii) / sizeof(radii[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return radii[static_cast<int>(el)];
}

typedef const char elname_t[3];

inline const char* element_name(El el) {
  static constexpr elname_t names[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
    "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    "D", ""
  };
  static_assert(static_cast<int>(El::Og) == 118, "Hmm");
  static_assert(names[118][0] == 'O', "Hmm");
  static_assert(sizeof(names) / sizeof(names[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return names[static_cast<int>(el)];
}

inline elname_t& element_uppercase_name(El el) {
  static constexpr elname_t names[] = {
    "X",  "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O", "F", "NE",
    "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR",
    "K",  "CA", "SC", "TI", "V",  "CR", "MN", "FE", "CO",
    "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
    "RB", "SR", "Y",  "ZR", "NB", "MO", "TC", "RU", "RH",
    "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
    "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU",
    "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
    "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN",
    "FR", "RA", "AC", "TH", "PA", "U",  "NP", "PU", "AM",
    "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR",
    "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN",
    "NH", "FL", "MC", "LV", "TS", "OG",
    "D", "", ""
  };
  static_assert(sizeof(names) / sizeof(names[0]) == 122, "not 122");
  return names[static_cast<int>(el)];
}


namespace impl {
// all the most common elements in the PDB are single-letter
inline El find_single_letter_element(char c) {
  switch (c) {
    case 'H': return El::H;
    case 'B': return El::B;
    case 'C': return El::C;
    case 'N': return El::N;
    case 'O': return El::O;
    case 'F': return El::F;
    case 'P': return El::P;
    case 'S': return El::S;
    case 'K': return El::K;
    case 'V': return El::V;
    case 'Y': return El::Y;
    case 'I': return El::I;
    case 'W': return El::W;
    case 'U': return El::U;
    case 'D': return El::D;
    default: return El::X;
  }
}
} // namespace impl

inline El find_element(const char* symbol) {
  if (symbol == nullptr || symbol[0] == '\0')
    return El::X;
  char first = symbol[0] & ~0x20;  // lower -> upper, space -> NUL
  char second = symbol[1] & ~0x20;
  if (first == '\0')
    return impl::find_single_letter_element(second);
  // To handle symbol being "S\n" we have the condition below.
  // In addition to \t, \v, \r and \n it catches also !"3$%&'()*+,- and
  // some control characters - inconsistent but not necessarily bad.
  if (second < 14)
    return impl::find_single_letter_element(first);
  elname_t* names = &element_uppercase_name(El::X);
  for (int i = 0; i != 120; ++i) {
    if (names[i][0] == first && names[i][1] == second)
      return static_cast<El>(i);
  }
  return El::X;
}

struct Element {
  El elem;

  /*implicit*/ Element(El e) noexcept : elem(e) {}
  explicit Element(const char* str) noexcept : elem(find_element(str)) {}
  explicit Element(const std::string& s) noexcept : Element(s.c_str()) {}
  explicit Element(int number) noexcept
    : elem(static_cast<El>(number > 0 && number <= 118 ? number : 0)) {}
  /*implicit*/ operator El() const { return elem; }
  bool operator==(El e) const { return elem == e; }
  bool operator!=(El e) const { return elem != e; }

  int atomic_number() const {
    return elem == El::D ? 1 : static_cast<int>(elem);
  }
  double weight() const { return molecular_weight(elem); }
  float covalent_r() const { return covalent_radius(elem); }
  // return name such as Mg (not MG)
  const char* name() const { return element_name(elem); }
  // return uppercase name such as MG
  const char* uname() const { return element_uppercase_name(elem); }
};

} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
