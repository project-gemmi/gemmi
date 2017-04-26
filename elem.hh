// Copyright 2017 Global Phasing Ltd.
//
// Elements from the periodic table.

#ifndef GEMMI_ELEM_HH
#define GEMMI_ELEM_HH
namespace gemmi {
namespace mol {

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

inline double molecular_weight(El el) {
  constexpr double weights[] = {
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

inline const char* element_name(El el) {
  constexpr const char* names[] = {
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
    "D", nullptr
  };
  static_assert(static_cast<int>(El::Og) == 118, "Hmm");
  static_assert(names[118][0] == 'O', "Hmm");
  static_assert(sizeof(names) / sizeof(names[0]) ==
                static_cast<int>(El::END) + 1, "Hmm");
  return names[static_cast<int>(el)];
}


inline El find_element(const char* symbol) {
  if (symbol == nullptr || symbol[0] == '\0')
    return El::X;
  char first = symbol[0] & ~0x20;  // lower -> upper, space -> NUL
  if (symbol[1] == '\0')
		// all the most common elements in the PDB are single-letter
    switch (first) {
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

  if (symbol[2] != '\0')
    return El::X;

#define EL(s) uint16_t(#s[0] << 8 | #s[1])
  const uint16_t ptable[119] = {
    EL(H),  EL(HE), EL(LI), EL(BE), EL(B),  EL(C),  EL(N),  EL(O),  EL(F),
    EL(NE), EL(NA), EL(MG), EL(AL), EL(SI), EL(P),  EL(S),  EL(CL), EL(AR),
    EL(K),  EL(CA), EL(SC), EL(TI), EL(V),  EL(CR), EL(MN), EL(FE), EL(CO),
    EL(NI), EL(CU), EL(ZN), EL(GA), EL(GE), EL(AS), EL(SE), EL(BR), EL(KR),
    EL(RB), EL(SR), EL(Y),  EL(ZR), EL(NB), EL(MO), EL(TC), EL(RU), EL(RH),
    EL(PD), EL(AG), EL(CD), EL(IN), EL(SN), EL(SB), EL(TE), EL(I),  EL(XE),
    EL(CS), EL(BA), EL(LA), EL(CE), EL(PR), EL(ND), EL(PM), EL(SM), EL(EU),
    EL(GD), EL(TB), EL(DY), EL(HO), EL(ER), EL(TM), EL(YB), EL(LU),
    EL(HF), EL(TA), EL(W),  EL(RE), EL(OS), EL(IR), EL(PT), EL(AU), EL(HG),
    EL(TL), EL(PB), EL(BI), EL(PO), EL(AT), EL(RN),
    EL(FR), EL(RA), EL(AC), EL(TH), EL(PA), EL(U),  EL(NP), EL(PU), EL(AM),
    EL(CM), EL(BK), EL(CF), EL(ES), EL(FM), EL(MD), EL(NO), EL(LR),
    EL(RF), EL(DB), EL(SG), EL(BH), EL(HS), EL(MT), EL(D)
  };
#undef EL
  uint16_t sym16 = (first << 8) | (symbol[1] & ~0x20);
  unsigned char offset = std::find(ptable, ptable + 119, sym16) - ptable;
  return offset != 119 ? static_cast<El>(offset + 1) : El::X;
}

struct Element {
  El elem;

  /*implicit*/ Element(El e) noexcept : elem(e) {}
  explicit Element(const char* str) noexcept : elem(find_element(str)) {}
  explicit Element(const std::string& s) noexcept : Element(s.c_str()) {}
  explicit Element(int number) noexcept
    : elem(static_cast<El>(number > 0 && number <= 118 ? number : 0)) {}

  int atomic_number() const {
    return elem == El::D ? 1 : static_cast<int>(elem);
  }
  double weight() const { return molecular_weight(elem); }
  // return name such as Mg (not MG)
  const char* name() const { return element_name(elem); }
};

} // namespace mol
} // namespace gemmi
#endif
// vim:sw=2:ts=2:et
