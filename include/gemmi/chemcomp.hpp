// Copyright 2018 Global Phasing Ltd.
//
// ChemComp - chemical component that represents a monomer from Refmac
// monomer library, or from PDB CCD.

#ifndef GEMMI_CHEMCOMP_HPP_
#define GEMMI_CHEMCOMP_HPP_

#include <string>
#include <vector>
#include "cifdoc.hpp"
#include "elem.hpp"  // for Element
#include "fail.hpp"  // for fail, unreachable
#include "numb.hpp"  // for as_number
#include "util.hpp"  // for istarts_with, join_str
#include "model.hpp" // for Residue, Atom

namespace gemmi {

enum class BondType {
  Unspec, Single, Double, Triple, Aromatic, Deloc, Metal
};
enum class ChiralityType { Positive, Negative, Both };

struct Restraints {
  struct AtomId {
    int comp;
    std::string atom;

    bool operator==(const AtomId& o) const {
      return comp == o.comp && atom == o.atom;
    }
    bool operator!=(const AtomId& o) const { return !operator==(o); }

    bool operator==(const std::string& name) const { return atom == name; }
    bool operator!=(const std::string& name) const { return atom != name; }

    bool operator<(const AtomId& o) const {
      return comp == o.comp ? atom < o.atom : comp < o.comp;
    }

    Atom* get_from(Residue& res1, Residue* res2, char altloc) const {
      Residue* residue = (comp == 1 || res2 == nullptr ? &res1 : res2);
      Atom* a = residue->find_atom(atom, altloc);
      // Special case: microheterogeneity may have shared atoms only in
      // the first residue. Example: in 1ejg N is shared between PRO and SER.
      if (a == nullptr && altloc != '\0' && residue->group_idx > 0)
        a = (residue - residue->group_idx)->find_atom(atom, altloc);
      return a;
    }
    const Atom* get_from(const Residue& res1, const Residue* res2, char altloc) const {
      return get_from(const_cast<Residue&>(res1), const_cast<Residue*>(res2), altloc);
    }
  };

  static std::string lexicographic_str(const std::string& name1,
                                       const std::string& name2) {
    return name1 < name2 ? name1 + "-" + name2 : name2 + "-" + name1;
  }

  enum class DistanceOf { ElectronCloud, Nucleus };

  struct Bond {
    AtomId id1, id2;
    BondType type;
    bool aromatic;
    double value;
    double esd;
    double value_nucleus;
    double esd_nucleus;
    std::string str() const { return id1.atom + "-" + id2.atom; }
    std::string lexicographic_str() const {
      return Restraints::lexicographic_str(id1.atom, id2.atom);
    }
    double distance(DistanceOf of) const {
      return of == DistanceOf::ElectronCloud ? value : value_nucleus;
    }
    template<typename T> const AtomId* other(const T& a) const {
      if (id1 == a) return &id2;
      if (id2 == a) return &id1;
      return nullptr;
    }
  };

  struct Angle {
    AtomId id1, id2, id3;
    double value;  // degrees
    double esd;
    double radians() const { return rad(value); }
    std::string str() const {
      return id1.atom + "-" + id2.atom + "-" + id3.atom;
    }
  };

  struct Torsion {
    std::string label;
    AtomId id1, id2, id3, id4;
    double value;
    double esd;
    int period;
    std::string str() const {
      return id1.atom + "-" + id2.atom + "-" + id3.atom + "-" + id4.atom;
    }
  };

  struct Chirality {
    AtomId id_ctr, id1, id2, id3;
    ChiralityType sign;

    bool is_wrong(double volume) const {
      return (sign == ChiralityType::Positive && volume < 0) ||
             (sign == ChiralityType::Negative && volume > 0);
    }
    std::string str() const {
      return id_ctr.atom + "," + id1.atom + "," + id2.atom + "," + id3.atom;
    }
  };

  struct Plane {
    std::string label;
    std::vector<AtomId> ids;
    double esd;
    std::string str() const {
      return join_str(ids, ',', [](const AtomId& a) { return a.atom; });
    }
  };

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Torsion> torsions;
  std::vector<Chirality> chirs;
  std::vector<Plane> planes;

  bool empty() const {
    return bonds.empty() && angles.empty() && torsions.empty() &&
           chirs.empty() && planes.empty();
  }

  template<typename T>
  std::vector<Bond>::iterator find_bond(const T& a1, const T& a2) {
    return std::find_if(bonds.begin(), bonds.end(), [&](const Bond& b) {
        return (b.id1 == a1 && b.id2 == a2) || (b.id1 == a2 && b.id2 == a1);
    });
  }
  template<typename T>
  std::vector<Bond>::const_iterator find_bond(const T& a1, const T& a2) const {
    return const_cast<Restraints*>(this)->find_bond(a1, a2);
  }
  const Bond& get_bond(const AtomId& a1, const AtomId& a2) const {
    auto it = find_bond(a1, a2);
    if (it == bonds.end())
      fail("Bond restraint not found: ", a1.atom, '-', a2.atom);
    return *it;
  }

  template<typename T>
  bool are_bonded(const T& a1, const T& a2) const {
    return find_bond(a1, a2) != bonds.end();
  }

  template<typename T>
  const AtomId* first_bonded_atom(const T& a) const {
    for (const Bond& bond : bonds)
      if (const AtomId* other = bond.other(a))
        return other;
    return nullptr;
  }

  // BFS
  std::vector<AtomId> find_shortest_path(const AtomId& a, const AtomId& b,
                                         std::vector<AtomId> visited) const {
    int start = (int) visited.size();
    int end = -1;
    visited.push_back(b);
    std::vector<int> parent(visited.size(), -1);
    for (int n = start; end == -1 && n != (int) visited.size(); ++n) {
      for (const Bond& bond : bonds)
        if (const AtomId* id = bond.other(visited[n])) {
          if (*id == a)
            end = (int) visited.size();
          if (!in_vector(*id, visited)) {
            visited.push_back(*id);
            parent.push_back(n);
          }
        }
    }
    std::vector<AtomId> path;
    for (int n = end; n != -1; n = parent[n])
      path.push_back(visited[n]);
    return path;
  }

  template<typename T>
  std::vector<Angle>::iterator find_angle(const T& a, const T& b, const T& c) {
    return std::find_if(angles.begin(), angles.end(), [&](const Angle& ang) {
        return ang.id2 == b && ((ang.id1 == a && ang.id3 == c) ||
                                (ang.id1 == c && ang.id3 == a));
    });
  }
  template<typename T> std::vector<Angle>::const_iterator
  find_angle(const T& a, const T& b, const T& c) const {
    return const_cast<Restraints*>(this)->find_angle(a, b, c);
  }
  const Angle& get_angle(const AtomId& a, const AtomId& b, const AtomId& c) const {
    auto it = const_cast<Restraints*>(this)->find_angle(a, b, c);
    if (it == angles.end())
      fail("Angle restraint not found: ", a.atom, '-', b.atom, '-', c.atom);
    return *it;
  }

  template<typename T>
  std::vector<Torsion>::iterator find_torsion(const T& a, const T& b,
                                              const T& c, const T& d) {
    return std::find_if(torsions.begin(), torsions.end(),
                        [&](const Torsion& t) {
        return (t.id1 == a && t.id2 == b && t.id3 == c && t.id4 == d) ||
               (t.id1 == d && t.id2 == c && t.id3 == b && t.id4 == a);
    });
  }
  template<typename T> std::vector<Torsion>::const_iterator
  find_torsion(const T& a, const T& b, const T& c, const T& d) const {
    return const_cast<Restraints*>(this)->find_torsion(a, b, c, d);
  }

  template<typename T>
  std::vector<Chirality>::iterator find_chir(const T& ctr, const T& a,
                                             const T& b, const T& c) {
    return std::find_if(chirs.begin(), chirs.end(), [&](const Chirality& t) {
        return t.id_ctr == ctr && ((t.id1 == a && t.id2 == b && t.id3 == c) ||
                                   (t.id1 == b && t.id2 == c && t.id3 == a) ||
                                   (t.id1 == c && t.id2 == a && t.id3 == b));
    });
  }
  template<typename T> std::vector<Chirality>::const_iterator
  find_chir(const T& ctr, const T& a, const T& b, const T& c) const {
    return const_cast<Restraints*>(this)->find_chir(ctr, a, b, c);
  }

  double chiral_abs_volume(const Restraints::Chirality& ch) const;

  std::vector<Plane>::iterator get_plane(const std::string& label) {
    return std::find_if(planes.begin(), planes.end(),
                        [&label](const Plane& p) { return p.label == label; });
  }

  Plane& get_or_add_plane(const std::string& label) {
    std::vector<Plane>::iterator it = get_plane(label);
    if (it != planes.end())
      return *it;
    planes.push_back(Plane{label, {}, 0.0});
    return planes.back();
  }

  void rename_atom(const AtomId& atom_id, const std::string& new_name) {
    auto rename_atom = [&](AtomId& id) {
      if (id == atom_id)
        id.atom = new_name;
    };
    for (Bond& bond : bonds) {
      rename_atom(bond.id1);
      rename_atom(bond.id2);
    }
    for (Angle& angle : angles) {
      rename_atom(angle.id1);
      rename_atom(angle.id2);
      rename_atom(angle.id3);
    }
    for (Torsion& tor : torsions) {
      rename_atom(tor.id1);
      rename_atom(tor.id2);
      rename_atom(tor.id3);
      rename_atom(tor.id4);
    }
    for (Chirality& chir : chirs) {
      rename_atom(chir.id_ctr);
      rename_atom(chir.id1);
      rename_atom(chir.id2);
      rename_atom(chir.id3);
    }
    for (Plane& plane : planes)
      for (AtomId& id : plane.ids)
        rename_atom(id);
  }
};

template<typename Restr>
double angle_z(double value_rad, const Restr& restr, double full=360.) {
  return angle_abs_diff(deg(value_rad), restr.value, full) / restr.esd;
}

inline double chiral_abs_volume(double bond1, double bond2, double bond3,
                                double angle1, double angle2, double angle3) {
  double mult = bond1 * bond2 * bond3;
  double x = 1;
  double y = 2;
  for (double a : {angle1, angle2, angle3}) {
    double cosine = a == 90. ? 0. : std::cos(rad(a));
    x -= cosine * cosine;
    y *= cosine;
  }
  return mult * std::sqrt(x + y);
}

inline double Restraints::chiral_abs_volume(const Restraints::Chirality& ch) const {
  return gemmi::chiral_abs_volume(get_bond(ch.id_ctr, ch.id1).value,
                                  get_bond(ch.id_ctr, ch.id2).value,
                                  get_bond(ch.id_ctr, ch.id3).value,
                                  get_angle(ch.id1, ch.id_ctr, ch.id2).value,
                                  get_angle(ch.id2, ch.id_ctr, ch.id3).value,
                                  get_angle(ch.id3, ch.id_ctr, ch.id1).value);
}

struct ChemComp {
  // Items used in _chem_comp.group and _chem_link.group_comp_N in CCP4.
  enum class Group {
    Peptide,      // "peptide"
    PPeptide,     // "P-peptide"
    MPeptide,     // "M-peptide"
    Dna,          // "DNA" - used in _chem_comp.group
    Rna,          // "RNA" - used in _chem_comp.group
    DnaRna,       // "DNA/RNA" - used in _chem_link.group_comp_N
    Pyranose,     // "pyranose"
    Ketopyranose, // "ketopyranose"
    Furanose,     // "furanose"
    NonPolymer,   // "non-polymer"
    Null
  };

  struct Atom {
    std::string id;
    Element el;
    // _chem_comp_atom.partial_charge can be non-integer,
    // _chem_comp_atom.charge is always integer (but sometimes has format
    //  '0.000' which is not correct but we ignore it).
    float charge;
    std::string chem_type;

    bool is_hydrogen() const { return gemmi::is_hydrogen(el); }
  };

  struct Aliasing {
    Group group;
    // pairs of (name in chem_comp, usual name in this group)
    std::vector<std::pair<std::string, std::string>> related;

    const std::string* name_from_alias(const std::string& atom_id) const {
      for (const auto& item : related)
        if (item.second == atom_id)
          return &item.first;
      return nullptr;
    }
  };

  std::string name;
  std::string type_or_group;  // _chem_comp.type or _chem_comp.group
  Group group = Group::Null;
  std::vector<Atom> atoms;
  std::vector<Aliasing> aliases;
  Restraints rt;

  const Aliasing& get_aliasing(Group g) const {
    for (const Aliasing& aliasing : aliases)
      if (aliasing.group == g)
        return aliasing;
    fail("aliasing not found");
  }

  static Group read_group(const std::string& str) {
    if (str.size() >= 3) {
      const char* cstr = str.c_str();
      if ((str[0] == '\'' || str[0] == '"') && str.size() >= 5)
        ++cstr;
      switch (ialpha4_id(cstr)) {
        case ialpha4_id("non-"): return Group::NonPolymer;
        case ialpha4_id("pept"): return Group::Peptide;
        case ialpha4_id("l-pe"): return Group::Peptide;
        case ialpha4_id("p-pe"): return Group::PPeptide;
        case ialpha4_id("m-pe"): return Group::MPeptide;
        case ialpha4_id("dna"):  return Group::Dna;
        case ialpha4_id("rna"):  return Group::Rna;
        case ialpha4_id("dna/"): return Group::DnaRna;
        case ialpha4_id("pyra"): return Group::Pyranose;
        case ialpha4_id("keto"): return Group::Ketopyranose;
        case ialpha4_id("fura"): return Group::Furanose;
      }
    }
    return Group::Null;
  }

  static const char* group_str(Group g) {
    switch (g) {
      case Group::Peptide: return "peptide";
      case Group::PPeptide: return "P-peptide";
      case Group::MPeptide: return "M-peptide";
      case Group::Dna: return "DNA";
      case Group::Rna: return "RNA";
      case Group::DnaRna: return "DNA/RNA";
      case Group::Pyranose: return "pyranose";
      case Group::Ketopyranose: return "ketopyranose";
      case Group::Furanose: return "furanose";
      case Group::NonPolymer: return "non-polymer";
      case Group::Null: return ".";
    }
    unreachable();
  }

  void set_group(const std::string& s) {
    type_or_group = s;
    group = read_group(s);
  }

  std::vector<Atom>::iterator find_atom(const std::string& atom_id) {
    return std::find_if(atoms.begin(), atoms.end(),
                        [&](const Atom& a) { return a.id == atom_id; });
  }
  std::vector<Atom>::const_iterator find_atom(const std::string& atom_id) const {
    return const_cast<ChemComp*>(this)->find_atom(atom_id);
  }
  bool has_atom(const std::string& atom_id) const {
    return find_atom(atom_id) == atoms.end();
  }

  int get_atom_index(const std::string& atom_id) const {
    auto it = find_atom(atom_id);
    if (it == atoms.end())
      fail("Chemical component ", name, " has no atom ", atom_id);
    return int(it - atoms.begin());
  }

  const Atom& get_atom(const std::string& atom_id) const {
    return atoms[get_atom_index(atom_id)];
  }

  /// Check if the group (M-|P-)peptide
  static bool is_peptide_group(Group g) {
    return g == Group::Peptide || g == Group::PPeptide || g == Group::MPeptide;
  }

  /// Check if the group is DNA/RNA
  static bool is_nucleotide_group(Group g) {
    return g == Group::Dna || g == Group::Rna || g == Group::DnaRna;
  }

  void remove_nonmatching_restraints() {
    vector_remove_if(rt.bonds, [&](const Restraints::Bond& x) {
      return !has_atom(x.id1.atom) ||
             !has_atom(x.id2.atom);
    });
    vector_remove_if(rt.angles, [&](const Restraints::Angle& x) {
      return !has_atom(x.id1.atom) ||
             !has_atom(x.id2.atom) ||
             !has_atom(x.id3.atom);
    });
    vector_remove_if(rt.torsions, [&](const Restraints::Torsion& x) {
      return !has_atom(x.id1.atom) ||
             !has_atom(x.id2.atom) ||
             !has_atom(x.id3.atom) ||
             !has_atom(x.id4.atom);
    });
    vector_remove_if(rt.chirs, [&](const Restraints::Chirality& x) {
      return !has_atom(x.id_ctr.atom) ||
             !has_atom(x.id1.atom) ||
             !has_atom(x.id2.atom) ||
             !has_atom(x.id3.atom);
    });
    for (Restraints::Plane& plane : rt.planes)
      vector_remove_if(plane.ids, [&](const Restraints::AtomId& x) {
        return !has_atom(x.atom);
      });
  }

  ChemComp& remove_hydrogens() {
    vector_remove_if(atoms, [](const ChemComp::Atom& a) {
      return a.is_hydrogen();
    });
    remove_nonmatching_restraints();
    return *this;
  }
};

inline BondType bond_type_from_string(const std::string& s) {
  if (istarts_with(s, "sing"))
    return BondType::Single;
  if (istarts_with(s, "doub"))
    return BondType::Double;
  if (istarts_with(s, "trip"))
    return BondType::Triple;
  if (istarts_with(s, "arom"))
    return BondType::Aromatic;
  if (istarts_with(s, "metal"))
    return BondType::Metal;
  if (istarts_with(s, "delo") || s == "1.5")
    return BondType::Deloc;
  if (cif::is_null(s))
    return BondType::Unspec;
  // program PDB2TNT produces a restraint file with bond type 'coval'
  if (s == "coval")
    return BondType::Unspec;
  throw std::out_of_range("Unexpected bond type: " + s);
}

inline const char* bond_type_to_string(BondType btype) {
  switch (btype) {
    case BondType::Unspec: return ".";
    case BondType::Single: return "single";
    case BondType::Double: return "double";
    case BondType::Triple: return "triple";
    case BondType::Aromatic: return "aromatic";
    case BondType::Deloc: return "deloc";
    case BondType::Metal: return "metal";
  }
  unreachable();
}

inline float order_of_bond_type(BondType btype) {
  switch (btype) {
    case BondType::Single: return 1.0f;
    case BondType::Double: return 2.0f;
    case BondType::Triple: return 3.0f;
    case BondType::Aromatic: return 1.5f;
    case BondType::Deloc: return 1.5f;
    case BondType::Metal: return 1.0f;
    case BondType::Unspec: return 0.0f;
  }
  unreachable();
}

// it doesn't handle crossN types from the monomer library
inline ChiralityType chirality_from_string(const std::string& s) {
  switch (s[0] | 0x20) {
    case 'p': return ChiralityType::Positive;
    case 'n': return ChiralityType::Negative;
    case 'b': return ChiralityType::Both;
    case '.': return ChiralityType::Both;
    default: throw std::out_of_range("Unexpected chirality: " + s);
  }
}

inline ChiralityType chirality_from_flag_and_volume(const std::string& s,
                                                    double volume) {
  switch (s[0] | 0x20) {
    case 's': return volume > 0 ? ChiralityType::Positive
                                : ChiralityType::Negative;
    case 'n': return ChiralityType::Both;
    default: throw std::out_of_range("Unexpected volume_flag: " + s);
  }
}

inline const char* chirality_to_string(ChiralityType chir_type) {
  switch (chir_type) {
    case ChiralityType::Positive: return "positive";
    case ChiralityType::Negative: return "negative";
    case ChiralityType::Both: return "both";
  }
  unreachable();
}

inline ChemComp make_chemcomp_from_block(const cif::Block& block_) {
  ChemComp cc;
  cc.name = block_.name.substr(starts_with(block_.name, "comp_") ? 5 : 0);
  cif::Block& block = const_cast<cif::Block&>(block_);
  // CCD uses _chem_comp.type, monomer libraries use .group in a separate block
  // named comp_list, but it'd be a better idea to have _chem_comp.group in the
  // same block, so we try read it.
  if (cif::Column group_col = block.find_values("_chem_comp.group"))
    cc.set_group(group_col.str(0));
  else if (cif::Column type_col = block.find_values("_chem_comp.type"))
    cc.type_or_group = type_col.str(0);
  for (auto row : block.find("_chem_comp_atom.",
                             {"atom_id", "type_symbol", "?type_energy",
                             "?charge", "?partial_charge"}))
    cc.atoms.push_back({row.str(0), Element(row.str(1)),
                        (float) cif::as_number(row.one_of(3, 4), 0.0),
                        row.has(2) ? row.str(2) : ""});
  for (auto row : block.find("_chem_comp_bond.",
                             {"atom_id_1", "atom_id_2",              // 0, 1
                              "?type", "?value_order",               // 2, 3
                              "?aromatic", "?pdbx_aromatic_flag",    // 4, 5
                              "?value_dist", "?value_dist_esd",      // 6, 7
                              "?value_dist_nucleus", "?value_dist_nucleus_esd"})) { // 8, 9
    bool aromatic_flag = (row.one_of(4, 5)[0] | 0x20) == 'y';
    double dist = row.has(6) ? cif::as_number(row[6]) : NAN;
    double esd = row.has(7) ? cif::as_number(row[7]) : NAN;
    double dist_nucl = row.has(8) ? cif::as_number(row[8]) : NAN;
    double esd_nucl = row.has(9) ? cif::as_number(row[9]) : NAN;
    cc.rt.bonds.push_back({{1, row.str(0)}, {1, row.str(1)},
                          bond_type_from_string(row.one_of(2, 3)),
                          aromatic_flag, dist, esd, dist_nucl, esd_nucl});
  }
  for (auto row : block.find("_chem_comp_angle.",
                             {"atom_id_1", "atom_id_2", "atom_id_3",
                              "value_angle", "value_angle_esd"}))
    cc.rt.angles.push_back({{1, row.str(0)}, {1, row.str(1)}, {1, row.str(2)},
                            cif::as_number(row[3]), cif::as_number(row[4])});
  for (auto row : block.find("_chem_comp_tor.",
                             {"id",
                              "atom_id_1", "atom_id_2",
                              "atom_id_3", "atom_id_4",
                              "value_angle", "value_angle_esd",
                              "period"}))
    cc.rt.torsions.push_back({row.str(0),
                              {1, row.str(1)}, {1, row.str(2)},
                              {1, row.str(3)}, {1, row.str(4)},
                              cif::as_number(row[5]), cif::as_number(row[6]),
                              cif::as_int(row[7])});
  for (auto row : block.find("_chem_comp_chir.",
                             {"atom_id_centre",
                              "atom_id_1", "atom_id_2", "atom_id_3",
                              "volume_sign"}))
    if (row[4][0] != 'c') // ignore crossN
      cc.rt.chirs.push_back({{1, row.str(0)},
                             {1, row.str(1)}, {1, row.str(2)}, {1, row.str(3)},
                             chirality_from_string(row[4])});
  // mmCIF compliant
  cif::Table chir_tab = block.find("_chem_comp_chir_atom.",
                                   {"chir_id", "atom_id"});
  if (chir_tab.ok()) {
    std::map<std::string, std::vector<std::string>> chir_atoms;
    for (auto chir : chir_tab)
      chir_atoms[chir[0]].push_back(chir[1]);
    for (auto row : block.find("_chem_comp_chir.",
                               {"id", "atom_id", "volume_flag", "volume_three"})) {
        auto atoms = chir_atoms.find(row.str(0));
        if (atoms != chir_atoms.end() && atoms->second.size() == 3)
          cc.rt.chirs.push_back({{1, row.str(1)},
                                 {1, atoms->second[0]}, {1, atoms->second[1]},
                                 {1, atoms->second[2]},
                                 chirality_from_flag_and_volume(row[2],
                                                                cif::as_number(row[3]))});
    }
  }
  for (auto row : block.find("_chem_comp_plane_atom.",
                             {"plane_id", "atom_id" , "dist_esd"})) {
    Restraints::Plane& plane = cc.rt.get_or_add_plane(row.str(0));
    if (plane.esd == 0.0)
      plane.esd = cif::as_number(row[2]);
    plane.ids.push_back({1, row.str(1)});
  }
  for (auto row : block.find("_chem_comp_alias.",
                             {"group", "atom_id", "atom_id_standard"})) {
    ChemComp::Group group = ChemComp::read_group(row.str(0));
    if (cc.aliases.empty() || cc.aliases.back().group != group) {
      cc.aliases.emplace_back();
      cc.aliases.back().group = group;
    }
    cc.aliases.back().related.emplace_back(row.str(1), row.str(2));
  }

  return cc;
}

} // namespace gemmi
#endif
