//! @file
//! @brief Chemical component (monomer) definitions with restraints.
//!
//! ChemComp - chemical component that represents a monomer from Refmac
//! monomer library, or from PDB CCD.

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

//! @brief Type of chemical bond.
//!
//! Used in restraints and chemical component definitions.
enum class BondType {
  Unspec, Single, Double, Triple, Aromatic, Deloc, Metal
};

//! @brief Chirality type for tetrahedral centers.
enum class ChiralityType { Positive, Negative, Both };

//! @brief Geometric restraints for a chemical component or link.
//!
//! Contains bond distances, angles, torsions, chirality, and planarity restraints
//! used in refinement. Can represent monomer restraints from Refmac library or
//! PDB Chemical Component Dictionary (CCD).
struct Restraints {
  //! @brief Identifier for an atom in restraints.
  //!
  //! The comp field identifies the residue (1 for current, 2 for next in a link).
  //! The atom field is the atom name.
  struct AtomId {
    int comp;  //!< Component number (1=current residue, 2=next residue in link)
    std::string atom;  //!< Atom name

    bool operator==(const AtomId& o) const {
      return comp == o.comp && atom == o.atom;
    }
    bool operator!=(const AtomId& o) const { return !operator==(o); }

    bool operator==(const std::string& name) const { return atom == name; }
    bool operator!=(const std::string& name) const { return atom != name; }

    bool operator<(const AtomId& o) const {
      return comp == o.comp ? atom < o.atom : comp < o.comp;
    }

    //! @brief Get atom pointer from residue(s).
    //! @param res1 First residue
    //! @param res2 Second residue (for links), or nullptr
    //! @param alt Alternate location indicator
    //! @param altloc2 Second altloc (rare case for links between different altlocs, e.g., 2e7z)
    //! @return Pointer to atom, or nullptr if not found
    //!
    //! altloc2 is needed only in rare case when we have a link between
    //! atoms with different altloc (example: 2e7z).
    Atom* get_from(Residue& res1, Residue* res2, char alt, char altloc2) const {
      Residue* residue = &res1;
      if (comp == 2 && res2 != nullptr) {
        residue = res2;
        if (altloc2 != '\0')
          alt = altloc2;
      }
      Atom* a = residue->find_atom(atom, alt, El::X, false);
      // Special case: microheterogeneity may have shared atoms only in
      // the first residue. Example: in 1ejg N is shared between PRO and SER.
      if (a == nullptr && alt != '\0' && residue->group_idx > 0)
        a = (residue - residue->group_idx)->find_atom(atom, alt, El::X, false);
      return a;
    }
    const Atom* get_from(const Residue& res1, const Residue* res2,
                         char alt, char alt2) const {
      return get_from(const_cast<Residue&>(res1), const_cast<Residue*>(res2), alt, alt2);
    }
  };

  static std::string lexicographic_str(const std::string& name1,
                                       const std::string& name2) {
    return name1 < name2 ? cat(name1, '-', name2) : cat(name2, '-', name1);
  }

  //! @brief What type of distance (electron cloud or nucleus).
  enum class DistanceOf { ElectronCloud, Nucleus };

  //! @brief Bond distance restraint.
  //!
  //! Stores target bond distance and its estimated standard deviation (esd).
  //! Can store both electron cloud and nucleus distances.
  struct Bond {
    static const char* what() { return "bond"; }
    AtomId id1, id2;  //!< Bonded atoms
    BondType type;  //!< Bond type (single, double, etc.)
    bool aromatic;  //!< Is aromatic bond
    double value;  //!< Target distance (electron cloud), in Angstroms
    double esd;  //!< Estimated standard deviation for electron cloud distance
    double value_nucleus;  //!< Target distance (nucleus), in Angstroms
    double esd_nucleus;  //!< Estimated standard deviation for nucleus distance
    std::string str() const { return cat(id1.atom, '-', id2.atom); }
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

  //! @brief Angle restraint.
  //!
  //! Restrains the angle formed by three atoms.
  struct Angle {
    static const char* what() { return "angle"; }
    AtomId id1, id2, id3;  //!< Three atoms forming the angle (id2 is center)
    double value;  //!< Target angle value in degrees
    double esd;  //!< Estimated standard deviation in degrees
    double radians() const { return rad(value); }
    std::string str() const {
      return cat(id1.atom, '-', id2.atom, '-', id3.atom);
    }
  };

  //! @brief Torsion angle restraint.
  //!
  //! Restrains the dihedral angle formed by four atoms.
  struct Torsion {
    static const char* what() { return "torsion"; }
    std::string label;  //!< Torsion label/identifier
    AtomId id1, id2, id3, id4;  //!< Four atoms defining the torsion
    double value;  //!< Target torsion angle in degrees
    double esd;  //!< Estimated standard deviation in degrees
    int period;  //!< Periodicity of the torsion
    std::string str() const {
      return cat(id1.atom, '-', id2.atom, '-', id3.atom, '-', id4.atom);
    }
  };

  //! @brief Chirality restraint.
  //!
  //! Restrains the handedness of a tetrahedral center.
  struct Chirality {
    static const char* what() { return "chirality"; }
    AtomId id_ctr, id1, id2, id3;  //!< Central atom and three neighbors
    ChiralityType sign;  //!< Expected chirality (positive/negative/both)

    bool is_wrong(double volume) const {
      return (sign == ChiralityType::Positive && volume < 0) ||
             (sign == ChiralityType::Negative && volume > 0);
    }
    std::string str() const {
      return cat(id_ctr.atom, ',', id1.atom, ',', id2.atom, ',', id3.atom);
    }
  };

  //! @brief Planarity restraint.
  //!
  //! Restrains a set of atoms to be coplanar.
  struct Plane {
    static const char* what() { return "plane"; }
    std::string label;  //!< Plane label/identifier
    std::vector<AtomId> ids;  //!< Atoms that should be coplanar
    double esd;  //!< Estimated standard deviation for planarity
    std::string str() const {
      return join_str(ids, ',', [](const AtomId& a) { return a.atom; });
    }
  };

  std::vector<Bond> bonds;  //!< Bond distance restraints
  std::vector<Angle> angles;  //!< Angle restraints
  std::vector<Torsion> torsions;  //!< Torsion angle restraints
  std::vector<Chirality> chirs;  //!< Chirality restraints
  std::vector<Plane> planes;  //!< Planarity restraints

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
                                         std::vector<AtomId> visited,
                                         int min_length=1) const {
    int start = (int) visited.size();
    int end = -1;
    visited.push_back(b);
    std::vector<int> parent(visited.size(), -1);
    std::vector<int> depth(visited.size(), -1);
    depth[start] = 0;
    for (int n = start; end == -1 && n != (int) visited.size(); ++n) {
      for (const Bond& bond : bonds)
        if (const AtomId* id = bond.other(visited[n])) {
          if (*id == a) {
            int cand_len = depth[n] + 1;
            if (cand_len >= min_length) {
              end = (int) visited.size();
              if (!in_vector(*id, visited)) {
                visited.push_back(*id);
                parent.push_back(n);
                depth.push_back(cand_len);
              }
            }
          } else if (!in_vector(*id, visited)) {
            visited.push_back(*id);
            parent.push_back(n);
            depth.push_back(depth[n] + 1);
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
  return mult * std::sqrt(std::max(0., x + y));
}

inline double Restraints::chiral_abs_volume(const Restraints::Chirality& ch) const {
  return gemmi::chiral_abs_volume(get_bond(ch.id_ctr, ch.id1).value,
                                  get_bond(ch.id_ctr, ch.id2).value,
                                  get_bond(ch.id_ctr, ch.id3).value,
                                  get_angle(ch.id1, ch.id_ctr, ch.id2).value,
                                  get_angle(ch.id2, ch.id_ctr, ch.id3).value,
                                  get_angle(ch.id3, ch.id_ctr, ch.id1).value);
}

//! @brief Chemical component (monomer) definition.
//!
//! Represents a monomer from Refmac monomer library or PDB Chemical Component
//! Dictionary (CCD). Contains atom list, restraints, and aliasing information.
struct ChemComp {
  //! @brief Chemical component group classification.
  //!
  //! Items used in _chem_comp.group and _chem_link.group_comp_N in CCP4.
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

  //! @brief Atom in a chemical component definition.
  struct Atom {
    std::string id;  //!< Atom identifier/name
    std::string old_id;  //!< Alternative atom ID (from _chem_comp_atom.alt_atom_id)
    Element el = El::X;  //!< Element type
    //! Partial charge. _chem_comp_atom.partial_charge can be non-integer,
    //! _chem_comp_atom.charge is always integer (but sometimes has format
    //! '0.000' which is not correct but we ignore it).
    float charge = 0;
    std::string chem_type;  //!< Chemical type descriptor
    std::string acedrg_type;  // read from _chem_comp_atom.atom_type
    Position xyz;  //!< 3D coordinates (if available)

    bool is_hydrogen() const { return gemmi::is_hydrogen(el); }
  };

  //! @brief Atom name aliasing for different chemical groups.
  //!
  //! Maps atom names in this component to standard names in specific groups
  //! (e.g., peptide, DNA, RNA).
  struct Aliasing {
    Group group;  //!< Chemical group this aliasing applies to
    //! Pairs of (name in chem_comp, usual name in this group).
    //! For example, maps local atom names to standard peptide atom names.
    std::vector<std::pair<std::string, std::string>> related;

    const std::string* name_from_alias(const std::string& atom_id) const {
      for (const auto& item : related)
        if (item.second == atom_id)
          return &item.first;
      return nullptr;
    }
  };

  std::string name;  //!< Component name/ID (e.g., "ALA", "GLY")
  std::string type_or_group;  //!< Value from _chem_comp.type or _chem_comp.group
  Group group = Group::Null;  //!< Component group classification
  bool has_coordinates = false;  //!< Whether 3D coordinates are available
  std::vector<Atom> atoms;  //!< Atoms in this component
  std::vector<Aliasing> aliases;  //!< Atom name aliases for different groups
  Restraints rt;  //!< Geometric restraints for this component

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
    return find_atom(atom_id) != atoms.end();
  }

  std::vector<Atom>::iterator find_atom_by_old_name(const std::string& old_id) {
    return std::find_if(atoms.begin(), atoms.end(),
                        [&](const Atom& a) { return a.old_id == old_id; });
  }
  std::vector<Atom>::const_iterator find_atom_by_old_name(const std::string& old_id) const {
    return const_cast<ChemComp*>(this)->find_atom_by_old_name(old_id);
  }
  bool has_old_names() const {
    return std::any_of(atoms.begin(), atoms.end(),
                       [&](const Atom& a) { return !a.old_id.empty() && a.old_id != a.id; });
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
  if (s.size() >= 3)
    switch (ialpha4_id(s.c_str())) {
      case ialpha4_id("sing"): return BondType::Single;
      case ialpha4_id("doub"): return BondType::Double;
      case ialpha4_id("trip"): return BondType::Triple;
      case ialpha4_id("arom"): return BondType::Aromatic;
      case ialpha4_id("meta"): return BondType::Metal;
      case ialpha4_id("delo"): return BondType::Deloc;
      case ialpha4_id("1.5"):  return BondType::Deloc; // rarely used
      // program PDB2TNT produces a restraint file with bond type 'coval'
      case ialpha4_id("cova"):  return BondType::Unspec;
    }
  if (cif::is_null(s))
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
                             "?charge", "?partial_charge", "?alt_atom_id",
                             "?atom_type"})) {
    ChemComp::Atom atom;
    atom.id = row.str(0);
    atom.old_id = row.has(5) ? row.str(5) : "";
    atom.el = Element(row.str(1));
    atom.charge = (float) cif::as_number(row.one_of(3, 4), 0.0);
    atom.chem_type = row.has(2) ? row.str(2) : "";
    atom.acedrg_type = row.has(6) ? row.str(6) : "";
    cc.atoms.push_back(std::move(atom));
  }
  // Also check _chem_comp_acedrg table for atom types (used by acedrg output)
  for (auto row : block.find("_chem_comp_acedrg.", {"atom_id", "atom_type"})) {
    std::string atom_id = row.str(0);
    std::string atom_type = row.str(1);
    for (auto& atom : cc.atoms)
      if (atom.id == atom_id && atom.acedrg_type.empty()) {
        atom.acedrg_type = atom_type;
        break;
      }
  }
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
