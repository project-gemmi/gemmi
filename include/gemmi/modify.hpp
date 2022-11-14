// Copyright 2017-2021 Global Phasing Ltd.
//
// Modify various properties of the model.

// For modifications that depend on entities or connectivity see polyheur.hpp.

#ifndef GEMMI_MODIFY_HPP_
#define GEMMI_MODIFY_HPP_

#include "model.hpp"
#include "util.hpp"      // for vector_remove_if
#include <set>

namespace gemmi {

/// Remove alternative conformations.
template<class T> void remove_alternative_conformations(T& obj) {
  for (auto& child : obj.children())
    remove_alternative_conformations(child);
}
template<> inline void remove_alternative_conformations(Chain& chain) {
  std::set<SeqId> seqids;
  for (size_t i = 0; i < chain.residues.size(); ) {
    if (seqids.insert(chain.residues[i].seqid).second)
      ++i;
    else
      chain.residues.erase(chain.residues.begin() + i);
  }
  for (Residue& residue : chain.residues) {
    std::set<std::string> names;
    for (size_t i = 0; i < residue.atoms.size(); ) {
      Atom& atom = residue.atoms[i];
      atom.altloc = '\0';
      if (names.insert(atom.name).second)
        ++i;
      else
        residue.atoms.erase(residue.atoms.begin() + i);
    }
  }
}

/// Remove hydrogens.
template<class T> void remove_hydrogens(T& obj) {
  for (auto& child : obj.children())
    remove_hydrogens(child);
}
template<> inline void remove_hydrogens(Residue& res) {
  vector_remove_if(res.atoms, [](const Atom& a) {
    return a.element == El::H || a.element == El::D;
  });
}

/// Set isotropic ADP to the range (b_min, b_max). Values smaller than
/// b_min are changed to b_min, values larger than b_max to b_max.
/// Anisotropic ADP is left unchanged.
template<class T> void assign_b_iso(T& obj, float b_min, float b_max) {
  for (auto& child : obj.children())
    assign_b_iso(child, b_min, b_max);
}
template<> inline void assign_b_iso(Atom& atom, float b_min, float b_max) {
  atom.b_iso = std::min(std::max(atom.b_iso, b_min), b_max);
}

/// Remove anisotropic ADP
template<class T> void remove_anisou(T& obj) {
  for (auto& child : obj.children())
    remove_anisou(child);
}
template<> inline void remove_anisou(Atom& atom) {
  atom.aniso = {0, 0, 0, 0, 0, 0};
}

/// Set absent ANISOU to value from B_iso
template<class T> void ensure_anisou(T& obj) {
  for (auto& child : obj.children())
    ensure_anisou(child);
}
template<> inline void ensure_anisou(Atom& atom) {
  if (!atom.aniso.nonzero()) {
    float u = float(1. / gemmi::u_to_b() * atom.b_iso);
    atom.aniso = {u, u, u, 0.f, 0.f, 0.f};
  }
}

/// apply Transform to both atom's position and ADP
template<class T> void transform_pos_and_adp(T& obj, const Transform& tr) {
  for (auto& child : obj.children())
    transform_pos_and_adp(child, tr);
}
template<> inline void transform_pos_and_adp(Atom& atom, const Transform& tr) {
  atom.pos = Position(tr.apply(atom.pos));
  if (atom.aniso.nonzero())
    atom.aniso = atom.aniso.transformed_by<float>(tr.mat);
}

/// set atom site serial numbers to 1, 2, ...
inline void assign_serial_numbers(Model& model) {
  int serial = 0;
  for (CRA cra : model.all())
    cra.atom->serial = ++serial;
}
inline void assign_serial_numbers(Structure& st) {
  for (Model& model : st.models)
    assign_serial_numbers(model);
}

/// Hydrogens modelled as H/D mixture (altlocs H and D with the same position
/// and ADP, but with refined fraction of D), it can be stored in a file either
/// as two atoms (H and D) or, using CCP4/Refmac extension, as H atoms with
/// the ccp4_deuterium_fraction parameter.
/// This function switches fraction -> altlocs
inline void expand_hd_mixture(Structure& st) {
  if (!st.has_d_fraction)
    return;
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        for (size_t i = res.atoms.size(); i-- != 0; ) {
          Atom& atom = res.atoms[i];
          float d_fraction = atom.fraction;
          if (atom.element == El::H && d_fraction > 0) {
            if (d_fraction >= 1) {
              atom.element = El::D;
              if (atom.name[0] == 'H')
                atom.name[0] = 'D';
            } else {
              int alt_offset = atom.altloc;
              if (alt_offset) {
                alt_offset -= 'A';
                // we don't expect 4+ altlocs - ignore such cases
                if (alt_offset < 0 || alt_offset >= 3)
                  continue;
              }
              atom.altloc = 'A' + alt_offset;
              float d_occ = atom.occ * d_fraction;
              atom.occ *= (1 - d_fraction);
              auto deut = res.atoms.insert(res.atoms.begin() + i + 1, atom);
              deut->altloc = 'D' + alt_offset;
              deut->element = El::D;
              deut->occ = d_occ;
              if (deut->name[0] == 'H')
                deut->name[0] = 'D';
            }
          }
        }
  st.has_d_fraction = false;
}

inline bool replace_deuterium_with_fraction(Residue& res) {
  bool found = false;
  for (auto d = res.atoms.end(); d-- != res.atoms.begin(); )
    if (d->element == El::D) {
      found = true;
      auto h = res.atoms.begin();
      for (; h != res.atoms.end(); ++h)
        if (h->element == El::H && h->pos.approx(d->pos, 1e-9))
          break;
      if (h != res.atoms.end()) {
        h->occ += d->occ;
        h->fraction = h->occ > 0.f ? d->occ / h->occ : 0.f;
        if (h->altloc) {
          bool keep_altloc = false;
          for (auto i = res.atoms.begin(); i != res.atoms.end(); ++i)
            if (i != d && i != h && (i->name == h->name || i->name == d->name))
              keep_altloc = true;
          if (!keep_altloc)
            h->altloc = '\0';
        }
        res.atoms.erase(d);
      } else {
        d->element = El::H;
        d->fraction = 1;
        // Atom name is left unchanged. prepare_topology() first calls this
        // function and then conditionally changes the name (Dxx -> Hxx).
      }
    }
  return found;
}

/// Switch H/D altlocs at the same position to H w/ ccp4_deuterium_fraction.
inline void collapse_hd_mixture(Structure& st) {
  if (st.has_d_fraction)
    return;
  for (Model& model : st.models)
    for (Chain& chain : model.chains)
      for (Residue& res : chain.residues)
        if (replace_deuterium_with_fraction(res))
          st.has_d_fraction = true;
}

} // namespace gemmi
#endif
