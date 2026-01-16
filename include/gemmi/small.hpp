//! @file
//! @brief Representation of small molecule or inorganic crystal.
//!
//! Flat list of atom sites with minimal functionality.

// Copyright 2018 Global Phasing Ltd.
//
// Representation of a small molecule or inorganic crystal.
// Flat list of atom sites. Minimal functionality.

#ifndef GEMMI_SMALL_HPP_
#define GEMMI_SMALL_HPP_

#include <cctype>        // for isalpha
#include <algorithm>     // for any_of
#include <bitset>
#include <string>
#include <vector>
#include "elem.hpp"      // Element
#include "math.hpp"      // SMat33
#include "symmetry.hpp"  // find_spacegroup_by_name
#include "unitcell.hpp"  // UnitCell, Fractional
#include "util.hpp"      // vector_remove_if

namespace gemmi {

//! @brief Check if group operations are complete (closed under multiplication).
//! @param gops Group operations
//! @return True if gops is complete
inline bool is_complete(const GroupOps& gops) {
  for (Op op1 : gops.sym_ops)
    for (Op op2 : gops.sym_ops)
      if (gops.find_by_rotation((op1 * op2).rot) == nullptr)
        return false;
  return true;
}

//! @brief Convert symmetry operation strings to Op objects.
//! @param symops Vector of symmetry triplets (e.g., "x,y,z")
//! @return Vector of Op objects
inline std::vector<Op> triplets_to_ops(const std::vector<std::string>& symops) {
  std::vector<Op> ops;
  ops.reserve(symops.size());
  for (const std::string& xyz : symops)
    ops.push_back(parse_triplet(xyz));
  return ops;
}

//! @brief Small molecule or inorganic crystal structure.
struct SmallStructure {
  //! @brief Atom site in small structure.
  struct Site {
    std::string label;  //!< Atom label
    std::string type_symbol;  //!< Atom type symbol
    Fractional fract;  //!< Fractional coordinates
    double occ = 1.0;  //!< Occupancy
    double u_iso = 0.;  //!< Isotropic displacement parameter
    SMat33<double> aniso = {0, 0, 0, 0, 0, 0};  //!< Anisotropic displacement parameters
    int disorder_group = 0;  //!< Disorder group ID
    Element element = El::X;  //!< Element
    signed char charge = 0;  //!< Charge [-8, +8]

    //! @brief Get orthogonal coordinates.
    //! @param cell_ Unit cell
    //! @return Orthogonalized position
    Position orth(const gemmi::UnitCell& cell_) const {
      return cell_.orthogonalize(fract);
    }

    //! @brief Get element and charge as string.
    //! @return String like "Fe3+" or "O2-"
    std::string element_and_charge_symbol() const {
      std::string s = element.name();
      if (charge != 0) {
        s += std::to_string(std::abs(charge));
        s += charge > 0 ? '+' : '-';
      }
      return s;
    }
  };

  //! @brief Atom type with dispersion corrections.
  struct AtomType {
    std::string symbol;  //!< Atom type symbol
    Element element = El::X;  //!< Element
    signed char charge = 0;  //!< Charge [-8, +8]
    double dispersion_real;  //!< Real part of anomalous dispersion (f')
    double dispersion_imag;  //!< Imaginary part of anomalous dispersion (f'')
  };

  std::string name;  //!< Structure name
  UnitCell cell;  //!< Unit cell
  const SpaceGroup* spacegroup = nullptr;  //!< Space group pointer
  std::string spacegroup_hm;  //!< Hermann-Mauguin symbol
  std::string spacegroup_hall;  //!< Hall symbol
  int spacegroup_number = 0;  //!< Space group number
  std::vector<std::string> symops;  //!< Symmetry operations as triplets
  std::vector<Site> sites;  //!< Atom sites
  std::vector<AtomType> atom_types;  //!< Atom types with dispersion
  double wavelength = 0.;  //!< Wavelength (first if multiple)

  //! @brief Get all symmetry-related sites in unit cell.
  //! @return Vector of all sites including symmetry mates
  std::vector<Site> get_all_unit_cell_sites() const;

  //! @brief Determine and set space group from available data.
  //! @param order Order of preference (s=symops, h=hall, 1/2=H-M, n=number)
  void determine_and_set_spacegroup(const char* order) {
    spacegroup = nullptr;
    if (order)
      for (const char* c = order; *c != '\0' && spacegroup == nullptr; ++c) {
        try {
          GroupOps gops;
          spacegroup = determine_spacegroup_from(*c, gops);
          if (!spacegroup && *(c+1) == '.') {
            // If symops don't correspond to tabulated settings,
            // we can't set spacegroup, but we can set UnitCell::images.
            if (gops.order() == (int) symops.size() && is_complete(gops)) {
              cell.set_cell_images_from_groupops(gops);
              return;
            }
            ++c;
          }
        } catch (std::exception&) {}
      }
    setup_cell_images();
  }

  //! @brief Determine space group from specific source.
  //! @param c Source character (s=symops, h=hall, 1/2=H-M, n=number)
  //! @param gops Group operations (output for 's' option)
  //! @return Space group pointer or nullptr
  const SpaceGroup* determine_spacegroup_from(char c, GroupOps& gops) const {
    switch (lower(c)) {
      case 's':
        if (symops.empty())
          return nullptr;
        gops = split_centering_vectors(triplets_to_ops(symops));
        return find_spacegroup_by_ops(gops);
      case 'h':
        if (spacegroup_hall.empty())
          return nullptr;
        return find_spacegroup_by_ops(symops_from_hall(spacegroup_hall.c_str()));
      case '1':
      case '2': {
        if (spacegroup_hm.empty())
          return nullptr;
        char prefer[] = {c, '\0'};
        return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma, prefer);
      }
      case 'n':
        if (spacegroup_number == 0)
          return nullptr;
        return find_spacegroup_by_number(spacegroup_number);
      default:
        throw std::invalid_argument("determine_and_set_spacegroup(): wrong character in 'order'");
    }
  }

  //! @brief Check consistency of space group information.
  //! @return Error message string (empty if OK)
  std::string check_spacegroup() const {
    std::string err;
    if (!symops.empty())
      try {
        std::vector<Op> ops = triplets_to_ops(symops);
        for (Op& op : ops)
          op.wrap();
        std::sort(ops.begin(), ops.end());
        GroupOps gops = split_centering_vectors(ops);
        if (!is_complete(gops))
          cat_to(err, "symops list is incomplete or incorrect\n");
        else if (gops.all_ops_sorted() != ops)
          cat_to(err, "symops list is incorrect or incomplete or redundant\n");
        const SpaceGroup* sg = find_spacegroup_by_ops(gops);
        if (!sg)
          cat_to(err, "space group from symops not found in the table\n");
        else if (sg != spacegroup)
          cat_to(err, "space group from symops differs: ", sg->xhm(), '\n');
      } catch (std::exception& e) {
        cat_to(err, "error while processing symops: ", e.what(), '\n');
      }
    if (!spacegroup_hall.empty())
      try {
        const SpaceGroup* sg = find_spacegroup_by_ops(symops_from_hall(spacegroup_hall.c_str()));
        if (!sg)
          cat_to(err, "space group from Hall symbol (", spacegroup_hall,
                 ") not found in the table\n");
        else if (spacegroup != sg)
          cat_to(err, "space group from Hall symbol (", spacegroup_hall,
                 ") differs: ", sg->xhm(), '\n');
      } catch (std::exception& e) {
        cat_to(err, "error while processing Hall symbol: ", e.what(), '\n');
      }
    if (!spacegroup_hm.empty()) {
      const SpaceGroup* sg = find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma, "2");
      if (!sg)
        cat_to(err, "H-M symbol (", spacegroup_hm, ") not found in the table\n");
      else if (!spacegroup || strcmp(spacegroup->hm, sg->hm) != 0)
        cat_to(err, "space group from H-M symbol (", spacegroup_hm,
               ") differs: ", sg->hm, '\n');
    }
    if (spacegroup_number != 0 && spacegroup && spacegroup->number != spacegroup_number)
      cat_to(err, "space group number (", spacegroup_number, ") differs\n");
    return err;
  }

  //! @brief Find atom type by symbol.
  //! @param symbol Atom type symbol
  //! @return Pointer to AtomType or nullptr if not found
  const AtomType* get_atom_type(const std::string& symbol) const {
    for (const AtomType& at : atom_types)
      if (at.symbol == symbol)
        return &at;
    return nullptr;
  }

  //! @brief Get set of elements present in structure.
  //! @return Bitset with elements present
  //!
  //! Similar to Model::present_elements() from model.hpp.
  std::bitset<(size_t)El::END> present_elements() const {
    std::bitset<(size_t)El::END> table;
    for (const Site& atom : sites)
      table.set((size_t)atom.element.elem);
    return table;
  }

  //! @brief Remove all hydrogen atoms.
  void remove_hydrogens() {
    vector_remove_if(sites, [](const Site& a) { return a.element.is_hydrogen(); });
  }

  //! @brief Convert occupancies from chemical to crystallographic.
  //! @param max_dist Maximum distance for special position check
  //!
  //! Precondition: atoms on special positions have "chemical" occupancy (i.e. not divided
  //! by n for n-fold symmetry).
  void change_occupancies_to_crystallographic(double max_dist=0.4) {
    for (Site& site : sites) {
      int n_mates = cell.is_special_position(site.fract, max_dist);
      if (n_mates != 0)
        site.occ /= (n_mates + 1);
    }
  }

  //! @brief Setup cell images from space group.
  void setup_cell_images() {
    cell.set_cell_images_from_spacegroup(spacegroup);
  }
};

//! @brief Parse element and charge from label string.
//! @tparam T Destination type (must have element and charge fields)
//! @param label Label string (e.g., "Fe3+", "O2-")
//! @param dest Destination object to set element and charge
template<typename T>
inline void split_element_and_charge(const std::string& label, T* dest) {
  int len = label.size() > 1 && std::isalpha(label[1]) ? 2 : 1;
  dest->element = len == 1 ? impl::find_single_letter_element(label[0] & ~0x20)
                           : find_element(label.c_str());
  if (dest->element != El::X && (label.back() == '+' || label.back() == '-')) {
    int sign = label.back() == '+' ? 1 : -1;
    if (label.size() - len == 1)
      dest->charge = sign;
    else if (label.size() - len == 2 && label[len] >= '0' && label[len] <= '9')
      dest->charge = sign * (label[len] - '0');
  }
}

inline std::vector<SmallStructure::Site>
SmallStructure::get_all_unit_cell_sites() const {
  const double SPECIAL_POS_TOL = 0.4;
  std::vector<Site> all;
  for (const Site& site : sites) {
    size_t start = all.size();
    all.push_back(site);
    for (const FTransform& image : cell.images) {
      Fractional fpos = image.apply(site.fract);
      if (std::any_of(all.begin() + start, all.end(), [&](const Site& other) {
            return cell.distance_sq(fpos, other.fract) < sq(SPECIAL_POS_TOL);
          }))
        continue;
      all.push_back(site);
      all.back().fract = fpos;
    }
  }
  return all;
}

} // namespace gemmi
#endif
