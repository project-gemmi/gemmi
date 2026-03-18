#include "doctest.h"

#include <vector>

#include <gemmi/ace_graph.hpp>

namespace {

void add_bond(gemmi::AceBondAdjacency& adj, size_t a, size_t b) {
  adj[a].push_back({b, gemmi::BondType::Single, false});
  adj[b].push_back({a, gemmi::BondType::Single, false});
}

std::vector<gemmi::AceAromaticAtom> make_benzene_like_atoms() {
  std::vector<gemmi::AceAromaticAtom> atoms(12);
  for (int i = 0; i < 6; ++i) {
    atoms[i].id = "C" + std::to_string(i + 1);
    atoms[i].el = gemmi::El::C;
    atoms[i].bonding_idx = 2;
    atoms[i].conn_atoms_no_metal = {
        (i + 5) % 6,
        (i + 1) % 6,
        i + 6,
    };
  }
  for (int i = 0; i < 6; ++i) {
    atoms[i + 6].id = "H" + std::to_string(i + 1);
    atoms[i + 6].el = gemmi::El::H;
    atoms[i + 6].bonding_idx = 1;
    atoms[i + 6].conn_atoms_no_metal = {i};
  }
  return atoms;
}

gemmi::AceBondAdjacency make_benzene_like_adjacency() {
  gemmi::AceBondAdjacency adj(12);
  for (size_t i = 0; i < 6; ++i) {
    add_bond(adj, i, (i + 1) % 6);
    add_bond(adj, i, i + 6);
  }
  return adj;
}

}  // namespace

TEST_CASE("Ace ring aromaticity can be computed without CodAtomInfo") {
  std::vector<gemmi::AceAromaticAtom> atoms = make_benzene_like_atoms();
  gemmi::AceBondAdjacency adj = make_benzene_like_adjacency();
  std::vector<gemmi::RingInfo> rings(1);
  rings[0].atoms = {0, 1, 2, 3, 4, 5};

  gemmi::set_ring_aromaticity_from_bonds(adj, atoms, rings);

  CHECK(rings[0].is_aromatic);
  CHECK(rings[0].is_aromatic_permissive);
}

TEST_CASE("Ace ring aromaticity still requires planar ring atoms") {
  std::vector<gemmi::AceAromaticAtom> atoms = make_benzene_like_atoms();
  gemmi::AceBondAdjacency adj = make_benzene_like_adjacency();
  std::vector<gemmi::RingInfo> rings(1);
  rings[0].atoms = {0, 1, 2, 3, 4, 5};
  atoms[0].bonding_idx = 3;

  gemmi::set_ring_aromaticity_from_bonds(adj, atoms, rings);

  CHECK(!rings[0].is_aromatic);
  CHECK(!rings[0].is_aromatic_permissive);
}
