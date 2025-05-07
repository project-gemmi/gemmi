// Copyright Global Phasing Ltd.
//
// Functions for working with sequences (other than alignment).

#ifndef GEMMI_DSSP_HPP_
#define GEMMI_DSSP_HPP_

#include "topo.hpp"
#include "neighbor.hpp"

namespace gemmi {

struct HBond {
  Topo::ResInfo *donor = nullptr, *acceptor = nullptr;
  char alt1 = '\0';
  char alt2 = '\0';
  double energy = 0;
};

std::vector<HBond>
dssp_determine_hydrogen_bonds(NeighborSearch& ns, Topo::ChainInfo& cinfo);

} // namespace gemmi
#endif
