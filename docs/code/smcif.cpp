#include <cassert>
#include <gemmi/cif.hpp>
#include <gemmi/smcif.hpp>

int main() {
  auto block = gemmi::cif::read_file("1011031.cif").sole_block();
  gemmi::SmallStructure SiC = gemmi::make_small_structure_from_block(block);
  assert(SiC.cell.a == 4.358);
  assert(SiC.spacegroup_hm == "F -4 3 m");
  assert(SiC.sites.size() == 2);
  assert(SiC.get_all_unit_cell_sites().size() == 8);
}
