#include <gemmi/resinfo.hpp>

gemmi::ResidueInfo info = gemmi::find_tabulated_residue("ALA");
bool is_it_aminoacid = info.is_amino_acid();
int approximate_number_of_h_atoms = info.hydrogen_count;
