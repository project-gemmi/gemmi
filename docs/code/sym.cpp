#include <assert.h>
#include <iostream>

#include <gemmi/symmetry.hpp>

int main() {
  // operators
  gemmi::Op op = gemmi::parse_triplet("-y,x-y,z"); // one operation from P 3
  assert(op * op == gemmi::parse_triplet("-x+y,-x,z")); // and the other one
  assert(op * op == op.inverse());

  // iteration over all tabulated settings
  for (const gemmi::SpaceGroup& sg : gemmi::spacegroup_tables::main)
    std::cout << sg.number << ' ' << sg.xhm() << "   " << sg.hall << '\n';

  // selecting a space group
  const gemmi::SpaceGroup* c2 = gemmi::find_spacegroup_by_number(5);
  assert(c2->xhm() == "C 1 2 1");
  assert(c2->ccp4 == 5);

  const gemmi::SpaceGroup* i2 = gemmi::find_spacegroup_by_name("I2");
  assert(i2->number == 5);
  assert(i2->ccp4 == 4005);
  assert(i2->xhm() == "I 1 2 1");

  gemmi::GroupOps ops = i2->operations();
  for (gemmi::Op operation : ops)
    std::cout << "   " << operation.triplet();
  // output:   x,y,z   -x,y,-z   x+1/2,y+1/2,z+1/2   -x+1/2,y+1/2,-z+1/2
  ops.change_basis_forward(gemmi::parse_triplet("x,y,x+z"));
  assert(gemmi::find_spacegroup_by_ops(ops) == c2);

  assert(gemmi::find_spacegroup_by_name("C m m e") ==  // new name with 'e'
         gemmi::find_spacegroup_by_name("C m m a"));
}
