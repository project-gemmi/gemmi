#include <assert.h>
#include <iostream>

#include <gemmi/symmetry.hpp>
namespace sym = gemmi::sym;

int main() {
  // operators
  sym::Op op = sym::parse_triplet("-y,x-y,z");        // one operation from P 3
  assert(op * op == sym::parse_triplet("-x+y,-x,z")); // and the other one
  assert(op * op == op.inverse());

  // iteration over all tabulated settings
  for (const sym::SpaceGroup& sg : sym::tables::main)
    std::cout << sg.number << ' ' << sg.xhm() << "   " << sg.hall << '\n';

  // selecting a space group
  const sym::SpaceGroup* c2 = sym::find_spacegroup_by_number(5);
  assert(c2->xhm() == "C 1 2 1");
  assert(c2->ccp4 == 5);

  const sym::SpaceGroup* i2 = sym::find_spacegroup_by_name("I2");
  assert(i2->number == 5);
  assert(i2->ccp4 == 4005);
  assert(i2->xhm() == "I 1 2 1");

  sym::GroupOps ops = i2->operations();
  for (const sym::Op& op : ops)
    std::cout << "   " << op.triplet();
  // output:   x,y,z   -x,y,-z   x+1/2,y+1/2,z+1/2   -x+1/2,y+1/2,-z+1/2
  ops.change_basis(sym::parse_triplet("x,y,x+z"));
  assert(sym::find_spacegroup_by_ops(ops) == c2);

  assert(sym::find_spacegroup_by_name("C m m e") ==  // new names may have 'e'
         sym::find_spacegroup_by_name("C m m a"));
}
