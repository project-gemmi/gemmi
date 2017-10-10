#include <iostream>
#include <gemmi/model.hpp>

using namespace std;

void iterate_over_everything(gemmi::Structure& st) {
  cout << "Unit cell: " << st.cell.a << " " << st.cell.b /*etc.*/ << endl;
  for (const gemmi::Model& model : st.models) {
    cout << "Model #" << model.name << endl;
    for (const gemmi::Chain& chain : model.chains) {
      cout << " Chain " << chain.name << endl;
      for (const gemmi::Residue& res : chain.residues) {
        cout << "  Residue " << res.name        // string
                      << " " << res.seq_id      // int
                      << " " << res.auth_seq_id // int
                      << " " << res.ins_code    // char
                      << " " << res.segment     // string
             << endl;
        for (const gemmi::Atom& atom : res.atoms) {
          // atom.group is 'A' for ATOM, 'H' for HETATM, or 0 if not set
          cout << "   " << (atom.group == 'H' ? "HETA" : "ATOM")
                 << " " << atom.name
                 << " " << atom.element.name() // Element
                 << " " << atom.altloc         // char
                 << " " << (int) atom.charge   // signed char
                 << " " << atom.pos.x << " " << atom.pos.y // and .z
                 << " " << atom.occ            // occupancy
                 << " " << atom.b_iso          // isotropic B
                 // anisotropic ADP: u11 u22 u33 u12 u13 u23
                 << " " << atom.u11 // (0 if not set)
                 << endl;
        }
      }
    }
  }
}
// vim:sw=2:ts=2:et:path^=../include,../third_party
