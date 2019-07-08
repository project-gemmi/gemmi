#include <iostream>
#include <gemmi/model.hpp>

using std::cout;

void iterate_over_everything(gemmi::Structure& st) {
  cout << "Unit cell: " << st.cell.a << " " << st.cell.b /*etc.*/ << '\n';
  for (const gemmi::Model& model : st.models) {
    cout << "Model #" << model.name << '\n';
    for (const gemmi::Chain& chain : model.chains) {
      cout << " Chain " << chain.name << '\n';
      for (const gemmi::Residue& res : chain.residues) {
        cout << "  Residue " << res.name             // string
                      << " " << res.seqid.str()      // sequence id
                      << " " << res.segment          // string
                      // het_flag: 'A' = ATOM, 'H' = HETATM, 0 = not set
                      << " " << (res.het_flag == 'H' ? "HETA" : "ATOM")
             << '\n';
        for (const gemmi::Atom& atom : res.atoms) {
          cout << "   " << atom.name
                 << " " << atom.element.name() // Element
                 << " " << atom.altloc         // char
                 << " " << (int) atom.charge   // signed char
                 << " " << atom.pos.x << " " << atom.pos.y // and .z
                 << " " << atom.occ            // occupancy
                 << " " << atom.b_iso          // isotropic B
                 // anisotropic ADP: u11 u22 u33 u12 u13 u23
                 << " " << atom.u11 // (0 if not set)
                 << '\n';
        }
      }
    }
  }
}
