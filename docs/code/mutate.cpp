#include <gemmi/model.hpp>

void met_to_mse(gemmi::Structure& st) {
  for (gemmi::Model& model : st.models)
    for (gemmi::Chain& chain : model.chains)
      for (gemmi::Residue& res : chain.residues)
        if (res.name == "MET") {
          res.name = "MSE";
          for (gemmi::Atom& atom : res.atoms)
            if (atom.name == "SD") {
              atom.name = "SE";
              atom.element = gemmi::El::Se;
            }
        }
}
