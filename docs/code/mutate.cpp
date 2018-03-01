#include <gemmi/model.hpp>

void met_to_mse(gemmi::Structure& st) {
  for (gemmi::Model& model : st.models)
    for (gemmi::Chain& chain : model.chains)
      for (gemmi::Residue& res : chain.residues)
        if (res.name == "MET")
          if (gemmi::Atom* s_atom = res.find_atom("SD")) {
            res.name = "MSE";
            s_atom->name = "SE";
            s_atom->element = gemmi::El::Se;
          }
}
