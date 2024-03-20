
const fs = require('node:fs');
const Gemmi = require('./gemmi.js')

function count_atoms_with_indices(st) {
  let occ_sum = 0.0;
  for (let model_idx = 0; model_idx < st.length; ++model_idx) {
    const model = st.at(model_idx);
    for (let chain_idx = 0; chain_idx < model.length; ++chain_idx) {
      const chain = model.at(chain_idx);
      for (let res_idx = 0; res_idx < chain.length; ++res_idx) {
        const res = chain.at(res_idx);
        for (let atom_idx = 0; atom_idx < res.length; ++atom_idx) {
          const atom = res.at(atom_idx);
          occ_sum += atom.occ;
        }
      }
    }
  }
  return occ_sum;
}

function count_atoms_with_iterators(st) {
  let occ_sum = 0.0;
  for (let model of st) {
    for (let chain of model) {
      for (let res of chain) {
        for (let atom of res) {
          occ_sum += atom.occ;
        }
      }
    }
  }
  return occ_sum;
}


test('counts atom occupancies', () => {
  Gemmi().then((gemmi) => {
    const path = '../tests/3dg1_final.cif';
    const result = 40.5;
    const buffer = fs.readFileSync(path);
    const st = gemmi.read_structure(buffer, path);
    // no WASM memory allocations in this function
    const heap_length = gemmi.HEAPU8.length;
    expect(count_atoms_with_indices(st)).toBe(result);
    expect(count_atoms_with_iterators(st)).toBe(result);
    expect(st.at(0).count_occupancies(null)).toBe(result);
    expect(gemmi.HEAPU8.length).toBe(heap_length);
    st.delete();
  });
});
