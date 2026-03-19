
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


test('counts atom occupancies', async () => {
  const gemmi = await Gemmi();
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

test('exposes residue secondary structure from file', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/1orc.pdb';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);
  const chain = st.at(0).at(0);
  const bySeqid = new Map();
  for (let i = 0; i < chain.length; ++i) {
    const res = chain.at(i);
    bySeqid.set(res.seqid_string, res);
  }

  expect(bySeqid.get('7').ss_from_file).toBe(gemmi.ResidueSs.Helix);
  expect(bySeqid.get('7').ss_from_file_string).toBe('Helix');
  expect(bySeqid.get('7').strand_sense_from_file)
    .toBe(gemmi.ResidueStrandSense.NotStrand);
  expect(bySeqid.get('15').ss_from_file).toBe(gemmi.ResidueSs.Coil);
  expect(bySeqid.get('39').ss_from_file).toBe(gemmi.ResidueSs.Strand);
  expect(bySeqid.get('39').strand_sense_from_file)
    .toBe(gemmi.ResidueStrandSense.First);
  expect(bySeqid.get('50').strand_sense_from_file)
    .toBe(gemmi.ResidueStrandSense.Antiparallel);
  expect(bySeqid.get('56C').strand_sense_from_file_string)
    .toBe('Antiparallel');

  st.delete();
});
