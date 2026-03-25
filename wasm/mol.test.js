
const fs = require('node:fs');
const Gemmi = require('./gemmi.js')

const SCREENING_HIT_PDB = `\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
HETATM    1  C1  LIG A   1       3.200   0.500   0.000  1.00 20.00           C
HETATM    2  C2  LIG A   1       0.000   0.500   0.000  1.00 20.00           C
END
`;

const SCREENING_MISS_PDB = `\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
HETATM    1  C1  LIG A   1       3.200   0.500   0.000  1.00 20.00           C
HETATM    2  C2  LIG A   1       3.100   0.500   0.000  1.00 20.00           C
END
`;

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

function get_bonds(gemmi, bondInfo) {
  const ptr = bondInfo.bond_data_ptr();
  const len = bondInfo.bond_data_size();
  return new Int32Array(gemmi.HEAPU8.buffer, ptr, len).slice();
}

function find_atom_index(st, seqid, atom_name, chain_name = 'A', resname = null) {
  let atom_index = 0;
  for (let model of st) {
    for (let chain of model) {
      for (let residue of chain) {
        for (let atom of residue) {
          if (chain.name === chain_name && residue.seqid_string === seqid &&
              atom.name === atom_name && (resname === null || residue.name === resname))
            return atom_index;
          ++atom_index;
        }
      }
    }
  }
  throw new Error(`atom not found: ${seqid} ${atom_name}`);
}

function expectBondBetween(bonds, idx1, idx2, bondType) {
  for (let i = 0; i < bonds.length; i += 3) {
    const a1 = bonds[i];
    const a2 = bonds[i + 1];
    if ((a1 === idx1 && a2 === idx2) || (a1 === idx2 && a2 === idx1)) {
      expect(bonds[i + 2]).toBe(bondType);
      return;
    }
  }
  throw new Error(`bond not found between ${idx1} and ${idx2}`);
}

function expectNoBondBetween(bonds, idx1, idx2) {
  for (let i = 0; i < bonds.length; i += 3) {
    const a1 = bonds[i];
    const a2 = bonds[i + 1];
    if ((a1 === idx1 && a2 === idx2) || (a1 === idx2 && a2 === idx1))
      throw new Error(`unexpected bond found between ${idx1} and ${idx2}`);
  }
}

function get_atom_position(st, seqid, atom_name, chain_name = 'A', resname = null) {
  for (let model of st) {
    for (let chain of model) {
      for (let residue of chain) {
        for (let atom of residue) {
          if (chain.name === chain_name && residue.seqid_string === seqid &&
              atom.name === atom_name && (resname === null || residue.name === resname))
            return [atom.pos[0], atom.pos[1], atom.pos[2]];
        }
      }
    }
  }
  throw new Error(`atom not found: ${seqid} ${atom_name}`);
}

function expectPositionClose(actual, expected, epsilon = 1e-6) {
  expect(Math.abs(actual[0] - expected[0])).toBeLessThan(epsilon);
  expect(Math.abs(actual[1] - expected[1])).toBeLessThan(epsilon);
  expect(Math.abs(actual[2] - expected[2])).toBeLessThan(epsilon);
}

function get_image_codes(images) {
  const codes = [];
  for (let i = 0; i < images.size(); ++i) {
    const image = images.get(i);
    codes.push(image.symmetry_code(true));
    if (typeof image.delete === 'function')
      image.delete();
  }
  return codes;
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

test('exports inferred peptide bonds', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/3dg1_final.cif';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);
  const bondInfo = new gemmi.BondInfo();
  bondInfo.get_bond_lines(st);
  const bonds = get_bonds(gemmi, bondInfo);
  const c_idx = find_atom_index(st, '1', 'C');
  const n_idx = find_atom_index(st, '2', 'N');
  expectBondBetween(bonds, c_idx, n_idx, 1);

  bondInfo.delete();
  st.delete();
});

test('exports inferred peptide bonds for pdb input', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/1orc.pdb';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);
  const bondInfo = new gemmi.BondInfo();
  bondInfo.get_bond_lines(st);
  const bonds = get_bonds(gemmi, bondInfo);
  const c_idx = find_atom_index(st, '14', 'C');
  const n_idx = find_atom_index(st, '15', 'N');
  expectBondBetween(bonds, c_idx, n_idx, 1);

  bondInfo.delete();
  st.delete();
});

test('uses embedded chemcomp bonds before monlib fallback', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/5i55.cif';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);
  const bondInfo = new gemmi.BondInfo();
  bondInfo.get_bond_lines(st);
  const bonds = get_bonds(gemmi, bondInfo);
  const c_idx = find_atom_index(st, '102', 'C', 'A', 'ACT');
  const o_idx = find_atom_index(st, '102', 'O', 'A', 'ACT');
  const ch3_idx = find_atom_index(st, '102', 'CH3', 'A', 'ACT');
  expectBondBetween(bonds, c_idx, o_idx, 2);
  expectBondBetween(bonds, c_idx, ch3_idx, 1);

  bondInfo.delete();
  st.delete();
});

test('ignores explicit bonds to symmetry images', async () => {
  const gemmi = await Gemmi();
  const cifText = `data_symm_link
_cell.length_a 10
_cell.length_b 10
_cell.length_c 10
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . GLY A 1 1 0 0 0 1 20 A 1 GLY 1
ATOM 2 C CB . ALA A 1 2 1 0 0 1 20 A 2 ALA 1

loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
sym1 covale A GLY 1 CA A ALA 2 CB A 1 A 2 1_555 18_555
`;
  const st = gemmi.read_structure(Buffer.from(cifText), 'symm-link.cif');
  const bondInfo = new gemmi.BondInfo();
  bondInfo.get_bond_lines(st);
  const bonds = get_bonds(gemmi, bondInfo);
  const caIdx = find_atom_index(st, '1', 'CA');
  const cbIdx = find_atom_index(st, '2', 'CB');
  expectNoBondBetween(bonds, caIdx, cbIdx);

  bondInfo.delete();
  st.delete();
});

test('matches cross-symmetry bonds by nearest image, not reported_sym', async () => {
  const gemmi = await Gemmi();
  const cifText = `data_symm_link
_cell.length_a 10
_cell.length_b 10
_cell.length_c 10
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . GLY A 1 1 0 9.8 0 1 20 A 1 GLY 1
ATOM 2 C CB . ALA A 1 2 0 0.2 0 1 20 A 2 ALA 1

loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
sym1 covale A GLY 1 CA A ALA 2 CB A 1 A 2 18_555 1_555
`;
  const st = gemmi.read_structure(Buffer.from(cifText), 'symm-link-nearest.cif');
  const point = get_atom_position(st, '1', 'CA');
  const images = gemmi.get_nearby_sym_ops(st, point, 1.0);
  expect(get_image_codes(images)).toEqual(['1_565']);

  const crossSymBonds = new gemmi.CrossSymBonds();
  const image0 = images.get(0);
  crossSymBonds.find(st, image0);
  const bonds = get_bonds(gemmi, crossSymBonds);
  const caIdx = find_atom_index(st, '1', 'CA');
  const cbIdx = find_atom_index(st, '2', 'CB');
  expectBondBetween(bonds, caIdx, cbIdx, 1);

  if (typeof image0.delete === 'function')
    image0.delete();
  crossSymBonds.delete();
  images.delete();
  st.delete();
});

test('lists monomer names missing chemcomp data', async () => {
  const gemmi = await Gemmi();
  const cif_path = '../tests/5i55.cif';
  const cif_st = gemmi.read_structure(fs.readFileSync(cif_path), cif_path);
  expect(gemmi.get_missing_monomer_names(cif_st)).toBe('');
  cif_st.delete();

  const pdb_path = '../tests/1orc.pdb';
  const pdb_st = gemmi.read_structure(fs.readFileSync(pdb_path), pdb_path);
  expect(gemmi.get_missing_monomer_names(pdb_st))
    .toBe(gemmi.get_residue_names(pdb_st));
  pdb_st.delete();
});

test('finds nearby symmetry images in a real structure', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/4oz7.pdb';
  const st = gemmi.read_structure(fs.readFileSync(path), path);
  const point = get_atom_position(st, '208', 'O', 'B', 'HOH');
  const images = gemmi.get_nearby_sym_ops(st, point, 3.0);
  expect(get_image_codes(images)).toEqual(['4_355', '3_545']);

  const image0 = images.get(0);
  expect(image0.symmetry_code(true)).toBe('4_355');
  const moved = gemmi.get_sym_image(st, image0);
  expect(count_atoms_with_iterators(moved)).toBe(count_atoms_with_iterators(st));

  moved.delete();
  if (typeof image0.delete === 'function')
    image0.delete();
  images.delete();
  st.delete();
});

test('uses residue screening for translated symmetry images', async () => {
  const gemmi = await Gemmi();
  const st = gemmi.read_structure(Buffer.from(SCREENING_HIT_PDB), 'screening-hit.pdb');
  const point = [0, 10.9, 0];
  const images = gemmi.get_nearby_sym_ops(st, point, 1.0);
  expect(get_image_codes(images)).toEqual(['1_565']);

  const image0 = images.get(0);
  expect(image0.sym_idx).toBe(0);
  expect(image0.same_asu()).toBe(false);
  const moved = gemmi.get_sym_image(st, image0);
  const origPos = get_atom_position(st, '1', 'C2', 'A', 'LIG');
  const movedPos = get_atom_position(moved, '1', 'C2', 'A', 'LIG');
  expectPositionClose(movedPos, [origPos[0], origPos[1] + st.cell.b, origPos[2]]);

  moved.delete();
  if (typeof image0.delete === 'function')
    image0.delete();
  images.delete();
  st.delete();
});

test('rejects screening-only residue misses', async () => {
  const gemmi = await Gemmi();
  const st = gemmi.read_structure(Buffer.from(SCREENING_MISS_PDB), 'screening-miss.pdb');
  const point = [0, 10.9, 0];
  const images = gemmi.get_nearby_sym_ops(st, point, 1.0);
  expect(images.size()).toBe(0);

  images.delete();
  st.delete();
});
