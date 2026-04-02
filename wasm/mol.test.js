
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

const SELECTION_PDB = `\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.200   0.000   0.000  1.00 20.00           C
ATOM      3  N   GLY A   2       2.400   0.000   0.000  1.00 20.00           N
ATOM      4  CA  GLY A   2       3.600   0.000   0.000  1.00 20.00           C
ATOM      5  N   SER B   1       0.000   3.000   0.000  1.00 20.00           N
ATOM      6  CA  SER B   1       1.200   3.000   0.000  1.00 20.00           C
TER
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

function find_atom(st, predicate) {
  for (let model of st) {
    for (let chain of model) {
      for (let residue of chain) {
        for (let atom of residue) {
          if (predicate(atom, residue, chain, model))
            return atom;
        }
      }
    }
  }
  throw new Error('atom not found');
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

test('exposes atom metal classification', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/4oz7.pdb';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);

  const copper = find_atom(st, (atom) => atom.element_uname === 'CU');
  const carbon = find_atom(st, (atom, residue) =>
    residue.name === 'ALA' && atom.name === 'CA');

  expect(copper.is_metal).toBe(true);
  expect(carbon.is_metal).toBe(false);

  st.delete();
});

test('exposes struct_conn metadata and append helpers', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/4oz7.pdb';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);

  const connections = st.connections;
  expect(connections.size()).toBeGreaterThan(0);
  const connection = connections.get(0);
  expect(connection.name).toBe('disulf1');
  expect(connection.type).toBe(gemmi.ConnectionType.Disulf);
  expect(connection.asu).toBe(gemmi.Asu.Same);
  expect(connection.partner1.chain_name).toBe('A');
  expect(connection.partner1.res_id.name).toBe('CYS');
  expect(connection.partner1.res_id.seqid_string).toBe('4');
  expect(connection.partner1.atom_name).toBe('SG');
  expect(connection.partner2.chain_name).toBe('A');
  expect(connection.partner2.res_id.seqid_string).toBe('10');
  expect(connection.partner2.atom_name).toBe('SG');

  const added = new gemmi.Connection();
  const p1 = new gemmi.AtomAddress();
  const p1res = new gemmi.ResidueId();
  p1res.name = 'ALA';
  p1res.seqid_string = '2';
  p1.chain_name = 'A';
  p1.res_id = p1res;
  p1.atom_name = 'N';

  const p2 = new gemmi.AtomAddress();
  const p2res = new gemmi.ResidueId();
  p2res.name = 'SER';
  p2res.seqid_string = '3';
  p2.chain_name = 'A';
  p2.res_id = p2res;
  p2.atom_name = 'OG';
  p2.altloc = 'A';

  added.name = 'custom1';
  added.link_id = 'link-demo';
  added.type = gemmi.ConnectionType.Hydrog;
  added.asu = gemmi.Asu.Different;
  added.partner1 = p1;
  added.partner2 = p2;
  added.reported_distance = 2.75;
  st.add_connection(added);

  const appended = st.connections.get(st.connections.size() - 1);
  expect(appended.name).toBe('custom1');
  expect(appended.link_id).toBe('link-demo');
  expect(appended.type).toBe(gemmi.ConnectionType.Hydrog);
  expect(appended.asu).toBe(gemmi.Asu.Different);
  expect(appended.partner2.altloc).toBe('A');
  expect(appended.reported_distance).toBeCloseTo(2.75);
  const cifText = gemmi.make_mmcif_string(st);
  expect(cifText).toMatch(/_struct_conn\.id/);
  expect(cifText).toMatch(/custom1/);
  expect(cifText).toMatch(/link-demo/);

  appended.delete();
  added.delete();
  p1.delete();
  p1res.delete();
  p2.delete();
  p2res.delete();
  connection.delete();
  connections.delete();
  st.delete();
});

test('removes selected atoms, residues and chains by CID', async () => {
  const gemmi = await Gemmi();

  const atomStructure = gemmi.read_structure(Buffer.from(SELECTION_PDB), 'selection-atom.pdb');
  const atomSelection = new gemmi.Selection('//A/1/N');
  atomSelection.remove_selected(atomStructure);
  expect(count_atoms_with_iterators(atomStructure)).toBe(5);
  expect(atomStructure.at(0).at(0).at(0).length).toBe(1);
  expect(atomStructure.at(0).at(0).at(0).at(0).name).toBe('CA');
  atomSelection.delete();
  atomStructure.delete();

  const residueStructure = gemmi.read_structure(Buffer.from(SELECTION_PDB), 'selection-residue.pdb');
  const residueSelection = new gemmi.Selection('//A/2');
  residueSelection.remove_selected(residueStructure);
  expect(count_atoms_with_iterators(residueStructure)).toBe(4);
  expect(residueStructure.at(0).at(0).length).toBe(1);
  expect(residueStructure.at(0).at(0).at(0).name).toBe('ALA');
  residueSelection.delete();
  residueStructure.delete();

  const chainStructure = gemmi.read_structure(Buffer.from(SELECTION_PDB), 'selection-chain.pdb');
  const chainSelection = new gemmi.Selection('//A');
  chainSelection.remove_selected(chainStructure);
  expect(chainStructure.at(0).length).toBe(1);
  expect(chainStructure.at(0).at(0).name).toBe('B');
  expect(count_atoms_with_iterators(chainStructure)).toBe(2);
  chainSelection.delete();
  chainStructure.delete();
});

test('keeps only selected chains with remove_not_selected', async () => {
  const gemmi = await Gemmi();
  const st = gemmi.read_structure(Buffer.from(SELECTION_PDB), 'selection-keep.pdb');
  const sel = new gemmi.Selection('//B');

  sel.remove_not_selected(st);
  expect(st.at(0).length).toBe(1);
  expect(st.at(0).at(0).name).toBe('B');
  expect(st.at(0).at(0).length).toBe(1);
  expect(count_atoms_with_iterators(st)).toBe(2);

  sel.delete();
  st.delete();
});

test('builds structure hierarchy with append methods and writable setters', async () => {
  const gemmi = await Gemmi();

  const st = new gemmi.Structure();
  const model = new gemmi.Model();
  const chain = new gemmi.Chain();
  const water = new gemmi.Residue();
  const oxygen = new gemmi.Atom();
  const metalResidue = new gemmi.Residue();
  const zinc = new gemmi.Atom();

  model.num = 7;
  chain.name = 'Z';

  water.name = 'HOH';
  water.seqid_string = '101B';
  oxygen.name = 'O';
  oxygen.element_uname = 'O';
  oxygen.pos = [1.0, 2.0, 3.0];
  water.add_atom(oxygen);

  metalResidue.name = 'ZN';
  metalResidue.set_seqid(102, 'A');
  zinc.name = 'ZN';
  zinc.set_element('zn');
  zinc.pos = [4.0, 5.0, 6.0];
  zinc.charge = 2;
  metalResidue.add_atom(zinc);

  chain.add_residue(water);
  chain.add_residue(metalResidue);
  model.add_chain(chain);
  st.add_model(model);

  oxygen.delete();
  zinc.delete();
  water.delete();
  metalResidue.delete();
  chain.delete();
  model.delete();

  expect(st.length).toBe(1);
  expect(st.at(0).num).toBe(7);
  expect(st.at(0).length).toBe(1);
  expect(st.at(0).at(0).name).toBe('Z');
  expect(st.at(0).at(0).length).toBe(2);
  expect(st.at(0).at(0).at(0).seqid_string).toBe('101B');
  expect(st.at(0).at(0).at(0).at(0).element_uname).toBe('O');
  expectPositionClose(st.at(0).at(0).at(0).at(0).pos, [1.0, 2.0, 3.0]);
  expect(st.at(0).at(0).at(1).seqid_string).toBe('102A');
  expect(st.at(0).at(0).at(1).at(0).element_uname).toBe('ZN');
  expectPositionClose(st.at(0).at(0).at(1).at(0).pos, [4.0, 5.0, 6.0]);
  expect(st.at(0).at(0).at(1).at(0).charge).toBe(2);
  expect(st.at(0).at(0).at(1).at(0).is_metal).toBe(true);
  expect(count_atoms_with_iterators(st)).toBe(2);

  st.delete();
});

test('exposes structure sites metadata and append helpers', async () => {
  const gemmi = await Gemmi();
  const path = '../tests/5moo_header.pdb';
  const buffer = fs.readFileSync(path);
  const st = gemmi.read_structure(buffer, path);

  const sites = st.sites;
  expect(sites.size()).toBe(3);
  const site = sites.get(0);
  expect(site.name).toBe('AC1');
  expect(site.evidence_code).toBe('SOFTWARE');
  expect(site.residue.chain_name).toBe('A');
  expect(site.residue.res_id.name).toBe('CA');
  expect(site.residue.res_id.seqid_string).toBe('301');
  expect(site.members.size()).toBe(6);
  const member = site.members.get(0);
  expect(member.residue_num).toBe(1);
  expect(member.auth.chain_name).toBe('A');
  expect(member.auth.res_id.name).toBe('GLU');
  expect(member.auth.res_id.seqid_string).toBe('70');
  expect(member.label_seq_string).toBe('');
  expect(member.label_alt_id).toBe('');

  const addedSite = new gemmi.StructSite('ZZ1');
  const siteResidue = new gemmi.AtomAddress();
  const siteResidueId = new gemmi.ResidueId();
  siteResidueId.name = 'CA';
  siteResidueId.seqid_string = '301';
  siteResidue.chain_name = 'A';
  siteResidue.res_id = siteResidueId;
  addedSite.residue = siteResidue;
  addedSite.evidence_code = 'Software';
  addedSite.residue_count = 1;
  addedSite.details = 'demo site';

  const addedMember = new gemmi.StructSiteMember();
  const memberAddr = new gemmi.AtomAddress();
  const memberResidueId = new gemmi.ResidueId();
  memberResidueId.name = 'GLU';
  memberResidueId.seqid_string = '70';
  memberAddr.chain_name = 'A';
  memberAddr.res_id = memberResidueId;
  addedMember.residue_num = 1;
  addedMember.auth = memberAddr;
  addedMember.label_comp_id = 'GLU';
  addedMember.label_asym_id = 'A';
  addedMember.label_seq_string = '70';
  addedMember.label_atom_id = 'CA';
  addedMember.label_alt_id = 'B';
  addedMember.symmetry = '1_555';
  addedSite.add_member(addedMember);
  st.add_site(addedSite);

  const added = st.sites.get(st.sites.size() - 1);
  expect(added.name).toBe('ZZ1');
  expect(added.members.size()).toBe(1);
  expect(added.members.get(0).label_seq_string).toBe('70');
  expect(added.members.get(0).label_alt_id).toBe('B');
  expect(gemmi.make_mmcif_string(st)).toMatch(/_struct_site\.id/);
  expect(gemmi.make_mmcif_string(st)).toMatch(/ZZ1/);

  member.delete();
  site.delete();
  sites.delete();
  added.delete();
  addedMember.delete();
  memberAddr.delete();
  memberResidueId.delete();
  siteResidue.delete();
  siteResidueId.delete();
  addedSite.delete();
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
