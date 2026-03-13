
export interface UnitCell {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
  delete(): void;
}

export interface Structure {
  readonly name: string;
  readonly cell: UnitCell;
  readonly length: number;
  at(index: number): Model;
  delete(): void;
}

export interface Model {
  readonly num: number;
  readonly length: number;
  at(index: number): Chain;
  count_occupancies(_0: any): number;
  delete(): void;
}

export interface Chain {
  readonly name: string;
  readonly length: number;
  at(index: number): Residue;
  delete(): void;
}

export interface Residue {
  readonly seqid_string: string;
  readonly segment: string;
  readonly name: string;
  readonly subchain: string;
  readonly entity_type_string: string;
  readonly length: number;
  at(index: number): Atom;
  delete(): void;
}

export type Position = [number, number, number];

export interface Atom {
  readonly name: string;
  readonly altloc: number;
  readonly charge: number;
  readonly element_uname: string;
  readonly serial: number;
  readonly pos: Position;
  readonly occ: number;
  readonly b_iso: number;
  delete(): void;
}

export interface Selection {
  delete(): void;
}

export interface BondInfo {
  add_monomer_cif(cif_text: string): void;
  get_bond_lines(st: Structure): void;
  bond_data_ptr(): number;
  bond_data_size(): number;
  delete(): void;
}

export interface SelectionResult {
  set_atom_indices(st: Structure, cid: string, model_index: number): void;
  atom_data_ptr(): number;
  atom_data_size(): number;
  delete(): void;
}

export interface Ccp4Map {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly mean: number;
  readonly rms: number;
  readonly last_error: string;
  read(_0: boolean): boolean;
  data(): Float32Array;
  delete(): void;
}

export interface Mtz {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly rmsd: number;
  readonly last_error: string;
  read(_0: number, _1: number): boolean;
  calculate_map(_0: boolean): any;
  calculate_map_from_labels(_0: string, _1: string): any;
  delete(): void;
}

export interface Module {
  read_structure(buf: string|ArrayBuffer, name: string, format?: string): Structure;
  get_residue_names(st: Structure): string;
  Selection: {
    new (): Selection;
  };
  BondInfo: {
    new (): BondInfo;
  };
  SelectionResult: {
    new (): SelectionResult;
  };
  readCcp4Map(map_buf: string|ArrayBuffer, expand_symmetry?: boolean): Ccp4Map;
  readMtz(mtz_buf: string|ArrayBuffer): Mtz;
  HEAPU8: Uint8Array;
}
