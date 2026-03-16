// TypeScript bindings for emscripten-generated code.  Automatically generated at compile time.
declare namespace RuntimeExports {
    function writeArrayToMemory(array: any, buffer: any): void;
    let HEAPU8: any;
}
interface WasmModule {
}

type EmbindString = ArrayBuffer|Uint8Array|Uint8ClampedArray|Int8Array|string;
export interface ClassHandle {
  isAliasOf(other: ClassHandle): boolean;
  delete(): void;
  deleteLater(): this;
  isDeleted(): boolean;
  // @ts-ignore - If targeting lower than ESNext, this symbol might not exist.
  [Symbol.dispose](): void;
  clone(): this;
}
export interface UnitCellParameters extends ClassHandle {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

export interface UnitCell extends UnitCellParameters {
  volume: number;
  is_crystal(): boolean;
  fractionalize(_0: Position): Fractional;
  orthogonalize(_0: Fractional): Position;
}

export interface Isosurface extends ClassHandle {
  readonly last_error: string;
  resize_input(_0: number): void;
  set_size(_0: number, _1: number, _2: number): void;
  calculate(_0: number, _1: EmbindString): boolean;
  input_points(): any;
  input_values(): any;
  vertices(): any;
  segments(): any;
}

export interface Structure extends ClassHandle {
  cell: UnitCell;
  readonly length: number;
  get name(): string;
  set name(value: EmbindString);
  at(_0: number): Model | null;
}

export interface Model extends ClassHandle {
  num: number;
  readonly length: number;
  at(_0: number): Chain | null;
  count_occupancies(_0: Selection | null): number;
}

export interface Chain extends ClassHandle {
  readonly length: number;
  get name(): string;
  set name(value: EmbindString);
  at(_0: number): Residue | null;
}

export interface ResidueId extends ClassHandle {
  readonly seqid_string: string;
  get segment(): string;
  set segment(value: EmbindString);
  get name(): string;
  set name(value: EmbindString);
}

export interface Residue extends ResidueId {
  readonly length: number;
  get subchain(): string;
  set subchain(value: EmbindString);
  readonly entity_type_string: string;
  at(_0: number): Atom | null;
}

export interface Atom extends ClassHandle {
  altloc: number;
  charge: number;
  serial: number;
  occ: number;
  b_iso: number;
  pos: Position;
  get name(): string;
  set name(value: EmbindString);
  readonly element_uname: string;
}

export interface Selection extends ClassHandle {
}

export interface BondInfo extends ClassHandle {
  get_bond_lines(_0: Structure): void;
  bond_data_ptr(): number;
  bond_data_size(): number;
  add_monomer_cif(_0: EmbindString): void;
}

export interface SelectionResult extends ClassHandle {
  atom_data_ptr(): number;
  atom_data_size(): number;
  set_atom_indices(_0: Structure, _1: EmbindString, _2: number): void;
}

export interface Ccp4Map extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly mean: number;
  readonly rms: number;
  readonly last_error: string;
  read(_0: boolean): boolean;
  extract_isosurface(_0: number, _1: number, _2: number, _3: number, _4: number, _5: EmbindString): boolean;
  data(): any;
  isosurface_vertices(): any;
  isosurface_segments(): any;
}

export interface Dsn6Map extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly mean: number;
  readonly rms: number;
  readonly last_error: string;
  read(): boolean;
  extract_isosurface(_0: number, _1: number, _2: number, _3: number, _4: number, _5: EmbindString): boolean;
  data(): any;
  isosurface_vertices(): any;
  isosurface_segments(): any;
}

export interface Mtz extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly rmsd: number;
  readonly last_error: string;
  read(): boolean;
  calculate_map(_0: boolean): any;
  calculate_map_from_labels(_0: EmbindString, _1: EmbindString): any;
}

export type Fractional = [ number, number, number ];

export type Position = [ number, number, number ];

interface EmbindModule {
  UnitCellParameters: {};
  UnitCell: {
    new(_0: number, _1: number, _2: number, _3: number, _4: number, _5: number): UnitCell;
  };
  Isosurface: {
    new(): Isosurface;
  };
  Structure: {
    new(): Structure;
  };
  Model: {
    new(): Model;
  };
  Chain: {
    new(): Chain;
  };
  ResidueId: {};
  Residue: {
    new(): Residue;
  };
  Atom: {
    new(): Atom;
  };
  Selection: {
    new(): Selection;
  };
  BondInfo: {
    new(): BondInfo;
  };
  SelectionResult: {
    new(): SelectionResult;
  };
  Ccp4Map: {
    new(_0: EmbindString): Ccp4Map;
  };
  Dsn6Map: {
    new(_0: EmbindString): Dsn6Map;
  };
  Mtz: {
    new(_0: EmbindString): Mtz;
  };
  get_residue_names(_0: Structure): string;
  _read_structure(_0: EmbindString, _1: EmbindString, _2: EmbindString): Structure;
}

export type MainModule = WasmModule & typeof RuntimeExports & EmbindModule;
export default function MainModuleFactory (options?: unknown): Promise<MainModule>;
