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
export interface NearestImage extends ClassHandle {
  sym_idx: number;
  readonly pbc_shift_x: number;
  readonly pbc_shift_y: number;
  readonly pbc_shift_z: number;
  same_asu(): boolean;
  dist(): number;
  symmetry_code(_0: boolean): string;
}

export interface NearestImageVector extends ClassHandle, Iterable<NearestImage> {
  push_back(_0: NearestImage): void;
  resize(_0: number, _1: NearestImage): void;
  size(): number;
  get(_0: number): NearestImage | undefined;
  set(_0: number, _1: NearestImage): boolean;
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
  triangles(): any;
}

export interface ResidueSsValue<T extends number> {
  value: T;
}
export type ResidueSs = ResidueSsValue<0>|ResidueSsValue<1>|ResidueSsValue<2>;

export interface ResidueStrandSenseValue<T extends number> {
  value: T;
}
export type ResidueStrandSense = ResidueStrandSenseValue<0>|ResidueStrandSenseValue<1>|ResidueStrandSenseValue<2>|ResidueStrandSenseValue<-1>;

export interface AsuValue<T extends number> {
  value: T;
}
export type Asu = AsuValue<0>|AsuValue<1>|AsuValue<2>;

export interface ConnectionTypeValue<T extends number> {
  value: T;
}
export type ConnectionType = ConnectionTypeValue<0>|ConnectionTypeValue<1>|ConnectionTypeValue<2>|ConnectionTypeValue<3>|ConnectionTypeValue<4>;

export interface AtomAddress extends ClassHandle {
  res_id: ResidueId;
  get chain_name(): string;
  set chain_name(value: EmbindString);
  get atom_name(): string;
  set atom_name(value: EmbindString);
  get altloc(): string;
  set altloc(value: EmbindString);
  str(): string;
}

export interface Connection extends ClassHandle {
  type: ConnectionType;
  asu: Asu;
  partner1: AtomAddress;
  partner2: AtomAddress;
  reported_distance: number;
  get name(): string;
  set name(value: EmbindString);
  get link_id(): string;
  set link_id(value: EmbindString);
}

export interface ConnectionVector extends ClassHandle, Iterable<Connection> {
  push_back(_0: Connection): void;
  resize(_0: number, _1: Connection): void;
  size(): number;
  get(_0: number): Connection | undefined;
  set(_0: number, _1: Connection): boolean;
}

export interface StructSiteMember extends ClassHandle {
  auth: AtomAddress;
  residue_num: number;
  get label_comp_id(): string;
  set label_comp_id(value: EmbindString);
  get label_asym_id(): string;
  set label_asym_id(value: EmbindString);
  get label_seq_string(): string;
  set label_seq_string(value: EmbindString);
  get label_atom_id(): string;
  set label_atom_id(value: EmbindString);
  get label_alt_id(): string;
  set label_alt_id(value: EmbindString);
  get symmetry(): string;
  set symmetry(value: EmbindString);
  get details(): string;
  set details(value: EmbindString);
}

export interface StructSiteMemberVector extends ClassHandle, Iterable<StructSiteMember> {
  push_back(_0: StructSiteMember): void;
  resize(_0: number, _1: StructSiteMember): void;
  size(): number;
  get(_0: number): StructSiteMember | undefined;
  set(_0: number, _1: StructSiteMember): boolean;
}

export interface StructSite extends ClassHandle {
  residue: AtomAddress;
  members: StructSiteMemberVector;
  residue_count: number;
  get name(): string;
  set name(value: EmbindString);
  get evidence_code(): string;
  set evidence_code(value: EmbindString);
  get details(): string;
  set details(value: EmbindString);
  add_member(_0: StructSiteMember): void;
}

export interface StructSiteVector extends ClassHandle, Iterable<StructSite> {
  push_back(_0: StructSite): void;
  resize(_0: number, _1: StructSite): void;
  size(): number;
  get(_0: number): StructSite | undefined;
  set(_0: number, _1: StructSite): boolean;
}

export interface Structure extends ClassHandle {
  cell: UnitCell;
  connections: ConnectionVector;
  sites: StructSiteVector;
  readonly length: number;
  get name(): string;
  set name(value: EmbindString);
  add_model(_0: Model): void;
  add_connection(_0: Connection): void;
  add_site(_0: StructSite): void;
  at(_0: number): Model | null;
}

export interface Model extends ClassHandle {
  num: number;
  readonly length: number;
  add_chain(_0: Chain): void;
  at(_0: number): Chain | null;
  count_occupancies(_0: Selection | null): number;
}

export interface Chain extends ClassHandle {
  readonly length: number;
  get name(): string;
  set name(value: EmbindString);
  add_residue(_0: Residue): void;
  at(_0: number): Residue | null;
}

export interface ResidueId extends ClassHandle {
  get seqid_string(): string;
  set seqid_string(value: EmbindString);
  get segment(): string;
  set segment(value: EmbindString);
  get name(): string;
  set name(value: EmbindString);
  set_seqid(_0: number, _1: EmbindString): void;
  set_seqid_string(_0: EmbindString): void;
}

export interface Residue extends ResidueId {
  ss_from_file: ResidueSs;
  strand_sense_from_file: ResidueStrandSense;
  readonly length: number;
  get subchain(): string;
  set subchain(value: EmbindString);
  readonly ss_from_file_string: string;
  readonly strand_sense_from_file_string: string;
  readonly entity_type_string: string;
  add_atom(_0: Atom): void;
  at(_0: number): Atom | null;
}

export interface Atom extends ClassHandle {
  readonly is_metal: boolean;
  altloc: number;
  charge: number;
  serial: number;
  occ: number;
  b_iso: number;
  pos: Position;
  get name(): string;
  set name(value: EmbindString);
  get element_uname(): string;
  set element_uname(value: EmbindString);
  set_element(_0: EmbindString): void;
}

export interface Selection extends ClassHandle {
  remove_selected(_0: Structure): void;
  remove_not_selected(_0: Structure): void;
}

export interface BondInfo extends ClassHandle {
  get_bond_lines(_0: Structure): void;
  bond_data_ptr(): number;
  bond_data_size(): number;
  add_monomer_cif(_0: EmbindString): void;
}

export interface CrossSymBonds extends ClassHandle {
  find(_0: Structure, _1: NearestImage): void;
  bond_data_ptr(): number;
  bond_data_size(): number;
}

export interface SelectionResult extends ClassHandle {
  atom_data_ptr(): number;
  atom_data_size(): number;
  set_atom_indices(_0: Structure, _1: EmbindString, _2: number): void;
}

export interface BlobSearchResult extends ClassHandle {
  size(): number;
  centroids(): any;
  peak_positions(): any;
  scores(): any;
  volumes(): any;
  peak_values(): any;
}

export interface MapData extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly mean: number;
  readonly rms: number;
  readonly last_error: string;
  find_blobs(_0: number, _1: number, _2: number, _3: number, _4: boolean, _5: Structure | null, _6: number, _7: number, _8: boolean): BlobSearchResult | null;
  extract_isosurface(_0: number, _1: number, _2: number, _3: number, _4: number, _5: EmbindString): boolean;
  data(): any;
  isosurface_vertices(): any;
  isosurface_triangles(): any;
}

export interface Ccp4Map extends MapData {
  read(_0: boolean): boolean;
}

export interface Dsn6Map extends MapData {
  read(): boolean;
}

export interface MtzMap extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly mean: number;
  readonly rms: number;
  readonly last_error: string;
  find_blobs(_0: number, _1: number, _2: number, _3: number, _4: boolean, _5: Structure | null, _6: number, _7: number, _8: boolean): BlobSearchResult | null;
  extract_isosurface(_0: number, _1: number, _2: number, _3: number, _4: number, _5: EmbindString): boolean;
  data(): any;
  isosurface_vertices(): any;
  isosurface_triangles(): any;
}

export interface Mtz extends ClassHandle {
  readonly cell: UnitCell;
  readonly nx: number;
  readonly ny: number;
  readonly nz: number;
  readonly rmsd: number;
  readonly last_error: string;
  read(): boolean;
  calculate_wasm_map(_0: boolean): MtzMap | null;
  calculate_wasm_map_from_labels(_0: EmbindString, _1: EmbindString): MtzMap | null;
  calculate_map(_0: boolean): any;
  calculate_map_from_labels(_0: EmbindString, _1: EmbindString): any;
}

export type Fractional = [ number, number, number ];

export type Position = [ number, number, number ];

interface EmbindModule {
  NearestImage: {
    new(): NearestImage;
  };
  NearestImageVector: {
    new(): NearestImageVector;
  };
  UnitCellParameters: {};
  UnitCell: {
    new(_0: number, _1: number, _2: number, _3: number, _4: number, _5: number): UnitCell;
  };
  Isosurface: {
    new(): Isosurface;
  };
  ResidueSs: {Coil: ResidueSsValue<0>, Helix: ResidueSsValue<1>, Strand: ResidueSsValue<2>};
  ResidueStrandSense: {NotStrand: ResidueStrandSenseValue<0>, Parallel: ResidueStrandSenseValue<1>, First: ResidueStrandSenseValue<2>, Antiparallel: ResidueStrandSenseValue<-1>};
  Asu: {Same: AsuValue<0>, Different: AsuValue<1>, Any: AsuValue<2>};
  ConnectionType: {Covale: ConnectionTypeValue<0>, Disulf: ConnectionTypeValue<1>, Hydrog: ConnectionTypeValue<2>, MetalC: ConnectionTypeValue<3>, Unknown: ConnectionTypeValue<4>};
  AtomAddress: {
    new(): AtomAddress;
  };
  Connection: {
    new(): Connection;
  };
  ConnectionVector: {
    new(): ConnectionVector;
  };
  StructSiteMember: {
    new(): StructSiteMember;
  };
  StructSiteMemberVector: {
    new(): StructSiteMemberVector;
  };
  StructSite: {
    new(): StructSite;
    new(_0: EmbindString): StructSite;
  };
  StructSiteVector: {
    new(): StructSiteVector;
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
  ResidueId: {
    new(): ResidueId;
  };
  Residue: {
    new(): Residue;
  };
  Atom: {
    new(): Atom;
  };
  Selection: {
    new(): Selection;
    new(_0: EmbindString): Selection;
  };
  BondInfo: {
    new(): BondInfo;
  };
  CrossSymBonds: {
    new(): CrossSymBonds;
  };
  SelectionResult: {
    new(): SelectionResult;
  };
  get_sym_image(_0: Structure, _1: NearestImage): Structure;
  BlobSearchResult: {
    new(): BlobSearchResult;
  };
  MapData: {};
  Ccp4Map: {
    new(_0: EmbindString): Ccp4Map;
  };
  Dsn6Map: {
    new(_0: EmbindString): Dsn6Map;
  };
  MtzMap: {};
  Mtz: {
    new(_0: EmbindString): Mtz;
  };
  readMtz(_0: EmbindString): Mtz;
  readCcp4Map(_0: EmbindString, _1?: boolean): Ccp4Map;
  readDsn6Map(_0: EmbindString): Dsn6Map;
  read_structure(_0: EmbindString, _1: EmbindString, _2?: EmbindString): Structure;
  get_nearby_sym_ops(_0: Structure, _1: Position, _2: number): NearestImageVector;
  get_residue_names(_0: Structure): string;
  get_missing_monomer_names(_0: Structure): string;
  get_monomer_names_in_cif(_0: EmbindString): string;
  extract_monomer_cifs(_0: EmbindString, _1: EmbindString): string;
  make_pdb_string(_0: Structure): string;
  make_mmcif_string(_0: Structure): string;
  _read_structure(_0: EmbindString, _1: EmbindString, _2: EmbindString): Structure;
}

export type MainModule = WasmModule & typeof RuntimeExports & EmbindModule;
export type GemmiModule = MainModule;
export default function MainModuleFactory (options?: unknown): Promise<MainModule>;
