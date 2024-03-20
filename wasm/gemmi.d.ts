
export interface UnitCell {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
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
  readMtz(mtz_buf: string|ArrayBuffer): Mtz;
}
