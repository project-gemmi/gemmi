const fs = require('node:fs');
const Gemmi = require('./gemmi.js');

function nearestAtomDistance(st, pos) {
  let best = Infinity;
  for (let model of st)
    for (let chain of model)
      for (let residue of chain)
        for (let atom of residue) {
          const dx = atom.pos[0] - pos[0];
          const dy = atom.pos[1] - pos[1];
          const dz = atom.pos[2] - pos[2];
          best = Math.min(best, Math.hypot(dx, dy, dz));
        }
  return best;
}

test('finds blobs in CCP4 maps and supports model masking', async () => {
  const gemmi = await Gemmi();
  const mapBuf = fs.readFileSync('../tests/5i55_tiny.ccp4');
  const modelBuf = fs.readFileSync('../tests/5i55.cif');
  const map = gemmi.readCcp4Map(mapBuf, true);
  const st = gemmi.read_structure(modelBuf, '5i55.cif');

  const unmasked = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, null, 0, 2.0, false);
  expect(unmasked.size()).toBe(1);
  expect(unmasked.centroids().length).toBe(3);
  expect(unmasked.peak_positions().length).toBe(3);
  expect(unmasked.scores().length).toBe(1);
  expect(unmasked.volumes().length).toBe(1);
  expect(unmasked.peak_values().length).toBe(1);
  expect(unmasked.scores()[0]).toBeGreaterThan(15.0);

  const masked = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, st, 0, 2.0, false);
  expect(masked.size()).toBe(0);
  expect(masked.centroids().length).toBe(0);

  masked.delete();
  unmasked.delete();
  st.delete();
  map.delete();
});

test('finds blobs in MTZ-derived maps', async () => {
  const gemmi = await Gemmi();
  const mtzData = fs.readFileSync('../dim/final.mtz');
  const modelData = fs.readFileSync('../dim/final.pdb');
  const mtzBuf = mtzData.buffer.slice(mtzData.byteOffset, mtzData.byteOffset + mtzData.byteLength);
  const modelBuf = modelData.buffer.slice(modelData.byteOffset,
                                          modelData.byteOffset + modelData.byteLength);
  const mtz = gemmi.readMtz(mtzBuf);
  const st = gemmi.read_structure(modelBuf, 'final.pdb');
  const map = mtz.calculate_wasm_map(true);

  expect(map).not.toBeNull();
  const unmasked = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, null, 0, 2.0, true);
  expect(unmasked.size()).toBeGreaterThan(0);

  const masked = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, st, 0, 2.0, true);
  expect(masked.size()).toBeGreaterThan(0);

  masked.delete();
  unmasked.delete();
  map.delete();
  st.delete();
  mtz.delete();
});

test('moves wasm blob positions to the symmetry image nearest the model', async () => {
  const gemmi = await Gemmi();
  const mapBuf = fs.readFileSync('../tests/5i55_tiny.ccp4');
  const modelBuf = fs.readFileSync('../tests/5i55.cif');
  const map = gemmi.readCcp4Map(mapBuf, true);
  const st = gemmi.read_structure(modelBuf, '5i55.cif');

  const raw = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, null, 0, 0.0, false);
  expect(raw.size()).toBe(1);
  const rawCentroid = Array.from(raw.centroids());

  const shift = st.cell.orthogonalize([1, 0, 0]);
  for (let model of st)
    for (let chain of model)
      for (let residue of chain)
        for (let atom of residue)
          atom.pos = [atom.pos[0] + shift[0],
                      atom.pos[1] + shift[1],
                      atom.pos[2] + shift[2]];

  const moved = map.find_blobs(map.rms, 10.0, 15.0, 0.0, false, st, 0, 0.0, false);
  expect(moved.size()).toBe(1);
  const movedCentroid = Array.from(moved.centroids());

  expect(nearestAtomDistance(st, rawCentroid)).toBeGreaterThan(20.0);
  expect(nearestAtomDistance(st, movedCentroid)).toBeLessThan(10.0);
  expect(Math.abs(movedCentroid[0] - rawCentroid[0])).toBeGreaterThan(20.0);

  moved.delete();
  raw.delete();
  st.delete();
  map.delete();
});
