const fs = require('node:fs');
const Gemmi = require('./gemmi.js');

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
