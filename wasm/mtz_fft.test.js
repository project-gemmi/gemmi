const fs = require('node:fs');
const zlib = require('node:zlib');
const Gemmi = require('./gemmi.js')


test('counts atom occupancies', () => {
  Gemmi().then((gemmi) => {
    const path = '../tests/5wkd_phases.mtz.gz';
    const buffer_gz = fs.readFileSync(path);
    const buffer = zlib.gunzipSync(new Buffer.from(buffer_gz));
    const mtz = gemmi.readMtz(buffer);
    expect(mtz.cell.a).toBe(50.347);
    expect(mtz.cell.gamma).toBe(90.);
    mtz.delete();
  });
});
