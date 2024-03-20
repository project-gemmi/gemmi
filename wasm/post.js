
Module['read_structure'] = function (buf, name, format) {
  return Module._read_structure(buf, name, format || 'unknown');
};

function gemmi_iterate() {
  let i = 0;
  return {
    next: () => {
      if (i >= this.length) return { done: true };
      return { done: false, value: this.at(i++) };
    }
  }
}

function finalize_gemmi() {
  for (let c of ['Structure', 'Model', 'Chain', 'Residue']) {
    Module[c]['prototype'][Symbol.iterator] = gemmi_iterate;
  }
}
