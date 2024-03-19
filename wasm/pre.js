
Module['read_structure'] = function (buf, name, format) {
  return Module._read_structure(buf, name, format || 'unknown');
};
