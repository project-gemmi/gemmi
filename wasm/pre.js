
Module['readMtz'] = function (mtz_buf/*:ArrayBuffer*/) {
  var mtz = new Module.Mtz(mtz_buf);
  if (!mtz.read()) {
    var last_error = mtz.last_error;
    mtz.delete();
    throw Error(last_error);
  }
  return mtz;
}

Module['readCcp4Map'] = function (map_buf/*:ArrayBuffer*/, expand_symmetry/*:?boolean*/) {
  var map = new Module.Ccp4Map(map_buf);
  if (!map.read(expand_symmetry !== false)) {
    var last_error = map.last_error;
    map.delete();
    throw Error(last_error);
  }
  return map;
}
