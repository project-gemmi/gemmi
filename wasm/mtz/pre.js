
Module['readMtz'] = function (mtz_buf/*:ArrayBuffer*/) {
  var arr = new Uint8Array(mtz_buf);
  var buffer = Module._malloc(arr.length);
  Module.writeArrayToMemory(arr, buffer);
  var mtz = new Module.Mtz;
  if (!mtz.read(buffer, arr.length)) {
    throw Error(mtz.last_error);
  }
  return mtz;
}
