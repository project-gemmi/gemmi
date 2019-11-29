
Module['readMtz'] = function (mtz_buf/*:ArrayBuffer*/) {
  var arr = new Uint8Array(mtz_buf);
  var buffer = Module._malloc(arr.length);
  Module.writeArrayToMemory(arr, buffer);
  var mtz = new Module.Mtz;
  if (!mtz.read(buffer, arr.length)) {
    var last_error = mtz.last_error;
    mtz.delete();
    throw Error(last_error);
  }
  return mtz;
}
