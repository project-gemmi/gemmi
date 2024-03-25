
Module['readMtz'] = function (mtz_buf/*:ArrayBuffer*/) {
  var mtz = new Module.Mtz(mtz_buf);
  if (!mtz.read()) {
    var last_error = mtz.last_error;
    mtz.delete();
    throw Error(last_error);
  }
  return mtz;
}
