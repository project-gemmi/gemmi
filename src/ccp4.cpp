// Copyright 2021 Global Phasing Ltd.

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped

namespace gemmi {

Ccp4<float> read_ccp4_map(const std::string& path, bool setup) {
  Ccp4<float> ccp4;
  ccp4.read_ccp4(MaybeGzipped(path));
  if (setup)
    ccp4.setup(NAN);
  return ccp4;
}

Ccp4<int8_t> read_ccp4_mask(const std::string& path, bool setup) {
  Ccp4<int8_t> ccp4;
  ccp4.read_ccp4(MaybeGzipped(path));
  if (setup)
    ccp4.setup(-1);
  return ccp4;
}

Ccp4Base read_ccp4_header(const std::string& path) {
  Ccp4Base ccp4;
  ccp4.read_ccp4_header_(nullptr, *MaybeGzipped(path).create_stream(), path);
  return ccp4;
}

} // namespace gemmi
