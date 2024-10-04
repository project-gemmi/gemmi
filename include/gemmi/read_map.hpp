// Copyright 2021 Global Phasing Ltd.
//
// Functions for reading possibly gzipped CCP4 map files.

#ifndef GEMMI_READ_MAP_HPP_
#define GEMMI_READ_MAP_HPP_

#include "ccp4.hpp"  // for Ccp4

namespace gemmi {

GEMMI_DLL Ccp4<float> read_ccp4_map(const std::string& path, bool setup);
GEMMI_DLL Ccp4<int8_t> read_ccp4_mask(const std::string& path, bool setup);

} // namespace gemmi

#endif
