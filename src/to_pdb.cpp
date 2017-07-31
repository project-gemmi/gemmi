// Copyright 2017 Global Phasing Ltd.

// This file exists to make compilation of gemmi-convert faster.

// set this before only one of stb_sprintf.h includes
#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/to_pdb.hpp"

void write_pdb(const gemmi::mol::Structure& st, std::ostream& os,
               bool iotbx_compat) {
  gemmi::mol::write_pdb(st, os, iotbx_compat);
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
