// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.
#pragma once
#include "gemmi/model.hpp"  // for Structure
#include "gemmi/cifdoc.hpp"  // for Block

void write_pdb(const gemmi::mol::Structure& st, std::ostream& os,
               bool iotbx_compat);

void update_cif_block(const gemmi::mol::Structure& st,
                      gemmi::cif::Block& block);

// vim:sw=2:ts=2:et:path^=../include,../third_party
