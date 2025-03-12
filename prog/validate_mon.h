// Copyright Global Phasing Ltd.
#pragma once

#include "gemmi/cifdoc.hpp"

void check_monomer(const gemmi::cif::Block& block, double z_score);
void compare_monomer_with_ccd(const gemmi::cif::Block& lib_block,
                              const gemmi::cif::Block& ccd_block,
                              bool verbose);
