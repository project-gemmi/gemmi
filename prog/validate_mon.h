// Copyright Global Phasing Ltd.
#pragma once

#include <map>
#include "gemmi/cifdoc.hpp"

void check_monomer_doc(const gemmi::cif::Document& doc, bool normal_checks, double z_score,
                       const std::map<std::string, gemmi::cif::Block>& ccd_map, bool verbose);
