// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.
#pragma once
#include "gemmi/cifdoc.hpp" // for Document
#include "gemmi/model.hpp"  // for Structure
#include <cstdlib>          // for exit
#include <cstdio>           // for fprintf

gemmi::cif::Document cif_read_any(const std::string& path);

gemmi::Structure mmcif_read_atoms(const gemmi::cif::Document& doc);

gemmi::Structure read_structure(const std::string& path,
                    gemmi::CoorFormat format=gemmi::CoorFormat::Unknown);

gemmi::CoorFormat coordinate_format_from_extension(const std::string& path);

// vim:sw=2:ts=2:et:path^=../include,../third_party
