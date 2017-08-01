// Copyright 2017 Global Phasing Ltd.

// This file exists only to make compilation faster.
#pragma once
#include "gemmi/cifdoc.hpp" // for Document
#include "gemmi/model.hpp"  // for Structure

gemmi::cif::Document cif_read_any(const std::string& path);

gemmi::mol::Structure mmcif_read_atoms(const gemmi::cif::Document& doc);

gemmi::mol::Structure pdb_read_any(const std::string& path);

// vim:sw=2:ts=2:et:path^=../include,../third_party
