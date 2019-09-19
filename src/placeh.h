// Copyright 2018 Global Phasing Ltd.
#pragma once

#include <gemmi/model.hpp> // for Atom
#include <gemmi/topo.hpp>  // for Topo

void place_hydrogens(const gemmi::Atom& atom, gemmi::Topo::ResInfo& ri,
                     const gemmi::Topo& topo);
