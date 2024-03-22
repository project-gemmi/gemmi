// Copyright Global Phasing Ltd.
#pragma once

#include <emscripten/bind.h>
namespace em = emscripten;

void add_cell();    // cell.cpp
void add_mol();     // mol.cpp
void add_mtz_fft(); // mtz_fft.cpp
