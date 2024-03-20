// Copyright Global Phasing Ltd.

#include "common.h"
#include <emscripten/em_asm.h>

EMSCRIPTEN_BINDINGS(Gemmi) {
  add_mol();
  add_mtz_fft();
  EM_ASM( finalize_gemmi(); );
}
