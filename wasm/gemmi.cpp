// Copyright Global Phasing Ltd.

#include "common.h"
#include <emscripten/em_asm.h>

EMSCRIPTEN_BINDINGS(Gemmi) {
  add_cell();
  add_iso();
  add_mol();
  add_map();
  add_mtz_fft();
  EM_ASM( if (typeof finalize_gemmi === 'function') finalize_gemmi(); );
}
