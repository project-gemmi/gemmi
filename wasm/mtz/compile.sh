#!/bin/bash -eu

em++ --bind -O3 -Wall -Wextra -std=c++17 -I../../include \
    -s WASM=1 \
    -s DISABLE_EXCEPTION_CATCHING=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["writeArrayToMemory"]' \
    mtz_fft.cpp -o mtz.js \
    #--emrun
