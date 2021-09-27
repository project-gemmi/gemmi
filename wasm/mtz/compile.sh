#!/bin/bash -eu

# last time compiled with Emscripten 1.39.3

em++ --version | head -1

em++ --bind -O3 -Wall -Wextra -std=c++17 -I../../include \
    -DEMSCRIPTEN_HAS_UNBOUND_TYPE_NAMES=0 -fno-rtti \
    -s WASM=1 \
    -s STRICT=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME=GemmiMtz \
    -s DISABLE_EXCEPTION_CATCHING=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s 'EXPORTED_RUNTIME_METHODS=["writeArrayToMemory"]' \
    --pre-js pre.js \
    --no-entry \
    mtz_fft.cpp -o mtz.js \
    #--emrun
