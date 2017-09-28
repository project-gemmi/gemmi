// Copyright 2017 Global Phasing Ltd.

#include "gemmi/symmetry.hpp"
#include <cstdio>

namespace sym = gemmi::sym;

void process_arg(const char* arg) {
    const sym::SpaceGroup* sg = sym::find_spacegroup_by_name(arg);
    if (sg == nullptr) {
        std::printf("H-M name not found: %s\n", arg);
    } else {
        std::printf("SG #%d: %s", sg->number, sg->hm);
        if (sg->ext)
            std::printf(":%c", sg->ext);
        std::printf("\n");
    }
}

int main(int argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        process_arg(argv[i]);
    return 0;
}
