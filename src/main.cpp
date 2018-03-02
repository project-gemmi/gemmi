// Copyright 2018 Global Phasing Ltd.
// Entry point for gemmi utility built as a single program with subcommands.

#include <stdio.h>
#include <cstring>
#include "gemmi/version.hpp"

int contents_main(int argc, char** argv);
int convert_main(int argc, char** argv);
int grep_main(int argc, char** argv);
int map_main(int argc, char** argv);
int mask_main(int argc, char** argv);
int sg_main(int argc, char** argv);
int validate_main(int argc, char** argv);

typedef int (*main_type)(int argc, char** argv);

struct SubCmd {
  const char* cmd;
  main_type func;
  const char* desc;
};

#define CMD(s, desc) { #s, &s##_main, desc }
static SubCmd subcommands[] = {
  CMD(contents, "info about content of a coordinate file (pdb, mmCIF, ...)"),
  CMD(convert, "convert file (CIF - JSON, mmCIF - PDB) or modify structure"),
  CMD(grep, "search for tags in CIF file(s)"),
  CMD(map, "print info or modify a CCP4 map"),
  CMD(mask, "make mask in the CCP4 format"),
  CMD(sg, "info about space groups"),
  CMD(validate, "validate CIF 1.1 syntax"),
};

static void print_usage() {
  printf("Usage: gemmi [--version] [--help] <command> [<args>]\n\n"
         "Commands:\n");
  for (SubCmd& sub : subcommands)
    printf(" %-13s %s\n", sub.cmd, sub.desc);
}

static bool eq(const char* a, const char* b) { return std::strcmp(a, b) == 0; }

static main_type get_subcommand_function(const char* cmd) {
  for (SubCmd& sub : subcommands)
    if (eq(cmd, sub.cmd))
      return sub.func;
  return nullptr;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage();
    return 1;
  }
  if (eq(argv[1], "--version") || eq(argv[1], "-V")) {
    printf("gemmi %s\n", GEMMI_VERSION);
    return 0;
  }
  if (eq(argv[1], "--help") || eq(argv[1], "-h") || eq(argv[1], "help")) {
    if (argc == 2) {
      print_usage();
      return 0;
    }
    main_type func = get_subcommand_function(argv[2]);
    if (!func) {
      printf("'%s' is not a gemmi command. See 'gemmi --help'.\n", argv[2]);
      return 1;
    }
    char help_str[] = "--help";
    char* args[] = { argv[0], argv[2], help_str };
    return (*func)(3, args);
  }
  if (argv[1][0] == '-') {
    printf("Invalid option '%s'. See 'gemmi --help'.\n", argv[1]);
    return 1;
  }
  // call function
  main_type func = get_subcommand_function(argv[1]);
  if (!func) {
    printf("'%s' is not a gemmi command. See 'gemmi --help'.\n", argv[1]);
    return 1;
  }
  return (*func)(argc - 1, argv + 1);
}

// vim:sw=2:ts=2:et
