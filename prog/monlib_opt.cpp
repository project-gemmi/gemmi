// Copyright Global Phasing Ltd.
//
// Options for reading a monomer library + user-provided CIF files.

#include "monlib_opt.h"
#include <cstdlib>  // for getenv, exit
#include <gemmi/read_cif.hpp>  // for read_cif_gz

#define GEMMI_PROG n/a
#include "options.h"

const option::Descriptor MonLibUsage[] = {
  { Monomers, 0, "", "monomers", Arg::Required,
    "  --monomers=DIR  \tMonomer library directory (default: $CLIBD_MON)." },
  { Libin, 0, "L", "lib", Arg::Required,
    "  -L CIF, --lib=CIF  \tUser's restraint file(s). See more info below." },
  { NoOp, 0, "", "", Arg::None,
    "\nOption -L/--lib can be used multiple times, in order of priority."
    "\nIts argument is either a file path or one of the two special values:"
    "\n    '+' = monomer blocks in mmCIF INPUT_FILE (ignored by default)"
    "\n    '@' = the priority of the monomer library (ML, default: lowest)"
    "\nExample 1:   -L file.cif -L+    order: file.cif, input file, ML"
    "\nExample 2:   -L@ -L file.cif    order: ML, file.cif" }
};

MonArguments get_monomer_args(const std::vector<option::Option>& options) {
  MonArguments args;
  const option::Option* mon = options[Monomers];
  args.monomer_dir = mon ? mon->arg : std::getenv("CLIBD_MON");
  if (args.monomer_dir == nullptr || args.monomer_dir[0] == '\0') {
    fprintf(stderr, "Set $CLIBD_MON or use option --monomers.\n");
    std::exit(1);
  }
  args.libin = options[Libin];
  args.verbose = options[Verbose].count();
  return args;
}

void read_monomer_lib_and_user_files(gemmi::MonLib& monlib,
                                     std::vector<std::string>& wanted,
                                     const MonArguments& args,
                                     const gemmi::cif::Document* st_doc) {
  auto read_user_file = [&](const char* path) {
    if (path[0] == '+' && path[1] == '\0') {
      if (args.verbose)
        fprintf(stderr, "Checking coordinate file for restraint library...\n");
      if (st_doc && !st_doc->blocks.empty())
        monlib.read_monomer_doc(*st_doc);
      else
        fprintf(stderr, "Input file is not mmCIF, ignoring option --lib=+\n");
    } else {
      if (args.verbose)
        fprintf(stderr, "Reading user's library %s ...\n", path);
      monlib.read_monomer_cif(path, gemmi::read_cif_gz);
    }
  };

  auto libin = args.libin;
  for (; libin; libin = libin->next()) {
    if (libin->arg[0] == '@' && libin->arg[1] == '\0')
      break;
    read_user_file(libin->arg);
  }
  if (!monlib.monomers.empty()) {
    if (args.verbose >= 2) {
      fprintf(stderr, "Monomers read first:");
      for (auto& mpair : monlib.monomers)
        fprintf(stderr, " %s", mpair.first.c_str());
      putc('\n', stderr);
    }
  }
  if (args.verbose)
    fprintf(stderr, "Reading monomer library...\n");
  std::string error;
  monlib.read_monomer_lib(args.monomer_dir, wanted, gemmi::read_cif_gz, &error);
  auto is_found = [&](const std::string& s) { return monlib.monomers.count(s); };
  gemmi::vector_remove_if(wanted, is_found);
  if (!error.empty())
    fprintf(stderr, "%s", error.c_str());
  if (libin) {
    for (libin = libin->next(); libin; libin = libin->next()) {
      if (libin->arg[0] == '@' && libin->arg[1] == '\0')
        gemmi::fail("option --lib=/ can be given only once");
      read_user_file(libin->arg);
    }
    gemmi::vector_remove_if(wanted, is_found);
  }
}
