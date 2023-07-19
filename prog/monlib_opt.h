// Copyright Global Phasing Ltd.
#pragma once

#include <optionparser.h>    // for option::Option
#include <gemmi/monlib.hpp>  // for MonLib

// used by prep, h, rmsz
enum MonLibOptions { Monomers=4, Libin, AfterMonLibOptions };

extern const option::Descriptor MonLibUsage[];

struct MonArguments {
  const char* monomer_dir;
  const option::Option* libin;
  int verbose;
};

// Exits w/ error if monomer dir is not specified (--monomers or $CLIBD_MON)
MonArguments get_monomer_args(const std::vector<option::Option>& options);

// Read monomer library and user-provided cif files according to args.
// st_doc is used if we read monomer blocks from coordinate mmCIF.
// Argument 'wanted' contains wanted residues names as input
// and not found residues as output.
void read_monomer_lib_and_user_files(gemmi::MonLib& monlib,
                                     std::vector<std::string>& wanted,
                                     const MonArguments& args,
                                     const gemmi::cif::Document* st_doc);

// avoid including both mapcoef.h and monlib_opt.h (they define different options as 4)
#ifdef GEMMI_OPTIONS_4
#error Conflicting headers
#endif
#define GEMMI_OPTIONS_4
