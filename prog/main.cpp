// Copyright 2018 Global Phasing Ltd.
// Entry point for gemmi utility built as a single program with subcommands.

#include <stdio.h>
#include <cstring>

void print_version(const char* program_name, bool verbose);  // in options.h

int align_main(int argc, char** argv);
int blobs_main(int argc, char** argv);
int cif2json_main(int argc, char** argv);
int cif2mtz_main(int argc, char** argv);
int cifdiff_main(int argc, char** argv);
int contact_main(int argc, char** argv);
int contents_main(int argc, char** argv);
int convert_main(int argc, char** argv);
int crd_main(int argc, char** argv);
int ecalc_main(int argc, char** argv);
int fprime_main(int argc, char** argv);
int grep_main(int argc, char** argv);
int h_main(int argc, char** argv);
int json2cif_main(int argc, char** argv);
int map2sf_main(int argc, char** argv);
int map_main(int argc, char** argv);
int mask_main(int argc, char** argv);
int merge_main(int argc, char** argv);
int mondiff_main(int argc, char** argv);
int mtz2cif_main(int argc, char** argv);
int mtz_main(int argc, char** argv);
int reindex_main(int argc, char** argv);
int residues_main(int argc, char** argv);
int rmsz_main(int argc, char** argv);
int set_main(int argc, char** argv);
int sf2map_main(int argc, char** argv);
int sfcalc_main(int argc, char** argv);
int sg_main(int argc, char** argv);
int tags_main(int argc, char** argv);
int validate_main(int argc, char** argv);
int wcn_main(int argc, char** argv);
int xds2mtz_main(int argc, char** argv);

namespace {

typedef int (*main_type)(int argc, char** argv);

struct SubCmd {
  const char* cmd;
  main_type func;
  const char* desc;
};

#define CMD(s, desc) { #s, &s##_main, desc }
SubCmd subcommands[] = {
  CMD(align, "sequence alignment (global, pairwise, affine gap penalty)"),
  CMD(blobs, "list unmodelled electron density blobs"),
  CMD(cif2json, "translate (mm)CIF to (mm)JSON"),
  CMD(cif2mtz, "convert structure factor mmCIF to MTZ"),
  CMD(cifdiff, "compare tags in two (mm)CIF files"),
  CMD(contact, "searches for contacts (neighbouring atoms)"),
  CMD(contents, "info about content of a coordinate file (pdb, mmCIF, ...)"),
  CMD(convert, "convert file (CIF - JSON, mmCIF - PDB) or modify structure"),
  CMD(crd, "prepare topology file (.crd) for Refmac"),
  CMD(ecalc, "calculate normalized amplitudes E"),
  CMD(fprime, "calculate anomalous scattering factors f' and f\""),
  CMD(grep, "search for tags in CIF file(s)"),
  CMD(h, "add or remove hydrogen atoms"),
  CMD(json2cif, "translate mmJSON to mmCIF"),
  CMD(map, "print info or modify a CCP4 map"),
  CMD(map2sf, "transform CCP4 map to map coefficients (in MTZ or mmCIF)"),
  CMD(mask, "make a bulk-solvent mask in the CCP4 format"),
  CMD(merge, "merge or compare intensities from reflection file(s)"),
  CMD(mondiff, "compare two monomer CIF files"),
  CMD(mtz, "print info about MTZ reflection file"),
  CMD(mtz2cif, "convert MTZ to structure factor mmCIF"),
  CMD(reindex, "reindex MTZ file"),
  CMD(residues, "list residues from a coordinate file"),
  CMD(rmsz, "validate geometry using monomer library"),
  CMD(set, "modify coordinate files (think pdbset)"),
  CMD(sf2map, "transform map coefficients (from MTZ or mmCIF) to map"),
  CMD(sfcalc, "calculate structure factors from a model"),
  CMD(sg, "info about space groups"),
  CMD(tags, "list tags from CIF file(s)"),
  CMD(validate, "validate CIF 1.1 syntax"),
  CMD(wcn, "calculate local density / contact numbers (WCN, CN, ACN, LDM)"),
  CMD(xds2mtz, "convert XDS_ASCII to MTZ"),
};

void print_usage() {
  print_version("gemmi", /*verbose=*/false);
  printf("Command-line utility that accompanies the GEMMI library,\n"
         "which is a joint project of CCP4 and Global Phasing Ltd.\n"
         "Licence: Mozilla Public License 2.0. Copyright Global Phasing Ltd.\n"
         "https://github.com/project-gemmi/gemmi\n\n"
         "Usage: gemmi [--version] [--help] <command> [<args>]\n\n"
         "Commands:\n");
  for (SubCmd& sub : subcommands)
    printf(" %-13s %s\n", sub.cmd, sub.desc);
}

bool eq(const char* a, const char* b) { return std::strcmp(a, b) == 0; }

main_type get_subcommand_function(const char* cmd) {
  for (SubCmd& sub : subcommands)
    if (eq(cmd, sub.cmd))
      return sub.func;
  return nullptr;
}

} // anonymous namespace

#if defined(_WIN32) && defined(_UNICODE)
#include <vector>
#include "gemmi/utf.hpp"

extern "C"
int wmain(int argc, wchar_t** argv_)
#else
int main(int argc, char** argv)
#endif
{
  if (argc < 2) {
    print_usage();
    return 1;
  }
#if defined(_WIN32) && defined(_UNICODE)
  std::vector<std::string> utf8_args(argc);
  std::vector<char*> argv(argc);
  for (int i = 0; i < argc; ++i) {
    utf8_args[i] = gemmi::wchar_to_UTF8(argv_[i]);
    argv[i] = &utf8_args[i][0];
  }
#endif
  bool verbose = false;
  bool version = false;
  bool help = false;
  int command = 0;
  int wrong_option = 0;
  for (int i = 1; i < argc && wrong_option == 0; ++i) {
    const char* arg = argv[i];
    if (arg[0] == '-' && arg[1] == '-') {          // long options
      if (eq(arg+2, "version"))
        version = true;
      else if (eq(arg+2, "help"))
        help = true;
      else if (eq(arg+2, "verbose"))
        verbose = true;
      else
        wrong_option = i;
    } else if (arg[0] == '-' && arg[1] != '-') {   // short options
      for (int j = 1; arg[j] != '\0'; ++j) {
        if (arg[j] == 'V')
          version = true;
        else if (arg[j] == 'h')
          help = true;
        else if (arg[j] == 'v')
          verbose = true;
        else
          wrong_option = i;
      }
    } else {                                       // not options
      if (eq(arg, "help")) {
        help = true;
      } else if (eq(arg, "version")) {
        version = true;
      } else {
        command = i;
        break;
      }
    }
  }
  if (wrong_option != 0) {
    printf("Invalid option '%s'. See 'gemmi --help'.\n", argv[wrong_option]);
    return 1;
  }
  if (version) {
    print_version("gemmi", verbose);
    return 0;
  }
  if (command == 0) {  // handles both "gemmi" and "gemmi --help"
    print_usage();
    return 0;
  }
  // call function
  main_type func = get_subcommand_function(argv[command]);
  if (!func) {
    printf("'%s' is not a gemmi command. See 'gemmi --help'.\n", argv[command]);
    return 1;
  }
  if (verbose)
    printf("Note: Option -v/--verbose before subcommand has no effect.\n");
  if (help) {  // if we are here, command != 0
    char help_str[] = "--help";
    char* args[] = { argv[0], argv[command], help_str };
    return (*func)(3, args);
  }
  return (*func)(argc - command, &argv[command]);
}
