// Copyright 2018 Global Phasing Ltd.

#include <stdio.h>
#include <string>
#include "gemmi/select.hpp"
#include "gemmi/polyheur.hpp"  // for setup_entities
#include "gemmi/align.hpp"     // for assign_label_seq_id
#include "gemmi/read_coor.hpp" // for read_structure_gz

#define GEMMI_PROG residues
#include "options.h"

namespace {

enum OptionIndex { FormatIn=3, Match, Label, CheckSeqId, NoAlt };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nPrints one residue per line, with atom names." },
  CommonUsage[Help],
  CommonUsage[Version],
  { FormatIn, 0, "", "format", Arg::CoorFormat,
    "  --format=FORMAT  \tInput format (default: from the file extension)." },
  { Match, 0, "m", "match", Arg::Required,
    "  -mSEL, --match=SEL  \tPrint residues/atoms matching the selection." },
  { Label, 0, "l", "label", Arg::None,
    "  -l, --label  \tPrint 'label' chain ID and seq ID in brackets." },
  { CheckSeqId, 0, "", "check-seqid", Arg::None,
    "  --check-seqid  \tCheck if sequence IDs are unique and exit." },
  { NoAlt, 0, "", "no-alt", Arg::None,
    "  --no-alt  \tDo not print altlocs." },
  { NoOp, 0, "", "", Arg::None,
    "INPUT is a coordinate file (mmCIF, PDB, etc)."
    "\nThe selection SEL has MMDB syntax:"
    "\n/mdl/chn/s1.i1(res)-s2.i2/at[el]:aloc (all fields are optional)\n" },
  { 0, 0, 0, 0, 0, 0 }
};

static char get_primary_altloc(const gemmi::Residue& res) {
  for (const gemmi::Atom& atom : res.atoms)
    if (atom.altloc == '\0')
      return ' ';
  return res.atoms.at(0).altloc;
}

static bool check_if_atoms_are_unique(const gemmi::Residue& res) {
  // Unoptimized - O(N^2).
  for (const gemmi::Atom& atom : res.atoms)
    for (const gemmi::Atom* a = res.atoms.data(); a < &atom; ++a)
      if (a->name == atom.name && a->altloc == atom.altloc)
        return false;
  return true;
}

static bool check_sequence_id(const gemmi::Structure& st) {
  bool error = false;
  std::string model_num;
  for (const gemmi::Model& model : st.models) {
    if (st.models.size() > 1)
      model_num = " (model " + model.name + ")";
    for (const gemmi::Chain& chain : model.chains) {
      const gemmi::Residue* prev_res = nullptr;
      for (const gemmi::Residue& res : chain.residues) {
        if (prev_res) {
          if (prev_res->seqid == res.seqid) {
            printf("Microheterogeneity in %s%s %s %s:  ",
                   st.name.c_str(), model_num.c_str(), chain.name.c_str(),
                   res.seqid.str().c_str());
            char alt = get_primary_altloc(res);
            bool no_alt = (alt == ' ');
            bool dup_alt = false;
            for (const gemmi::Residue* r = prev_res; r < &res; ++r) {
              char alt2 = get_primary_altloc(*r);
              if (alt2 == ' ')
                no_alt = true;
              else if (alt2 == alt)
                dup_alt = true;
              printf("%s (%c) = ", r->name.c_str(), alt2);
            }
            printf("%s (%c)\n", res.name.c_str(), alt);
            if (no_alt || dup_alt) {
              error = true;
              printf(" ERROR: microheterogeneity %s altloc\n",
                     no_alt ? "without" : "with duplicated");
            }
          } else if (res.seqid < prev_res->seqid) {
            printf("Unordered sequence ID in %s%s %s: %s > %s\n",
                   st.name.c_str(), model_num.c_str(), chain.name.c_str(),
                   prev_res->seqid.str().c_str(), res.seqid.str().c_str());
            for (const gemmi::Residue* r = chain.residues.data(); r < prev_res; ++r)
              if (r->seqid == res.seqid) {
                error = true;
                printf("ERROR: duplicated sequence ID in %s%s %s: %s\n",
                       st.name.c_str(), model_num.c_str(), chain.name.c_str(),
                       res.seqid.str().c_str());
                break;
              }
          }
        }
        // Check if atoms in the residue (name + altloc) are unique.
        // If not, it may be actually a problem with seqid (two residues
        // had the same seqid and got read as one residue).
        if (!check_if_atoms_are_unique(res)) {
          error = true;
          printf("ERROR: duplicated atoms in %s%s %s %s\n",
                 st.name.c_str(), model_num.c_str(), chain.name.c_str(),
                 res.seqid.str().c_str());
        }
        prev_res = &res;
      }
    }
  }
  return !error;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  const char* cid = p.options[Match] ? p.options[Match].arg : "*";
  gemmi::CoorFormat format = coor_format_as_enum(p.options[FormatIn]);
  bool print_alt = !p.options[NoAlt];
  int status = 0;
  try {
    gemmi::Selection sel(cid);
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      gemmi::Structure st = gemmi::read_structure_gz(input, format);
      if (p.options[Label] && st.input_format == gemmi::CoorFormat::Pdb) {
        gemmi::setup_entities(st);
        // hidden feature: -ll generates label_seq even if SEQRES is missing
        bool force = p.options[Label].count() > 1;
        gemmi::assign_label_seq_id(st, force);
      }
      if (p.options[CheckSeqId]) {
        bool ok = check_sequence_id(st);
        if (!ok)
          ++status;
        continue;
      }
      if (i != 0)
        putchar('\n');
      printf("%s\n", input.c_str());
      for (gemmi::Model& model : sel.models(st)) {
        if (st.models.size() != 1)
          printf("Model %s\n", model.name.c_str());
        for (gemmi::Chain& chain : sel.chains(model)) {
          int line_count = 0;
          for (gemmi::Residue& res : sel.residues(chain)) {
            auto sel_atoms = sel.atoms(res);
            auto begin = sel_atoms.begin();
            auto end = sel_atoms.end();
            if (begin != end) {
              if (p.options[Label])
                printf("%s (%-3s %4s%c (%-4s %s ",
                       chain.name.c_str(), (res.subchain + ")").c_str(),
                       res.seqid.num.str().c_str(), res.seqid.icode,
                       (res.label_seq.str('.') + ")").c_str(),
                       res.name.c_str());
              else
                printf("%s %4s%c %s ",
                       chain.name.c_str(),
                       res.seqid.num.str().c_str(), res.seqid.icode,
                       res.name.c_str());
              const std::string* prev = nullptr;
              for (auto at = begin; at != end; ++at)
                if (!prev || *prev != at->name) {
                  printf(" %s", at->name.c_str());
                  if (print_alt && at->altloc)
                    printf(":%c", at->altloc);
                  prev = &at->name;
                } else {
                  if (print_alt) {
                    putchar(',');
                    if (at->altloc)
                      putchar(at->altloc);
                  }
                }
              putchar('\n');
              line_count++;
            }
          }
          if (line_count != 0)
            putchar('\n');
        }
      }
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return status;
}
