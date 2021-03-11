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

enum OptionIndex { FormatIn=3, Match, Label };

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
  { NoOp, 0, "", "", Arg::None,
    "INPUT is a coordinate file (mmCIF, PDB, etc)."
    "\nThe selection SEL has MMDB syntax:"
    "\n/mdl/chn/s1.i1(res)-s2.i2/at[el]:aloc (all fields are optional)\n" },
  { 0, 0, 0, 0, 0, 0 }
};

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  const char* cid = p.options[Match] ? p.options[Match].arg : "*";
  gemmi::CoorFormat format = coor_format_as_enum(p.options[FormatIn]);
  try {
    gemmi::Selection sel = gemmi::parse_cid(cid);
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      gemmi::Structure st = gemmi::read_structure_gz(input, format);
      if (p.options[Label] && st.input_format == gemmi::CoorFormat::Pdb) {
        gemmi::setup_entities(st);
        // hidden feature: -ll generates label_seq even if SEQRES is missing
        bool force = p.options[Label].count() > 1;
        gemmi::assign_label_seq_id(st, force);
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
                printf("%s (%-3s %4s%c (%-4s %s:",
                       chain.name.c_str(), (res.subchain + ")").c_str(),
                       res.seqid.num.str().c_str(), res.seqid.icode,
                       (res.label_seq.str('.') + ")").c_str(),
                       res.name.c_str());
              else
                printf("%s %4s%c %s:",
                       chain.name.c_str(),
                       res.seqid.num.str().c_str(), res.seqid.icode,
                       res.name.c_str());
              for (auto at = begin; at != end; ++at)
                printf(" %s", at->name.c_str());
              printf("\n");
              line_count++;
            }
          }
          if (line_count != 0)
            printf("\n");
        }
      }
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return 0;
}
// vim:sw=2:ts=2:et
