// Copyright 2019 Global Phasing Ltd.
//
// This program compares sequence from SEQRES and from the model.

#include <gemmi/model.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/seqalign.hpp>  // for align_sequences
#include <gemmi/labelseq.hpp>  // for align_polymer

#include <cstdio>   // for printf, fprintf, putchar
#include <cstring>  // for strlen
#define GEMMI_PROG seq
#include "options.h"

using std::printf;

enum OptionIndex { Match=4, Mismatch, GapOpen, GapExt,
                   CheckMmcif, PrintOneLetter, TextAlign };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nCompares sequence (SEQRES) with model from a PDB or mmCIF file."
    "\nFor testing, it can also compare strings with option --text-align."
    "\nPerforms global alignment with scoring matrix and affine gap penalty."
    "\nCurrently, we only use match/mismatch scoring matrix.\n" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { CheckMmcif, 0, "", "check-mmcif", Arg::None,
    "  --check-mmcif  \tCompare alignment with mmCIF _atom_site.label_seq_id" },
  { PrintOneLetter, 0, "p", "", Arg::None,
    "  -p  \tprint formatted alignment with one-letter codes." },
  { TextAlign, 0, "", "text-align", Arg::None,
    "  --text-align  \tAlign characters in two strings (for testing)." },

  { NoOp, 0, "", "", Arg::None, "\nScoring (absolute values):" },
  { Match, 0, "", "match", Arg::Int,
    "  --match=INT  \tMatch score (default: 1)." },
  { Mismatch, 0, "", "mism", Arg::Int,
    "  --mism=INT  \tMismatch penalty (default: -1)." },
  { GapOpen, 0, "", "gapo", Arg::Int,
    "  --gapo=INT  \tGap opening penalty (default: -1)." },
  { GapExt, 0, "", "gape", Arg::Int,
    "  --gape=INT  \tGap extension penalty (default: -1)." },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_text_alignment(const char* text1, const char* text2,
                                 const gemmi::AlignmentScoring& scoring,
                                 bool print_one_letter) {
  std::vector<std::string> v1(std::strlen(text1));
  for (size_t i = 0; i != v1.size(); ++i)
    v1[i] = text1[i];
  std::vector<std::string> v2(std::strlen(text2));
  for (size_t i = 0; i != v2.size(); ++i)
    v2[i] = text2[i];
  std::vector<bool> free_gapo(1, 1);
  gemmi::AlignmentResult result
      = gemmi::align_string_sequences(v1, v2, free_gapo, scoring);
  printf("Score: %d   CIGAR: %s\n", result.score, result.cigar_str().c_str());
  if (print_one_letter) {
    printf("%s\n", result.add_gaps(text1, 1).c_str());
    printf("%s\n", result.match_string.c_str());
    printf("%s\n", result.add_gaps(text2, 2).c_str());
  }
}

static void print_alignment_details(const gemmi::AlignmentResult& result,
                                    const std::string& chain_name,
                                    const gemmi::ConstResidueSpan& polymer,
                                    const gemmi::Entity& ent) {
  std::vector<bool> gaps = prepare_free_gapo(polymer, ent.polymer_type);
  auto gap = gaps.begin();
  int seq_pos = 0;
  auto model_residues = polymer.first_conformer();
  auto res = model_residues.begin();
  for (gemmi::AlignmentResult::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      std::string fmon = gemmi::Entity::first_mon(ent.full_sequence[seq_pos]);
      printf("  %s ", chain_name.c_str());
      if (op == 'I' || op == 'M') {
        seq_pos++;
        printf("%4d %-3s -", seq_pos, fmon.c_str());
      } else {
        printf("         -");
      }
      if (op == 'D' || op == 'M') {
        printf(" %-3s %4d%c",
               res->name.c_str(), *res->seqid.num, res->seqid.icode);
        if (res->label_seq)
          printf("   id:%4d %c",
                 *res->label_seq, *res->label_seq == seq_pos ? ' ' : '!');
        std::putchar(*gap++ ? '^' : ' ');
        if (op == 'D' || fmon != res->name)
          printf("    <-- BAD");
        res++;
      }
      printf("\n");
    }
  }
}

static void check_label_seq_id(const gemmi::AlignmentResult& result,
                               const gemmi::ConstResidueSpan& polymer) {
  int seq_pos = 1;
  auto residues = polymer.first_conformer();
  auto res = residues.begin();
  for (gemmi::AlignmentResult::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      if (op == 'D' || op == 'M') {
        if (*res->label_seq != seq_pos)
          printf("NOTE: %s with label_seq_id %d , expected %d.\n",
                 res->name.c_str(), *res->label_seq, seq_pos);
        res++;
      }
      if (op == 'I' || op == 'M')
        seq_pos++;
    }
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  gemmi::AlignmentScoring scoring;
  if (p.options[Match])
    scoring.match = std::atoi(p.options[Match].arg);
  if (p.options[Mismatch])
    scoring.mismatch = std::atoi(p.options[Mismatch].arg);
  if (p.options[GapOpen])
    scoring.gapo = std::atoi(p.options[GapOpen].arg);
  if (p.options[GapExt])
    scoring.gape = std::atoi(p.options[GapExt].arg);
  if (p.options[TextAlign]) {
    p.require_positional_args(2);
    print_text_alignment(p.nonOption(0), p.nonOption(1), scoring,
                         p.options[PrintOneLetter]);
    return 0;
  }
  p.require_input_files_as_args();
  bool verbose = p.options[Verbose];
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (i > 0)
        printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        printf("File: %s\n", input.c_str());
      gemmi::Structure st = gemmi::read_structure_gz(input);
      gemmi::setup_entities(st);
      if (st.models.empty())
        gemmi::fail("No atoms found. Wrong input file?");
      if (st.models.size() > 1)
        printf("Warning: using only model 1 of %zu.\n", st.models.size());
      const gemmi::Model& model = st.models[0];
      for (const gemmi::Chain& chain : model.chains) {
        gemmi::ConstResidueSpan polymer = chain.get_polymer();
        if (!polymer)
          continue;
        const gemmi::Entity* ent = st.get_entity_of(polymer);
        if (!ent)
          gemmi::fail("No sequence (SEQRES) for chain " + chain.name);
        if (gemmi::seqid_matches_seqres(polymer, *ent))
          printf("Sequence numbers are wrt the full sequence (SEQRES).\n");
        gemmi::AlignmentResult result =
            gemmi::align_sequence_to_polymer(ent->full_sequence, polymer,
                                             ent->polymer_type, scoring);
        printf("%s chain %s CIGAR: %s\n",
               st.name.c_str(), chain.name.c_str(), result.cigar_str().c_str());
        if (p.options[CheckMmcif])
          check_label_seq_id(result, polymer);
        if (p.options[PrintOneLetter]) {
          std::string text1 = gemmi::one_letter_code(ent->full_sequence);
          std::string text2 = gemmi::one_letter_code(polymer);
          printf("%s\n", result.add_gaps(text1, 1).c_str());
          printf("%s\n", result.match_string.c_str());
          printf("%s\n", result.add_gaps(text2, 2).c_str());
        }
        if (verbose)
          print_alignment_details(result, chain.name, polymer, *ent);
      }
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
