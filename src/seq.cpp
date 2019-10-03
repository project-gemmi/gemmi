// Copyright 2019 Global Phasing Ltd.
//
// This program compares sequence from SEQRES and from the model.

#include <gemmi/model.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/seqalign.hpp>  // for align_sequences
#include <gemmi/labelseq.hpp>  // for align_polymer

#define GEMMI_PROG seq
#include "options.h"
#include <stdio.h>
#include <cstring>  // for strlen

enum OptionIndex { Verbose=3, Match, Mismatch, GapOpen, GapExt,
                   CheckMmcif, TextAlign };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nCompares sequence (SEQRES) with model from a PDB or mmCIF file."
    "\nFor testing, it can also compare strings with option --text-align."
    "\nPerforms global alignment with scoring matrix and affine gap penalty."
    "\nCurrently, we only use match/mismatch scoring matrix.\n" },
  CommonUsage[Help],
  CommonUsage[Version],
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { CheckMmcif, 0, "", "check-mmcif", Arg::None,
    "  --check-mmcif  \tCompare alignment with mmCIF _atom_site.label_seq_id" },
  { TextAlign, 0, "", "text-align", Arg::None,
    "  --text-align  \tAlign characters in two strings (for testing)." },

  { NoOp, 0, "", "", Arg::None, "\nScoring (absolute values):" },
  { Match, 0, "", "match", Arg::Int,
    "  --match=INT  \tMatch score (default: 1)." },
  { Mismatch, 0, "", "mism", Arg::Int,
    "  --mism=INT  \tMismatch penalty (default: 1)." },
  { GapOpen, 0, "", "gapo", Arg::Int,
    "  --gapo=INT  \tGap opening penalty (default: 1)." },
  { GapExt, 0, "", "gape", Arg::Int,
    "  --gape=INT  \tGap extension penalty (default: 1)." },
  { 0, 0, 0, 0, 0, 0 }
};

static void print_text_alignment(const char* text1, const char* text2,
                                 const gemmi::AlignmentScoring& scoring) {
  int len1 = (int) std::strlen(text1);
  int len2 = (int) std::strlen(text2);
  std::vector<uint8_t> v1(len1);
  for (int i = 0; i != len1; ++i)
    v1[i] = text1[i] >= 32 && text1[i] < 127 ? text1[i] - 32 : 0;
  std::vector<uint8_t> v2(len2);
  for (int i = 0; i != len2; ++i)
    v2[i] = text2[i] >= 32 && text2[i] < 127 ? text2[i] - 32 : 0;
  int m = 127 - 32;
  std::vector<int8_t> score_matrix(m * m, scoring.mismatch);
  for (int i = 0; i < m; ++i)
    score_matrix[i * m + i] = scoring.match;
  std::vector<bool> free_gapo(1, 1);
  gemmi::Alignment result = gemmi::align_sequences(
      len1, v1.data(),
      len2, v2.data(),
      free_gapo, m, score_matrix.data(),
      scoring.gapo, scoring.gape);
  printf("Score: %d   CIGAR: %s\n", result.score, result.cigar_str().c_str());
  size_t pos1 = 0;
  size_t pos2 = 0;
  std::string match, out1, out2;
  for (gemmi::Alignment::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      if (op == 'I') {
        match += 'I';
        out1 += text1[pos1++];
        out2 += '-';
      } else if (op == 'D') {
        match += 'D';
        out1 += '-';
        out2 += text2[pos2++];
      } else /* if (op == 'M') */ {
        match += v1[pos1] == v2[pos2] ? '=' : 'X';
        out1 += text1[pos1++];
        out2 += text2[pos2++];
      }
    }
  }
  printf("%s\n%s\n%s\n", match.c_str(), out1.c_str(), out2.c_str());
}

static void print_alignment_details(const gemmi::Alignment& result,
                                    const std::string& chain_name,
                                    const gemmi::ConstResidueSpan& polymer,
                                    const gemmi::Entity& ent) {
  std::vector<bool> gaps = prepare_free_gapo(polymer, ent.polymer_type);
  auto gap = gaps.begin();
  int seq_pos = 0;
  auto model_residues = polymer.first_conformer();
  auto res = model_residues.begin();
  for (gemmi::Alignment::Item item : result.cigar) {
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
        putchar(*gap++ ? '^' : ' ');
        if (op == 'D' || fmon != res->name)
          printf("    <-- BAD");
        res++;
      }
      printf("\n");
    }
  }
}

static void check_label_seq_id(const gemmi::Alignment& result,
                               const gemmi::ConstResidueSpan& polymer) {
  int seq_pos = 1;
  auto residues = polymer.first_conformer();
  auto res = residues.begin();
  for (gemmi::Alignment::Item item : result.cigar) {
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
    scoring.match = std::abs(std::atoi(p.options[Match].arg));
  if (p.options[Mismatch])
    scoring.mismatch = -std::abs(std::atoi(p.options[Mismatch].arg));
  if (p.options[GapOpen])
    scoring.gapo = std::abs(std::atoi(p.options[GapOpen].arg));
  if (p.options[GapExt])
    scoring.gape = std::abs(std::atoi(p.options[GapExt].arg));
  if (p.options[TextAlign]) {
    p.require_positional_args(2);
    print_text_alignment(p.nonOption(0), p.nonOption(1), scoring);
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
        gemmi::Alignment result = align_polymer(polymer, *ent, scoring);
        printf("%s chain %s CIGAR: %s\n",
               st.name.c_str(), chain.name.c_str(), result.cigar_str().c_str());
        if (p.options[CheckMmcif])
          check_label_seq_id(result, polymer);
        if (verbose)
          print_alignment_details(result, chain.name, polymer, *ent);
      }
    }
  } catch (std::runtime_error& e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
