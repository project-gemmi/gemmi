// Copyright 2019 Global Phasing Ltd.
//
// This program compares sequence from SEQRES and from the model.

#include <gemmi/model.hpp>
#include <gemmi/gzread.hpp>
#include <gemmi/polyheur.hpp>  // for add_entity_types, setup_entities
#include <gemmi/seqalign.hpp>  // for align_sequences

#define GEMMI_PROG seq
#include "options.h"
#include <stdio.h>
#include <cstring>  // for strlen

enum OptionIndex { Verbose=3, Match, Mismatch, GapOpen, GapExt, TextAlign };

static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nCompares sequence (SEQRES) with model from a PDB or mmCIF file."
    "\nFor testing, it can also compare strings with option --text-align."
    "\nPerforms global alignment with scoring matrix and affine gap penalty."
    "\nCurrently, we only use match/mismatch scoring matrix.\n" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "v", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { TextAlign, 0, "", "text-align", Arg::None,
    "  --text-align  \tAlign characters in two strings (for testing)." },

  { NoOp, 0, "", "", Arg::None, "\nScoring (absolute values):" },
  { Match, 0, "", "match", Arg::Int,
    "  --match=INT  \tMatch score (default: 0)." },
  { Mismatch, 0, "", "mism", Arg::Int,
    "  --mism=INT  \tMismatch penalty (default: 1)." },
  { GapOpen, 0, "", "gapo", Arg::Int,
    "  --gapo=INT  \tGap opening penalty (default: 1)." },
  { GapExt, 0, "", "gape", Arg::Int,
    "  --gape=INT  \tGap extension penalty (default: 1)." },
  { 0, 0, 0, 0, 0, 0 }
};

struct Scoring {
  int match = 0;
  int mismatch = -1;
  int gapo = 1;
  int gape = 1;
};

static void print_text_alignment(const char* text1, const char* text2,
                                 const Scoring& scoring) {
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
  gemmi::Alignment result = gemmi::align_sequences(
      len1, v1.data(),
      len2, v2.data(),
      m, score_matrix.data(),
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


static void print_sequence_alignment(const gemmi::Structure& st,
                                     const Scoring& scoring, bool verbose) {
  if (st.models.size() > 1)
    std::fprintf(stderr, "Warning: using only the first model out of %zu.\n",
                 st.models.size());
  const gemmi::Model& model = st.models.at(0);
  for (const gemmi::Chain& chain : model.chains) {
    gemmi::ConstResidueSpan polymer = chain.get_polymer();
    if (!polymer)
      continue;
    const gemmi::Entity* ent = st.get_entity_of(polymer);
    if (!ent)
      gemmi::fail("Full sequence (SEQRES) not found for chain " + chain.name);
    std::map<std::string, uint8_t> encoding;
    for (const gemmi::Residue& res : polymer)
      encoding.emplace(res.name, 0);
    for (const std::string& mon_list : ent->full_sequence)
      encoding.emplace(gemmi::Entity::first_mon(mon_list), 0);
    std::vector<std::string> all_monomers;
    size_t n_mon = encoding.size();
    all_monomers.reserve(n_mon);
    for (auto& item : encoding) {
      item.second = (uint8_t) all_monomers.size();
      all_monomers.push_back(item.first);
    }
    std::vector<uint8_t> model_seq;
    model_seq.reserve(polymer.size());
    // TODO include gaps based on distances
    for (const gemmi::Residue& res : polymer)
      model_seq.push_back(encoding.at(res.name));
    std::vector<uint8_t> full_seq;
    full_seq.reserve(ent->full_sequence.size());
    for (const std::string& mon_list : ent->full_sequence)
      full_seq.push_back(encoding.at(gemmi::Entity::first_mon(mon_list)));
    std::vector<int8_t> score_matrix(n_mon * n_mon, scoring.mismatch);
    for (size_t i = 0; i != n_mon; ++i)
      score_matrix[i * n_mon + i] = scoring.match;
    gemmi::Alignment result = gemmi::align_sequences(
        full_seq.size(), full_seq.data(),
        model_seq.size(), model_seq.data(),
        n_mon, score_matrix.data(),
        scoring.gapo, scoring.gape);
    printf("%s chain %s CIGAR: %s\n",
           st.name.c_str(), chain.name.c_str(), result.cigar_str().c_str());
    if (verbose) {
      size_t model_pos = 0;
      size_t full_pos = 0;
      for (gemmi::Alignment::Item item : result.cigar) {
        char op = item.op();
        for (uint32_t i = 0; i < item.len(); ++i) {
          bool ok = op == 'I' ||
                    (op == 'M' && full_seq[full_pos] == model_seq[model_pos]);
          printf("  %s ", chain.name.c_str());
          if (op == 'I' || op == 'M') {
            uint8_t enc = full_seq.at(full_pos++);
            printf("%4zu %-3s -", full_pos, all_monomers.at(enc).c_str());
          } else {
            printf("         -");
          }
          if (op == 'D' || op == 'M') {
            const gemmi::Residue& res = polymer.at(model_pos++);
            printf(" %-3s %4d%c",
                   res.name.c_str(), (int) res.seqid.num, res.seqid.icode);
          }
          if (!ok)
            printf("    <-- BAD");
          printf("\n");
        }
      }
    }
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  Scoring scoring;
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
        std::printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        std::printf("File: %s\n", input.c_str());
      gemmi::Structure st = gemmi::read_structure_gz(input);
      gemmi::setup_entities(st);
      print_sequence_alignment(st, scoring, verbose);
    }
  } catch (std::runtime_error& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
