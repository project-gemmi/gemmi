// Copyright 2019 Global Phasing Ltd.
//
// This program compares sequence from SEQRES and from the model.

#include <gemmi/model.hpp>
#include <gemmi/polyheur.hpp>  // for setup_entities
#include <gemmi/seqtools.hpp>  // for one_letter_code
#include <gemmi/align.hpp>     // for align_sequence_to_polymer
#include <gemmi/seqalign.hpp>  // for align_string_sequences
#include <gemmi/mmread_gz.hpp>

#include <cstdio>   // for printf, fprintf, putchar
#include <cstdlib>  // for atoi
#define GEMMI_PROG align
#include "options.h"

namespace {

using std::printf;

enum OptionIndex { Match=4, Mismatch, GapOpen, GapExt, Blosum62, Partial,
                   CheckMmcif, PrintOneLetter, Query, Target, TextAlign, Rmsd };

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Pairwise sequence alignment with scoring matrix and affine gap penalty."
    "\n\nUsage:"
    "\n\n" EXE_NAME " [options] FILE[...]"
    "\n    Aligns sequence from the model to the full sequence (SEQRES)."
    "\n    Both are from the same FILE - either in the PDB or mmCIF format."
    "\n    If the mmCIF format is used, option --check-mmcif can be used."
    "\n\n" EXE_NAME " [options] --query=CHAIN1 --target=CHAIN2 FILE [FILE2]"
    "\n    Aligns CHAIN1 from FILE to CHAIN2 from FILE2 (if given) or FILE."
    "\n    By default, the sequence of residues in the model is used."
    "\n    To use SEQRES prepend '+' to the chain name (e.g. --query=+A),"
    "\n    or, when using mmCIF, prepend '@' to entity name (--query=@1)."
    "\n\n" EXE_NAME " [options] --text-align STRING1 STRING2"
    "\n    Aligns two ASCII strings (used for testing)."
    "\n\nOptions:" },
  CommonUsage[Help],
  CommonUsage[Version],

  { CheckMmcif, 0, "", "check-mmcif", Arg::None,
    "  --check-mmcif  \tchecks alignment against _atom_site.label_seq_id" },
  { Query, 0, "", "query", Arg::Required,
    "  --query=[+|@]CHAIN  \tAlign CHAIN from file INPUT1." },
  { Target, 0, "", "target", Arg::Required,
    "  --target=[+|@]CHAIN  \tAlign CHAIN from file INPUT2." },
  { TextAlign, 0, "", "text-align", Arg::None,
    "  --text-align  \tAlign characters in two strings (for testing)." },

  { NoOp, 0, "", "", Arg::None, "\nScoring (absolute values):" },
  { Blosum62, 0, "", "blosum62", Arg::None,
    "  --blosum62  \tUse BLOSUM62 score matrix." },
  { Partial, 0, "", "partial", Arg::YesNo,
    "  --partial=y|n  \tUse scoring meant to align partially-modelled polymer"
    " to its full sequence (default in 1st mode)." },
  { Match, 0, "", "match", Arg::Int,
    "  --match=INT  \tMatch score (default: 1)." },
  { Mismatch, 0, "", "mism", Arg::Int,
    "  --mism=INT  \tMismatch penalty (default: -1)." },
  { GapOpen, 0, "", "gapo", Arg::Int,
    "  --gapo=INT  \tGap opening penalty (default: -1)." },
  { GapExt, 0, "", "gape", Arg::Int,
    "  --gape=INT  \tGap extension penalty (default: -1)." },

  { NoOp, 0, "", "", Arg::None, "\nOutput options:" },
  { PrintOneLetter, 0, "p", "", Arg::None,
    "  -p  \tPrint formatted alignment with one-letter codes." },
  { Rmsd, 0, "", "rmsd", Arg::None,
    "  --rmsd  \tIn addition to aligning two CHAINs (--query and --target), "
                 "superpose structures and print RMSD." },
  CommonUsage[Verbose],
  { 0, 0, 0, 0, 0, 0 }
};

void print_alignment_details(const gemmi::AlignmentResult& result,
                             const std::string& chain_name,
                             const gemmi::ConstResidueSpan& polymer,
                             const gemmi::Entity& ent) {
  std::vector<int> gaps = prepare_target_gapo(polymer, ent.polymer_type);
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
        std::putchar(gap < gaps.end() && *(gap++) == 0 ? '^' : ' ');
        if (op == 'D' || fmon != res->name)
          printf("    <-- BAD");
        res++;
      }
      printf("\n");
    }
  }
}

// similar to print_alignment_details(), but for --target/--query
void print_alignment_details_tq(const gemmi::AlignmentResult& result,
                                const std::vector<std::string>& seq1,
                                const std::vector<std::string>& seq2) {
  size_t n1 = 0;
  size_t n2 = 0;
  printf("  query  -  target\n");
  for (gemmi::AlignmentResult::Item item : result.cigar) {
    char op = item.op();
    for (uint32_t i = 0; i < item.len(); ++i) {
      if (op == 'I' || op == 'M') {
        printf("%4zu %-3s -", n1+1, seq1[n1].c_str());
        n1++;
      } else {
        printf("         -");
      }
      if (op == 'D' || op == 'M') {
        printf(" %-3s %4zu", seq2[n2].c_str(), n2+1);
        if (op == 'M' && seq1[n1-1] != seq2[n2])
          printf("    <-- DIFFERS");
        n2++;
      }
      printf("\n");
    }
  }
}

void check_label_seq_id(const gemmi::AlignmentResult& result,
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

const gemmi::Model& get_first_model(gemmi::Structure& st) {
  gemmi::setup_entities(st);
  if (st.models.empty())
    gemmi::fail("No atoms found. Wrong input file?");
  return st.models[0];
}

gemmi::ConstResidueSpan get_polymer(const gemmi::Model& model,
                                    const std::string& chain_name) {
  const gemmi::Chain* ch = model.find_chain(chain_name);
  if (!ch)
    gemmi::fail("No such chain: " + chain_name);
  gemmi::ConstResidueSpan polymer = ch->get_polymer();
  if (!polymer)
    gemmi::fail("Polymer not found in chain " + chain_name);
  return polymer;
}

const gemmi::Entity* take_entity(gemmi::Structure& st, const char* arg) {
  if (arg[0] == '+') {
    auto polymer = get_polymer(get_first_model(st), arg+1);
    if (const gemmi::Entity* ent = st.get_entity_of(polymer))
      return ent;
    gemmi::fail("No sequence (SEQRES) for chain ", arg+1);
  } else if (arg[0] == '@') {
    if (const gemmi::Entity* ent = st.get_entity(arg+1))
      return ent;
    gemmi::fail("No such entity: ", arg+1);
  }
  return nullptr;
}

void print_one_letter_alignment(const gemmi::AlignmentResult& result,
                                const std::string& a, const std::string& b) {
  std::fputs(result.formatted(a, b).c_str(), stdout);
}

void print_result_summary(const gemmi::AlignmentResult& result) {
  printf("identity: %.2f%% / %.2f%%  CIGAR: %s\n",
         result.calculate_identity(1), result.calculate_identity(2),
         result.cigar_str().c_str());
}

std::vector<std::string> string_to_vector(const std::string& s) {
  std::vector<std::string> v(s.size());
  for (size_t i = 0; i != v.size(); ++i)
    v[i] = s[i];
  return v;
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.check_exclusive_pair(Blosum62, Partial);
  if ((bool)p.options[Query] != (bool)p.options[Target]) {
    std::fputs("Options --query and --target must be used together.\n", stderr);
    return 1;
  }
  p.check_exclusive_pair(TextAlign, Query);

  gemmi::AlignmentScoring scoring;
  bool first_mode = !(p.options[TextAlign] || p.options[Query]);
  if (p.options[Blosum62])
    scoring = *gemmi::AlignmentScoring::blosum62();
  else if (p.is_yes(Partial, first_mode))
    scoring = *gemmi::AlignmentScoring::partial_model();
  else
    scoring = *gemmi::AlignmentScoring::simple();

  if (p.options[Match])
    scoring.match = std::atoi(p.options[Match].arg);
  if (p.options[Mismatch])
    scoring.mismatch = std::atoi(p.options[Mismatch].arg);
  if (p.options[GapOpen])
    scoring.gapo = std::atoi(p.options[GapOpen].arg);
  if (p.options[GapExt])
    scoring.gape = std::atoi(p.options[GapExt].arg);
  bool verbose = p.options[Verbose];

  if (p.options[TextAlign]) {
    p.require_positional_args(2);
    if (p.nonOptionsCount() != 2) {
      std::fputs("two input strings are expected with --text-align\n", stderr);
      return 1;
    }
    std::string text1 = p.nonOption(0);
    std::string text2 = p.nonOption(1);
    std::vector<int> target_gapo(1, scoring.good_gapo);
    auto result = gemmi::align_string_sequences(string_to_vector(text1),
                                                string_to_vector(text2),
                                                target_gapo, &scoring);
    printf("Score: %d   CIGAR: %s\n", result.score, result.cigar_str().c_str());
    if (p.options[PrintOneLetter])
      print_one_letter_alignment(result, text1, text2);
    return 0;
  }

  try {
    if (p.options[Query]) {
      int n_files = p.nonOptionsCount();
      if (n_files != 2 && n_files != 1)
        gemmi::fail("one or two input files are expected with --query/--target");
      std::vector<std::string> query;
      gemmi::PolymerType ptype = gemmi::PolymerType::Unknown;
      gemmi::Structure st1 = gemmi::read_structure_gz(p.coordinate_input_file(0));
      if (const gemmi::Entity* ent = take_entity(st1, p.options[Query].arg)) {
        query = ent->full_sequence;
        ptype = ent->polymer_type;
      } else {
        query = get_polymer(get_first_model(st1), p.options[Query].arg)
                .extract_sequence();
      }
      gemmi::AlignmentResult result;
      gemmi::Structure st2_;
      if (n_files == 2)
        st2_ = gemmi::read_structure_gz(p.coordinate_input_file(1));
      gemmi::Structure& st2 = n_files == 2 ? st2_ : st1;
      if (const gemmi::Entity* ent = take_entity(st2, p.options[Target].arg)) {
        std::vector<int> target_gapo(1, scoring.good_gapo);
        result = gemmi::align_string_sequences(query, ent->full_sequence,
                                               target_gapo, &scoring);
        print_result_summary(result);
        if (p.options[PrintOneLetter])
          print_one_letter_alignment(result, gemmi::one_letter_code(query),
                                   gemmi::one_letter_code(ent->full_sequence));
        if (verbose)
          print_alignment_details_tq(result, query, ent->full_sequence);
      } else {
        auto polymer = get_polymer(get_first_model(st2), p.options[Target].arg);
        if (ptype == gemmi::PolymerType::Unknown)
          if (const gemmi::Entity* entity = st2.get_entity_of(polymer))
            ptype = entity->polymer_type;
        result = gemmi::align_sequence_to_polymer(query, polymer, ptype, &scoring);
        print_result_summary(result);
        if (p.options[PrintOneLetter])
          print_one_letter_alignment(result, gemmi::one_letter_code(query),
                                     gemmi::one_letter_code(polymer.extract_sequence()));
        if (verbose)
          print_alignment_details_tq(result, query, polymer.extract_sequence());
      }
      if (p.options[Rmsd]) {
        auto poly1 = get_polymer(st1.models.at(0), p.options[Query].arg);
        auto poly2 = get_polymer(st2.models.at(0), p.options[Target].arg);
        printf("We superpose polymers using only matching residues and atoms w/o altloc.\n");
        gemmi::SupResult r;
        r = gemmi::calculate_superposition(poly1, poly2, ptype, gemmi::SupSelect::All);
        printf("RMSD of %zu atoms: %g\n", r.count, r.rmsd);
        r = gemmi::calculate_superposition(poly1, poly2, ptype, gemmi::SupSelect::CaP);
        printf("RMSD of %zu CA/P atoms: %g\n", r.count, r.rmsd);

        // this last part is not particularly useful
        auto mpoly2 = st2.models[0].find_chain(p.options[Target].arg)->get_polymer();
        for (gemmi::Residue& res : mpoly2)
          for (gemmi::Atom& atom : res.atoms)
            atom.pos = gemmi::Position(r.transform.apply(atom.pos));
        r = gemmi::calculate_current_rmsd(poly1, mpoly2, ptype, gemmi::SupSelect::All);
        printf("   the same rotation+shift applied to %zu atoms: %g\n", r.count, r.rmsd);
      }
      return 0;
    }

    p.require_input_files_as_args();
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      std::string input = p.coordinate_input_file(i);
      if (i > 0)
        printf("\n");
      if (verbose || p.nonOptionsCount() > 1)
        printf("File: %s\n", input.c_str());
      gemmi::Structure st = gemmi::read_structure_gz(input);
      const gemmi::Model& model = get_first_model(st);
      if (st.models.size() > 1)
        printf("Warning: using only model 1 of %zu.\n", st.models.size());
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
                                             ent->polymer_type, &scoring);
        printf("%s chain %s  ", st.name.c_str(), chain.name.c_str());
        print_result_summary(result);
        if (p.options[CheckMmcif]) {
          if (st.input_format == gemmi::CoorFormat::Pdb)
            printf("Option --check-mmcif ignored for PDB file: %s\n", input.c_str());
          else
            check_label_seq_id(result, polymer);
        }
        if (p.options[PrintOneLetter])
          print_one_letter_alignment(result,
                                     gemmi::one_letter_code(ent->full_sequence),
                                     gemmi::one_letter_code(polymer.extract_sequence()));
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
