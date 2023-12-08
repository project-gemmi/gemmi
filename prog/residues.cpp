// Copyright 2018 Global Phasing Ltd.

#include <stdio.h>
#include <string>
#include "gemmi/select.hpp"
#include "gemmi/polyheur.hpp"  // for setup_entities
#include "gemmi/align.hpp"     // for assign_label_seq_id
#include "gemmi/mmread_gz.hpp" // for read_structure_gz
#include "gemmi/enumstr.hpp"   // for polymer_type_to_string, entity_type_to_string
#include "histogram.h"         // for terminal_columns

#define GEMMI_PROG residues
#include "options.h"

namespace {

enum OptionIndex {
  FormatIn=4, Match, Label, CheckSeqId, NoAlt, Short, Chains, Ent
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:\n " EXE_NAME " [options] INPUT[...]"
    "\nPrints chains, residues and atoms. Or without atoms. Or only entities." },
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
  { Short, 0, "s", "--short", Arg::None,
    "  -s, --short  \tShorter output (no atom info). Can be given 2x or 3x." },
  { Ent, 0, "e", "entities", Arg::None,
    "  -e, --entities  \tList (so-called, in mmCIF speak) entities." },
  { Chains, 0, "c", "chains", Arg::None,
    "  -c, --chains  \tList chain IDs." },
  { NoOp, 0, "", "", Arg::None,
    "INPUT is a coordinate file (mmCIF, PDB, etc)."
    "\nThe optional selection SEL has MMDB syntax:"
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

void print_long_info(const gemmi::Model& model, OptParser& p) {
  bool print_alt = !p.options[NoAlt];
  for (const gemmi::Chain& chain : model.chains) {
    int line_count = 0;
    for (const gemmi::Residue& res : chain.residues) {
      if (res.atoms.empty())
        continue;
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
      for (const gemmi::Atom& at : res.atoms)
        if (!prev || *prev != at.name) {
          printf(" %s", at.name.c_str());
          if (print_alt && at.altloc)
            printf(":%c", at.altloc);
          prev = &at.name;
        } else {
          if (print_alt) {
            putchar(',');
            if (at.altloc)
              putchar(at.altloc);
          }
        }
      putchar('\n');
      line_count++;

    }
    if (line_count != 0)
      putchar('\n');
  }
}

void print_short_info(const gemmi::Model& model, OptParser& p) {
  int short_level = p.options[Short].count();
  const int kWrap = 5;    // for level 1
  const int kLimit = terminal_columns() - 23;  // for levels 2 and 3
  for (const gemmi::Chain& chain : model.chains) {
    int col = 0;
    int counter = 0;
    auto prev = gemmi::EntityType::Unknown;
    for (const gemmi::Residue& res : chain.residues) {
      if (!chain.is_first_in_group(res)) {  // microheterogeneity
        if (counter < 8 && short_level < 3)
          col += printf("/%s", res.name.c_str());
        continue;
      }
      if (res.entity_type != prev) {
        if (counter != 0) {
          if (short_level > 1 && col >= kLimit)
            printf("...  (%d residues)", counter);
          putchar('\n');
        }
        const char* etype = entity_type_to_string(res.entity_type);
        if (p.options[Label])
          col = printf("%s (%s) %-11s", chain.name.c_str(), res.subchain.c_str(), etype);
        else
          col = printf("%-3s %-11s ", chain.name.c_str(), etype);
        counter = 0;
        prev = res.entity_type;
      }
      if (short_level == 1) {
        if (counter == kWrap) {
          printf("\n                  ");
          counter = 0;
        }
        printf(" %5s%c %-3s",
               res.seqid.num.str().c_str(), res.seqid.icode, res.name.c_str());
      } else {
        if (col < kLimit) {
          if (short_level == 2) {
            col += printf(" %-3s", res.name.c_str());
          } else { // short_level > 2
            // cf. pdbx_one_letter_code()
            char c = gemmi::find_tabulated_residue(res.name).fasta_code();
            if (res.entity_type == gemmi::EntityType::Polymer && c != 'X')
              col += printf("%c", c);
            else
              col += printf("(%s)", res.name.c_str());
          }
        }
      }
      ++counter;
    }
    if (short_level > 1 && col >= kLimit)
      printf("...  (%d residues)", counter);
    putchar('\n');
  }
}

void print_chain_info(const gemmi::Model& model) {
  for (const gemmi::Chain& chain : model.chains) {
    printf("%s  length/count:", chain.name.c_str());
    gemmi::EntityType prev_et = gemmi::EntityType::Unknown;
    int counter = 0;
    for (const gemmi::Residue& res :  chain.first_conformer()) {
      if (res.entity_type != prev_et) {
        if (counter != 0) {
          printf("  %s %d", gemmi::entity_type_to_string(prev_et), counter);
        }
        counter = 0;
        prev_et = res.entity_type;
      }
      ++counter;
    }
    printf("  %s %d", gemmi::entity_type_to_string(prev_et), counter);
    putchar('\n');
  }
}

void print_entity_info(const gemmi::Structure& st) {
  if (st.models.size() > 1)
    printf("Checking only the first model.\n");
  const gemmi::Model& model = st.models.at(0);
  printf("Polymers\n");
  std::map<std::string, std::string> sub_to_strand = model.subchain_to_chain();
  for (const gemmi::Entity& ent : st.entities)
    if (ent.entity_type == gemmi::EntityType::Polymer) {
      printf("  entity %s, %s, length %zu, subchains:\n",
             ent.name.c_str(),
             gemmi::polymer_type_to_string(ent.polymer_type),
             ent.full_sequence.size());
      for (const std::string& sub : ent.subchains) {
        auto strand = sub_to_strand.find(sub);
        if (strand == sub_to_strand.end())
          throw std::runtime_error("error in subchain_to_chain()");
        auto polymer = model.get_subchain(sub);
        auto conformer = polymer.first_conformer();
        int length = 0;
        std::vector<std::pair<int,int>> gaps;
        int prev = gemmi::SeqId::OptionalNum::None;
        for (const gemmi::Residue& res : conformer) {
          if (prev != gemmi::SeqId::OptionalNum::None && prev != *res.label_seq - 1)
            gaps.emplace_back(prev+1, *res.label_seq-1);
          prev = *res.label_seq;
          ++length;
        }
        printf("    - %s from strand %s, %d residues",
               sub.c_str(), strand->second.c_str(), length);
        if (!polymer.empty()) {
          printf(": %s-%s", polymer.front().label_seq.str().c_str(),
                              polymer.back().label_seq.str().c_str());
          if (!gaps.empty()) {
            printf(" except");
            for (std::pair<int, int> gap : gaps) {
              printf(" %d", gap.first);
              if (gap.second != gap.first)
                printf("-%d", gap.second);
            }
          }
        }
        putchar('\n');
      }
    }
  printf("Others\n");
  for (const gemmi::Entity& ent : st.entities)
    if (ent.entity_type != gemmi::EntityType::Polymer) {
      printf("  entity %s, %s",
             ent.name.c_str(), gemmi::entity_type_to_string(ent.entity_type));
      if (ent.entity_type != gemmi::EntityType::Branched) {
        // one residue is expected
        std::string name;
        for (const std::string& sub : ent.subchains)
          for (const gemmi::Residue& res : model.get_subchain(sub)) {
            if (name.empty()) {
              name = res.name;
            } else if (name != res.name) {
              name += ", ...";
              break;
            }
          }
        printf(" (%s)", name.c_str());
      }
      printf(", subchains: %s\n", gemmi::join_str(ent.subchains, ' ').c_str());
    }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_input_files_as_args();
  gemmi::CoorFormat format = coor_format_as_enum(p.options[FormatIn]);
  int status = 0;
  try {
    for (int i = 0; i < p.nonOptionsCount(); ++i) {
      if (i != 0)
        putchar('\n');
      std::string input = p.coordinate_input_file(i);
      printf("%s\n", input.c_str());
      gemmi::Structure st = gemmi::read_structure_gz(input, format);
      if (p.options[Match]) {
        gemmi::Selection sel(p.options[Match].arg);
        sel.remove_not_selected(st);
      }
      if (p.options[Label] || p.options[Ent]) {
        gemmi::setup_entities(st);
        // hidden feature: -ll generates label_seq even if SEQRES is missing
        bool force = p.options[Label].count() > 1;
        gemmi::assign_label_seq_id(st, force);
      } else if (p.options[Short]) {
        gemmi::add_entity_types(st, false);
      }
      if (p.options[CheckSeqId]) {
        bool ok = check_sequence_id(st);
        if (!ok)
          ++status;
        continue;
      }
      if (p.options[Ent]) {
        print_entity_info(st);
        continue;
      }
      if (p.options[Chains])
        st.merge_chain_parts();
      for (gemmi::Model& model : st.models) {
        if (st.models.size() != 1)
          printf("Model %s\n", model.name.c_str());
        if (p.options[Chains])
          print_chain_info(model);
        else if (p.options[Short])
          print_short_info(model, p);
        else
          print_long_info(model, p);
      }
    }
  } catch (std::exception& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  }
  return status;
}
