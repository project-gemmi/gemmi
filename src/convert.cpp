// Copyright 2017 Global Phasing Ltd.

#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/polyheur.hpp"  // for setup_entities, remove_hydrogens, ...
#include "gemmi/align.hpp"     // for assign_label_seq_id
#include "gemmi/to_pdb.hpp"    // for write_pdb, ...
#include "gemmi/fstream.hpp"   // for Ofstream, Ifstream
#include "gemmi/to_mmcif.hpp"  // for update_mmcif_block
#include "gemmi/remarks.hpp"   // for read_metadata_from_remarks
#include "gemmi/assembly.hpp"  // for ChainNameGenerator, change_to_assembly
#include "gemmi/pirfasta.hpp"  // for read_pir_or_fasta
#include "gemmi/resinfo.hpp"   // for expand_protein_one_letter
#include "gemmi/read_coor.hpp" // for read_structure_gz
#include "gemmi/select.hpp"    // for parse_cid

#include <cstring>
#include <iostream>
#include <algorithm>           // for sort
#include <map>

#define GEMMI_PROG convert
#include "options.h"
#include "cifmod.h"

namespace cif = gemmi::cif;
using gemmi::CoorFormat;
using gemmi::HowToNameCopiedChain;

namespace {

struct ConvArg: public Arg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"mmjson", "pdb", "mmcif", "ccd"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }

  static option::ArgStatus AnisouChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"yes", "no", "heavy"});
  }

  static option::ArgStatus NcsChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"dup", "num", "x"});
  }
};

enum OptionIndex {
  FormatIn=AfterCifModOptions, FormatOut, PdbxStyle, BlockName,
  ExpandNcs, AsAssembly,
  RemoveH, RemoveWaters, RemoveLigWat, TrimAla, Select, Remove,
  ShortTer, Linkr, CopyRemarks, Minimal, ShortenCN, RenameChain, SetSeq, Anisou,
  SegmentAsChain, OldPdb, ForceLabel
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions between PDB, mmCIF and mmJSON."
    "\nFORMAT can be specified as one of: mmcif, mmjson, pdb, ccd (read-only)."
    "\nccd = coordinates of a chemical component from CCD or monomer library."
    "\n\nGeneral options:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { FormatIn, 0, "", "from", ConvArg::FileFormat,
    "  --from=FORMAT  \tInput format (default: from the file extension)." },
  { FormatOut, 0, "", "to", ConvArg::FileFormat,
    "  --to=FORMAT  \tOutput format (default: from the file extension)." },

  { NoOp, 0, "", "", Arg::None, "\nCIF output options:" },
  { PdbxStyle, 0, "", "pdbx-style", Arg::None,
    "  --pdbx-style  \tSimilar styling (formatting) as in wwPDB." },
  { BlockName, 0, "b", "block", Arg::Required,
    "  -b NAME, --block=NAME  \tSet block name and default _entry.id" },
  CifModUsage[SortCif],
  CifModUsage[SkipCat],

  { NoOp, 0, "", "", Arg::None, "\nPDB input options:" },
  { SegmentAsChain, 0, "", "segment-as-chain", Arg::None,
    "  --segment-as-chain \tAppend segment id to label_asym_id (chain name)." },
  { OldPdb, 0, "", "old-pdb", Arg::None,
    "  --old-pdb \tRead only the first 72 characters in line." },
  { ForceLabel, 0, "L", "force-label", Arg::None,
    "  -L, --force-label  \tAdd label_seq_id even if SEQRES is missing" },

  { NoOp, 0, "", "", Arg::None, "\nPDB output options:" },
  { ShortTer, 0, "", "short-ter", Arg::None,
    "  --short-ter  \tWrite PDB TER records without numbers (iotbx compat.)." },
  { Linkr, 0, "", "linkr", Arg::None,
    "  --linkr  \tWrite LINKR record (for Refmac) if link_id is known." },
  { CopyRemarks, 0, "", "copy-remarks", Arg::None,
    "  --copy-remarks  \t(pdb->pdb only) Copy REMARK records." },

  { NoOp, 0, "", "", Arg::None, "\nAny output options:" },
  { Minimal, 0, "", "minimal", Arg::None,
    "  --minimal  \tWrite only the most essential records." },
  { ShortenCN, 0, "", "shorten", Arg::None,
    "  --shorten  \tShorten chain names to 1 (if # < 63) or 2 characters." },
  { RenameChain, 0, "", "rename-chain", Arg::ColonPair,
    "  --rename-chain=OLD:NEW  \tRename chain OLD to NEW "
    "(--rename-chain=:A adds missing chain IDs)." },
  { SetSeq, 0, "s", "", Arg::Required,
    "  -s FILE  \tUse sequence from FILE (PIR or FASTA format), "
    "which must contain either one sequence (for all chains) "
    "or as many sequences as there are chains." },
  { Anisou, 0, "", "anisou", ConvArg::AnisouChoice,
    "  --anisou=yes|no|heavy  \tAdd or remove ANISOU records." },

  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { ExpandNcs, 0, "", "expand-ncs", ConvArg::NcsChoice,
    "  --expand-ncs=dup|num|x  \tExpand strict NCS from in MTRIXn or"
    " _struct_ncs_oper. New chain names are the same, have added numbers,"
    " or the shortest unused names are picked."},
  { AsAssembly, 0, "", "assembly", Arg::Required,
    "  --assembly=ID  \tOutput bioassembly with given ID (1, 2, ...)." },
  { RemoveH, 0, "", "remove-h", Arg::None,
    "  --remove-h  \tRemove hydrogens." },
  { RemoveWaters, 0, "", "remove-waters", Arg::None,
    "  --remove-waters  \tRemove waters." },
  { RemoveLigWat, 0, "", "remove-lig-wat", Arg::None,
    "  --remove-lig-wat  \tRemove ligands and waters." },
  { TrimAla, 0, "", "trim-to-ala", Arg::None,
    "  --trim-to-ala  \tTrim aminoacids to alanine." },
  { Select, 0, "", "select", Arg::Required,
    "  --select=CID  \tOutput only the selection." },
  { Remove, 0, "", "remove", Arg::Required,
    "  --remove=CID  \tRemove the selection." },

  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

std::string format_as_string(CoorFormat format) {
  switch (format) {
    case CoorFormat::Unknown: return "unknown";
    case CoorFormat::UnknownAny: return "unknown";
    case CoorFormat::Pdb: return "pdb";
    case CoorFormat::Mmcif: return "mmcif";
    case CoorFormat::Mmjson: return "mmjson";
    case CoorFormat::ChemComp: return "chemcomp";
  }
  gemmi::unreachable();
}

void convert(gemmi::Structure& st,
             const std::string& output, CoorFormat output_type,
             const std::vector<option::Option>& options) {
  if (st.models.empty())
    gemmi::fail("No atoms in the input file. Wrong file format?");

  if (st.input_format == CoorFormat::Pdb) {
    gemmi::read_metadata_from_remarks(st);
    gemmi::setup_entities(st);
    if (!options[SetSeq])
      gemmi::assign_label_seq_id(st, options[ForceLabel]);
    if (!options[CopyRemarks])
      st.raw_remarks.clear();
  }

  if (options[Anisou]) {
    char anisou_opt = options[Anisou].arg[0];
    if (anisou_opt == 'n') {
      gemmi::remove_anisou(st);
    } else if (anisou_opt == 'y') {
      gemmi::ensure_anisou(st);
    } else if (anisou_opt == 'h') {
      for (gemmi::Model& model : st.models)
        for (gemmi::Chain& chain : model.chains)
          for (gemmi::Residue& res : chain.residues)
            for (gemmi::Atom& atom : res.atoms)
              if (!atom.is_hydrogen())
                gemmi::ensure_anisou(atom);
    }
  }

  if (options[Select])
    gemmi::parse_cid(options[Select].arg).remove_not_selected(st);
  if (options[Remove])
    gemmi::parse_cid(options[Remove].arg).remove_selected(st);
  if (st.models.empty())
    gemmi::fail("all models got removed");


  for (const option::Option* opt = options[RenameChain]; opt; opt = opt->next()) {
    const char* sep = std::strchr(opt->arg, ':');
    std::string old_name(opt->arg, sep);
    std::string new_name(sep+1);
    if (gemmi::Chain* chain = st.first_model().find_chain(old_name))
      gemmi::rename_chain(st, *chain, new_name);
  }
  if (options[ShortenCN]) {
    shorten_chain_names(st);
  } else if (output_type == CoorFormat::Pdb) {
    for (const gemmi::Chain& chain : st.models[0].chains)
      if (chain.name.size() > 2)
        gemmi::fail("long chain name cannot be written in the PDB format: " +
                    chain.name + "\nTry option --shorten");
  }
  if (options[SetSeq]) {
    gemmi::Ifstream stream(options[SetSeq].arg);
    std::string seq = gemmi::read_pir_or_fasta(stream.ref());
    for (gemmi::Entity& ent : st.entities)
      if (ent.entity_type == gemmi::EntityType::Polymer) {
        ent.full_sequence.clear();
        for (char letter : seq) {
          const char* str = gemmi::expand_protein_one_letter(letter);
          if (!str)
            gemmi::fail("unexpected letter in protein sequence: ", letter);
          ent.full_sequence.emplace_back(str);
        }
      }
    gemmi::deduplicate_entities(st);
    gemmi::assign_label_seq_id(st, options[ForceLabel]);
  }

  HowToNameCopiedChain how = HowToNameCopiedChain::AddNumber;
  if (output_type == CoorFormat::Pdb)
    how = HowToNameCopiedChain::Short;
  if (options[AsAssembly]) {
    gemmi::change_to_assembly(st, options[AsAssembly].arg, how,
                              options[Verbose] ? &std::cerr : nullptr);
    // After this change Assembly instructions can be outdated.
    // Should they be preserved anyway or removed? Currently - removing.
    st.assemblies.clear();
  }

  if (options[ExpandNcs]) {
    HowToNameCopiedChain how_ncs;
    switch (options[ExpandNcs].arg[0]) {
      case 'd': how_ncs = HowToNameCopiedChain::Dup; break;
      case 'x': how_ncs = HowToNameCopiedChain::Short; break;
      default: how_ncs = HowToNameCopiedChain::AddNumber; break;
    }
    gemmi::expand_ncs(st, how_ncs);
  }

  if (options[RemoveH])
    remove_hydrogens(st);

  if (options[RemoveWaters]) {
    remove_waters(st);
    remove_empty_chains(st);
  }

  if (options[RemoveLigWat]) {
    remove_ligands_and_waters(st);
    remove_empty_chains(st);
  }

  if (options[TrimAla])
    for (gemmi::Model& model : st.models)
      for (gemmi::Chain& chain : model.chains)
        trim_to_alanine(chain);


  if (options[SegmentAsChain])
    for (gemmi::Model& model : st.models)
      split_chains_by_segments(model, gemmi::HowToNameCopiedChain::Dup);

  gemmi::Ofstream os(output, &std::cout);

  if (output_type == CoorFormat::Mmcif || output_type == CoorFormat::Mmjson) {
    if (options[BlockName])
      st.name = options[BlockName].arg;
    cif::Document doc;
    doc.blocks.resize(1);
    if (options[Minimal])
      gemmi::add_minimal_mmcif_data(st, doc.blocks[0]);
    else
      gemmi::update_mmcif_block(st, doc.blocks[0], gemmi::MmcifOutputGroups(true));
    apply_cif_doc_modifications(doc, options);

    if (output_type == CoorFormat::Mmcif) {
      auto style = options[PdbxStyle] ? cif::Style::Pdbx
                                      : cif::Style::PreferPairs;
      write_cif_to_stream(os.ref(), doc, style);
    } else /*output_type == CoorFormat::Mmjson*/ {
      cif::JsonWriter writer(os.ref());
      writer.set_mmjson();
      writer.write_json(doc);
    }
  } else if (output_type == CoorFormat::Pdb) {
    // call wrapper from output.cpp - to make building faster
    gemmi::PdbWriteOptions opt;
    if (options[ShortTer])
      opt.numbered_ter = false;
    if (options[Linkr])
      opt.use_linkr = true;
    if (options[Minimal])
      gemmi::write_minimal_pdb(st, os.ref(), opt);
    else
      gemmi::write_pdb(st, os.ref(), opt);
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  std::map<std::string, CoorFormat> filetypes{{"pdb", CoorFormat::Pdb},
                                              {"mmcif", CoorFormat::Mmcif},
                                              {"mmjson", CoorFormat::Mmjson},
                                              {"ccd", CoorFormat::ChemComp}};
  CoorFormat in_type = p.options[FormatIn] ? filetypes[p.options[FormatIn].arg]
                                           : CoorFormat::UnknownAny;

  char pdb_code_type = in_type == CoorFormat::Pdb ? 'P' : 'M';
  std::string input = p.coordinate_input_file(0, pdb_code_type);
  const char* output = p.nonOption(1);

  CoorFormat out_type = p.options[FormatOut]
    ? filetypes[p.options[FormatOut].arg]
    : gemmi::coor_format_from_ext_gz(output);
  if (out_type == CoorFormat::ChemComp) {
    std::cerr << "The output format cannot be ccd.\n";
    return 1;
  }
  if (out_type == CoorFormat::Unknown) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }
  if (p.options[Verbose])
    std::cerr << "Converting " << input << " to " << format_as_string(out_type)
              << "..." << std::endl;
  try {
    gemmi::Structure st;
    if (p.options[OldPdb]) {
      gemmi::PdbReadOptions options;
      options.max_line_length = 72;
      st = gemmi::read_pdb_gz(input, options);
    } else {
      st = gemmi::read_structure_gz(input, in_type);
    }
    convert(st, output, out_type, p.options);
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  } catch (std::invalid_argument& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 3;
  }
  if (p.options[Verbose])
    std::cerr << "Done." << std::endl;
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
