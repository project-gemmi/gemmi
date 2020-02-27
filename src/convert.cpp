// Copyright 2017 Global Phasing Ltd.

#include "gemmi/gzread.hpp"
#include "gemmi/mmread.hpp"    // for CoorFormat
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/polyheur.hpp"  // for remove_hydrogens, ...
#include "gemmi/to_pdb.hpp"    // for write_pdb, ...
#include "gemmi/ofstream.hpp"  // for Ofstream
#include "gemmi/to_mmcif.hpp"  // for update_cif_block
#include "gemmi/remarks.hpp"   // for read_metadata_from_remarks
#include "gemmi/labelseq.hpp"  // for assign_label_seq_id

#include <cstring>
#include <iostream>
#include <algorithm>           // for sort
#include <map>

#define GEMMI_PROG convert
#include "options.h"
#include "cifmod.h"

namespace cif = gemmi::cif;
using gemmi::CoorFormat;

struct ConvArg: public Arg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"mmjson", "pdb", "mmcif", "ccd"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }

  static option::ArgStatus NcsChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"dup", "new"});
  }
};

enum OptionIndex { FormatIn=AfterCifModOptions, FormatOut, PdbxStyle, BlockName,
                   ExpandNcs, RemoveH, RemoveWaters, RemoveLigWat, TrimAla,
                   ShortTer, Linkr, SegmentAsChain, };
static const option::Descriptor Usage[] = {
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
  { NoOp, 0, "", "", Arg::None, "\nPDB output options:" },
  { ShortTer, 0, "", "short-ter", Arg::None,
    "  --short-ter  \tWrite PDB TER records without numbers (iotbx compat.)." },
  { Linkr, 0, "", "linkr", Arg::None,
    "  --linkr  \tWrite LINKR record (for Refmac) if link_id is known." },

  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { ExpandNcs, 0, "", "expand-ncs", ConvArg::NcsChoice,
    "  --expand-ncs=dup|new  \tExpand strict NCS specified in MTRIXn or"
    " equivalent. New chain names are the same or have added numbers." },
  { RemoveH, 0, "", "remove-h", Arg::None,
    "  --remove-h  \tRemove hydrogens." },
  { RemoveWaters, 0, "", "remove-waters", Arg::None,
    "  --remove-waters  \tRemove waters." },
  { RemoveLigWat, 0, "", "remove-lig-wat", Arg::None,
    "  --remove-lig-wat  \tRemove ligands and waters." },
  { TrimAla, 0, "", "trim-to-ala", Arg::None,
    "  --trim-to-ala  \tTrim aminoacids to alanine." },

  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

static const char symbols[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "abcdefghijklmnopqrstuvwxyz0123456789";

enum class ChainNaming { Short, AddNum, Dup };

static void expand_ncs(gemmi::Structure& st, ChainNaming ch_naming) {
  int n_ops = std::count_if(st.ncs.begin(), st.ncs.end(),
                            [](gemmi::NcsOp& op){ return !op.given; });
  if (n_ops == 0)
    return;
  for (gemmi::Model& model : st.models) {
    size_t new_length = model.chains.size() * (n_ops + 1);
    if (new_length >= 63 && ch_naming == ChainNaming::Short)
      ch_naming = ChainNaming::AddNum;
    model.chains.reserve(new_length);
    if (ch_naming == ChainNaming::Dup)
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
          res.segment = "0";
    auto orig_end = model.chains.cend();
    std::vector<std::string> used_names;
    if (ch_naming != ChainNaming::Dup)
      for (const gemmi::Chain& chain : model.chains)
        used_names.push_back(chain.name);
    for (const gemmi::NcsOp& op : st.ncs)
      if (!op.given) {
        for (auto ch = model.chains.cbegin(); ch != orig_end; ++ch) {
          // the most difficult part - choosing names for new chains
          std::string name;
          if (ch_naming == ChainNaming::Short) {
            for (char symbol : symbols) {
              name = std::string(1, symbol);
              if (!gemmi::in_vector(name, used_names))
                break;
            }
            used_names.push_back(name);
          } else if (ch_naming == ChainNaming::AddNum) {
            name = ch->name + op.id;
            while (gemmi::in_vector(name, used_names))
              name += "a";
            used_names.push_back(name);
          } else { // ChainNaming::Dup
            name = ch->name;
          }

          model.chains.push_back(*ch);
          gemmi::Chain& new_chain = model.chains.back();
          new_chain.name = name;

          for (gemmi::Residue& res : new_chain.residues) {
            for (gemmi::Atom& a : res.atoms)
              a.pos = op.apply(a.pos);
            if (!res.subchain.empty())
              res.subchain = name + ":" + res.subchain;
            if (ch_naming == ChainNaming::Dup)
              res.segment = op.id;
          }
        }
      }
  }
  for (gemmi::NcsOp& op : st.ncs)
    op.given = true;
}


std::vector<gemmi::Chain> split_by_segments(gemmi::Chain& orig) {
  std::vector<gemmi::Chain> chains;
  std::vector<gemmi::Residue> orig_res;
  orig_res.swap(orig.residues);
  for (auto seg_start = orig_res.begin(); seg_start != orig_res.end(); ) {
    const std::string& sn = seg_start->segment;
    auto ch = std::find_if(chains.begin(), chains.end(), [&](gemmi::Chain& c) {
                return !c.residues.empty() && c.residues[0].segment == sn; });
    if (ch == chains.end()) {
      chains.push_back(orig);
      ch = chains.end() - 1;
      // it's not clear how chain naming should be handled
      ch->name += sn;
    }
    auto seg_end = std::find_if(seg_start, orig_res.end(),
                            [&](gemmi::Residue& r) { return r.segment != sn; });
    ch->residues.insert(ch->residues.end(), std::make_move_iterator(seg_start),
                                            std::make_move_iterator(seg_end));
    seg_start = seg_end;
  }
  return chains;
}

static std::string format_as_string(CoorFormat format) {
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

static void convert(gemmi::Structure& st,
                    const std::string& output, CoorFormat output_type,
                    const std::vector<option::Option>& options) {
  if (st.models.empty())
    gemmi::fail("No atoms in the input file. Wrong file format?");

  if (st.input_format == CoorFormat::Pdb) {
    gemmi::read_metadata_from_remarks(st);
    setup_entities(st);
    assign_label_seq_id(st);
  }

  if (options[ExpandNcs]) {
    ChainNaming ch_naming = ChainNaming::AddNum;
    if (options[ExpandNcs].arg[0] == 'd')
     ch_naming = ChainNaming::Dup;
    else if (output_type == CoorFormat::Pdb)
      ch_naming = ChainNaming::Short;
    expand_ncs(st, ch_naming);
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
    for (gemmi::Model& model : st.models) {
      std::vector<gemmi::Chain> new_chains;
      for (gemmi::Chain& chain : model.chains)
        gemmi::vector_move_extend(new_chains, split_by_segments(chain));
      model.chains = std::move(new_chains);
    }

  gemmi::Ofstream os(output, &std::cout);

  if (output_type == CoorFormat::Mmcif || output_type == CoorFormat::Mmjson) {
    if (options[BlockName])
      st.name = options[BlockName].arg;
    cif::Document doc;
    doc.blocks.resize(1);
    update_cif_block(st, doc.blocks[0], /*with_atoms=*/true);
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
    gemmi::write_pdb(st, os.ref(), opt);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  std::string input = p.coordinate_input_file(0);
  const char* output = p.nonOption(1);

  std::map<std::string, CoorFormat> filetypes{{"pdb", CoorFormat::Pdb},
                                              {"mmcif", CoorFormat::Mmcif},
                                              {"mmjson", CoorFormat::Mmjson},
                                              {"ccd", CoorFormat::ChemComp}};

  CoorFormat in_type = p.options[FormatIn] ? filetypes[p.options[FormatIn].arg]
                                           : CoorFormat::UnknownAny;
  if (in_type == CoorFormat::Unknown) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

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
    gemmi::Structure st = gemmi::read_structure_gz(input, in_type);
    convert(st, output, out_type, p.options);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
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
