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

enum OptionIndex {
  FormatIn=AfterCifModOptions, FormatOut, PdbxStyle, BlockName,
  ExpandNcs, ExpandAssembly, RemoveH, RemoveWaters, RemoveLigWat, TrimAla,
  ShortTer, Linkr, SegmentAsChain,
};

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
  { ExpandAssembly, 0, "", "assembly", Arg::Required,
    "  --assembly=ID  \tOutput bioassembly with given ID (1, 2, ...)." },
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

struct ChainNameGenerator {
  enum class How { Short, AddNum, Dup };
  How how;
  std::vector<std::string> used_names;

  ChainNameGenerator(const gemmi::Model& model, How how_) : how(how_) {
    if (how != How::Dup)
      for (const gemmi::Chain& chain : model.chains)
        used_names.push_back(chain.name);
  }
  bool has(const std::string& name) const {
    return gemmi::in_vector(name, used_names);
  }
  const std::string& added(const std::string& name) {
    used_names.push_back(name);
    return name;
  }

  std::string make_short_name() {
    static const char symbols[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                  "abcdefghijklmnopqrstuvwxyz0123456789";
    std::string name(1, 'A');
    for (char symbol : symbols) {
      name[0] = symbol;
      if (!has(name))
        return added(name);
    }
    name += 'A';
    for (char symbol1 : symbols) {
      name[0] = symbol1;
      for (char symbol2 : symbols) {
        name[1] = symbol2;
        if (!has(name))
          return added(name);
      }
    }
    gemmi::fail("run out of 1- and 2-letter chain names");
  }

  std::string make_name_with_numeric_postfix(const std::string& base, int n) {
    std::string name = base;
    name += std::to_string(n);
    while (has(name)) {
      name.resize(base.size());
      name += std::to_string(++n);
    }
    return added(name);
  }

  std::string make_new_name(const std::string& old, int n) {
    switch (how) {
      case How::Short: return make_short_name();
      case How::AddNum: return make_name_with_numeric_postfix(old, n);
      case How::Dup: return old;
    }
    gemmi::unreachable();
  }
};

static void expand_ncs(gemmi::Structure& st, ChainNameGenerator::How how) {
  for (gemmi::Model& model : st.models) {
    size_t orig_size = model.chains.size();
    ChainNameGenerator namegen(model, how);
    for (const gemmi::NcsOp& op : st.ncs)
      if (!op.given) {
        for (size_t i = 0; i != orig_size; ++i) {
          if (how == ChainNameGenerator::How::Dup)
            for (gemmi::Residue& res : model.chains[i].residues)
              res.segment = "0";

          model.chains.push_back(model.chains[i]);
          gemmi::Chain& new_chain = model.chains.back();
          new_chain.name = namegen.make_new_name(new_chain.name, (int)i+1);

          for (gemmi::Residue& res : new_chain.residues) {
            for (gemmi::Atom& a : res.atoms)
              a.pos = op.apply(a.pos);
            if (!res.subchain.empty())
              res.subchain = new_chain.name + ":" + res.subchain;
            if (how == ChainNameGenerator::How::Dup)
              res.segment = op.id;
          }
        }
      }
  }
  for (gemmi::NcsOp& op : st.ncs)
    op.given = true;
}

void expand_assembly(gemmi::Model& model, const gemmi::Assembly& assembly,
                     ChainNameGenerator::How how, bool verbose) {
  using namespace gemmi;
  ChainNameGenerator namegen(model, how);
  std::map<std::string, std::string> subs = model.subchain_to_chain();
  for (const Assembly::Gen& gen : assembly.generators)
    for (const Assembly::Oper& oper : gen.opers) {
      if (verbose) {
        std::cerr << "Applying " << oper.name << " to";
        if (!gen.chains.empty())
          std::cerr << " chains: " << join_str(gen.chains, ',');
        else if (!gen.subchains.empty())
          std::cerr << " subchains: " << join_str(gen.subchains, ',');
        std::cerr << std::endl;
        for (const std::string& chain_name : gen.chains)
          if (!model.find_chain(chain_name))
            std::cerr << "Warning: no chain " << chain_name << std::endl;
        for (const std::string& subchain_name : gen.subchains)
          if (subs.find(subchain_name) == subs.end())
            std::cerr << "Warning: no subchain " << subchain_name << std::endl;
      }
      if (!gen.chains.empty()) {
        // chains are not merged here, multiple chains may have the same name
        std::map<std::string, std::string> new_names;
        size_t orig_size = model.chains.size();
        for (size_t i = 0; i != orig_size; ++i) {
          if (in_vector(model.chains[i].name, gen.chains)) {
            model.chains.push_back(model.chains[i]);
            Chain& new_chain = model.chains.back();
            auto name_iter = new_names.find(model.chains[i].name);
            if (name_iter == new_names.end()) {
              new_chain.name = namegen.make_new_name(new_chain.name, 1);
              new_names.emplace(model.chains[i].name, new_chain.name);
            } else {
              new_chain.name = name_iter->second;
            }
            for (Residue& res : new_chain.residues) {
              for (Atom& a : res.atoms)
                a.pos = Position(oper.transform.apply(a.pos));
              if (!res.subchain.empty())
                res.subchain = new_chain.name + ":" + res.subchain;
            }
          }
        }
      } else if (!gen.subchains.empty()) {
        std::map<std::string, std::string> new_names;
        for (const std::string& subchain_name : gen.subchains) {
          auto sub_iter = subs.find(subchain_name);
          if (sub_iter == subs.end())
            continue;
          const Chain& old_chain = *model.find_chain(sub_iter->second);
          auto name_iter = new_names.find(old_chain.name);
          Chain* new_chain;
          if (name_iter == new_names.end()) {
            std::string new_name = namegen.make_new_name(old_chain.name, 1);
            new_names.emplace(old_chain.name, new_name);
            model.chains.emplace_back(new_name);
            new_chain = &model.chains.back();
          } else {
            new_chain = model.find_chain(name_iter->second);
          }
          for (const Residue& res : old_chain.get_subchain(subchain_name)) {
            new_chain->residues.push_back(res);
            Residue& new_res = new_chain->residues.back();
            new_res.subchain = new_chain->name + ":" + res.subchain;
            for (Atom& a : new_res.atoms)
              a.pos = Position(oper.transform.apply(a.pos));
          }
        }
      }
    }
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

  ChainNameGenerator::How how = ChainNameGenerator::How::AddNum;
  if (output_type == CoorFormat::Pdb)
    how = ChainNameGenerator::How::Short;
  if (options[ExpandAssembly]) {
    gemmi::Assembly* assembly = st.find_assembly(options[ExpandAssembly].arg);
    if (!assembly) {
      if (st.assemblies.empty())
        gemmi::fail("no bioassemblies are listed for this structure");
      gemmi::fail("wrong assembly name, use one of: " +
                  gemmi::join_str(st.assemblies, ' ',
                           [](const gemmi::Assembly& a) { return a.name; }));
    }
    for (gemmi::Model& model : st.models)
      expand_assembly(model, *assembly, how, options[Verbose]);
  }

  if (options[ExpandNcs]) {
    if (options[ExpandNcs].arg[0] == 'd')
     how = ChainNameGenerator::How::Dup;
    expand_ncs(st, how);
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
