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
#include "gemmi/assembly.hpp"  // for ChainNameGenerator, change_to_assembly

#include <cstring>
#include <iostream>
#include <algorithm>           // for sort
#include <map>

#define GEMMI_PROG convert
#include "options.h"
#include "cifmod.h"

namespace cif = gemmi::cif;
using gemmi::CoorFormat;
using gemmi::HowToNameCopiedChains;

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
  ExpandNcs, AsAssembly, RemoveH, RemoveWaters, RemoveLigWat, TrimAla,
  ShortTer, Linkr, MinimalPdb, ShortenCN, SegmentAsChain, OldPdb
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
  { OldPdb, 0, "", "old-pdb", Arg::None,
    "  --old-pdb \tRead only the first 72 characters in line." },

  { NoOp, 0, "", "", Arg::None, "\nPDB output options:" },
  { ShortTer, 0, "", "short-ter", Arg::None,
    "  --short-ter  \tWrite PDB TER records without numbers (iotbx compat.)." },
  { Linkr, 0, "", "linkr", Arg::None,
    "  --linkr  \tWrite LINKR record (for Refmac) if link_id is known." },
  { MinimalPdb, 0, "", "minimal-pdb", Arg::None,
    "  --minimal-pdb  \tWrite only the most essential records." },

  { NoOp, 0, "", "", Arg::None, "\nAny output options:" },
  { ShortenCN, 0, "", "shorten", Arg::None,
    "  --shorten  \tShorten chain names to 1 (if # < 63) or 2 characters." },

  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { ExpandNcs, 0, "", "expand-ncs", ConvArg::NcsChoice,
    "  --expand-ncs=dup|new  \tExpand strict NCS specified in MTRIXn or"
    " equivalent. New chain names are the same or have added numbers." },
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

  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

static void expand_ncs(gemmi::Structure& st, HowToNameCopiedChains how) {
  for (gemmi::Model& model : st.models) {
    size_t orig_size = model.chains.size();
    gemmi::ChainNameGenerator namegen(model, how);
    for (const gemmi::NcsOp& op : st.ncs)
      if (!op.given) {
        for (size_t i = 0; i != orig_size; ++i) {
          if (how == HowToNameCopiedChains::Dup)
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
            if (how == HowToNameCopiedChains::Dup)
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

static void rename_chain_in_address(gemmi::AtomAddress& aa,
                                    const std::string& old_name,
                                    const std::string& new_name) {
  if (aa.chain_name == old_name)
    aa.chain_name = new_name;
}

// chain is assumed to be from st.models[0]
static void rename_chain(gemmi::Structure& st, gemmi::Chain& chain,
                         const std::string& new_name) {
  for (gemmi::Connection& con : st.connections) {
    rename_chain_in_address(con.partner1, chain.name, new_name);
    rename_chain_in_address(con.partner2, chain.name, new_name);
  }
  for (gemmi::Helix& helix : st.helices) {
    rename_chain_in_address(helix.start, chain.name, new_name);
    rename_chain_in_address(helix.end, chain.name, new_name);
  }
  for (gemmi::Sheet& sheet : st.sheets)
    for (gemmi::Sheet::Strand& strand : sheet.strands) {
      rename_chain_in_address(strand.start, chain.name, new_name);
      rename_chain_in_address(strand.end, chain.name, new_name);
      rename_chain_in_address(strand.hbond_atom2, chain.name, new_name);
      rename_chain_in_address(strand.hbond_atom1, chain.name, new_name);
    }
  for (gemmi::RefinementInfo& ri : st.meta.refinement)
    for (gemmi::TlsGroup& tls : ri.tls_groups)
      for (gemmi::TlsGroup::Selection& sel : tls.selections)
        if (sel.chain == chain.name)
          sel.chain = new_name;
  for (auto it = st.models.begin() + 1; it != st.models.end(); ++it)
    if (gemmi::Chain* ch = it->find_chain(chain.name))
      ch->name = new_name;
  chain.name = new_name;
}

static void convert(gemmi::Structure& st,
                    const std::string& output, CoorFormat output_type,
                    const std::vector<option::Option>& options) {
  if (st.models.empty())
    gemmi::fail("No atoms in the input file. Wrong file format?");

  if (st.input_format == CoorFormat::Pdb) {
    gemmi::read_metadata_from_remarks(st);
    setup_for_mmcif(st);
  }

  if (options[ShortenCN]) {
    gemmi::ChainNameGenerator namegen(HowToNameCopiedChains::Short);
    gemmi::Model& model0 = st.models[0];
    size_t max_len = model0.chains.size() < 63 ? 1 : 2;
    for (const gemmi::Chain& chain : model0.chains)
      if (chain.name.length() <= max_len)
        namegen.used_names.push_back(chain.name);
    for (gemmi::Chain& chain : model0.chains)
      if (chain.name.length() > max_len)
        rename_chain(st, chain, namegen.make_short_name(
                                                chain.name.substr(0, max_len)));
  } else if (output_type == CoorFormat::Pdb) {
    for (const gemmi::Chain& chain : st.models[0].chains)
      if (chain.name.size() > 2)
        gemmi::fail("long chain name cannot be written in the PDB format: " +
                    chain.name + "\nTry option --shorten");
  }

  HowToNameCopiedChains how = HowToNameCopiedChains::AddNumber;
  if (output_type == CoorFormat::Pdb)
    how = HowToNameCopiedChains::Short;
  if (options[AsAssembly])
    gemmi::change_to_assembly(st, options[AsAssembly].arg, how,
                              options[Verbose] ? &std::cerr : nullptr);

  if (options[ExpandNcs]) {
    if (options[ExpandNcs].arg[0] == 'd')
     how = HowToNameCopiedChains::Dup;
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
    if (options[MinimalPdb])
      gemmi::write_minimal_pdb(st, os.ref(), opt);
    else
      gemmi::write_pdb(st, os.ref(), opt);
  }
}

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
    if (p.options[OldPdb])
      st = gemmi::read_pdb_gz(input, 72);
    else
      st = gemmi::read_structure_gz(input, in_type);
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
