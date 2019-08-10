// Copyright 2017 Global Phasing Ltd.

#include "gemmi/gzread.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/polyheur.hpp"  // for remove_hydrogens, ...
#include "gemmi/to_pdb.hpp"    // for write_pdb, ...
#include "gemmi/to_mmcif.hpp"  // for update_cif_block
#include "gemmi/chemcomp_xyz.hpp" // for make_structure_from_chemcomp_block
#include "gemmi/remarks.hpp"   // for read_metadata_from_remarks

#include <cstring>
#include <iostream>
#include <map>

#define GEMMI_PROG convert
#include "options.h"

namespace cif = gemmi::cif;
using gemmi::CoorFormat;

struct ConvArg: public Arg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"json", "pdb", "cif"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }

  static option::ArgStatus NcsChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"dup", "addnum"});
  }
};

enum OptionIndex { Verbose=3, FormatIn, FormatOut,
                   Comcifs, Mmjson, Bare, Numb, CifDot, PdbxStyle, SkipCat,
                   ExpandNcs, RemoveH, RemoveWaters, RemoveLigWat, TrimAla,
                   ShortTer, SegmentAsChain, Translate };
static const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions CIF-JSON, and mmCIF-PDB-mmJSON."
    "\nFORMAT can be specified as one of: cif, json, pdb."
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { FormatIn, 0, "", "from", ConvArg::FileFormat,
    "  --from=FORMAT  \tInput format (default: from the file extension)." },
  { FormatOut, 0, "", "to", ConvArg::FileFormat,
    "  --to=FORMAT  \tOutput format (default: from the file extension)." },
  { NoOp, 0, "", "", Arg::None, "\nJSON output options:" },
  { Comcifs, 0, "c", "comcifs", Arg::None,
    "  -c, --comcifs  \tConform to the COMCIFS CIF-JSON standard draft." },
  { Mmjson, 0, "m", "mmjson", Arg::None,
    "  -m, --mmjson   \tCompatible with mmJSON from PDBj." },
  { Bare, 0, "b", "bare-tags", Arg::None,
    "  -b, --bare-tags  \tOutput tags without the first underscore." },
  { Numb, 0, "", "numb", ConvArg::NumbChoice,
    "  --numb=quote|nosu|mix  \tConvert the CIF numb type to one of:"
                             "\v  quote - string in quotes,"
                             "\v  nosu - number without s.u.,"
                             "\v  mix (default) - quote only numbs with s.u." },
  { CifDot, 0, "", "dot", Arg::Required,
    "  --dot=STRING  \tJSON representation of CIF's '.' (default: null)." },
  { NoOp, 0, "", "", Arg::None, "\nCIF output options:" },
  { PdbxStyle, 0, "", "pdbx-style", Arg::None,
    "  --pdbx-style  \tSimilar styling (formatting) as in wwPDB." },
  { SkipCat, 0, "", "skip-category", Arg::Required,
    "  --skip-category=CAT  \tDo not output tags starting with _CAT" },
  { NoOp, 0, "", "", Arg::None, "\nPDB input/output options:" },
  { SegmentAsChain, 0, "", "segment-as-chain", Arg::None,
    "  --segment-as-chain \tAppend segment id to label_asym_id (chain name)." },
  { ShortTer, 0, "", "short-ter", Arg::None,
    "  --short-ter  \tWrite PDB TER records without numbers (iotbx compat.)." },
  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { ExpandNcs, 0, "", "expand-ncs", ConvArg::NcsChoice,
    "  --expand-ncs=dup|addn  \tExpand strict NCS specified in MTRIXn or"
    " equivalent. New chain names are the same or have added numbers." },
  { RemoveH, 0, "", "remove-h", Arg::None,
    "  --remove-h  \tRemove hydrogens." },
  { RemoveWaters, 0, "", "remove-waters", Arg::None,
    "  --remove-waters  \tRemove waters." },
  { RemoveLigWat, 0, "", "remove-lig-wat", Arg::None,
    "  --remove-lig-wat  \tRemove ligands and waters." },
  { TrimAla, 0, "", "trim-to-ala", Arg::None,
    "  --trim-to-ala  \tTrim aminoacids to alanine." },
  { Translate, 0, "", "translate", Arg::None, 0 },
  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

static const char symbols[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "abcdefghijklmnopqrstuvwxyz0123456789";

enum class ChainNaming { Short, AddNum, Dup };

static bool is_mmcif_compatible(CoorFormat format) {
  return format == CoorFormat::Mmcif || format == CoorFormat::Mmjson;
}

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

static void convert(const std::string& input, CoorFormat input_type,
                    const std::string& output, CoorFormat output_type,
                    const std::vector<option::Option>& options,
                    bool transcribe) {
  cif::Document doc;
  gemmi::Structure st;
  // for cif->cif we do either cif->DOM->Structure->DOM->cif or cif->DOM->cif
  if (input_type == CoorFormat::Mmcif || input_type == CoorFormat::Mmjson) {
    doc = input_type == CoorFormat::Mmcif ? gemmi::read_cif_gz(input)
                                          : gemmi::read_mmjson_gz(input);
    if (!transcribe) {
      int n = gemmi::check_chemcomp_block_number(doc);
      // first handle special case - refmac dictionary or CCD file
      if (n != -1)
        st = gemmi::make_structure_from_chemcomp_block(doc.blocks[n]);
      else
        st = gemmi::make_structure(doc);
      if (st.models.empty())
        gemmi::fail("No atoms in the input file. Is it mmCIF?");
    }
  } else if (input_type == CoorFormat::Pdb) {
    st = gemmi::read_pdb_gz(input);
    gemmi::read_metadata_from_remarks(st);
    setup_entities(st);
  } else {
    gemmi::fail("Unexpected input format.");
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

  std::ostream* os;
  std::unique_ptr<std::ostream> os_deleter;
  if (output != "-") {
    os_deleter.reset(new std::ofstream(output));
    os = os_deleter.get();
    if (!os || !*os)
      gemmi::fail("Failed to open for writing: " + output);
  } else {
    os = &std::cout;
  }

  if (output_type == CoorFormat::Mmcif || output_type == CoorFormat::Mmjson) {
    if (!transcribe) {
      doc.blocks.clear();  // temporary, for testing
      doc.blocks.resize(1);
      update_cif_block(st, doc.blocks[0]);
    }
    for (const option::Option* opt = options[SkipCat]; opt; opt = opt->next()) {
      std::string category = opt->arg;
      if (category[0] != '_')
        category.insert(0, 1, '_');
      for (cif::Block& block : doc.blocks)
        for (cif::Item& item : block.items)
          if (item.has_prefix(category))
            item.erase();
    }
    if (output_type == CoorFormat::Mmcif) {
      auto style = options[PdbxStyle] ? cif::Style::Pdbx
                                      : cif::Style::PreferPairs;
      write_cif_to_stream(*os, doc, style);
    } else /*output_type == CoorFormat::Mmjson*/ {
      cif::JsonWriter writer(*os);
      if (options[Comcifs])
        writer.set_comcifs();
      if (options[Mmjson])
        writer.set_mmjson();
      if (options[Bare])
        writer.bare_tags = true;
      if (options[Numb]) {
        char first_letter = options[Numb].arg[0];
        if (first_letter == 'q')
          writer.quote_numbers = 2;
        else if (first_letter == 'n')
          writer.quote_numbers = 0;
      }
      if (options[CifDot])
        writer.cif_dot = options[CifDot].arg;
      writer.write_json(doc);
    }
  } else if (output_type == CoorFormat::Pdb) {
    // call wrapper from output.cpp - to make building faster
    gemmi::PdbWriteOptions opt;
    if (options[ShortTer])
      opt.numbered_ter = false;
    gemmi::write_pdb(st, *os, opt);
  }
}

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  std::string input = p.coordinate_input_file(0);
  const char* output = p.nonOption(1);

  // CoorFormat::Mmcif here stands for any CIF files,
  // CoorFormat::Mmjson may not be strictly mmJSON, but also CIF-JSON.
  std::map<std::string, CoorFormat> filetypes {{"json", CoorFormat::Mmjson},
                                               {"pdb", CoorFormat::Pdb},
                                               {"cif", CoorFormat::Mmcif}};

  CoorFormat in_type = p.options[FormatIn]
    ? filetypes[p.options[FormatIn].arg]
    : gemmi::coordinate_format_from_extension_gz(input);
  if (in_type == CoorFormat::Unknown) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

  CoorFormat out_type = p.options[FormatOut]
    ? filetypes[p.options[FormatOut].arg]
    : gemmi::coordinate_format_from_extension_gz(output);
  if (out_type == CoorFormat::Unknown) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }
  bool transcribe = is_mmcif_compatible(in_type) &&
                    is_mmcif_compatible(out_type) &&
                    !(p.options[Translate] || p.options[ExpandNcs] ||
                      p.options[RemoveH] || p.options[RemoveWaters] ||
                      p.options[RemoveLigWat] || p.options[TrimAla] ||
                      p.options[SegmentAsChain]);
  if (p.options[Verbose])
    std::cerr << (transcribe ? "Transcribing " : "Converting ")
              << input << " to " << format_as_string(out_type) << "..."
              << std::endl;
  try {
    convert(input, in_type, output, out_type, p.options, transcribe);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  }
  if (p.options[Verbose])
    std::cerr << "Done." << std::endl;
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
