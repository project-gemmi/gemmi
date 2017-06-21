// Copyright 2017 Global Phasing Ltd.

#include <iostream> // temporary, for debugging
#include "gemmi/cifgz.hpp"
#include "gemmi/mmcif.hpp"
#include "gemmi/pdbgz.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/version.hpp"
// set this before only one of stb_sprintf.h includes
#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/to_mmcif.hpp"

#include <cstring>
#include <iostream>
#include <map>
#include <optionparser.h>

#define EXE_NAME "gemmi-convert"

namespace mol = gemmi::mol;
namespace cif = gemmi::cif;

enum class FileType : char { Json, Pdb, Cif, Crd, Null, Unknown };

struct Arg: public option::Arg {
  static option::ArgStatus Required(const option::Option& option, bool msg) {
    if (option.arg != nullptr)
      return option::ARG_OK;
    if (msg)
      std::cerr << "Option '" << option.name << "' requires an argument\n";
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Choice(const option::Option& option, bool msg,
                                  std::vector<const char*> choices) {
    if (Required(option, msg) == option::ARG_ILLEGAL)
      return option::ARG_ILLEGAL;
    for (const char* a : choices)
      if (strcmp(option.arg, a) == 0)
        return option::ARG_OK;
    if (msg)
      std::cerr << "Invalid argument for "
                << std::string(option.name, option.namelen) << ": "
                << option.arg << "\n";
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    // the hidden option "none" is for testing only
    return Arg::Choice(option, msg, {"json", "pdb", "cif", "none"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }
};

enum OptionIndex { Unknown, Help, Version, Verbose, FormatIn, FormatOut,
                   Comcifs, Bare, Numb, CifDot,
                   ExpandNcs, IotbxCompat, SegmentAsChain };
static const option::Descriptor Usage[] = {
  { Unknown, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions: cif->json and cif<->pdb."
    "\n\nGeneral options:" },
  { Help, 0, "h", "help", Arg::None, "  -h, --help  \tPrint usage and exit." },
  { Version, 0, "V", "version", Arg::None,
    "  -V, --version  \tPrint version and exit." },
  { Verbose, 0, "", "verbose", Arg::None, "  --verbose  \tVerbose output." },
  { FormatIn, 0, "", "from", Arg::FileFormat,
    "  --from=pdb|cif  \tInput format (default: from the file extension)." },
  { FormatOut, 0, "", "to", Arg::FileFormat,
    "  --to=json|pdb  \tOutput format (default: from the file extension)." },
  { Unknown, 0, "", "", Arg::None, "\nCIF output options:" },
  { Comcifs, 0, "c", "comcifs", Arg::None,
    "  -c, --comcifs  \tConform to the COMCIFS CIF-JSON standard draft." },
  { Bare, 0, "b", "bare-tags", Arg::None,
    "  -b, --bare-tags  \tOutput tags without the first underscore." },
  { Numb, 0, "", "numb", Arg::NumbChoice,
    "  --numb=quote|nosu|mix  \tConvert the CIF numb type to one of:"
                             "\v  quote - string in quotes,"
                             "\v  nosu - number without s.u.,"
                             "\v  mix (default) - quote only numbs with s.u." },
  { CifDot, 0, "", "dot", Arg::Required,
    "  --dot=STRING  \tJSON representation of CIF's '.' (default: null)." },
  { Unknown, 0, "", "", Arg::None, "\nMacromolecular options:" },
  { ExpandNcs, 0, "", "expand-ncs", Arg::None,
    "  --expand-ncs  \tExpand strict NCS specified in MTRIXn or equivalent." },
  { IotbxCompat, 0, "", "iotbx-compat", Arg::None,
    "  --iotbx-compat  \tLimited compatibility with iotbx (details in docs)." },
  { SegmentAsChain, 0, "", "segment-as-chain", Arg::None,
    "  --segment-as-chain \tAppend segment id to label_asym_id (chain name)." },
  { Unknown, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

FileType get_format_from_extension(const std::string& path) {
  using gemmi::iends_with;
  if (iends_with(path, ".pdb") || iends_with(path, ".ent") ||
      iends_with(path, ".pdb.gz") || iends_with(path, ".ent.gz"))
    return FileType::Pdb;
  if (iends_with(path, ".js") || iends_with(path, ".json"))
    return FileType::Json;
  if (iends_with(path, ".cif") || iends_with(path, ".cif.gz"))
    return FileType::Cif;
  if (iends_with(path, ".crd"))
    return FileType::Crd;
  if (path == "/dev/null")
    return FileType::Null;
  return FileType::Unknown;
}

[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

static const char symbols[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "abcdefghijklmnopqrstuvwxyz01234567890";

enum class ChainNaming { Short, AddNum, Dup };

static void expand_ncs(mol::Structure& st, ChainNaming ch_naming) {
  int n_ops = std::count_if(st.ncs.begin(), st.ncs.end(),
                            [](mol::NcsOp& op){ return !op.given; });
  if (n_ops == 0)
    return;
  for (mol::Model& model : st.models) {
    size_t new_length = model.chains.size() * (n_ops + 1);
    if (new_length >= 63 && ch_naming == ChainNaming::Short)
      ch_naming = ChainNaming::AddNum;
    model.chains.reserve(new_length);
    // for compatibility with iotbx we reset serial in the 'Dup' mode
    if (ch_naming == ChainNaming::Dup && !model.chains.empty()) {
      model.chains[0].force_pdb_serial = 1;
      for (mol::Chain& chain : model.chains)
        for (mol::Residue& res : chain.residues)
          res.segment = "0";
    }
    auto orig_end = model.chains.cend();
    for (const mol::NcsOp& op : st.ncs)
      if (!op.given)
        for (auto ch = model.chains.cbegin(); ch != orig_end; ++ch) {
          // the most difficult part - choosing names for new chains
          std::string name;
          if (ch_naming == ChainNaming::Short) {
            for (char symbol : symbols) {
              name = std::string(1, symbol);
              if (!model.find_chain(name))
                break;
            }
          } else if (ch_naming == ChainNaming::AddNum) {
            name = ch->name + op.id;
            while (model.find_chain(name))
              name += "a";
          } else { // ChainNaming::Dup
            name = ch->name;
          }

          model.chains.push_back(*ch);
          mol::Chain& new_chain = model.chains.back();
          new_chain.name = name;
          new_chain.auth_name = "";

          for (mol::Residue& res : new_chain.residues) {
            for (mol::Atom& a : res.atoms) {
              linalg::vec<double,4> pos = {a.pos.x, a.pos.y, a.pos.z, 1.0};
              pos = linalg::mul(op.transform, pos);
              a.pos = {pos.x, pos.y, pos.z};
            }
            if (ch_naming == ChainNaming::Dup)
              res.segment = op.id;
          }
        }
  }
  for (mol::NcsOp& op : st.ncs)
    op.given = true;
}


std::vector<mol::Chain> split_by_segments(mol::Chain& orig) {
  std::vector<mol::Chain> chains;
  std::vector<mol::Residue> orig_res;
  orig_res.swap(orig.residues);
  for (auto seg_start = orig_res.begin(); seg_start != orig_res.end(); ) {
    const std::string& sn = seg_start->segment;
    auto ch = std::find_if(chains.begin(), chains.end(), [&](mol::Chain& c) {
                return !c.residues.empty() && c.residues[0].segment == sn; });
    if (ch == chains.end()) {
      chains.push_back(orig);
      ch = chains.end() - 1;
      // it's not clear how chain naming should be handled
      ch->name += sn;
      ch->auth_name += sn;
    }
    auto seg_end = std::find_if(seg_start, orig_res.end(),
                            [&](mol::Residue& r) { return r.segment != sn; });
    ch->residues.insert(ch->residues.end(), std::make_move_iterator(seg_start),
                                            std::make_move_iterator(seg_end));
    seg_start = seg_end;
  }
  return chains;
}

// for Refmac, to be merged with update_cif_block()
cif::Document make_crd(const mol::Structure& st) {
  cif::Document crd;
  if (st.models.empty())
    return crd;
  auto e_id = st.info.find("_entry.id");
  std::string id = (e_id != st.info.end() ? e_id->second : st.name);
  crd.blocks.emplace_back("structure_" + id);
  cif::Block& block = crd.blocks[0];
  auto& items = block.items;

  items.emplace_back("_entry.id", id);
  items.emplace_back("_database_2.code_PDB", id);
  auto keywords = st.info.find("_struct_keywords.pdbx_keywords");
  if (keywords != st.info.end())
    items.emplace_back("_struct_keywords.text", cif::quote(keywords->second));
  auto title = st.info.find("_struct.title");
  if (title != st.info.end())
    items.emplace_back(title->first, cif::quote(title->second));
  auto initial_date =
         st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
  if (initial_date != st.info.end())
    items.emplace_back("_audit.creation_date", initial_date->second);

  items.emplace_back(cif::CommentArg{"############\n"
                                     "## ENTITY ##\n"
                                     "############"});
  cif::Loop& entity_loop = block.clear_or_add_loop("_entity.", {"id", "type"});
  for (const auto& ent : st.entities) {
    entity_loop.values.push_back(ent->id);
    entity_loop.values.push_back(ent->type_as_string());
  }
  items.emplace_back(cif::CommentArg{"#####################\n"
                                     "## ENTITY_POLY_SEQ ##\n"
                                     "#####################"});
  cif::Loop& poly_loop = block.clear_or_add_loop("_entity_poly_seq.", {
              "mon_id", "ccp4_auth_seq_id", "entity_id",
              "ccp4_back_connect_type", "ccp4_num_mon_back", "ccp4_mod_id"});
  for (const auto& ent : st.entities)
    if (ent->type == mol::EntityType::Polymer)
      for (const mol::SequenceItem& si : ent->sequence) {
        poly_loop.values.emplace_back(si.mon);
        // TODO: real auth_seq_id
        std::string auth_seq_id = si.num >= 0 ? std::to_string(si.num) : "?";
        poly_loop.values.emplace_back(auth_seq_id);
        poly_loop.values.emplace_back(ent->id);
        poly_loop.values.emplace_back("?"); // ccp4_back_connect_type
        poly_loop.values.emplace_back("?"); // ccp4_num_mon_back
        poly_loop.values.emplace_back("?"); // ccp4_mod_id
      }
  items.emplace_back(cif::CommentArg{"##########\n"
                                     "## CELL ##\n"
                                     "##########"});
  items.emplace_back("_cell.entry_id", id);
  items.emplace_back("_cell.length_a",    mol::to_str(st.cell.a));
  items.emplace_back("_cell.length_b",    mol::to_str(st.cell.b));
  items.emplace_back("_cell.length_c",    mol::to_str(st.cell.c));
  items.emplace_back("_cell.angle_alpha", mol::to_str(st.cell.alpha));
  items.emplace_back("_cell.angle_beta",  mol::to_str(st.cell.beta));
  items.emplace_back("_cell.angle_gamma", mol::to_str(st.cell.gamma));
  items.emplace_back(cif::CommentArg{"##############################\n"
                                     "## FRACTIONALISATION MATRIX ##\n"
                                     "##############################"});
  if (st.cell.explicit_matrices || st.cell.frac.a11 != 1.0) {
    std::string prefix = "_atom_sites.fract_transf_";
    items.emplace_back(prefix + "matrix[1][1]", mol::to_str(st.cell.frac.a11));
    items.emplace_back(prefix + "matrix[1][2]", mol::to_str(st.cell.frac.a12));
    items.emplace_back(prefix + "matrix[1][3]", mol::to_str(st.cell.frac.a13));
    items.emplace_back(prefix + "matrix[2][1]", mol::to_str(st.cell.frac.a21));
    items.emplace_back(prefix + "matrix[2][2]", mol::to_str(st.cell.frac.a22));
    items.emplace_back(prefix + "matrix[2][3]", mol::to_str(st.cell.frac.a23));
    items.emplace_back(prefix + "matrix[3][1]", mol::to_str(st.cell.frac.a31));
    items.emplace_back(prefix + "matrix[3][2]", mol::to_str(st.cell.frac.a32));
    items.emplace_back(prefix + "matrix[3][3]", mol::to_str(st.cell.frac.a33));
    items.emplace_back(prefix + "vector[1]",    mol::to_str(st.cell.shift.x));
    items.emplace_back(prefix + "vector[2]",    mol::to_str(st.cell.shift.y));
    items.emplace_back(prefix + "vector[3]",    mol::to_str(st.cell.shift.z));
  }
  items.emplace_back(cif::CommentArg{"##############\n"
                                     "## SYMMETRY ##\n"
                                     "##############"});
  items.emplace_back("_symmetry.entry_id", id);
  items.emplace_back("_symmetry.space_group_name_H-M", cif::quote(st.sg_hm));
  items.emplace_back(cif::CommentArg{"#################\n"
                                     "## STRUCT_ASYM ##\n"
                                     "#################"});
  // _struct_asym
  cif::Loop& asym_loop = block.clear_or_add_loop("_struct_asym.",
                                                 {"id", "entity_id"});
  for (const auto& ch : st.get_chains()) {
    asym_loop.values.push_back(ch.name);
    asym_loop.values.push_back(ch.entity ? ch.entity->id : "?");
  }
  items.emplace_back(cif::CommentArg{"###############\n"
                                     "## ATOM_SITE ##\n"
                                     "###############"});
  cif::Loop& atom_loop = block.clear_or_add_loop("_atom_site.", {
      "group_PDB",
      "id",
      "label_atom_id",
      "label_alt_id",
      "label_comp_id",
      "label_asym_id",
      "auth_seq_id",
      //"pdbx_PDB_ins_code",
      "Cartn_x",
      "Cartn_y",
      "Cartn_z",
      "occupancy",
      "B_iso_or_equiv",
      "type_symbol",
      "calc_flag",
      "label_seg_id",
      "auth_atom_id",
      "label_chem_id"});
  std::vector<std::string>& vv = atom_loop.values;
  vv.reserve(count_atom_sites(st) * atom_loop.tags.size());
  int serial = 0;
  for (const mol::Model& model : st.models) {
    for (const mol::Chain& chain : model.chains) {
      for (const mol::Residue& res : chain.residues) {
        //std::string seq_id = std::to_string(res.seq_id);
        std::string auth_seq_id = std::to_string(res.auth_seq_id);
        //std::string ins_code(1, res.ins_code ? res.ins_code : '?');
        for (const mol::Atom& a : res.atoms) {
          vv.emplace_back("ATOM");
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc ? a.altloc : '.');
          vv.emplace_back(res.name);
          vv.emplace_back(chain.name);
          vv.emplace_back(auth_seq_id);
          //vv.emplace_back(ins_code);
          vv.emplace_back(mol::to_str(a.pos.x));
          vv.emplace_back(mol::to_str(a.pos.y));
          vv.emplace_back(mol::to_str(a.pos.z));
          vv.emplace_back(mol::to_str(a.occ));
          vv.emplace_back(mol::to_str(a.b_iso));
          vv.emplace_back(a.element.uname());
          vv.emplace_back("."); // calc_flag
          vv.emplace_back("."); // label_seg_id
          vv.emplace_back(a.name); // again
          vv.emplace_back("?"); // label_chem_id
        }
      }
    }
  }
  return crd;
}

void convert(const char* input, FileType input_type,
             const char* output, FileType output_type,
             const std::vector<option::Option>& options) {
  cif::Document cif_in;
  mol::Structure st;
  // for cif->cif we do either cif->DOM->Structure->DOM->cif or cif->DOM->cif
  bool modify_structure = (options[ExpandNcs] || options[SegmentAsChain]);
  if (input_type == FileType::Cif) {
    cif_in = cif::read_any(input);
    if ((output_type == FileType::Json || output_type == FileType::Cif) &&
        !modify_structure) {
      // no need to interpret the structure
    } else {
      st = mol::read_atoms(cif_in);
      if (st.models.empty())
        fail("No atoms in the input file. Is it mmCIF?");
    }
  } else if (input_type == FileType::Pdb) {
    st = mol::read_pdb_any(input);
  } else {
    fail("Unexpected input format.");
  }

  if (options[ExpandNcs]) {
    ChainNaming ch_naming = ChainNaming::AddNum;
    if (options[IotbxCompat])
     ch_naming = ChainNaming::Dup;
    else if (output_type == FileType::Pdb)
      ch_naming = ChainNaming::Short;
    expand_ncs(st, ch_naming);
  }

  if (options[SegmentAsChain])
    for (mol::Model& model : st.models) {
      std::vector<mol::Chain> new_chains;
      for (mol::Chain& chain : model.chains)
        gemmi::vector_move_extend(new_chains, split_by_segments(chain));
      model.chains = std::move(new_chains);
      // We should also modify Entity instances, but for now we don't.
      add_backlinks(model); // fix backlinks from Residue/Atom
    }

  std::ostream* os;
  std::unique_ptr<std::ostream> os_deleter;
  if (output != std::string("-")) {
    os_deleter.reset(new std::ofstream(output));
    os = os_deleter.get();
    if (!os || !*os)
      fail("Failed to open for writing: " + std::string(output));
  } else {
    os = &std::cout;
  }

  if (output_type == FileType::Json) {
    if (input_type != FileType::Cif)
      fail("Conversion to JSON is possible only from CIF");
    cif::JsonWriter writer(*os);
    if (options[Comcifs])
      writer.set_comcifs();
    writer.use_bare_tags = options[Bare];
    if (options[Numb]) {
      char first_letter = options[Numb].arg[0];
      if (first_letter == 'q')
        writer.quote_numbers = 2;
      else if (first_letter == 'n')
        writer.quote_numbers = 0;
    }
    if (options[CifDot])
      writer.cif_dot = options[CifDot].arg;
    writer.write_json(cif_in);
  } else if (output_type == FileType::Pdb || output_type == FileType::Null) {
    if (output_type == FileType::Pdb)
      mol::write_pdb(st, *os, options[IotbxCompat]);
    else {
      *os << st.name << ": " << count_atom_sites(st) << " atom locations";
      if (st.models.size() > 1)
        *os << " (total in " << st.models.size() << " models)";
      *os << ".\n";
    }
  } else if (output_type == FileType::Cif) {
    // cif to cif round trip is for testing only
    if (input_type != FileType::Cif || modify_structure) {
      cif_in.blocks.clear();  // temporary, for testing
      cif_in.blocks.resize(1);
      mol::update_cif_block(st, cif_in.blocks[0]);
    }
    *os << cif_in;
  } else if (output_type == FileType::Crd) {
    *os << make_crd(st);
  }
}

int main(int argc, char **argv) {
  if (argc < 1)
    return 2;
  std::ios_base::sync_with_stdio(false);
  option::Stats stats(Usage, argc-1, argv+1);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(Usage, argc-1, argv+1, options.data(), buffer.data());
  if (parse.error()) {
    option::printUsage(std::cerr, Usage);
    return 1;
  }
  if (options[Help]) {
    option::printUsage(std::cout, Usage);
    return 0;
  }
  if (options[Version]) {
    std::cout << EXE_NAME " " GEMMI_VERSION "\n";
    return 0;
  }
  if (options[Unknown]) {
    std::cerr << "Invalid option.\n";
    option::printUsage(std::cerr, Usage);
    return 1;
  }
  if (parse.nonOptionsCount() != 2) {
    std::cerr << "This program requires 2 arguments (input and output), "
              << parse.nonOptionsCount() << " given.\n"
                 "Try '" EXE_NAME " --help' for more information.\n";
    return 1;
  }

  const char* input = parse.nonOption(0);
  const char* output = parse.nonOption(1);

  std::map<std::string, FileType> filetypes {{"json", FileType::Json},
                                             {"pdb", FileType::Pdb},
                                             {"cif", FileType::Cif},
                                             {"crd", FileType::Crd},
                                             {"none", FileType::Null}};

  FileType in_type = options[FormatIn] ? filetypes[options[FormatIn].arg]
                                       : get_format_from_extension(input);
  if (in_type == FileType::Unknown) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

  FileType out_type = options[FormatOut] ? filetypes[options[FormatOut].arg]
                                         : get_format_from_extension(output);
  if (out_type == FileType::Unknown) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }
  if (options[Verbose])
    std::cerr << "Converting " << input << " ..." << std::endl;

  try {
    convert(input, in_type, output, out_type, options);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=include,third_party
