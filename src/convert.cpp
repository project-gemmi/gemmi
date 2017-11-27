// Copyright 2017 Global Phasing Ltd.

#include "input.h"
#include "output.h"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/sprintf.hpp"

#include <cstring>
#include <iostream>
#include <map>

#define EXE_NAME "gemmi-convert"
#include "options.h"

namespace cif = gemmi::cif;

enum class FileType : char { Json, Pdb, Cif, Crd, Null, Unknown };

struct ConvArg: public Arg {
  static option::ArgStatus FileFormat(const option::Option& option, bool msg) {
    // the hidden option "none" is for testing only
    return Arg::Choice(option, msg, {"json", "pdb", "cif", "none"});
  }

  static option::ArgStatus NumbChoice(const option::Option& option, bool msg) {
    return Arg::Choice(option, msg, {"quote", "nosu", "mix"});
  }
};

enum OptionIndex { Verbose=3, FormatIn, FormatOut,
                   Comcifs, Mmjson, Bare, Numb, CifDot,
                   ExpandNcs, IotbxCompat, SegmentAsChain };
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
  { NoOp, 0, "", "", Arg::None, "\nCIF output options:" },
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
  { NoOp, 0, "", "", Arg::None, "\nMacromolecular options:" },
  { ExpandNcs, 0, "", "expand-ncs", Arg::None,
    "  --expand-ncs  \tExpand strict NCS specified in MTRIXn or equivalent." },
  { IotbxCompat, 0, "", "iotbx-compat", Arg::None,
    "  --iotbx-compat  \tLimited compatibility with iotbx (details in docs)." },
  { SegmentAsChain, 0, "", "segment-as-chain", Arg::None,
    "  --segment-as-chain \tAppend segment id to label_asym_id (chain name)." },
  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output." },
  { 0, 0, 0, 0, 0, 0 }
};

FileType get_format_from_extension(const std::string& path) {
  gemmi::CoorFormat format = coordinate_format_from_extension(path);
  if (format == gemmi::CoorFormat::Pdb)
    return FileType::Pdb;
  if (format == gemmi::CoorFormat::Json)
    return FileType::Json;
  if (format == gemmi::CoorFormat::Cif)
    return FileType::Cif;
  if (gemmi::iends_with(path, ".crd"))
    return FileType::Crd;
  if (path == "/dev/null")
    return FileType::Null;
  return FileType::Unknown;
}

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
    // for compatibility with iotbx we reset serial in the 'Dup' mode
    if (ch_naming == ChainNaming::Dup && !model.chains.empty()) {
      model.chains[0].force_pdb_serial = 1;
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
          res.segment = "0";
    }
    auto orig_end = model.chains.cend();
    for (const gemmi::NcsOp& op : st.ncs)
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
          gemmi::Chain& new_chain = model.chains.back();
          new_chain.name = name;
          new_chain.auth_name = "";

          for (gemmi::Residue& res : new_chain.residues) {
            for (gemmi::Atom& a : res.atoms) {
              linalg::vec<double,4> pos = {a.pos.x, a.pos.y, a.pos.z, 1.0};
              pos = linalg::mul(op.transform, pos);
              a.pos = {pos.x, pos.y, pos.z};
            }
            if (ch_naming == ChainNaming::Dup)
              res.segment = op.id;
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
      ch->auth_name += sn;
    }
    auto seg_end = std::find_if(seg_start, orig_res.end(),
                            [&](gemmi::Residue& r) { return r.segment != sn; });
    ch->residues.insert(ch->residues.end(), std::make_move_iterator(seg_start),
                                            std::make_move_iterator(seg_end));
    seg_start = seg_end;
  }
  return chains;
}

// for Refmac, to be merged with update_cif_block()
cif::Document make_crd(const gemmi::Structure& st) {
  using gemmi::to_str;
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
  cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
  for (const auto& ent : st.entities) {
    entity_loop.values.push_back(ent->id);
    entity_loop.values.push_back(ent->type_as_string());
  }
  items.emplace_back(cif::CommentArg{"#####################\n"
                                     "## ENTITY_POLY_SEQ ##\n"
                                     "#####################"});
  cif::Loop& poly_loop = block.init_mmcif_loop("_entity_poly_seq.", {
              "mon_id", "ccp4_auth_seq_id", "entity_id",
              "ccp4_back_connect_type", "ccp4_num_mon_back", "ccp4_mod_id"});
  for (const auto& ent : st.entities)
    if (ent->type == gemmi::EntityType::Polymer)
      for (const gemmi::SequenceItem& si : ent->sequence) {
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
  items.emplace_back("_cell.length_a",    to_str(st.cell.a));
  items.emplace_back("_cell.length_b",    to_str(st.cell.b));
  items.emplace_back("_cell.length_c",    to_str(st.cell.c));
  items.emplace_back("_cell.angle_alpha", to_str(st.cell.alpha));
  items.emplace_back("_cell.angle_beta",  to_str(st.cell.beta));
  items.emplace_back("_cell.angle_gamma", to_str(st.cell.gamma));
  items.emplace_back(cif::CommentArg{"##############################\n"
                                     "## FRACTIONALISATION MATRIX ##\n"
                                     "##############################"});
  if (st.cell.explicit_matrices || st.cell.frac.a11 != 1.0) {
    std::string prefix = "_atom_sites.fract_transf_";
    items.emplace_back(prefix + "matrix[1][1]", to_str(st.cell.frac.a11));
    items.emplace_back(prefix + "matrix[1][2]", to_str(st.cell.frac.a12));
    items.emplace_back(prefix + "matrix[1][3]", to_str(st.cell.frac.a13));
    items.emplace_back(prefix + "matrix[2][1]", to_str(st.cell.frac.a21));
    items.emplace_back(prefix + "matrix[2][2]", to_str(st.cell.frac.a22));
    items.emplace_back(prefix + "matrix[2][3]", to_str(st.cell.frac.a23));
    items.emplace_back(prefix + "matrix[3][1]", to_str(st.cell.frac.a31));
    items.emplace_back(prefix + "matrix[3][2]", to_str(st.cell.frac.a32));
    items.emplace_back(prefix + "matrix[3][3]", to_str(st.cell.frac.a33));
    items.emplace_back(prefix + "vector[1]",    to_str(st.cell.shift.x));
    items.emplace_back(prefix + "vector[2]",    to_str(st.cell.shift.y));
    items.emplace_back(prefix + "vector[3]",    to_str(st.cell.shift.z));
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
  cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
                                               {"id", "entity_id"});
  for (const auto& ch : st.get_chains()) {
    asym_loop.values.push_back(ch.name);
    asym_loop.values.push_back(ch.entity ? ch.entity->id : "?");
  }
  items.emplace_back(cif::CommentArg{"###############\n"
                                     "## ATOM_SITE ##\n"
                                     "###############"});
  cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
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
  for (const gemmi::Model& model : st.models) {
    for (const gemmi::Chain& chain : model.chains) {
      for (const gemmi::Residue& res : chain.residues) {
        //std::string seq_id = std::to_string(res.seq_id);
        std::string auth_seq_id = std::to_string(res.snic.seq_num);
        //std::string ins_code(1, res.snic.ins_code ? res.snic.ins_code : '?');
        for (const gemmi::Atom& a : res.atoms) {
          vv.emplace_back("ATOM");
          vv.emplace_back(std::to_string(++serial));
          vv.emplace_back(a.name);
          vv.emplace_back(1, a.altloc ? a.altloc : '.');
          vv.emplace_back(res.name);
          vv.emplace_back(chain.name);
          vv.emplace_back(auth_seq_id);
          //vv.emplace_back(ins_code);
          vv.emplace_back(to_str(a.pos.x));
          vv.emplace_back(to_str(a.pos.y));
          vv.emplace_back(to_str(a.pos.z));
          vv.emplace_back(to_str(a.occ));
          vv.emplace_back(to_str(a.b_iso));
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

void convert(const std::string& input, FileType input_type,
             const std::string& output, FileType output_type,
             const std::vector<option::Option>& options) {
  cif::Document cif_in;
  gemmi::Structure st;
  // for cif->cif we do either cif->DOM->Structure->DOM->cif or cif->DOM->cif
  bool modify_structure = (options[ExpandNcs] || options[SegmentAsChain]);
  if (input_type == FileType::Cif || input_type == FileType::Json) {
    cif_in = cif_read_any(input);
    if ((output_type == FileType::Json || output_type == FileType::Cif) &&
        !modify_structure) {
      // no need to interpret the structure
    } else {
      st = mmcif_read_atoms(cif_in);
      if (st.models.empty())
        gemmi::fail("No atoms in the input file. Is it mmCIF?");
    }
  } else if (input_type == FileType::Pdb) {
    st = read_structure(input, gemmi::CoorFormat::Pdb);
  } else {
    gemmi::fail("Unexpected input format.");
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
    for (gemmi::Model& model : st.models) {
      std::vector<gemmi::Chain> new_chains;
      for (gemmi::Chain& chain : model.chains)
        gemmi::vector_move_extend(new_chains, split_by_segments(chain));
      model.chains = std::move(new_chains);
      // We should also modify Entity instances, but for now we don't.
      add_backlinks(model); // fix backlinks from Residue/Atom
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

  if (output_type == FileType::Json) {
    if (input_type != FileType::Cif && input_type != FileType::Json)
      gemmi::fail("Conversion to JSON is possible only from CIF");
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
    writer.write_json(cif_in);
  } else if (output_type == FileType::Pdb) {
    // call wrapper from output.cpp - to make building faster
    write_pdb(st, *os, options[IotbxCompat]);
  } else if (output_type == FileType::Null) {
    *os << st.name << ": " << count_atom_sites(st) << " atom locations";
    if (st.models.size() > 1)
      *os << " (total in " << st.models.size() << " models)";
    *os << ".\n";
  } else if (output_type == FileType::Cif) {
    if ((input_type != FileType::Cif && input_type != FileType::Json)
        || modify_structure) {
      cif_in.blocks.clear();  // temporary, for testing
      cif_in.blocks.resize(1);
      update_cif_block(st, cif_in.blocks[0]);
    }
    *os << cif_in;
  } else if (output_type == FileType::Crd) {
    *os << make_crd(st);
  }
}

int main(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p;
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  std::string input = p.nonOption(0);
  const char* output = p.nonOption(1);

  std::map<std::string, FileType> filetypes {{"json", FileType::Json},
                                             {"pdb", FileType::Pdb},
                                             {"cif", FileType::Cif},
                                             {"crd", FileType::Crd},
                                             {"none", FileType::Null}};

  FileType in_type = p.options[FormatIn] ? filetypes[p.options[FormatIn].arg]
                                         : get_format_from_extension(input);
  if (in_type == FileType::Unknown && gemmi::is_pdb_code(input)) {
    input = expand_pdb_code_to_path_or_fail(input);
    in_type = FileType::Cif;
  }
  if (in_type == FileType::Unknown) {
    std::cerr << "The input format cannot be determined from input"
                 " filename. Use option --from.\n";
    return 1;
  }

  FileType out_type = p.options[FormatOut] ? filetypes[p.options[FormatOut].arg]
                                           : get_format_from_extension(output);
  if (out_type == FileType::Unknown) {
    std::cerr << "The output format cannot be determined from output"
                 " filename. Use option --to.\n";
    return 1;
  }
  if (p.options[Verbose])
    std::cerr << "Converting " << input << " ..." << std::endl;

  try {
    convert(input, in_type, output, out_type, p.options);
  } catch (tao::pegtl::parse_error& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 2;
  }
  return 0;
}

// vim:sw=2:ts=2:et:path^=../include,../third_party
