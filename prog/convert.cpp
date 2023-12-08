// Copyright 2017 Global Phasing Ltd.

#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"
#include "gemmi/polyheur.hpp"  // for setup_entities, remove_waters, ...
#include "gemmi/modify.hpp"    // for remove_hydrogens, remove_anisou
#include "gemmi/align.hpp"     // for assign_label_seq_id
#include "gemmi/to_pdb.hpp"    // for write_pdb, ...
#include "gemmi/fstream.hpp"   // for Ofstream, Ifstream
#include "gemmi/to_mmcif.hpp"  // for update_mmcif_block
#include "gemmi/assembly.hpp"  // for ChainNameGenerator, transform_to_assembly
#include "gemmi/pirfasta.hpp"  // for read_pir_or_fasta
#include "gemmi/mmread_gz.hpp" // for read_structure_gz
#include "gemmi/select.hpp"    // for Selection
#include "gemmi/enumstr.hpp"   // for polymer_type_to_string

#include <cstring>
#include <iostream>
#include <algorithm>           // for sort

#define GEMMI_PROG convert
#include "options.h"
#include "cifmod.h"

namespace cif = gemmi::cif;
using gemmi::CoorFormat;
using gemmi::HowToNameCopiedChain;

namespace {

struct ConvArg: public Arg {
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
  FormatIn=AfterCifModOptions, FormatOut, CifStyle, AllAuth, BlockName,
  ExpandNcs, AsAssembly,
  RemoveH, RemoveWaters, RemoveLigWat, TrimAla, Select, Remove, ApplySymop,
  Reframe, ShortTer, Linkr, CopyRemarks, Minimal, ShortenCN, RenameChain,
  ChangeCcdCode, SetSeq,
  SiftsNum, Biso, Anisou, SetCis, SegmentAsChain, OldPdb, ForceLabel
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nwith possible conversions between PDB, mmCIF and mmJSON."
    "\nFORMAT can be specified as one of: mmcif, mmjson, pdb, chemcomp (read-only)."
    "\nchemcomp = example coordinates of a component from CCD or monomer library."
    "\n\nGeneral options:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { FormatIn, 0, "", "from", Arg::CoorFormat,
    "  --from=FORMAT  \tInput format (default: from the file extension)." },
  { FormatOut, 0, "", "to", Arg::CoorFormat,
    "  --to=FORMAT  \tOutput format (default: from the file extension)." },

  { NoOp, 0, "", "", Arg::None, "\nmmCIF output options:" },
  { CifStyle, 0, "", "style", Arg::CifStyle,
    "  --style=STYLE  \tone of: default, pdbx (categories separated with #),"
                     " aligned (left-aligned columns)." },
  { AllAuth, 0, "", "all-auth", Arg::None,
    "  --all-auth  \tWrite _atom_site.auth_atom_id (same as label_atom_id)"
    " and auth_comp_id (same as label_comp_id)." },
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
  { ChangeCcdCode, 0, "", "monomer", Arg::ColonPair,
    "  --monomer=OLD:NEW  \tChange monomer name (CCD code) OLD to NEW." },
  { SetSeq, 0, "s", "", Arg::Required,
    "  -s FILE  \tUse sequence(s) from FILE in PIR or FASTA format. Each chain"
    " is assigned the best matching sequence, if any." },
  { SiftsNum, 0, "", "sifts-num", Arg::None,
    "  --sifts-num  \tUse SIFTS-mapped position in UniProt sequence as sequence ID." },
  { Biso, 0, "B", "", Arg::Required,
    "  -B MIN[:MAX]  \tSet isotropic B-factors to a single value or change values "
      "out of given range to MIN/MAX." },
  { Anisou, 0, "", "anisou", ConvArg::AnisouChoice,
    "  --anisou=yes|no|heavy  \tAdd or remove ANISOU records." },
  // disabled: probably not used and implementing it would require either
  // using Topo::set_cispeps_in_structure() or duplicating the code.
  //{ SetCis, 0, "", "set-cispep", Arg::None,
  //  "  --set-cispep  \tReset CISPEP records from omega angles." },

  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL  \tOutput only the selection." },
  { Remove, 0, "", "remove", Arg::Required,
    "  --remove=SEL  \tRemove the selection." },
  { ApplySymop, 0, "", "apply-symop", Arg::Required,
    "  --apply-symop=OP  \tApply symmetry operation (e.g. '-x,y+1/2,-z'." },
  { Reframe, 0, "", "reframe", Arg::None,
    "  --reframe  \tStandardize the coordinate system (frame)." },
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
  { NoOp, 0, "", "", Arg::None,
    "\nWhen output file is -, write to standard output (default format: pdb)." },
  { 0, 0, 0, 0, 0, 0 }
};

std::string format_as_string(CoorFormat format) {
  switch (format) {
    case CoorFormat::Unknown: return "unknown";
    case CoorFormat::Detect: return "unknown";
    case CoorFormat::Pdb: return "pdb";
    case CoorFormat::Mmcif: return "mmcif";
    case CoorFormat::Mmjson: return "mmjson";
    case CoorFormat::ChemComp: return "chemcomp";
  }
  gemmi::unreachable();
}

std::string read_whole_file(std::istream& stream) {
  constexpr std::size_t buf_size = 4096;
  char buffer[buf_size];
  stream.exceptions(std::ios_base::badbit);  // throw on fail/bad
  std::string out;
  while (stream) {
    stream.read(buffer, buf_size);
    out.append(buffer, stream.gcount());
  }
  return out;
}

void convert(gemmi::Structure& st,
             const std::string& output, CoorFormat output_type,
             const std::vector<option::Option>& options) {
  if (st.models.empty())
    gemmi::fail("No atoms in the input file. Wrong file format?");
  if (st.ter_status == 'e')
    std::cerr << "WARNING: TER records in the input PDB are clearly where they shouldn't be."
                 "\nWARNING: Ignoring all TER records." << std::endl;

  for (const option::Option* opt = options[ChangeCcdCode]; opt; opt = opt->next()) {
    const char* sep = std::strchr(opt->arg, ':');
    std::string old_name(opt->arg, sep);
    std::string new_name(sep+1);
    if (options[Verbose])
      std::cerr << "Renaming " << old_name << " to " << new_name << std::endl;
    gemmi::change_ccd_code(st, old_name, new_name);
  }

  gemmi::setup_entities(st);
  if (st.input_format == CoorFormat::Pdb) {
    if (!options[SetSeq])
      gemmi::assign_label_seq_id(st, options[ForceLabel]);
    if (!options[CopyRemarks])
      st.raw_remarks.clear();
  } else {
    // handles special tag from Refmac's mmCIF
    store_deuterium_as_fraction(st, false);
  }

  if (options[Select])
    gemmi::Selection(options[Select].arg).remove_not_selected(st);
  if (options[Remove])
    gemmi::Selection(options[Remove].arg).remove_selected(st);
  if (st.models.empty())
    gemmi::fail("all models got removed");
  if (options[ApplySymop]) {
    gemmi::Op op = gemmi::parse_triplet(options[ApplySymop].arg);
    transform_pos_and_adp(st, st.cell.op_as_transform(op));
  }
  if (options[Reframe])
    standardize_crystal_frame(st);

  if (options[Biso]) {
    const char* start = options[Biso].arg;
    char* endptr = nullptr;
    float value1 = std::strtof(start, &endptr);
    float value2 = value1;
    if (endptr != start && *endptr == ':') {
      start = endptr + 1;
      value2 = std::strtof(start, &endptr);
    }
    if (endptr != start && *endptr == '\0')
      assign_b_iso(st, value1, value2);
    else
      gemmi::fail("argument for -B should be a number or number:number");
  }

  for (const option::Option* opt = options[Anisou]; opt; opt = opt->next()) {
    char anisou_opt = opt->arg[0];
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

  for (const option::Option* opt = options[RenameChain]; opt; opt = opt->next()) {
    const char* sep = std::strchr(opt->arg, ':');
    std::string old_name(opt->arg, sep);
    std::string new_name(sep+1);
    gemmi::rename_chain(st, old_name, new_name);
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
    std::vector<std::string> fasta_sequences;
    for (const option::Option* opt = options[SetSeq]; opt; opt = opt->next()) {
      gemmi::Ifstream stream(opt->arg);
      std::string str = read_whole_file(stream.ref());
      for (gemmi::FastaSeq& fs : gemmi::read_pir_or_fasta(str))
        fasta_sequences.push_back(std::move(fs.seq));
    }
    if (options[Verbose])
      std::cerr << fasta_sequences.size() << " sequence(s) was read..." << std::endl;
    gemmi::clear_sequences(st);
    gemmi::assign_best_sequences(st, fasta_sequences);
    gemmi::deduplicate_entities(st);
    gemmi::assign_label_seq_id(st, options[ForceLabel]);
    for (gemmi::Entity& ent : st.entities) {
      if (ent.entity_type == gemmi::EntityType::Polymer && ent.full_sequence.empty())
        std::cerr << "No sequence found for "
                  << polymer_type_to_string(ent.polymer_type) << " entity " << ent.name
                  << " (" << gemmi::join_str(ent.subchains, ',') << ')' << std::endl;
    }
  }

  if (options[SiftsNum]) {
    // Currently, we change seqid only for residues _pdbx_sifts_xref_db.unp_num
    // set, and leave seqid for other residues.
    // A more robust method would be to re-assign the whole polymer always when
    // nup_num is set for any residue. This would mean:
    // using insertion codes for insertions (up to 26 residues),
    // renumbering ligands and waters if needed, and
    // failing if we get duplicated seqid (multiple mappings for one chain).
    bool changed = false;
    for (gemmi::Model& model: st.models)
      for (gemmi::Chain& chain : model.chains) {
        const gemmi::Residue* prev_res = nullptr;
        bool wrong_order = false;
        for (gemmi::Residue& res : chain.residues) {
          if (res.sifts_unp.res) {
            res.seqid = gemmi::SeqId(res.sifts_unp.num, ' ');
            if (prev_res && !(prev_res->seqid < res.seqid))
              wrong_order = true;
            changed = true;
          }
          prev_res = &res;
        }
        if (wrong_order) {
          std::cerr << "WARNING: new sequence IDs are in wrong order in chain "
                    << chain.name;
          if (st.models.size() > 1)
            std::cerr << " (model " << model.name << ')';
          std::cerr << std::endl;
        }
      }
    if (options[Verbose] && !changed)
      std::cerr << "Option --sifts-num had no effect, "
                   "it works only with _pdbx_sifts_xref_db.unp_num." << std::endl;
  }

  HowToNameCopiedChain how = HowToNameCopiedChain::AddNumber;
  if (output_type == CoorFormat::Pdb)
    how = HowToNameCopiedChain::Short;
  if (options[AsAssembly]) {
    std::ostream* out = options[Verbose] ? &std::cerr : nullptr;
    gemmi::transform_to_assembly(st, options[AsAssembly].arg, how, out);
  }

  if (options[ExpandNcs]) {
    HowToNameCopiedChain how_ncs;
    switch (options[ExpandNcs].arg[0]) {
      case 'd': how_ncs = HowToNameCopiedChain::Dup; break;
      case 'x': how_ncs = HowToNameCopiedChain::Short; break;
      default: how_ncs = HowToNameCopiedChain::AddNumber; break;
    }
    gemmi::expand_ncs(st, how_ncs);
    for (gemmi::Model& model : st.models)
      gemmi::merge_atoms_in_expanded_model(model, st.cell);
  }

  if (options[RemoveH])
    remove_hydrogens(st);

  if (options[RemoveWaters]) {
    remove_waters(st);
    st.remove_empty_chains();
  }

  if (options[RemoveLigWat]) {
    remove_ligands_and_waters(st);
    st.remove_empty_chains();
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
    if (options[Minimal]) {
      gemmi::add_minimal_mmcif_data(st, doc.blocks[0]);
    } else {
      gemmi::MmcifOutputGroups groups(true);
      groups.auth_all = options[AllAuth];
      gemmi::update_mmcif_block(st, doc.blocks[0], groups);
    }
    apply_cif_doc_modifications(doc, options);

    if (output_type == CoorFormat::Mmcif) {
      write_cif_to_stream(os.ref(), doc, cif_write_options(options[CifStyle]));
    } else /*output_type == CoorFormat::Mmjson*/ {
      cif::JsonWriter writer(os.ref());
      writer.set_mmjson();
      writer.write_json(doc);
    }
  } else if (output_type == CoorFormat::Pdb) {
    shorten_ccd_codes(st);
    gemmi::PdbWriteOptions opt;
    if (options[Minimal])
      opt = gemmi::PdbWriteOptions::minimal();
    if (options[ShortTer])
      opt.numbered_ter = false;
    if (options[Linkr])
      opt.use_linkr = true;
    gemmi::write_pdb(st, os.ref(), opt);
  }
}

} // anonymous namespace

int GEMMI_MAIN(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  OptParser p(EXE_NAME);
  p.simple_parse(argc, argv, Usage);
  p.require_positional_args(2);

  CoorFormat in_type = coor_format_as_enum(p.options[FormatIn]);
  if (in_type == CoorFormat::Unknown)
    in_type = CoorFormat::Detect;
  char pdb_code_type = in_type == CoorFormat::Pdb ? 'P' : 'M';
  std::string input = p.coordinate_input_file(0, pdb_code_type);
  const char* output = p.nonOption(1);

  CoorFormat out_type = coor_format_as_enum(p.options[FormatOut]);
  if (out_type == CoorFormat::Unknown) {
    if (output[0] == '-' && output[1] == '\0')
      out_type = CoorFormat::Pdb;
    else
      out_type = gemmi::coor_format_from_ext_gz(output);
  }
  if (out_type == CoorFormat::ChemComp) {
    std::cerr << "The output format cannot be chemcomp.\n";
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
