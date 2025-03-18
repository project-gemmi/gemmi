// Copyright 2017 Global Phasing Ltd.

#include "gemmi/to_cif.hpp"
#include "gemmi/to_json.hpp"   // for write_json_to_stream
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
#include "gemmi/calculate.hpp" // for parse_triplet_as_ftransform

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

  static option::ArgStatus CoorFormatIn(const option::Option& option, bool msg) {
    return Choice(option, msg, {"cif", "mmcif", "pdb", "json", "mmjson",
                                "chemcomp", "chemcomp:m", "chemcomp:i"});
  }

  static option::ArgStatus RecordChoice(const option::Option& option, bool msg) {
    auto status = Arg::Optional(option, msg);
    if (status == option::ARG_OK && option.arg[0] != 'A' && option.arg[0] != 'H') {
      if (msg)
        fprintf(stderr, "If option %.*s has argument, it must be ATOM or HETATM.\n",
                option.namelen, option.name);
      status = option::ARG_ILLEGAL;
    }
    return status;
  }
};

enum OptionIndex {
  FormatIn=AfterCifModOptions, FormatOut, CifStyle, AllAuth, BlockName,
  ExpandNcs, AsAssembly,
  RemoveH, RemoveWaters, RemoveLigWat, TrimAla, Select, Remove, ApplySymop,
  Reframe, ShortTer, Linkr, CopyRemarks, Minimal, ShortenCN, RenameChain,
  ShortenTLC, ChangeCcdCode, SetSeq, SiftsNum,
  Biso, BisoScale, AddTls, Anisou, AssignRecords,
  SetCis, SegmentAsChain, OldPdb, ForceLabel
};

const option::Descriptor Usage[] = {
  { NoOp, 0, "", "", Arg::None,
    "Usage:"
    "\n " EXE_NAME " [options] INPUT_FILE OUTPUT_FILE"
    "\n\nAllows conversion between PDB, mmCIF, and mmJSON formats."
    "\n\nGeneral options:" },
  CommonUsage[Help],
  CommonUsage[Version],
  CommonUsage[Verbose],
  { FormatIn, 0, "", "from", ConvArg::CoorFormatIn,
    "  --from=FORMAT  \tInput format (default: inferred from file extension)." },
  { FormatOut, 0, "", "to", Arg::CoorFormat,
    "  --to=FORMAT  \tOutput format (default: inferred from file extension)." },

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
    "  --old-pdb \tRead only the first 72 characters of each line." },
  { ForceLabel, 0, "L", "force-label", Arg::None,
    "  -L, --force-label  \tAdd label_seq_id even if SEQRES is missing" },

  { NoOp, 0, "", "", Arg::None, "\nPDB output options:" },
  { ShortTer, 0, "", "short-ter", Arg::None,
    "  --short-ter  \tWrite PDB TER records without numbers (iotbx compat.)." },
  { Linkr, 0, "", "linkr", Arg::None,
    "  --linkr  \tWrite LINKR record (for Refmac) if link_id is known." },
  { CopyRemarks, 0, "", "copy-remarks", Arg::None,
    "  --copy-remarks  \t(pdb->pdb only) Copy REMARK records." },

  { NoOp, 0, "", "", Arg::None, "\nGeneral output options:" },
  { Minimal, 0, "", "minimal", Arg::None,
    "  --minimal  \tWrite only the most essential records." },
  { ShortenCN, 0, "", "shorten", Arg::None,
    "  --shorten  \tShorten chain names to 1 (if # < 63) or 2 characters." },
  { RenameChain, 0, "", "rename-chain", Arg::ColonPair,
    "  --rename-chain=OLD:NEW  \tRename chain OLD to NEW "
    "(--rename-chain=:A adds missing chain IDs)." },
  { ShortenTLC, 0, "", "shorten-tlc", Arg::None,
    "  --shorten-tlc  \tChange 5-character monomer names to 3-char. aliases." },
  { ChangeCcdCode, 0, "", "monomer", Arg::ColonPair,
    "  --monomer=OLD:NEW  \tChange monomer name (CCD code) OLD to NEW." },
  { SetSeq, 0, "s", "", Arg::Required,
    "  -s FILE  \tUse sequence(s) from FILE in PIR or FASTA format. Each chain"
    " is assigned the best matching sequence, if any." },
  { SiftsNum, 0, "", "sifts-num", Arg::Optional,
    "  --sifts-num[=AC,...]  \tSet sequence ID to SIFTS-mapped UniProt positions,"
    " add 5000+ to non-mapped seqnums. See docs for details." },
  { Biso, 0, "B", "", Arg::Required,
    "  -B MIN[:MAX]  \tSet isotropic B-factors to a single value or constrain them to a range." },
  { BisoScale, 0, "", "scale-biso", Arg::Float,
    "  --scale-biso=MULT  \tMultiply isotropic B-factors by MULT." },
  { AddTls, 0, "", "add-tls", Arg::None,
    "  --add-tls  \tconvert from residual to full B-factors." },
  { Anisou, 0, "", "anisou", ConvArg::AnisouChoice,
    "  --anisou=yes|no|heavy  \tAdd or remove ANISOU records." },
  { AssignRecords, 0, "", "assign-records", ConvArg::RecordChoice,
    "  --assign-records[=A|H]  \tRe-assign ATOM/HETATM (w/o argument: auto)." },
  // disabled: probably not used and implementing it would require either
  // using Topo::set_cispeps_in_structure() or duplicating the code.
  //{ SetCis, 0, "", "set-cispep", Arg::None,
  //  "  --set-cispep  \tReset CISPEP records from omega angles." },

  { NoOp, 0, "", "", Arg::None, "\nMacromolecular operations:" },
  { Select, 0, "", "select", Arg::Required,
    "  --select=SEL  \tOutput only the specified selection." },
  { Remove, 0, "", "remove", Arg::Required,
    "  --remove=SEL  \tRemove the specified selection." },
  { ApplySymop, 0, "", "apply-symop", Arg::Required,
    "  --apply-symop=OP  \tApply operation, e.g. '-x,y+1/2,-z' or 'x,y,z+0.1'." },
  { Reframe, 0, "", "reframe", Arg::None,
    "  --reframe  \tStandardize the coordinate system (frame)." },
  { ExpandNcs, 0, "", "expand-ncs", ConvArg::NcsChoice,
    "  --expand-ncs=dup|num|x  \tExpand strict NCS from MTRIXn or _struct_ncs_oper. "
    "Argument controls naming of new chains; see docs." },
  { AsAssembly, 0, "", "assembly", Arg::Required,
    "  --assembly=ID  \tOutput bioassembly with specified ID (1, 2, ...)." },
  { RemoveH, 0, "", "remove-h", Arg::None,
    "  --remove-h  \tRemove hydrogens." },
  { RemoveWaters, 0, "", "remove-waters", Arg::None,
    "  --remove-waters  \tRemove waters." },
  { RemoveLigWat, 0, "", "remove-lig-wat", Arg::None,
    "  --remove-lig-wat  \tRemove ligands and waters." },
  { TrimAla, 0, "", "trim-to-ala", Arg::None,
    "  --trim-to-ala  \tTrim aminoacids to alanine." },
  { NoOp, 0, "", "", Arg::None,
    "\nFORMAT can be specified as one of: mmcif, mmjson, pdb. chemcomp (read-only)."
    "\nchemcomp = coordinates of a component from CCD or monomer library (see docs)."
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

int select_ac_index(const std::vector<std::string>& acs, // Entity::sifts_unp_acc
                    const gemmi::ConstResidueSpan& polymer,
                    const std::vector<std::string>& preferred_acs) {
  for (const std::string& pa : preferred_acs) {
    if (pa == "=")
      return -2;
    if (acs.empty())
      continue;
    if (pa == "*") {
      // return SIFTS mapping with most residues
      std::vector<int> counts(acs.size(), 0);
      for (const gemmi::Residue& res : polymer)
        if (res.sifts_unp.res && res.sifts_unp.acc_index < acs.size())
          ++counts[res.sifts_unp.acc_index];
      auto max_el = std::max_element(counts.begin(), counts.end());
      if (*max_el > 0)
        return int(max_el - counts.begin());
    } else {
      auto it = std::find(acs.begin(), acs.end(), pa);
      if (it != acs.end())
        return int(it - acs.begin());
    }
  }
  return -1;
}

// Set seqid corresponding to one UniProt sequence.
// Other residues have seqid changed to avoid conflicts.
void to_sifts_num(gemmi::Structure& st, const std::vector<std::string>& preferred_acs) {
  using Key = std::pair<std::string, gemmi::SeqId>;
  std::map<Key, gemmi::SeqId> seqid_map;

  // find new sequence IDs and store them in seqid_map
  std::map<std::string, uint16_t> chain_offsets;
  for (gemmi::Model& model: st.models) {
    for (gemmi::Chain& chain : model.chains) {
      auto polymer = chain.get_polymer();
      gemmi::Entity* ent = st.get_entity_of(polymer);
      int ac_index = -1;
      if (ent)
        ac_index = select_ac_index(ent->sifts_unp_acc, polymer, preferred_acs);
      std::uint16_t max_unp_num = 0;
      // first pass - set seqid_map for AC-corresponding residues
      if (ac_index >= 0)
        for (gemmi::Residue& res : chain.residues) {
          if (res.sifts_unp.res && (int) res.sifts_unp.acc_index == ac_index) {
            // assert(ent && res.entity_id == ent->name);
            gemmi::SeqId new_seqid(res.sifts_unp.num, ' ');
            seqid_map.emplace(std::make_pair(chain.name, res.seqid), new_seqid);
            res.seqid = new_seqid;
            max_unp_num = std::max(max_unp_num, res.sifts_unp.num);
          }
        }
      // retrieve or (for a new chain) store the offset for non-AC residues
      auto result = chain_offsets.emplace(chain.name, 0);
      std::uint16_t& offset = result.first->second;
      if (result.second && ac_index != -2) {
        if (max_unp_num <= 4950)
          offset = 5000;
        else
          offset = (max_unp_num + 1049) / 1000 * 1000;  // N * 1000, N > 5
      }
      // second pass - set seqid_map for non-AC residues
      for (gemmi::Residue& res : chain.residues)
        if (!res.sifts_unp.res || res.sifts_unp.acc_index != ac_index) {
          gemmi::SeqId orig_seqid = res.seqid;
          res.seqid.num += offset;
          seqid_map.emplace(std::make_pair(chain.name, orig_seqid), res.seqid);
        }
    }
  }

  auto update_seqid = [&](const std::string& chain_name, gemmi::SeqId& seqid) {
    auto it = seqid_map.find(std::make_pair(chain_name, seqid));
    if (it != seqid_map.end())
      seqid = it->second;
    else
      seqid.num = gemmi::SeqId::OptionalNum();
  };
  process_sequence_ids(st, update_seqid);
  // Just remove outdated DbRef::seq_*; it is needed when writing a PDB file,
  // but if it's absent, it's determined from DbRef::label_seq_*.
  for (gemmi::Entity& ent : st.entities)
    for (gemmi::Entity::DbRef& dbref : ent.dbrefs)
      dbref.seq_begin = dbref.seq_end = gemmi::SeqId();
}

// simplified check, can return false positives
bool is_uniprot_ac_format(const std::string& ac) {
  return ac.size() >= 6 && ac[0] >= 'A' && ac[0] <= 'Z' && ac[1] >= '0' && ac[1] <= '9';
}

const gemmi::TlsGroup*
get_tls_group_by_id(const std::vector<gemmi::TlsGroup>& tls_groups, short num_id) {
  if (size_t(num_id - 1) < tls_groups.size() && tls_groups[num_id - 1].num_id == num_id)
    return &tls_groups[num_id - 1];
  for (const gemmi::TlsGroup& tg : tls_groups)
    if (tg.num_id == num_id)
      return &tg;
  return nullptr;
}

void convert(gemmi::Structure& st,
             const std::string& output, CoorFormat output_type,
             const std::vector<option::Option>& options) {
  if (st.models.empty())
    gemmi::fail("No atoms in the input (", format_as_string(st.input_format), ") file. "
                "Wrong file format?");
  if (st.ter_status == 'e')
    std::cerr << "WARNING: TER records in the input PDB are clearly where they shouldn't be."
                 "\nWARNING: Ignoring all TER records." << std::endl;

  for (const option::Option* opt = options[ChangeCcdCode]; opt; opt = opt->next()) {
    const char* sep = std::strchr(opt->arg, ':');
    std::string old_name(opt->arg, sep);
    std::string new_name(sep+1);
    if (options[Verbose])
      std::cerr << "Renaming " << old_name << " to " << new_name << std::endl;
    gemmi::rename_residues(st, old_name, new_name);
  }

  if (options[AssignRecords])
    // avoid using initial ATOM/HETATM in setup_entities()
    gemmi::assign_het_flags(st, '\0');
  // gemmi::setup_entities(st) with postponed deduplicate_entities()
  gemmi::add_entity_types(st, /*overwrite=*/false);
  assign_subchains(st, /*force=*/false);
  ensure_entities(st);
  // call deduplicate_entities() after add_microhetero_to_sequences()
  if (st.input_format == CoorFormat::Pdb) {
    if (!options[SetSeq]) {
      gemmi::assign_label_seq_id(st, options[ForceLabel]);
      gemmi::add_microhetero_to_sequences(st);
    }
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
    gemmi::FTransform frac_tr = gemmi::parse_triplet_as_ftransform(options[ApplySymop].arg);
    transform_pos_and_adp(st, st.cell.orthogonalize_transform(frac_tr));
  }

  if (options[Reframe])
    standardize_crystal_frame(st);

  if (options[BisoScale]) {
    double mult = std::atof(options[BisoScale].arg);
    for (gemmi::Model& model: st.models)
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
          for (gemmi::Atom& atom : res.atoms)
            atom.b_iso = float(atom.b_iso * mult);
  }

  if (options[Biso]) {
    float value1, value2;
    bool ok = parse_number_or_range(options[Biso].arg, &value1, &value2);
    if (!ok)
      gemmi::fail("argument for -B should be a number or number:number");
    assign_b_iso(st, value1, value2);
  }

  if (options[AddTls]) {
    const std::vector<gemmi::TlsGroup>* tls_groups = st.meta.get_tls_groups();
    if (!tls_groups)
      return;
    gemmi::add_tls_group_ids(st);
    for (gemmi::Model& model: st.models)
      for (gemmi::Chain& chain : model.chains)
        for (gemmi::Residue& res : chain.residues)
          for (gemmi::Atom& atom : res.atoms) {
            if (atom.aniso.nonzero())
              continue;
            if (const gemmi::TlsGroup* tg = get_tls_group_by_id(*tls_groups, atom.tls_group_id)) {
              gemmi::SMat33<double> u = gemmi::calculate_u_from_tls(*tg, atom.pos)
                                        .added_kI(1. / gemmi::u_to_b() * atom.b_iso);
              atom.aniso = {(float)u.u11, (float)u.u22, (float)u.u33,
                            (float)u.u12, (float)u.u13, (float)u.u23};
              atom.b_iso = float(gemmi::u_to_b() / 3. * u.trace());
            }
          }
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

  if (const option::Option* opt = options[AssignRecords]) {
    char flag = opt->arg ? opt->arg[0] : '?';
    gemmi::assign_het_flags(st, flag);
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
    gemmi::assign_label_seq_id(st, options[ForceLabel]);
    gemmi::add_microhetero_to_sequences(st);
    for (gemmi::Entity& ent : st.entities) {
      if (ent.entity_type == gemmi::EntityType::Polymer && ent.full_sequence.empty())
        std::cerr << "No sequence found for "
                  << polymer_type_to_string(ent.polymer_type) << " entity " << ent.name
                  << " (" << gemmi::join_str(ent.subchains, ',') << ')' << std::endl;
    }
  }
  // call it after add_microhetero_to_sequences()
  gemmi::deduplicate_entities(st);

  if (const option::Option* opt = options[SiftsNum]) {
    if (std::all_of(st.entities.begin(), st.entities.end(),
                    [](const gemmi::Entity& ent) { return ent.sifts_unp_acc.empty(); }))
      gemmi::fail("--sifts-num failed: SIFTS annotation not found in the file.\n"
                  "Check if tags such as _pdbx_sifts_xref_db.unp_num are present.");
    std::vector<std::string> preferred_acs;
    if (opt->arg) {
      gemmi::split_str_into(opt->arg, ',', preferred_acs);
      for (const std::string& ac : preferred_acs)
        if (ac != "*" && ac != "=" && !is_uniprot_ac_format(ac))
          gemmi::fail(ac + " is not in UniProtKB AC format, from: " + opt->name);
    } else {
      preferred_acs.emplace_back("*");
    }
    to_sifts_num(st, preferred_acs);
  }

  HowToNameCopiedChain how = HowToNameCopiedChain::AddNumber;
  if (output_type == CoorFormat::Pdb)
    how = HowToNameCopiedChain::Short;
  if (options[AsAssembly]) {
    gemmi::Logger logger{&gemmi::Logger::to_stderr, options[Verbose] ? 6 : 3};
    gemmi::transform_to_assembly(st, options[AsAssembly].arg, how, logger);
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

  if (options[ShortenTLC] || output_type == CoorFormat::Pdb)
    shorten_ccd_codes(st);

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
      cif::write_mmjson_to_stream(os.ref(), doc);
    }
  } else if (output_type == CoorFormat::Pdb) {
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
      if (in_type == CoorFormat::ChemComp) {
        int which = 7;
        // chemcomp:m or chemcomp:i
        if (p.options[FormatIn].arg && std::strlen(p.options[FormatIn].arg) > 9)
          // cf. ChemCompModel
          which = p.options[FormatIn].arg[9] == 'i' ? 4 : 2;
        st = gemmi::read_structure_from_chemcomp_gz(input, nullptr, which);
      } else {
        st = gemmi::read_structure_gz(input, in_type);
      }
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
