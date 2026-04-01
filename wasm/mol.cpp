// Copyright Global Phasing Ltd.

#include "common.h"
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <gemmi/model.hpp>
#include <gemmi/assembly.hpp> // for get_nearby_sym_ops, get_sym_image
#include <gemmi/select.hpp>   // for Selection
#include <gemmi/mmread.hpp>   // for read_structure_from_memory
#include <gemmi/to_cif.hpp>   // for cif::write_cif_to_stream
#include <gemmi/to_mmcif.hpp> // for make_mmcif_document
#include <gemmi/to_pdb.hpp>   // for make_pdb_string
#include <gemmi/polyheur.hpp> // for setup_entities
#include <gemmi/enumstr.hpp>  // for entity_type_to_string
#include <gemmi/calculate.hpp> // for count_occupancies
#include <gemmi/monlib.hpp>   // for MonLib
#include <gemmi/read_cif.hpp> // for cif::read_string
#include <emscripten/val.h>

gemmi::CoorFormat format_to_enum(const std::string& format) {
  using gemmi::CoorFormat;
  if (format == "unknown")
    return CoorFormat::Unknown;
  if (format == "detect")
    return CoorFormat::Detect;
  if (format == "pdb")
    return CoorFormat::Pdb;
  if (format == "mmcif")
    return CoorFormat::Mmcif;
  if (format == "mmjson")
    return CoorFormat::Mmjson;
  if (format == "chemcomp")
    return CoorFormat::ChemComp;
  gemmi::fail("unknown file format: " + format);
}

// IIUC passing string by value is OK here, it's copied the JS side anyway
gemmi::Structure read_structure(std::string buf, std::string name, std::string format) {
  auto st = gemmi::read_structure_from_memory(buf.data(), buf.size(), name,
                                              format_to_enum(format));
  setup_entities(st);
  return st;
}

template <typename T>
size_t get_children_length(const T& t) { return t.children().size(); }

template <typename T>
typename T::child_type* get_child(T& t, int n) {
  size_t idx = n >= 0 ? (size_t) n : n + t.children().size();
  return &t.children().at(idx);
}

template <typename T, typename... Args >
decltype(auto) wrap_children() {
  return em::class_<T, Args...>(T::what())
    .template constructor<>()
    .property("length", &get_children_length<T>)
    .function("at", &get_child<T>, em::allow_raw_pointers())
    ;
}

std::string element_uname(const gemmi::Atom& atom) {
  return atom.element.uname();
}

void set_element(gemmi::Atom& atom, const std::string& element) {
  gemmi::Element parsed(element);
  if (parsed == gemmi::El::X && element != "X" && element != "x")
    gemmi::fail("unknown element: " + element);
  atom.element = parsed;
}

bool atom_is_metal(const gemmi::Atom& atom) {
  return atom.element.is_metal();
}

std::string get_seqid_string(const gemmi::ResidueId& res) {
  return res.seqid.str();
}

char get_seqid_icode(const std::string& icode) {
  if (icode.size() > 1)
    gemmi::fail("insertion code must be empty or one character");
  return icode.empty() ? ' ' : icode[0];
}

void set_seqid(gemmi::ResidueId& res, int num, const std::string& icode) {
  res.seqid = gemmi::SeqId(num, get_seqid_icode(icode));
}

void set_seqid_string(gemmi::ResidueId& res, const std::string& seqid) {
  char* endptr;
  long num = std::strtol(seqid.c_str(), &endptr, 10);
  if (endptr == seqid.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
    throw std::invalid_argument("Not a seqid: " + seqid);
  res.seqid.num = static_cast<int>(num);
  res.seqid.icode = *endptr == '\0' ? ' ' : *endptr;
}

void structure_add_model(gemmi::Structure& st, const gemmi::Model& model) {
  st.models.push_back(model);
}

void model_add_chain(gemmi::Model& model, const gemmi::Chain& chain) {
  model.chains.push_back(chain);
}

void chain_add_residue(gemmi::Chain& chain, const gemmi::Residue& residue) {
  chain.residues.push_back(residue);
}

void residue_add_atom(gemmi::Residue& residue, const gemmi::Atom& atom) {
  residue.atoms.push_back(atom);
}

std::string get_entity_type_string(const gemmi::Residue& res) {
  return entity_type_to_string(res.entity_type);
}

bool has_chemcomp_data(const gemmi::ChemComp& cc) {
  return !cc.atoms.empty() || !cc.rt.bonds.empty() || !cc.aliases.empty() ||
         !cc.type_or_group.empty();
}

const gemmi::ChemComp* find_chemcomp(const gemmi::Structure& st,
                                     const gemmi::MonLib& monlib,
                                     const gemmi::Residue& res) {
  auto st_it = st.chemcomps.find(res.name);
  if (st_it != st.chemcomps.end() && has_chemcomp_data(st_it->second))
    return &st_it->second;
  auto mon_it = monlib.monomers.find(res.name);
  return mon_it != monlib.monomers.end() ? &mon_it->second : nullptr;
}

const gemmi::ChemComp::Aliasing* find_polymer_aliasing(const gemmi::ChemComp& cc,
                                                       gemmi::PolymerType polymer_type) {
  for (const gemmi::ChemComp::Aliasing& aliasing : cc.aliases)
    if ((gemmi::is_polypeptide(polymer_type) &&
         gemmi::ChemComp::is_peptide_group(aliasing.group)) ||
        (gemmi::is_polynucleotide(polymer_type) &&
         gemmi::ChemComp::is_nucleotide_group(aliasing.group)))
      return &aliasing;
  return nullptr;
}

std::string remap_polymer_atom_name(const std::string& atom_name,
                                    const gemmi::ChemComp::Aliasing* aliasing) {
  if (aliasing)
    if (const std::string* ptr = aliasing->name_from_alias(atom_name))
      return *ptr;
  return atom_name;
}

template<typename It>
It residue_group_end(It it, It end) {
  It next = it + 1;
  while (next != end && next->group_key() == it->group_key())
    ++next;
  return next;
}

template<typename AddBond>
void add_inferred_polymer_bonds(gemmi::Structure& st, gemmi::Model& model,
                                const gemmi::MonLib& monlib, AddBond&& add_bond) {
  for (gemmi::Chain& chain : model.chains)
    for (gemmi::ResidueSpan sub : chain.subchains()) {
      const gemmi::Entity* ent = st.get_entity_of(sub);
      gemmi::PolymerType polymer_type = gemmi::get_or_check_polymer_type(ent, sub);
      if (!gemmi::is_polypeptide(polymer_type) &&
          !gemmi::is_polynucleotide(polymer_type))
        continue;
      if (sub.size() < 2)
        continue;
      auto prev_begin = sub.begin();
      auto prev_end = residue_group_end(prev_begin, sub.end());
      while (prev_end != sub.end()) {
        auto group_begin = prev_end;
        auto group_end = residue_group_end(group_begin, sub.end());
        for (auto ri = group_begin; ri != group_end; ++ri) {
          const gemmi::ChemComp* cc2 = find_chemcomp(st, monlib, *ri);
          const gemmi::ChemComp::Aliasing* alias2 = nullptr;
          if (cc2 && ((gemmi::is_polypeptide(polymer_type) &&
                       !gemmi::ChemComp::is_peptide_group(cc2->group)) ||
                      (gemmi::is_polynucleotide(polymer_type) &&
                       !gemmi::ChemComp::is_nucleotide_group(cc2->group))))
            alias2 = find_polymer_aliasing(*cc2, polymer_type);
          std::string atom2_name = gemmi::is_polypeptide(polymer_type)
                                 ? remap_polymer_atom_name("N", alias2)
                                 : remap_polymer_atom_name("P", alias2);
          gemmi::El atom2_el = gemmi::is_polypeptide(polymer_type) ? gemmi::El::N
                                                                   : gemmi::El::P;
          for (auto prev_ri = prev_begin; prev_ri != prev_end; ++prev_ri) {
            const gemmi::ChemComp* cc1 = find_chemcomp(st, monlib, *prev_ri);
            const gemmi::ChemComp::Aliasing* alias1 = nullptr;
            if (cc1 && ((gemmi::is_polypeptide(polymer_type) &&
                         !gemmi::ChemComp::is_peptide_group(cc1->group)) ||
                        (gemmi::is_polynucleotide(polymer_type) &&
                         !gemmi::ChemComp::is_nucleotide_group(cc1->group))))
              alias1 = find_polymer_aliasing(*cc1, polymer_type);
            std::string atom1_name = gemmi::is_polypeptide(polymer_type)
                                   ? remap_polymer_atom_name("C", alias1)
                                   : remap_polymer_atom_name("O3'", alias1);
            gemmi::El atom1_el = gemmi::is_polypeptide(polymer_type) ? gemmi::El::C
                                                                     : gemmi::El::O;
            for (const gemmi::Atom& a1 : prev_ri->atoms)
              if (a1.name == atom1_name && a1.element == atom1_el)
                for (const gemmi::Atom& a2 : ri->atoms)
                  if (a2.name == atom2_name && a2.element == atom2_el &&
                      (a2.altloc == a1.altloc || a2.altloc == '\0' || a1.altloc == '\0') &&
                      (gemmi::is_polypeptide(polymer_type)
                        ? gemmi::in_peptide_bond_distance(&a1, &a2)
                        : gemmi::in_nucleotide_bond_distance(&a1, &a2)))
                    add_bond(&a1, &a2, static_cast<int>(gemmi::BondType::Single));
          }
        }
        prev_begin = group_begin;
        prev_end = group_end;
      }
    }
}

std::string get_residue_ss_string(const gemmi::Residue& res) {
  switch (res.ss_from_file) {
    case gemmi::ResidueSs::Coil:
      return "Coil";
    case gemmi::ResidueSs::Helix:
      return "Helix";
    case gemmi::ResidueSs::Strand:
      return "Strand";
  }
  return "";
}

std::string get_residue_strand_sense_string(const gemmi::Residue& res) {
  switch (res.strand_sense_from_file) {
    case gemmi::ResidueStrandSense::NotStrand:
      return "NotStrand";
    case gemmi::ResidueStrandSense::Parallel:
      return "Parallel";
    case gemmi::ResidueStrandSense::First:
      return "First";
    case gemmi::ResidueStrandSense::Antiparallel:
      return "Antiparallel";
  }
  return "";
}

std::string get_residue_names(gemmi::Structure& st) {
  auto names = st.models.at(0).get_all_residue_names();
  std::string result;
  for (size_t i = 0; i < names.size(); ++i) {
    if (i > 0)
      result += ',';
    result += names[i];
  }
  return result;
}

std::string get_missing_monomer_names(gemmi::Structure& st) {
  auto names = st.models.at(0).get_all_residue_names();
  std::string result;
  for (size_t i = 0; i < names.size(); ++i) {
    auto it = st.chemcomps.find(names[i]);
    if (it != st.chemcomps.end() && has_chemcomp_data(it->second))
      continue;
    if (!result.empty())
      result += ',';
    result += names[i];
  }
  return result;
}

std::string make_mmcif_string(const gemmi::Structure& st) {
  gemmi::cif::Document doc = gemmi::make_mmcif_document(st);
  std::ostringstream os;
  gemmi::cif::write_cif_to_stream(os, doc);
  return os.str();
}

std::string make_pdb_string_default(const gemmi::Structure& st) {
  return gemmi::make_pdb_string(st);
}

// Wrapper that holds MonLib and the bond list result buffer.
// After calling get_bond_lines(), use bond_data_ptr/bond_data_size
// to read [idx1, idx2, bond_type, ...] triples from WASM memory.
struct BondInfo {
  gemmi::MonLib monlib;
  std::vector<int32_t> bond_data;

  void add_monomer_cif(const std::string& cif_text) {
    gemmi::cif::Document doc = gemmi::cif::read_string(cif_text);
    for (const gemmi::cif::Block& block : doc.blocks)
      monlib.add_monomer_if_present(block);
  }

  void get_bond_lines(gemmi::Structure& st) {
    gemmi::Model& model = st.models.at(0);

    // build atom* -> flat index map
    std::unordered_map<const gemmi::Atom*, int> atom_idx;
    int idx = 0;
    for (auto cra : model.all())
      atom_idx[cra.atom] = idx++;

    bond_data.clear();
    std::unordered_set<uint64_t> seen_pairs;

    auto add_bond = [&](const gemmi::Atom* a1, const gemmi::Atom* a2, int bond_type) {
      auto it1 = atom_idx.find(a1);
      auto it2 = atom_idx.find(a2);
      if (it1 == atom_idx.end() || it2 == atom_idx.end())
        return;
      uint32_t i1 = static_cast<uint32_t>(it1->second);
      uint32_t i2 = static_cast<uint32_t>(it2->second);
      if (i2 < i1)
        std::swap(i1, i2);
      uint64_t key = (uint64_t(i1) << 32) | i2;
      if (!seen_pairs.insert(key).second)
        return;
      bond_data.push_back(static_cast<int32_t>(i1));
      bond_data.push_back(static_cast<int32_t>(i2));
      bond_data.push_back(bond_type);
    };

    // 1. intra-residue bonds from embedded chem_comp data or monomer library
    for (gemmi::Chain& chain : model.chains)
      for (gemmi::Residue& res : chain.residues) {
        const gemmi::ChemComp* cc = find_chemcomp(st, monlib, res);
        if (!cc)
          continue;
        std::string altlocs;
        add_distinct_altlocs(res, altlocs);
        if (altlocs.empty())
          altlocs += '*';
        for (const gemmi::Restraints::Bond& bond : cc->rt.bonds)
          for (char alt : altlocs) {
            const gemmi::Atom* a1 = res.find_atom(bond.id1.atom, alt);
            const gemmi::Atom* a2 = res.find_atom(bond.id2.atom, alt);
            if (a1 && a2) {
              add_bond(a1, a2, static_cast<int>(bond.type));
              if (!a1->altloc && !a2->altloc)
                break;
            }
          }
      }

    // 2. inter-residue bonds from _struct_conn / LINK / SSBOND records
    for (const gemmi::Connection& conn : st.connections) {
      if (conn.asu == gemmi::Asu::Different)
        continue;  // explicit bonds to symmetry images are ignored for now
      const gemmi::Atom* a1 = model.find_atom(conn.partner1);
      const gemmi::Atom* a2 = model.find_atom(conn.partner2);
      if (a1 && a2) {
        int bt = (conn.type == gemmi::Connection::MetalC)
                     ? static_cast<int>(gemmi::BondType::Metal)
                     : static_cast<int>(gemmi::BondType::Single);
        add_bond(a1, a2, bt);
      }
    }

    // 3. inferred polymer links between residues (for example peptide bonds)
    add_inferred_polymer_bonds(st, model, monlib, add_bond);
  }

  uintptr_t bond_data_ptr() const {
    return reinterpret_cast<uintptr_t>(bond_data.data());
  }
  size_t bond_data_size() const { return bond_data.size(); }
};

// Returns bond lines for struct_conn records with asu=Different
// that match the given symmetry image.
// The returned data is [idx1, idx2, bond_type, ...] where idx1 is an atom
// index in model1 (the original) and idx2 is an atom index in model2 (sym image).
struct CrossSymBonds {
  std::vector<int32_t> bond_data;

  void find(const gemmi::Structure& st, const gemmi::NearestImage& image) {
    const gemmi::Model& model = st.models.at(0);
    // build atom index map
    std::unordered_map<const gemmi::Atom*, int> atom_idx;
    int idx = 0;
    for (auto cra : model.all())
      atom_idx[cra.atom] = idx++;

    bond_data.clear();
    for (const gemmi::Connection& conn : st.connections) {
      if (conn.asu != gemmi::Asu::Different)
        continue;
      const gemmi::Atom* a1 = model.find_atom(conn.partner1);
      const gemmi::Atom* a2 = model.find_atom(conn.partner2);
      if (!a1 || !a2)
        continue;
      gemmi::NearestImage conn_image =
          st.cell.find_nearest_image(a1->pos, a2->pos, conn.asu);
      if (conn_image.sym_idx != image.sym_idx ||
          conn_image.pbc_shift[0] != image.pbc_shift[0] ||
          conn_image.pbc_shift[1] != image.pbc_shift[1] ||
          conn_image.pbc_shift[2] != image.pbc_shift[2])
        continue;
      auto it1 = atom_idx.find(a1);
      auto it2 = atom_idx.find(a2);
      if (it1 == atom_idx.end() || it2 == atom_idx.end())
        continue;
      int bt = (conn.type == gemmi::Connection::MetalC)
                   ? static_cast<int>(gemmi::BondType::Metal)
                   : static_cast<int>(gemmi::BondType::Single);
      bond_data.push_back(it1->second);
      bond_data.push_back(it2->second);
      bond_data.push_back(bt);
    }
  }

  uintptr_t bond_data_ptr() const {
    return reinterpret_cast<uintptr_t>(bond_data.data());
  }
  size_t bond_data_size() const { return bond_data.size(); }
};

struct SelectionResult {
  std::vector<int32_t> atom_data;

  void set_atom_indices(gemmi::Structure& st, const std::string& cid,
                        int model_index) {
    gemmi::Selection sel(cid);
    atom_data.clear();
    for (size_t imodel = 0; imodel < st.models.size(); ++imodel) {
      gemmi::Model& model = st.models.at(imodel);
      if (model_index >= 0 && (int) imodel != model_index)
        continue;
      if (!sel.matches(model))
        continue;
      int atom_index = 0;
      for (auto cra : model.all()) {
        if (sel.matches(cra))
          atom_data.push_back(atom_index);
        ++atom_index;
      }
    }
  }

  uintptr_t atom_data_ptr() const {
    return reinterpret_cast<uintptr_t>(atom_data.data());
  }
  size_t atom_data_size() const { return atom_data.size(); }
};

void selection_remove_selected(gemmi::Selection& sel, gemmi::Structure& st) {
  sel.remove_selected(st);
}

void selection_remove_not_selected(gemmi::Selection& sel, gemmi::Structure& st) {
  sel.remove_not_selected(st);
}

void add_mol() {
  em::enum_<gemmi::ResidueSs>("ResidueSs")
    .value("Coil", gemmi::ResidueSs::Coil)
    .value("Helix", gemmi::ResidueSs::Helix)
    .value("Strand", gemmi::ResidueSs::Strand)
    ;

  em::enum_<gemmi::ResidueStrandSense>("ResidueStrandSense")
    .value("NotStrand", gemmi::ResidueStrandSense::NotStrand)
    .value("Parallel", gemmi::ResidueStrandSense::Parallel)
    .value("First", gemmi::ResidueStrandSense::First)
    .value("Antiparallel", gemmi::ResidueStrandSense::Antiparallel)
    ;

  wrap_children<gemmi::Structure>()
    .property("name", &gemmi::Structure::name)
    .property("cell", &gemmi::Structure::cell)
    .function("add_model", &structure_add_model)
    ;

  wrap_children<gemmi::Model>()
    .property("num", &gemmi::Model::num)
    .function("add_chain", &model_add_chain)
    .function("count_occupancies", &gemmi::count_occupancies<gemmi::Model>,
              em::allow_raw_pointers())
    ;

  wrap_children<gemmi::Chain>()
    .property("name", &gemmi::Chain::name)
    .function("add_residue", &chain_add_residue)
    ;

  em::class_<gemmi::ResidueId>("ResidueId")
    .property("seqid_string", &get_seqid_string, &set_seqid_string)
    .property("segment", &gemmi::ResidueId::segment)
    .property("name", &gemmi::ResidueId::name)
    .function("set_seqid", &set_seqid)
    .function("set_seqid_string", &set_seqid_string)
    ;

  wrap_children<gemmi::Residue, em::base<gemmi::ResidueId>>()
    .property("subchain", &gemmi::Residue::subchain)
    .property("ss_from_file", &gemmi::Residue::ss_from_file)
    .property("ss_from_file_string", &get_residue_ss_string)
    .property("strand_sense_from_file", &gemmi::Residue::strand_sense_from_file)
    .property("strand_sense_from_file_string", &get_residue_strand_sense_string)
    .property("entity_type_string", &get_entity_type_string)
    .function("add_atom", &residue_add_atom)
    ;

  em::value_array<gemmi::Position>("Position")
    .element(&gemmi::Position::x)
    .element(&gemmi::Position::y)
    .element(&gemmi::Position::z)
    ;

  em::class_<gemmi::Atom>("Atom")
    .constructor<>()
    .property("name", &gemmi::Atom::name)
    .property("altloc", &gemmi::Atom::altloc)
    .property("charge", &gemmi::Atom::charge)
    .property("element_uname", &element_uname, &set_element)
    .property("is_metal", &atom_is_metal)
    .property("serial", &gemmi::Atom::serial)
    .property("pos", &gemmi::Atom::pos)
    .property("occ", &gemmi::Atom::occ)
    .property("b_iso", &gemmi::Atom::b_iso)
    .function("set_element", &set_element)
    ;

  em::class_<gemmi::Selection>("Selection")
    .constructor<>()
    .constructor<const std::string&>()
    .function("remove_selected", &selection_remove_selected,
              em::allow_raw_pointers())
    .function("remove_not_selected", &selection_remove_not_selected,
              em::allow_raw_pointers())
    ;

  em::class_<BondInfo>("BondInfo")
    .constructor<>()
    .function("add_monomer_cif", &BondInfo::add_monomer_cif)
    .function("get_bond_lines", &BondInfo::get_bond_lines, em::allow_raw_pointers())
    .function("bond_data_ptr", &BondInfo::bond_data_ptr)
    .function("bond_data_size", &BondInfo::bond_data_size)
    ;

  em::class_<CrossSymBonds>("CrossSymBonds")
    .constructor<>()
    .function("find", &CrossSymBonds::find)
    .function("bond_data_ptr", &CrossSymBonds::bond_data_ptr)
    .function("bond_data_size", &CrossSymBonds::bond_data_size)
    ;

  em::class_<SelectionResult>("SelectionResult")
    .constructor<>()
    .function("set_atom_indices", &SelectionResult::set_atom_indices,
              em::allow_raw_pointers())
    .function("atom_data_ptr", &SelectionResult::atom_data_ptr)
    .function("atom_data_size", &SelectionResult::atom_data_size)
    ;

  em::function("get_residue_names", &get_residue_names, em::allow_raw_pointers());
  em::function("get_missing_monomer_names", &get_missing_monomer_names,
               em::allow_raw_pointers());
  em::function("get_nearby_sym_ops", &gemmi::get_nearby_sym_ops);
  em::function("get_sym_image", &gemmi::get_sym_image);
  em::function("make_pdb_string", &make_pdb_string_default);
  em::function("make_mmcif_string", &make_mmcif_string);

  // wrapped in post.js to add default value
  em::function("_read_structure", &read_structure);
}
