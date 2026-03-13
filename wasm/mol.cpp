// Copyright Global Phasing Ltd.

#include "common.h"
#include <unordered_map>
#include <gemmi/model.hpp>
#include <gemmi/select.hpp>   // for Selection
#include <gemmi/mmread.hpp>   // for read_structure_from_memory
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

std::string get_seqid_string(const gemmi::ResidueId& res) {
  return res.seqid.str();
}

std::string get_entity_type_string(const gemmi::Residue& res) {
  return entity_type_to_string(res.entity_type);
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

    // 1. intra-residue bonds from monomer library
    for (gemmi::Chain& chain : model.chains)
      for (gemmi::Residue& res : chain.residues) {
        auto it = monlib.monomers.find(res.name);
        if (it == monlib.monomers.end())
          continue;
        std::string altlocs;
        add_distinct_altlocs(res, altlocs);
        if (altlocs.empty())
          altlocs += '*';
        for (const gemmi::Restraints::Bond& bond : it->second.rt.bonds)
          for (char alt : altlocs) {
            const gemmi::Atom* a1 = res.find_atom(bond.id1.atom, alt);
            const gemmi::Atom* a2 = res.find_atom(bond.id2.atom, alt);
            if (a1 && a2) {
              bond_data.push_back(atom_idx[a1]);
              bond_data.push_back(atom_idx[a2]);
              bond_data.push_back(static_cast<int>(bond.type));
              if (!a1->altloc && !a2->altloc)
                break;
            }
          }
      }

    // 2. inter-residue bonds from _struct_conn / LINK / SSBOND records
    for (const gemmi::Connection& conn : st.connections) {
      const gemmi::Atom* a1 = model.find_atom(conn.partner1);
      const gemmi::Atom* a2 = model.find_atom(conn.partner2);
      if (a1 && a2) {
        auto it1 = atom_idx.find(a1);
        auto it2 = atom_idx.find(a2);
        if (it1 != atom_idx.end() && it2 != atom_idx.end()) {
          bond_data.push_back(it1->second);
          bond_data.push_back(it2->second);
          int bt = (conn.type == gemmi::Connection::MetalC)
                       ? static_cast<int>(gemmi::BondType::Metal)
                       : static_cast<int>(gemmi::BondType::Single);
          bond_data.push_back(bt);
        }
      }
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

void add_mol() {
  wrap_children<gemmi::Structure>()
    .property("name", &gemmi::Structure::name)
    .property("cell", &gemmi::Structure::cell)
    ;

  wrap_children<gemmi::Model>()
    .property("num", &gemmi::Model::num)
    .function("count_occupancies", &gemmi::count_occupancies<gemmi::Model>,
              em::allow_raw_pointers())
    ;

  wrap_children<gemmi::Chain>()
    .property("name", &gemmi::Chain::name)
    ;

  em::class_<gemmi::ResidueId>("ResidueId")
    .property("seqid_string", &get_seqid_string)
    .property("segment", &gemmi::ResidueId::segment)
    .property("name", &gemmi::ResidueId::name)
    ;

  wrap_children<gemmi::Residue, em::base<gemmi::ResidueId>>()
    .property("subchain", &gemmi::Residue::subchain)
    .property("entity_type_string", &get_entity_type_string)
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
    .property("element_uname", &element_uname)
    .property("serial", &gemmi::Atom::serial)
    .property("pos", &gemmi::Atom::pos)
    .property("occ", &gemmi::Atom::occ)
    .property("b_iso", &gemmi::Atom::b_iso)
    ;

  em::class_<gemmi::Selection>("Selection")
    .constructor<>()
    ;

  em::class_<BondInfo>("BondInfo")
    .constructor<>()
    .function("add_monomer_cif", &BondInfo::add_monomer_cif)
    .function("get_bond_lines", &BondInfo::get_bond_lines, em::allow_raw_pointers())
    .function("bond_data_ptr", &BondInfo::bond_data_ptr)
    .function("bond_data_size", &BondInfo::bond_data_size)
    ;

  em::class_<SelectionResult>("SelectionResult")
    .constructor<>()
    .function("set_atom_indices", &SelectionResult::set_atom_indices,
              em::allow_raw_pointers())
    .function("atom_data_ptr", &SelectionResult::atom_data_ptr)
    .function("atom_data_size", &SelectionResult::atom_data_size)
    ;

  em::function("get_residue_names", &get_residue_names, em::allow_raw_pointers());

  // wrapped in post.js to add default value
  em::function("_read_structure", &read_structure);
}
