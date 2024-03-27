// Copyright Global Phasing Ltd.

#include "common.h"
#include <gemmi/model.hpp>
#include <gemmi/select.hpp>   // for Selection
#include <gemmi/mmread.hpp>   // for read_structure_from_memory
#include <gemmi/polyheur.hpp> // for setup_entities
#include <gemmi/enumstr.hpp>  // for entity_type_to_string
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

void add_mol() {
  wrap_children<gemmi::Structure>()
    .property("name", &gemmi::Structure::name)
    .property("cell", &gemmi::Structure::cell)
    ;

  wrap_children<gemmi::Model>()
    .property("name", &gemmi::Model::name)
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

  // wrapped in post.js to add default value
  em::function("_read_structure", &read_structure);
}
