// Copyright 2017 Global Phasing Ltd.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <ostream>

#include "gemmi/seqid.hpp"
#include "gemmi/metadata.hpp"
#include "gemmi/entstr.hpp"
#include "gemmi/tostr.hpp"
#include "meta.h"

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Assembly::Gen>)
PYBIND11_MAKE_OPAQUE(std::vector<Assembly::Operator>)

namespace gemmi {
  // not inline
  std::ostream& operator<< (std::ostream& os, const Entity& ent) {
    os << "<gemmi.Entity '" << ent.name << "' "
       << entity_type_to_string(ent.entity_type);
    if (ent.polymer_type != PolymerType::Unknown)
      os << ' ' << polymer_type_to_qstring(ent.polymer_type);
    os << " object at " << (void*)&ent << '>';
    return os;
  }
}

void add_meta(py::module& m) {

  // seqid.hpp

  py::class_<SeqId>(m, "SeqId")
    .def(py::init<int, char>())
    .def(py::init<const std::string&>())
    .def_readwrite("num", &SeqId::num)
    .def_readwrite("icode", &SeqId::icode)
    .def("__str__", &SeqId::str)
    .def("__repr__", [](const SeqId& self) {
        return "<gemmi.SeqId " + self.str() + ">";
    });

  py::class_<ResidueId>(m, "ResidueId")
    .def(py::init<>())
    .def_readwrite("name", &ResidueId::name)
    .def_readwrite("seqid", &ResidueId::seqid)
    .def_readwrite("segment", &ResidueId::segment)
    .def("__str__", &ResidueId::str)
    .def("__repr__", [](const ResidueId& self) {
        return "<gemmi.ResidueId " + self.str() + ">";
    });

  py::class_<AtomAddress>(m, "AtomAddress")
    .def(py::init<>())
    .def(py::init<const std::string&, const SeqId&, const std::string&,
                  const std::string&, char>(),
         py::arg("chain"), py::arg("seqid"), py::arg("resname"),
         py::arg("atom"), py::arg("altloc")='\0')
    .def_readwrite("chain_name", &AtomAddress::chain_name)
    .def_readwrite("res_id", &AtomAddress::res_id)
    .def_readwrite("atom_name", &AtomAddress::atom_name)
    .def_readwrite("altloc", &AtomAddress::altloc)
    .def("__str__", &AtomAddress::str)
    .def("__repr__", [](const AtomAddress& self) {
        return tostr("<gemmi.AtomAddress ", self.str(), '>');
    });


  // metadata.hpp

  py::class_<NcsOp>(m, "NcsOp")
    .def(py::init<>())
    .def_readwrite("id", &NcsOp::id)
    .def_readwrite("given", &NcsOp::given)
    .def_readonly("tr", &NcsOp::tr)
    .def("apply", &NcsOp::apply)
    .def("__repr__", [](const NcsOp& self) {
        return tostr("<gemmi.NcsOp ", self.id,
                     " |shift|=", self.tr.vec.length(),
                     (self.given ? " (" : " (not "), "given)>");
    });

  py::enum_<EntityType>(m, "EntityType")
    .value("Unknown", EntityType::Unknown)
    .value("Polymer", EntityType::Polymer)
    .value("NonPolymer", EntityType::NonPolymer)
    .value("Water", EntityType::Water);

  py::enum_<PolymerType>(m, "PolymerType")
    .value("PeptideL", PolymerType::PeptideL)
    .value("PeptideD", PolymerType::PeptideD)
    .value("Dna", PolymerType::Dna)
    .value("Rna", PolymerType::Rna)
    .value("DnaRnaHybrid", PolymerType::DnaRnaHybrid)
    .value("SaccharideD", PolymerType::SaccharideD)
    .value("SaccharideL", PolymerType::SaccharideL)
    .value("Pna", PolymerType::Pna)
    .value("CyclicPseudoPeptide", PolymerType::CyclicPseudoPeptide)
    .value("Other", PolymerType::Other)
    .value("Unknown", PolymerType::Unknown);

  py::class_<Entity>(m, "Entity")
    .def(py::init<std::string>())
    .def_readwrite("name", &Entity::name)
    .def_readwrite("subchains", &Entity::subchains)
    .def_readwrite("entity_type", &Entity::entity_type)
    .def_readwrite("polymer_type", &Entity::polymer_type)
    .def_readwrite("full_sequence", &Entity::full_sequence)
    .def_static("first_mon", &Entity::first_mon)
    .def("__repr__", [](const Entity& self) { return tostr(self); });


  py::enum_<Connection::Type>(m, "ConnectionType")
    .value("Covale", Connection::Type::Covale)
    .value("Disulf", Connection::Type::Disulf)
    .value("Hydrog", Connection::Type::Hydrog)
    .value("MetalC", Connection::Type::MetalC)
    .value("Unknown", Connection::Type::Unknown);

  py::class_<Connection>(m, "Connection")
    .def(py::init<>())
    .def_readwrite("name", &Connection::name)
    .def_readwrite("link_id", &Connection::link_id)
    .def_readwrite("type", &Connection::type)
    .def_readwrite("asu", &Connection::asu)
    .def_readwrite("partner1", &Connection::partner1)
    .def_readwrite("partner2", &Connection::partner2)
    .def_readwrite("reported_distance", &Connection::reported_distance)
    .def("__repr__", [](const Connection& self) {
        return tostr("<gemmi.Connection ", self.name, "  ",
                     self.partner1.str(), " - ", self.partner2.str(), '>');
    });


  py::class_<Assembly> assembly(m, "Assembly");
  py::bind_vector<std::vector<Assembly::Gen>>(assembly, "GenList");
  py::bind_vector<std::vector<Assembly::Operator>>(assembly, "OperatorList");
  py::class_<Assembly::Operator>(assembly, "Operator")
    .def(py::init<>())
    .def_readonly("name", &Assembly::Operator::name)
    .def_readonly("type", &Assembly::Operator::type)
    .def_readonly("transform", &Assembly::Operator::transform);

  py::class_<Assembly::Gen>(assembly, "Gen")
    .def(py::init<>())
    .def_readonly("chains", &Assembly::Gen::chains)
    .def_readonly("subchains", &Assembly::Gen::subchains)
    .def_readonly("operators", &Assembly::Gen::operators);

  assembly
    .def_readonly("name", &Assembly::name)
    .def_readonly("oligomeric_details", &Assembly::oligomeric_details)
    .def_readonly("generators", &Assembly::generators)
    ;
}
