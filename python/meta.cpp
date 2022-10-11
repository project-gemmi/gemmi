// Copyright 2017 Global Phasing Ltd.

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <ostream>

#include "gemmi/seqid.hpp"
#include "gemmi/metadata.hpp"
#include "gemmi/enumstr.hpp"
#include "tostr.hpp"
#include "meta.h"

namespace py = pybind11;
using namespace gemmi;

namespace gemmi {
  // not inline
  std::ostream& operator<< (std::ostream& os, const Entity& ent) {
    os << "<gemmi.Entity '" << ent.name << "' "
       << entity_type_to_string(ent.entity_type);
    if (ent.polymer_type != PolymerType::Unknown)
      os << ' ' << polymer_type_to_string(ent.polymer_type);
    os << " object at 0x" << std::hex << (size_t)&ent << std::dec << '>';
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
    })
    .def(py::self == py::self)
    .def(py::self < py::self)
    .def(py::pickle(
        [](const SeqId &self) {
            return py::make_tuple(self.num, self.icode);
        },
        [](py::tuple t) {
            if (t.size() != 2)
                throw std::runtime_error("invalid tuple size");
            return SeqId(t[0].cast<int>(), t[1].cast<char>());
        }
    ));

  py::class_<ResidueId>(m, "ResidueId")
    .def(py::init<>())
    .def_readwrite("name", &ResidueId::name)
    .def_readwrite("seqid", &ResidueId::seqid)
    .def_readwrite("segment", &ResidueId::segment)
    .def("__str__", &ResidueId::str)
    .def("__repr__", [](const ResidueId& self) {
        return "<gemmi.ResidueId " + self.str() + ">";
    })
    .def(py::self == py::self)
    .def(py::pickle(
        [](const ResidueId &self) {
            return py::make_tuple(self.seqid, self.segment, self.name);
        },
        [](py::tuple t) {
            if (t.size() != 3)
                throw std::runtime_error("invalid tuple size");
            ResidueId r;
            r.seqid = t[0].cast<SeqId>();
            r.segment = t[1].cast<std::string>();
            r.name = t[2].cast<std::string>();
            return r;
        }
    ));


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
    })
    .def(py::self == py::self)
    .def(py::pickle(
        [](const AtomAddress &self) {
            return py::make_tuple(self.chain_name, self.res_id, self.atom_name, self.altloc);
        },
        [](py::tuple t) {
            if (t.size() != 4)
                throw std::runtime_error("invalid tuple size");
            return AtomAddress(t[0].cast<std::string>(), t[1].cast<ResidueId>(),
                               t[2].cast<std::string>(), t[3].cast<char>());
        }
    ));
  // metadata.hpp

  py::class_<NcsOp>(m, "NcsOp")
    .def(py::init<>())
    .def(py::init([](const Transform& tr, const std::string& id, bool given) {
      NcsOp* op = new NcsOp();
      op->tr = tr;
      op->id = id;
      op->given = given;
      return op;
    }), py::arg("tr"), py::arg("id")="", py::arg("given")=false)
    .def_readwrite("id", &NcsOp::id)
    .def_readwrite("given", &NcsOp::given)
    .def_readonly("tr", &NcsOp::tr)
    .def("apply", &NcsOp::apply)
    .def("__repr__", [](const NcsOp& self) {
        return tostr("<gemmi.NcsOp ", self.id,
                     " |shift|=", self.tr.vec.length(),
                     (self.given ? " (" : " (not "), "given)>");
    });
  py::bind_vector<std::vector<NcsOp>>(m, "NcsOpList");

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
    .def_readwrite("sifts_unp_acc", &Entity::sifts_unp_acc)
    .def_readwrite("full_sequence", &Entity::full_sequence)
    .def_static("first_mon", &Entity::first_mon)
    .def("__repr__", [](const Entity& self) { return tostr(self); });
  py::bind_vector<std::vector<Entity>>(m, "EntityList");

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
  py::bind_vector<std::vector<Connection>>(m, "ConnectionList");

  py::class_<Helix> helix(m, "Helix");
  py::enum_<Helix::HelixClass>(helix, "HelixClass")
    .value("UnknownHelix", Helix::HelixClass::UnknownHelix)
    .value("RAlpha", Helix::HelixClass::RAlpha)
    .value("ROmega", Helix::HelixClass::ROmega)
    .value("RPi", Helix::HelixClass::RPi)
    .value("RGamma", Helix::HelixClass::RGamma)
    .value("R310", Helix::HelixClass::R310)
    .value("LAlpha", Helix::HelixClass::LAlpha)
    .value("LOmega", Helix::HelixClass::LOmega)
    .value("LGamma", Helix::HelixClass::LGamma)
    .value("Helix27", Helix::HelixClass::Helix27)
    .value("HelixPolyProlineNone", Helix::HelixClass::HelixPolyProlineNone);

  helix
    .def(py::init<>())
    .def_readwrite("start", &Helix::start)
    .def_readwrite("end", &Helix::end)
    .def_readwrite("pdb_helix_class", &Helix::pdb_helix_class)
    .def_readwrite("length", &Helix::length);
  py::bind_vector<std::vector<Helix>>(m, "HelixList");

  py::class_<Sheet> sheet(m, "Sheet");
  py::class_<Sheet::Strand>(sheet, "Strand")
    .def(py::init<>())
    .def_readwrite("start", &Sheet::Strand::start)
    .def_readwrite("end", &Sheet::Strand::end)
    .def_readwrite("hbond_atom2", &Sheet::Strand::hbond_atom2)
    .def_readwrite("hbond_atom1", &Sheet::Strand::hbond_atom1)
    .def_readwrite("sense", &Sheet::Strand::sense)
    .def_readwrite("name", &Sheet::Strand::name);
  py::bind_vector<std::vector<Sheet::Strand>>(sheet, "StrandList");

  sheet
    .def(py::init<std::string>())
    .def_readwrite("name", &Sheet::name)
    .def_readwrite("strands", &Sheet::strands);
  py::bind_vector<std::vector<Sheet>>(m, "SheetList");

  py::class_<Assembly> assembly(m, "Assembly");
  py::class_<Assembly::Operator>(assembly, "Operator")
    .def(py::init<>())
    .def_readwrite("name", &Assembly::Operator::name)
    .def_readwrite("type", &Assembly::Operator::type)
    .def_readwrite("transform", &Assembly::Operator::transform);
  py::bind_vector<std::vector<Assembly::Operator>>(assembly, "OperatorList");

  py::class_<Assembly::Gen>(assembly, "Gen")
    .def(py::init<>())
    .def_readwrite("chains", &Assembly::Gen::chains)
    .def_readwrite("subchains", &Assembly::Gen::subchains)
    .def_readonly("operators", &Assembly::Gen::operators);
  py::bind_vector<std::vector<Assembly::Gen>>(assembly, "GenList");

  py::enum_<Assembly::SpecialKind>(m, "AssemblySpecialKind")
    .value("NA", Assembly::SpecialKind::NA)
    .value("CompleteIcosahedral", Assembly::SpecialKind::CompleteIcosahedral)
    .value("RepresentativeHelical", Assembly::SpecialKind::RepresentativeHelical)
    .value("CompletePoint", Assembly::SpecialKind::CompletePoint);

  assembly
    .def(py::init<const std::string&>())
    .def_readwrite("name", &Assembly::name)
    .def_readwrite("author_determined", &Assembly::author_determined)
    .def_readwrite("software_determined", &Assembly::software_determined)
    .def_readwrite("oligomeric_details", &Assembly::oligomeric_details)
    .def_readonly("generators", &Assembly::generators)
    .def_readwrite("special_kind", &Assembly::special_kind)
    ;
  py::bind_vector<std::vector<Assembly>>(m, "AssemblyList");

  py::class_<Metadata>(m, "Metadata")
    .def_readwrite("authors", &Metadata::authors)
    ;
}
