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
        return cat("<gemmi.AtomAddress ", self.str(), '>');
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
    .value("Branched", EntityType::Branched)
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
        return cat("<gemmi.Connection ", self.name, "  ",
                   self.partner1.str(), " - ", self.partner2.str(), '>');
    });
  py::bind_vector<std::vector<Connection>>(m, "ConnectionList");

  py::class_<CisPep>(m, "CisPep")
    .def(py::init<>())
    .def_readwrite("partner_c", &CisPep::partner_c)
    .def_readwrite("partner_n", &CisPep::partner_n)
    .def_readwrite("model_str", &CisPep::model_str)
    .def_readwrite("only_altloc", &CisPep::only_altloc)
    .def_readwrite("reported_angle", &CisPep::reported_angle)
    ;

  py::class_<ModRes>(m, "ModRes")
    .def(py::init<>())
    .def_readwrite("chain_name", &ModRes::chain_name)
    .def_readwrite("res_id", &ModRes::res_id)
    .def_readwrite("parent_comp_id", &ModRes::parent_comp_id)
    .def_readwrite("mod_id", &ModRes::mod_id)
    .def_readwrite("details", &ModRes::details)
    ;

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

  py::class_<SoftwareItem> softitem(m, "SoftwareItem");
  py::enum_<SoftwareItem::Classification>(softitem, "Classification")
    .value("DataCollection", SoftwareItem::Classification::DataCollection)
    .value("DataExtraction", SoftwareItem::Classification::DataExtraction)
    .value("DataProcessing", SoftwareItem::Classification::DataProcessing)
    .value("DataReduction", SoftwareItem::Classification::DataReduction)
    .value("DataScaling", SoftwareItem::Classification::DataScaling)
    .value("ModelBuilding", SoftwareItem::Classification::ModelBuilding)
    .value("Phasing", SoftwareItem::Classification::Phasing)
    .value("Refinement", SoftwareItem::Classification::Refinement)
    .value("Unspecified", SoftwareItem::Classification::Unspecified)
    ;
  softitem
    .def(py::init<>())
    .def_readwrite("name", &SoftwareItem::name)
    .def_readwrite("version", &SoftwareItem::version)
    .def_readwrite("date", &SoftwareItem::date)
    .def_readwrite("classification", &SoftwareItem::classification)
    ;
  py::class_<ReflectionsInfo>(m, "ReflectionsInfo")
    .def(py::init<>())
    .def_readwrite("resolution_high", &ReflectionsInfo::resolution_high)
    .def_readwrite("resolution_low", &ReflectionsInfo::resolution_low)
    .def_readwrite("completeness", &ReflectionsInfo::completeness)
    .def_readwrite("redundancy", &ReflectionsInfo::redundancy)
    .def_readwrite("r_merge", &ReflectionsInfo::r_merge)
    .def_readwrite("r_sym", &ReflectionsInfo::r_sym)
    .def_readwrite("mean_I_over_sigma", &ReflectionsInfo::mean_I_over_sigma)
    ;
  py::class_<ExperimentInfo>(m, "ExperimentInfo")
    .def(py::init<>())
    .def_readwrite("method", &ExperimentInfo::method)
    .def_readwrite("number_of_crystals", &ExperimentInfo::number_of_crystals)
    .def_readwrite("unique_reflections", &ExperimentInfo::unique_reflections)
    .def_readwrite("reflections", &ExperimentInfo::reflections)
    .def_readwrite("b_wilson", &ExperimentInfo::b_wilson)
    .def_readwrite("shells", &ExperimentInfo::shells)
    .def_readwrite("diffraction_ids", &ExperimentInfo::diffraction_ids)
    ;
  py::class_<DiffractionInfo>(m, "DiffractionInfo")
    .def(py::init<>())
    .def_readwrite("id", &DiffractionInfo::id)
    .def_readwrite("temperature", &DiffractionInfo::temperature)
    .def_readwrite("source", &DiffractionInfo::source)
    .def_readwrite("source_type", &DiffractionInfo::source_type)
    .def_readwrite("synchrotron", &DiffractionInfo::synchrotron)
    .def_readwrite("beamline", &DiffractionInfo::beamline)
    .def_readwrite("wavelengths", &DiffractionInfo::wavelengths)
    .def_readwrite("scattering_type", &DiffractionInfo::scattering_type)
    .def_readwrite("mono_or_laue", &DiffractionInfo::mono_or_laue)
    .def_readwrite("monochromator", &DiffractionInfo::monochromator)
    .def_readwrite("collection_date", &DiffractionInfo::collection_date)
    .def_readwrite("optics", &DiffractionInfo::optics)
    .def_readwrite("detector", &DiffractionInfo::detector)
    .def_readwrite("detector_make", &DiffractionInfo::detector_make)
    ;
  py::class_<CrystalInfo>(m, "CrystalInfo")
    .def(py::init<>())
    .def_readwrite("id", &CrystalInfo::id)
    .def_readwrite("description", &CrystalInfo::description)
    .def_readwrite("ph", &CrystalInfo::ph)
    .def_readwrite("ph_range", &CrystalInfo::ph_range)
    .def_readwrite("diffractions", &CrystalInfo::diffractions)
    ;
  py::class_<TlsGroup> tlsgroup(m, "TlsGroup");
  py::class_<TlsGroup::Selection>(tlsgroup, "Selection")
    .def(py::init<>())
    .def_readwrite("chain", &TlsGroup::Selection::chain)
    .def_readwrite("res_begin", &TlsGroup::Selection::res_begin)
    .def_readwrite("res_end", &TlsGroup::Selection::res_end)
    .def_readwrite("details", &TlsGroup::Selection::details)
    ;
  tlsgroup
    .def(py::init<>())
    .def_readwrite("id", &TlsGroup::id)
    .def_readwrite("selections", &TlsGroup::selections)
    .def_readwrite("origin", &TlsGroup::origin)
    .def_readwrite("T", &TlsGroup::T)
    .def_readwrite("L", &TlsGroup::L)
    .def_readwrite("S", &TlsGroup::S)
    ;
  py::class_<BasicRefinementInfo>(m, "BasicRefinementInfo")
    .def(py::init<>())
    .def_readwrite("resolution_high", &BasicRefinementInfo::resolution_high)
    .def_readwrite("resolution_low", &BasicRefinementInfo::resolution_low)
    .def_readwrite("completeness", &BasicRefinementInfo::completeness)
    .def_readwrite("reflection_count", &BasicRefinementInfo::reflection_count)
    .def_readwrite("rfree_set_count", &BasicRefinementInfo::rfree_set_count)
    .def_readwrite("r_all", &BasicRefinementInfo::r_all)
    .def_readwrite("r_work", &BasicRefinementInfo::r_work)
    .def_readwrite("r_free", &BasicRefinementInfo::r_free)
    ;
  py::class_<RefinementInfo, BasicRefinementInfo> refinfo(m, "RefinementInfo");
  py::class_<RefinementInfo::Restr>(refinfo, "Restr")
    .def(py::init<const std::string&>())
    .def_readwrite("name", &RefinementInfo::Restr::name)
    .def_readwrite("count", &RefinementInfo::Restr::count)
    .def_readwrite("weight", &RefinementInfo::Restr::weight)
    .def_readwrite("function", &RefinementInfo::Restr::function)
    .def_readwrite("dev_ideal", &RefinementInfo::Restr::dev_ideal)
    ;
  refinfo
    .def(py::init<>())
    .def_readwrite("id", &RefinementInfo::id)
    .def_readwrite("cross_validation_method", &RefinementInfo::cross_validation_method)
    .def_readwrite("rfree_selection_method", &RefinementInfo::rfree_selection_method)
    .def_readwrite("bin_count", &RefinementInfo::bin_count)
    .def_readwrite("bins", &RefinementInfo::bins)
    .def_readwrite("mean_b", &RefinementInfo::mean_b)
    .def_readwrite("aniso_b", &RefinementInfo::aniso_b)
    .def_readwrite("luzzati_error", &RefinementInfo::luzzati_error)
    .def_readwrite("dpi_blow_r", &RefinementInfo::dpi_blow_r)
    .def_readwrite("dpi_blow_rfree", &RefinementInfo::dpi_blow_rfree)
    .def_readwrite("dpi_cruickshank_r", &RefinementInfo::dpi_cruickshank_r)
    .def_readwrite("dpi_cruickshank_rfree", &RefinementInfo::dpi_cruickshank_rfree)
    .def_readwrite("cc_fo_fc", &RefinementInfo::cc_fo_fc)
    .def_readwrite("cc_fo_fc_free", &RefinementInfo::cc_fo_fc_free)
    .def_readwrite("restr_stats", &RefinementInfo::restr_stats)
    .def_readwrite("tls_groups", &RefinementInfo::tls_groups)
    .def_readwrite("remarks", &RefinementInfo::remarks)
    ;
  py::class_<Metadata>(m, "Metadata")
    .def_readwrite("authors", &Metadata::authors)
    //.def_readwrite("experiments", &Metadata::experiments)
    //.def_readwrite("crystals", &Metadata::crystals)
    .def_readwrite("refinement", &Metadata::refinement)
    .def_readwrite("software", &Metadata::software)
    //.def_readwrite("solved_by", &Metadata::solved_by)
    //.def_readwrite("starting_model", &Metadata::starting_model)
    //.def_readwrite("remark_300_detail", &Metadata::remark_300_detail)
    ;
}
