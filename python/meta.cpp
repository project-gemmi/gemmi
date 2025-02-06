// Copyright 2017 Global Phasing Ltd.

#include "common.h"
#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "gemmi/seqid.hpp"
#include "gemmi/metadata.hpp"
#include "gemmi/enumstr.hpp"
#include "gemmi/sprintf.hpp"
#include "meta.h"

using namespace gemmi;

static std::string repr_end(const void* ptr) {
  char buf[32];
  // avoiding %p, because stb_sprintf formats it without 0x
  snprintf_z(buf, 32, " object at %#zx>", (size_t)ptr);
  return buf;
}

void add_meta(nb::module_& m) {

  // seqid.hpp

  nb::class_<SeqId>(m, "SeqId")
    .def(nb::init<int, char>())
    .def(nb::init<const std::string&>())
    .def_rw("num", &SeqId::num, nb::arg().none())
    .def_rw("icode", &SeqId::icode)
    .def("__str__", &SeqId::str)
    .def("__repr__", [](const SeqId& self) {
        return "<gemmi.SeqId " + self.str() + ">";
    })
    // NOLINTBEGIN(misc-redundant-expression)
    .def(nb::self == nb::self, nb::sig("def __eq__(self, arg: object, /) -> bool"))
    .def(nb::self < nb::self)
    .def("__hash__", [](const SeqId& self) {
        // TODO: implement in hash<gemmi::SeqId> in seqid.hpp?
        return nb::hash(nb::make_tuple(*self.num, self.icode));
    })
    ;

  nb::class_<ResidueId>(m, "ResidueId")
    .def(nb::init<>())
    .def_rw("name", &ResidueId::name)
    .def_rw("seqid", &ResidueId::seqid)
    .def_rw("segment", &ResidueId::segment)
    .def("__str__", &ResidueId::str)
    .def("__repr__", [](const ResidueId& self) {
        return "<gemmi.ResidueId " + self.str() + ">";
    })
    .def(nb::self == nb::self, nb::sig("def __eq__(self, arg: object, /) -> bool"))
    ;

  nb::class_<AtomAddress>(m, "AtomAddress")
    .def(nb::init<>())
    .def(nb::init<const std::string&, const SeqId&, const std::string&,
                  const std::string&, char>(),
         nb::arg("chain"), nb::arg("seqid"), nb::arg("resname"),
         nb::arg("atom"), nb::arg("altloc")='\0')
    .def_rw("chain_name", &AtomAddress::chain_name)
    .def_rw("res_id", &AtomAddress::res_id)
    .def_rw("atom_name", &AtomAddress::atom_name)
    .def_rw("altloc", &AtomAddress::altloc)
    .def("__str__", &AtomAddress::str)
    .def("__repr__", [](const AtomAddress& self) {
        return cat("<gemmi.AtomAddress ", self.str(), '>');
    })
    .def(nb::self == nb::self, nb::sig("def __eq__(self, arg: object, /) -> bool"))
    // NOLINTEND(misc-redundant-expression)
    ;

  // metadata.hpp

  nb::class_<NcsOp>(m, "NcsOp")
    .def(nb::init<>())
    .def("__init__", [](NcsOp* op, const Transform& tr, const std::string& id, bool given) {
      new(op) NcsOp;
      op->tr = tr;
      op->id = id;
      op->given = given;
    }, nb::arg("tr"), nb::arg("id")="", nb::arg("given")=false)
    .def_rw("id", &NcsOp::id)
    .def_rw("given", &NcsOp::given)
    .def_ro("tr", &NcsOp::tr)
    .def("apply", &NcsOp::apply)
    .def("__repr__", [](const NcsOp& self) {
        return cat("<gemmi.NcsOp ", self.id, " |shift|=", std::to_string(self.tr.vec.length()),
                   "( ", (self.given ? "" : "not "), "given)>");
    });
  nb::bind_vector<std::vector<NcsOp>, rv_ri>(m, "NcsOpList");

  nb::enum_<EntityType>(m, "EntityType")
    .value("Unknown", EntityType::Unknown)
    .value("Polymer", EntityType::Polymer)
    .value("NonPolymer", EntityType::NonPolymer)
    .value("Branched", EntityType::Branched)
    .value("Water", EntityType::Water);

  nb::enum_<PolymerType>(m, "PolymerType")
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

  nb::class_<Entity>(m, "Entity")
    .def(nb::init<std::string>())
    .def_rw("name", &Entity::name)
    .def_rw("subchains", &Entity::subchains)
    .def_rw("entity_type", &Entity::entity_type)
    .def_rw("polymer_type", &Entity::polymer_type)
    .def_rw("sifts_unp_acc", &Entity::sifts_unp_acc)
    .def_rw("full_sequence", &Entity::full_sequence)
    .def_static("first_mon", &Entity::first_mon)
    .def("__repr__", [](const Entity& self) {
      std::string s = cat("<gemmi.Entity '", self.name, "' ");
      s += entity_type_to_string(self.entity_type);
      if (self.polymer_type != PolymerType::Unknown)
        cat_to(s, ' ', polymer_type_to_string(self.polymer_type));
      return s + repr_end(&self);
    });
  nb::bind_vector<std::vector<Entity>, rv_ri>(m, "EntityList");

  nb::enum_<Connection::Type>(m, "ConnectionType")
    .value("Covale", Connection::Type::Covale)
    .value("Disulf", Connection::Type::Disulf)
    .value("Hydrog", Connection::Type::Hydrog)
    .value("MetalC", Connection::Type::MetalC)
    .value("Unknown", Connection::Type::Unknown);

  nb::class_<Connection>(m, "Connection")
    .def(nb::init<>())
    .def_rw("name", &Connection::name)
    .def_rw("link_id", &Connection::link_id)
    .def_rw("type", &Connection::type)
    .def_rw("asu", &Connection::asu)
    .def_rw("partner1", &Connection::partner1)
    .def_rw("partner2", &Connection::partner2)
    .def_rw("reported_distance", &Connection::reported_distance)
    .def("__repr__", [](const Connection& self) {
        return cat("<gemmi.Connection ", self.name, "  ",
                   self.partner1.str(), " - ", self.partner2.str(), '>');
    });
  nb::bind_vector<std::vector<Connection>, rv_ri>(m, "ConnectionList");

  nb::class_<CisPep>(m, "CisPep")
    .def(nb::init<>())
    .def_rw("partner_c", &CisPep::partner_c)
    .def_rw("partner_n", &CisPep::partner_n)
    .def_rw("model_num", &CisPep::model_num)
    .def_rw("only_altloc", &CisPep::only_altloc)
    .def_rw("reported_angle", &CisPep::reported_angle)
    ;

  nb::class_<ModRes>(m, "ModRes")
    .def(nb::init<>())
    .def_rw("chain_name", &ModRes::chain_name)
    .def_rw("res_id", &ModRes::res_id)
    .def_rw("parent_comp_id", &ModRes::parent_comp_id)
    .def_rw("mod_id", &ModRes::mod_id)
    .def_rw("details", &ModRes::details)
    ;

  nb::class_<Helix> helix(m, "Helix");
  nb::enum_<Helix::HelixClass>(helix, "HelixClass")
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
    .def(nb::init<>())
    .def_rw("start", &Helix::start)
    .def_rw("end", &Helix::end)
    .def_rw("pdb_helix_class", &Helix::pdb_helix_class)
    .def_rw("length", &Helix::length);
  nb::bind_vector<std::vector<Helix>, rv_ri>(m, "HelixList");

  nb::class_<Sheet> sheet(m, "Sheet");
  nb::class_<Sheet::Strand>(sheet, "Strand")
    .def(nb::init<>())
    .def_rw("start", &Sheet::Strand::start)
    .def_rw("end", &Sheet::Strand::end)
    .def_rw("hbond_atom2", &Sheet::Strand::hbond_atom2)
    .def_rw("hbond_atom1", &Sheet::Strand::hbond_atom1)
    .def_rw("sense", &Sheet::Strand::sense)
    .def_rw("name", &Sheet::Strand::name);
  nb::bind_vector<std::vector<Sheet::Strand>, rv_ri>(sheet, "StrandList");

  sheet
    .def(nb::init<std::string>())
    .def_rw("name", &Sheet::name)
    .def_rw("strands", &Sheet::strands);
  nb::bind_vector<std::vector<Sheet>, rv_ri>(m, "SheetList");

  nb::class_<Assembly> assembly(m, "Assembly");
  nb::class_<Assembly::Operator>(assembly, "Operator")
    .def(nb::init<>())
    .def_rw("name", &Assembly::Operator::name)
    .def_rw("type", &Assembly::Operator::type)
    .def_rw("transform", &Assembly::Operator::transform);
  nb::bind_vector<std::vector<Assembly::Operator>, rv_ri>(assembly, "OperatorList");

  nb::class_<Assembly::Gen>(assembly, "Gen")
    .def(nb::init<>())
    .def_rw("chains", &Assembly::Gen::chains)
    .def_rw("subchains", &Assembly::Gen::subchains)
    .def_ro("operators", &Assembly::Gen::operators);
  nb::bind_vector<std::vector<Assembly::Gen>, rv_ri>(assembly, "GenList");

  nb::enum_<Assembly::SpecialKind>(m, "AssemblySpecialKind")
    .value("NA", Assembly::SpecialKind::NA)
    .value("CompleteIcosahedral", Assembly::SpecialKind::CompleteIcosahedral)
    .value("RepresentativeHelical", Assembly::SpecialKind::RepresentativeHelical)
    .value("CompletePoint", Assembly::SpecialKind::CompletePoint);

  assembly
    .def(nb::init<const std::string&>())
    .def_rw("name", &Assembly::name)
    .def_rw("author_determined", &Assembly::author_determined)
    .def_rw("software_determined", &Assembly::software_determined)
    .def_rw("oligomeric_details", &Assembly::oligomeric_details)
    .def_ro("generators", &Assembly::generators)
    .def_rw("special_kind", &Assembly::special_kind)
    .def("__repr__", [](const Assembly& self) {
        return cat("<gemmi.Assembly ", self.name, ">");
    });
    ;
  nb::bind_vector<std::vector<Assembly>, rv_ri>(m, "AssemblyList");

  nb::class_<SoftwareItem> softitem(m, "SoftwareItem");
  nb::enum_<SoftwareItem::Classification>(softitem, "Classification")
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
    .def(nb::init<>())
    .def_rw("name", &SoftwareItem::name)
    .def_rw("version", &SoftwareItem::version)
    .def_rw("date", &SoftwareItem::date)
    .def_rw("classification", &SoftwareItem::classification)
    .def_rw("contact_author", &SoftwareItem::contact_author)
    .def_rw("contact_author_email", &SoftwareItem::contact_author_email)
    ;
  nb::class_<ReflectionsInfo>(m, "ReflectionsInfo")
    .def(nb::init<>())
    .def_rw("resolution_high", &ReflectionsInfo::resolution_high)
    .def_rw("resolution_low", &ReflectionsInfo::resolution_low)
    .def_rw("completeness", &ReflectionsInfo::completeness)
    .def_rw("redundancy", &ReflectionsInfo::redundancy)
    .def_rw("r_merge", &ReflectionsInfo::r_merge)
    .def_rw("r_sym", &ReflectionsInfo::r_sym)
    .def_rw("mean_I_over_sigma", &ReflectionsInfo::mean_I_over_sigma)
    ;
  nb::class_<ExperimentInfo>(m, "ExperimentInfo")
    .def(nb::init<>())
    .def_rw("method", &ExperimentInfo::method)
    .def_rw("number_of_crystals", &ExperimentInfo::number_of_crystals)
    .def_rw("unique_reflections", &ExperimentInfo::unique_reflections)
    .def_rw("reflections", &ExperimentInfo::reflections)
    .def_rw("b_wilson", &ExperimentInfo::b_wilson)
    .def_rw("shells", &ExperimentInfo::shells)
    .def_rw("diffraction_ids", &ExperimentInfo::diffraction_ids)
    ;
  nb::class_<DiffractionInfo>(m, "DiffractionInfo")
    .def(nb::init<>())
    .def_rw("id", &DiffractionInfo::id)
    .def_rw("temperature", &DiffractionInfo::temperature)
    .def_rw("source", &DiffractionInfo::source)
    .def_rw("source_type", &DiffractionInfo::source_type)
    .def_rw("synchrotron", &DiffractionInfo::synchrotron)
    .def_rw("beamline", &DiffractionInfo::beamline)
    .def_rw("wavelengths", &DiffractionInfo::wavelengths)
    .def_rw("scattering_type", &DiffractionInfo::scattering_type)
    .def_rw("mono_or_laue", &DiffractionInfo::mono_or_laue)
    .def_rw("monochromator", &DiffractionInfo::monochromator)
    .def_rw("collection_date", &DiffractionInfo::collection_date)
    .def_rw("optics", &DiffractionInfo::optics)
    .def_rw("detector", &DiffractionInfo::detector)
    .def_rw("detector_make", &DiffractionInfo::detector_make)
    ;
  nb::class_<CrystalInfo>(m, "CrystalInfo")
    .def(nb::init<>())
    .def_rw("id", &CrystalInfo::id)
    .def_rw("description", &CrystalInfo::description)
    .def_rw("ph", &CrystalInfo::ph)
    .def_rw("ph_range", &CrystalInfo::ph_range)
    .def_rw("diffractions", &CrystalInfo::diffractions)
    ;
  nb::class_<TlsGroup> tlsgroup(m, "TlsGroup");
  nb::class_<TlsGroup::Selection>(tlsgroup, "Selection")
    .def(nb::init<>())
    .def_rw("chain", &TlsGroup::Selection::chain)
    .def_rw("res_begin", &TlsGroup::Selection::res_begin)
    .def_rw("res_end", &TlsGroup::Selection::res_end)
    .def_rw("details", &TlsGroup::Selection::details)
    ;
  tlsgroup
    .def(nb::init<>())
    .def_rw("id", &TlsGroup::id)
    .def_rw("selections", &TlsGroup::selections)
    .def_rw("origin", &TlsGroup::origin)
    .def_rw("T", &TlsGroup::T)
    .def_rw("L", &TlsGroup::L)
    .def_rw("S", &TlsGroup::S)
    ;
  nb::class_<BasicRefinementInfo>(m, "BasicRefinementInfo")
    .def(nb::init<>())
    .def_rw("resolution_high", &BasicRefinementInfo::resolution_high)
    .def_rw("resolution_low", &BasicRefinementInfo::resolution_low)
    .def_rw("completeness", &BasicRefinementInfo::completeness)
    .def_rw("reflection_count", &BasicRefinementInfo::reflection_count)
    .def_rw("work_set_count", &BasicRefinementInfo::work_set_count)
    .def_rw("rfree_set_count", &BasicRefinementInfo::rfree_set_count)
    .def_rw("r_all", &BasicRefinementInfo::r_all)
    .def_rw("r_work", &BasicRefinementInfo::r_work)
    .def_rw("r_free", &BasicRefinementInfo::r_free)
    .def_rw("cc_fo_fc_work", &RefinementInfo::cc_fo_fc_work)
    .def_rw("cc_fo_fc_free", &RefinementInfo::cc_fo_fc_free)
    .def_rw("fsc_work", &RefinementInfo::fsc_work)
    .def_rw("fsc_free", &RefinementInfo::fsc_free)
    .def_rw("cc_intensity_work", &RefinementInfo::cc_intensity_work)
    .def_rw("cc_intensity_free", &RefinementInfo::cc_intensity_free)
    ;
  nb::class_<RefinementInfo, BasicRefinementInfo> refinfo(m, "RefinementInfo");
  nb::class_<RefinementInfo::Restr>(refinfo, "Restr")
    .def(nb::init<const std::string&>())
    .def_rw("name", &RefinementInfo::Restr::name)
    .def_rw("count", &RefinementInfo::Restr::count)
    .def_rw("weight", &RefinementInfo::Restr::weight)
    .def_rw("function", &RefinementInfo::Restr::function)
    .def_rw("dev_ideal", &RefinementInfo::Restr::dev_ideal)
    ;
  refinfo
    .def(nb::init<>())
    .def_rw("id", &RefinementInfo::id)
    .def_rw("cross_validation_method", &RefinementInfo::cross_validation_method)
    .def_rw("rfree_selection_method", &RefinementInfo::rfree_selection_method)
    .def_rw("bin_count", &RefinementInfo::bin_count)
    .def_rw("bins", &RefinementInfo::bins)
    .def_rw("mean_b", &RefinementInfo::mean_b)
    .def_rw("aniso_b", &RefinementInfo::aniso_b)
    .def_rw("luzzati_error", &RefinementInfo::luzzati_error)
    .def_rw("dpi_blow_r", &RefinementInfo::dpi_blow_r)
    .def_rw("dpi_blow_rfree", &RefinementInfo::dpi_blow_rfree)
    .def_rw("dpi_cruickshank_r", &RefinementInfo::dpi_cruickshank_r)
    .def_rw("dpi_cruickshank_rfree", &RefinementInfo::dpi_cruickshank_rfree)
    .def_rw("restr_stats", &RefinementInfo::restr_stats)
    .def_rw("tls_groups", &RefinementInfo::tls_groups)
    .def_rw("remarks", &RefinementInfo::remarks)
    ;
  nb::class_<Metadata>(m, "Metadata")
    .def_rw("authors", &Metadata::authors)
    .def_rw("experiments", &Metadata::experiments)
    .def_rw("crystals", &Metadata::crystals)
    .def_rw("refinement", &Metadata::refinement)
    .def_rw("software", &Metadata::software)
    .def_rw("solved_by", &Metadata::solved_by)
    .def_rw("starting_model", &Metadata::starting_model)
    .def_rw("remark_300_detail", &Metadata::remark_300_detail)
    ;
}
