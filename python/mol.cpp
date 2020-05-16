// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/entstr.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/model.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/polyheur.hpp"
#include "gemmi/assembly.hpp"
#include "gemmi/labelseq.hpp"   // for setup_for_mmcif
#include "gemmi/to_pdb.hpp"
#include "gemmi/to_mmcif.hpp"
#include "gemmi/tostr.hpp"
#include "gemmi/ofstream.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Connection>)
PYBIND11_MAKE_OPAQUE(std::vector<NcsOp>)
PYBIND11_MAKE_OPAQUE(std::vector<Entity>)
PYBIND11_MAKE_OPAQUE(std::vector<Assembly>)
using info_map_type = std::map<std::string, std::string>;
PYBIND11_MAKE_OPAQUE(info_map_type)

namespace gemmi {
  inline std::ostream& operator<< (std::ostream& os, const Entity& ent) {
    os << "<gemmi.Entity '" << ent.name << "' "
       << entity_type_to_string(ent.entity_type);
    if (ent.polymer_type != PolymerType::Unknown)
      os << ' ' << polymer_type_to_qstring(ent.polymer_type);
    os << " object at " << (void*)&ent << '>';
    return os;
  }
}

namespace pybind11 { namespace detail {
  template<> struct type_caster<SeqId::OptionalNum>
    : optional_caster<SeqId::OptionalNum> {};
}} // namespace pybind11::detail

template<typename T, typename C>
C& add_item(T& container, C child, int pos) {
  if ((size_t) pos > container.size()) // true also for negative pos
    pos = (int) container.size();
  return *container.insert(container.begin() + pos, std::move(child));
}

template<typename P, typename C>
C& add_child(P& parent, C child, int pos) {
  return add_item(parent.children(), std::move(child), pos);
}

template<typename T> int normalize_index(int index, const T& container) {
  if (index < 0)
    index += (int) container.size();
  if ((size_t) index >= container.size())
    throw py::index_error();
  return index;
}

template<typename P, typename C> C& get_child(P& parent, int index) {
  auto& children = parent.children();
  return children[normalize_index(index, children)];
}

template<typename P> void remove_child(P& parent, int index) {
  auto& children = parent.children();
  children.erase(children.begin() + normalize_index(index, children));
}

template<typename P> void remove_children(P& parent, py::slice slice) {
  ssize_t start, stop, step, slength;
  auto& children = parent.children();
  if (!slice.compute(children.size(), &start, &stop, &step, &slength))
    throw py::error_already_set();
  for (int i = 0; i < slength; ++i) {
    ssize_t idx = start + (step > 0 ? slength - 1 - i : i) * step;
    children.erase(children.begin() + idx);
  }
}

static std::vector<std::string>
expand_protein_one_letter_string(const std::string& s) {
  std::vector<std::string> r;
  r.reserve(s.size());
  for (char c : s)
    r.push_back(expand_protein_one_letter(c));
  return r;
}

void add_mol(py::module& m) {

  // "Forward declaration" of python classes to avoid
  // C++ signatures in docstrings
  py::class_<AtomAddress> pyAtomAddress(m, "AtomAddress");
  py::class_<ResidueId> pyResidueId(m, "ResidueId");
  py::class_<Atom> pyAtom(m, "Atom");
  py::class_<Residue, ResidueId> pyResidue(m, "Residue");
  py::class_<Chain> pyChain(m, "Chain");
  py::class_<Model> pyModel(m, "Model");
  py::class_<SeqId> pySeqId(m, "SeqId");
  py::class_<ResidueSpan> pyResidueSpan(m, "ResidueSpan");
  py::class_<ResidueGroup, ResidueSpan> pyResidueGroup(m, "ResidueGroup");
  py::class_<CraProxy> pyCraGenerator(m, "CraGenerator");

  py::class_<ResidueInfo>(m, "ResidueInfo")
    .def_readonly("one_letter_code", &ResidueInfo::one_letter_code)
    .def_readonly("hydrogen_count", &ResidueInfo::hydrogen_count)
    .def_readonly("weight", &ResidueInfo::weight)
    .def("found", &ResidueInfo::found)
    .def("is_standard", &ResidueInfo::is_standard)
    .def("is_water", &ResidueInfo::is_water)
    .def("is_nucleic_acid", &ResidueInfo::is_nucleic_acid)
    .def("is_amino_acid", &ResidueInfo::is_amino_acid);

  m.def("find_tabulated_residue", &find_tabulated_residue, py::arg("name"),
        "Find chemical component information in the internal table.");
  m.def("expand_protein_one_letter", &expand_protein_one_letter);
  m.def("expand_protein_one_letter_string", &expand_protein_one_letter_string);
  m.def("one_letter_code",
        (std::string (*)(const std::vector<std::string>&)) &one_letter_code);
  m.def("one_letter_code",
        [](const ResidueSpan& span) { return one_letter_code(span); });

  py::enum_<CoorFormat>(m, "CoorFormat")
    .value("Unknown", CoorFormat::Unknown)
    .value("UnknownAny", CoorFormat::UnknownAny)
    .value("Pdb", CoorFormat::Pdb)
    .value("Mmcif", CoorFormat::Mmcif)
    .value("Mmjson", CoorFormat::Mmjson)
    .value("ChemComp", CoorFormat::ChemComp);

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

  py::enum_<HowToNameCopiedChains>(m, "HowToNameCopiedChains")
    .value("Short", HowToNameCopiedChains::Short)
    .value("AddNumber", HowToNameCopiedChains::AddNumber)
    .value("Dup", HowToNameCopiedChains::Dup)
    ;

  py::class_<Entity>(m, "Entity")
    .def(py::init<std::string>())
    .def_readwrite("name", &Entity::name)
    .def_readwrite("subchains", &Entity::subchains)
    .def_readwrite("entity_type", &Entity::entity_type)
    .def_readwrite("polymer_type", &Entity::polymer_type)
    .def_readwrite("full_sequence", &Entity::full_sequence)
    .def_static("first_mon", &Entity::first_mon)
    .def("__repr__", [](const Entity& self) { return tostr(self); });

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

  py::class_<Assembly>(m, "Assembly")
    .def_readonly("name", &Assembly::name)
    .def_readonly("oligomeric_details", &Assembly::oligomeric_details)
    .def("make_assembly", [](const Assembly& self, const Model& model,
                             HowToNameCopiedChains how) {
        return make_assembly(self, model, how, nullptr);
    })
    ;

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
  py::bind_vector<std::vector<NcsOp>>(m, "NcsOpList");
  py::bind_vector<std::vector<Entity>>(m, "EntityList");
  py::bind_vector<std::vector<Assembly>>(m, "AssemblyList");
  py::bind_map<info_map_type>(m, "InfoMap");

  py::class_<Structure>(m, "Structure")
    .def(py::init<>())
    .def_readwrite("name", &Structure::name)
    .def_readwrite("cell", &Structure::cell)
    .def_readwrite("spacegroup_hm", &Structure::spacegroup_hm)
    .def_readwrite("ncs", &Structure::ncs)
    .def_readwrite("resolution", &Structure::resolution)
    .def_readwrite("entities", &Structure::entities)
    .def_readwrite("connections", &Structure::connections)
    .def_readwrite("assemblies", &Structure::assemblies)
    .def_readwrite("info", &Structure::info)
    .def("get_entity",
         (Entity* (Structure::*)(const std::string&)) &Structure::get_entity,
         py::arg("subchain"), py::return_value_policy::reference_internal)
    .def("get_entity_of", [](Structure& st, ResidueSpan& span) {
        return st.get_entity_of(span);
    }, py::arg("subchain"), py::return_value_policy::reference_internal)
    .def("__len__", [](const Structure& st) { return st.models.size(); })
    .def("__iter__", [](const Structure& st) {
        return py::make_iterator(st.models);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Structure, Model>, py::arg("index"),
         py::return_value_policy::reference_internal)
    .def("__getitem__", [](Structure& st, const std::string& name) -> Model& {
        return *impl::find_iter(st.models, name);
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__delitem__", &remove_child<Structure>, py::arg("index"))
    .def("__delitem__", &remove_children<Structure>)
    .def("__delitem__", &Structure::remove_model, py::arg("name"))
    .def("find_connection", &Structure::find_connection,
         py::arg("partner1"), py::arg("partner2"),
         py::return_value_policy::reference_internal)
    .def("find_or_add_model", &Structure::find_or_add_model,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("add_model", add_child<Structure, Model>,
         py::arg("model"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("renumber_models", &Structure::renumber_models)
    .def("merge_chain_parts", &Structure::merge_chain_parts,
         py::arg("min_sep")=0)
    .def("setup_cell_images", &Structure::setup_cell_images)
    .def("make_pdb_headers", &make_pdb_headers)
    .def("write_pdb", [](const Structure& st, const std::string& path,
                         bool seqres_records, bool ssbond_records,
                         bool link_records, bool cispep_records,
                         bool ter_records, bool numbered_ter,
                         bool use_linkr) {
       PdbWriteOptions options;
       options.seqres_records = seqres_records;
       options.ssbond_records = ssbond_records;
       options.link_records = link_records;
       options.cispep_records = cispep_records;
       options.ter_records = ter_records;
       options.numbered_ter = numbered_ter;
       options.use_linkr = use_linkr;
       Ofstream f(path);
       write_pdb(st, f.ref(), options);
    }, py::arg("path"),
       py::arg("seqres_records")=true, py::arg("ssbond_records")=true,
       py::arg("link_records")=true, py::arg("cispep_records")=true,
       py::arg("ter_records")=true, py::arg("numbered_ter")=true,
       py::arg("use_linkr")=false)
    .def("write_minimal_pdb",
         [](const Structure& st, const std::string& path) {
       Ofstream f(path);
       write_minimal_pdb(st, f.ref());
    }, py::arg("path"))
    .def("make_minimal_pdb", [](const Structure& st) -> std::string {
       std::ostringstream os;
       write_minimal_pdb(st, os);
       return os.str();
    })
    .def("make_mmcif_document", &make_mmcif_document)
    .def("make_mmcif_headers", &make_mmcif_headers)
    .def("add_entity_types", (void (*)(Structure&, bool)) &add_entity_types,
         py::arg("overwrite")=false)
    .def("assign_subchains", (void (*)(Structure&, bool)) &assign_subchains,
         py::arg("overwrite")=false)
    .def("ensure_entities", &ensure_entities)
    .def("deduplicate_entities", &deduplicate_entities)
    .def("setup_entities", &setup_entities)
    .def("setup_for_mmcif", &setup_for_mmcif)
    .def("remove_alternative_conformations",
         remove_alternative_conformations<Structure>)
    .def("remove_hydrogens", remove_hydrogens<Structure>)
    .def("remove_waters", remove_waters<Structure>)
    .def("remove_ligands_and_waters",
         (void (*)(Structure&)) &remove_ligands_and_waters)
    .def("remove_empty_chains", (void (*)(Structure&)) &remove_empty_chains)
    .def("shorten_chain_names", &shorten_chain_names)
    .def("clone", [](const Structure& self) { return new Structure(self); })
    .def("__repr__", [](const Structure& self) {
        return tostr("<gemmi.Structure ", self.name, " with ",
                     self.models.size(), " model(s)>");
    });

  pyAtomAddress
    .def(py::init<>())
    .def(py::init<const Chain&, const Residue&, const Atom&>())
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

  py::class_<CRA>(m, "CRA")
    .def_readonly("chain", &CRA::chain)
    .def_readonly("residue", &CRA::residue)
    .def_readonly("atom", &CRA::atom)
    .def("__str__", [](const CRA& self) { return atom_str(self); })
    .def("__repr__", [](const CRA& self) {
        return tostr("<gemmi.CRA ", atom_str(self), '>');
    });

  pyModel
    .def(py::init<std::string>())
    .def_readwrite("name", &Model::name)
    .def("__len__", [](const Model& self) { return self.chains.size(); })
    .def("__iter__", [](Model& self) {
        return py::make_iterator(self.chains);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Model, Chain>, py::arg("index"),
         py::return_value_policy::reference_internal)
    .def("__getitem__", [](Model& self, const std::string& name) -> Chain& {
        return *impl::find_iter(self.chains, name);
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("all", (CraProxy (Model::*)()) &Model::all, py::keep_alive<0, 1>())
    .def("get_subchain",
         (ResidueSpan (Model::*)(const std::string&)) &Model::get_subchain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("subchains", (std::vector<ResidueSpan> (Model::*)()) &Model::subchains,
         py::return_value_policy::reference_internal)
    .def("find_residue_group", &Model::find_residue_group,
         py::arg("chain"), py::arg("seqid"),
         py::keep_alive<0, 1>())
    .def("sole_residue", &Model::sole_residue,
         py::arg("chain"), py::arg("seqid"),
         py::return_value_policy::reference_internal)
    .def("get_all_residue_names", &Model::get_all_residue_names)
    .def("find_cra", (CRA (Model::*)(const AtomAddress&)) &Model::find_cra)
    .def("find_chain",
         (Chain* (Model::*)(const std::string&))&Model::find_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("find_last_chain", &Model::find_last_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("add_chain", add_child<Model, Chain>,
         py::arg("chain"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    // for compatibility with older gemmi version
    .def("add_chain", [](Model& self, const std::string& name) -> Chain& {
        self.chains.emplace_back(name);
        return self.chains.back();
     }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("remove_chain", &Model::remove_chain, py::arg("name"))
    .def("__delitem__", &Model::remove_chain, py::arg("name"))
    .def("__delitem__", remove_child<Model>, py::arg("index"))
    .def("__delitem__", remove_children<Model>)
    .def("count_atom_sites", &count_atom_sites<Model>)
    .def("count_occupancies", &count_occupancies<Model>)
    .def("calculate_mass", &calculate_mass<Model>)
    .def("calculate_center_of_mass", [](const Model& self) {
        return calculate_center_of_mass(self).get();
    })
    .def("clone", [](const Model& self) { return new Model(self); })
    .def("__repr__", [](const Model& self) {
        return tostr("<gemmi.Model ", self.name, " with ",
                     self.chains.size(), " chain(s)>");
    });

  py::class_<UniqProxy<Residue>>(m, "FirstConformerRes")
    .def("__iter__", [](UniqProxy<Residue>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  py::class_<UniqProxy<Residue, ResidueSpan>>(m, "FirstConformerResSpan")
    .def("__iter__", [](UniqProxy<Residue, ResidueSpan>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  py::class_<UniqProxy<Atom>>(m, "FirstConformerAtoms")
    .def("__iter__", [](UniqProxy<Atom>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  py::class_<ResidueSpan::GroupingProxy>(m, "ResidueSpanGroups")
    .def("__iter__", [](ResidueSpan::GroupingProxy& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  pyCraGenerator
    .def("__iter__", [](CraProxy& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>());

  pyChain
    .def(py::init<std::string>())
    .def_readwrite("name", &Chain::name)
    .def("__len__", [](const Chain& ch) { return ch.residues.size(); })
    .def("__iter__", [](const Chain& ch) {
        return py::make_iterator(ch.residues);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Chain& ch, const std::string& seqid) {
        return ch.find_residue_group(SeqId(seqid));
    }, py::arg("pdb_seqid"), py::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Chain, Residue>, py::arg("index"),
         py::return_value_policy::reference_internal)
    .def("__getitem__", [](Chain &ch, py::slice slice) -> py::list {
        ssize_t start, stop, step, slength;
        if (!slice.compute(ch.residues.size(), &start, &stop, &step, &slength))
          throw py::error_already_set();
        py::list l;
        for (ssize_t i = 0; i < slength; ++i)
          l.append(py::cast(&ch.residues[start + i * step]));
        return l;
    }, py::return_value_policy::reference_internal)
    .def("__delitem__", remove_child<Chain>, py::arg("index"))
    .def("__delitem__", remove_children<Chain>)
    .def("add_residue", add_child<Chain, Residue>,
         py::arg("residue"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("subchains", (std::vector<ResidueSpan> (Chain::*)()) &Chain::subchains)
    .def("whole", (ResidueSpan (Chain::*)()) &Chain::whole)
    .def("get_polymer", (ResidueSpan (Chain::*)()) &Chain::get_polymer)
    .def("get_ligands", (ResidueSpan (Chain::*)()) &Chain::get_ligands)
    .def("get_waters", (ResidueSpan (Chain::*)()) &Chain::get_waters)
    .def("get_subchain",
         (ResidueSpan (Chain::*)(const std::string&)) &Chain::get_subchain)
    .def("has_subchains_assigned", &has_subchains_assigned)
    .def("previous_residue", &Chain::previous_residue,
         py::return_value_policy::reference_internal)
    .def("next_residue", &Chain::next_residue,
         py::return_value_policy::reference_internal)
    .def("append_residues", &Chain::append_residues,
         py::arg("new_residues"), py::arg("min_sep")=0)
    .def("count_atom_sites", &count_atom_sites<Chain>)
    .def("count_occupancies", &count_occupancies<Chain>)
    .def("calculate_mass", &calculate_mass<Chain>)
    .def("trim_to_alanine", (void (*)(Chain&)) &trim_to_alanine)
    .def("first_conformer",
         (UniqProxy<Residue> (Chain::*)()) &Chain::first_conformer,
         py::keep_alive<0, 1>())
    .def("clone", [](const Chain& self) { return new Chain(self); })
    .def("__repr__", [](const Chain& self) {
        return tostr("<gemmi.Chain ", self.name,
                     " with ", self.residues.size(), " res>");
    });

  pyResidueSpan
    .def("__len__", &ResidueSpan::size)
    .def("__iter__", [](ResidueSpan& g) { return py::make_iterator(g); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const ResidueSpan &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](ResidueSpan& g, int index) -> Residue& {
        return g[normalize_index(index, g)];
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](ResidueSpan& self, const std::string& seqid) {
        return self.find_residue_group(SeqId(seqid));
    }, py::arg("pdb_seqid"), py::keep_alive<0, 1>())
    .def("__delitem__", [](ResidueSpan& self, int index) {
        self.erase(self.begin() + normalize_index(index, self));
    }, py::arg("index"))
    .def("add_residue", add_item<ResidueSpan, Residue>,
         py::arg("residue"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("first_conformer",
         (UniqProxy<Residue, ResidueSpan> (ResidueSpan::*)())
                                                &ResidueSpan::first_conformer,
         py::keep_alive<0, 1>())
    .def("residue_groups", &ResidueSpan::residue_groups, py::keep_alive<0, 1>())
    .def("length", &ResidueSpan::length)
    .def("subchain_id", &ResidueSpan::subchain_id)
    .def("label_seq_id_to_auth", &ResidueSpan::label_seq_id_to_auth)
    .def("auth_seq_id_to_label", &ResidueSpan::auth_seq_id_to_label)
    .def("check_polymer_type", [](const ResidueSpan& span) {
        return check_polymer_type(span);
    })
    .def("make_one_letter_sequence", [](const ResidueSpan& span) {
        return make_one_letter_sequence(span);
    })
    .def("__repr__", [](const ResidueSpan& self) {
        int N = (int) self.size();
        std::string r = "<gemmi.ResidueSpan of " + std::to_string(N) + ": [";
        for (int i = 0; i < (N < 5 ? N : 3); ++i) {
          if (i != 0) r += ' ';
          r += self[i].str();
        }
        if (N >= 5)
          r += " ... " + self[N-1].str();
        return r + "]>";
    });

  pyResidueGroup
    // need to duplicate it so it is visible
    .def("__getitem__", [](ResidueGroup& g, int index) -> Residue& {
        return g[normalize_index(index, g)];
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", &ResidueGroup::by_resname,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("__delitem__", &ResidueGroup::remove_residue, py::arg("name"))
    .def("__repr__", [](const ResidueGroup& self) {
        return "<gemmi.ResidueGroup [" +
               join_str(self, ' ', [](const Residue& r) { return r.str(); }) +
               "]>";
    });

  py::class_<AtomGroup>(m, "AtomGroup")
    .def("__len__", &AtomGroup::size)
    .def("__iter__", [](AtomGroup& g) { return py::make_iterator(g); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const ResidueSpan &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](AtomGroup& g, int index) -> Atom& {
        return g[normalize_index(index, g)];
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", &AtomGroup::by_altloc,
         py::arg("altloc"), py::return_value_policy::reference_internal)
    .def("name", &AtomGroup::name)
    .def("__repr__", [](const AtomGroup& self) {
        return tostr("<gemmi.AtomGroup ", self.name(), ", sites: ",
                     self.size(), '>');
    });

  pySeqId
    .def(py::init<int, char>())
    .def(py::init<const std::string&>())
    .def_readwrite("num", &SeqId::num)
    .def_readwrite("icode", &SeqId::icode)
    .def("__str__", &SeqId::str)
    .def("__repr__", [](const SeqId& self) {
        return tostr("<gemmi.SeqId ", self.str(), '>');
    });

  pyResidueId
    .def(py::init<>())
    .def_readwrite("name", &ResidueId::name)
    .def_readwrite("seqid", &ResidueId::seqid)
    .def_readwrite("segment", &ResidueId::segment)
    .def("__str__", &ResidueId::str)
    .def("__repr__", [](const ResidueId& self) {
        return tostr("<gemmi.ResidueId ", self.str(), '>');
    });

  pyResidue
    .def(py::init<>())
    .def_readwrite("subchain", &Residue::subchain)
    .def_readwrite("entity_type", &Residue::entity_type)
    .def_readwrite("het_flag", &Residue::het_flag)
    .def_readwrite("flag", &Residue::flag)
    .def_readwrite("label_seq", &Residue::label_seq)
    .def("__len__", [](const Residue& res) { return res.atoms.size(); })
    .def("__contains__", [](const Residue& res, const std::string& name) {
        return res.find_atom(name, '*') != nullptr;
    })
    .def("__iter__", [](const Residue& res) {
        return py::make_iterator(res.atoms);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Residue, Atom>, py::arg("index"),
         py::return_value_policy::reference_internal)
    .def("__getitem__", &Residue::get,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("__delitem__", remove_child<Residue>, py::arg("index"))
    .def("__delitem__", remove_children<Residue>)
    .def("find_atom", [](Residue& self, const std::string& name, char altloc) {
           return self.find_atom(name, altloc);
         },
         //(Atom* (Residue::*)(const std::string&, char, El)) &Residue::find_atom,
         py::arg("name"), py::arg("altloc"), /*py::arg("el")=El::X,*/
         py::return_value_policy::reference_internal)
    .def("remove_atom",
         [](Residue& self, const std::string& name, char altloc/*, El el*/) {
           self.atoms.erase(self.find_atom_iter(name, altloc/*, el*/));
    }, py::arg("name"), py::arg("altloc") /*, py::arg("el")=El::X*/)
    .def("add_atom", add_child<Residue, Atom>,
         py::arg("atom"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("first_conformer",
         (UniqProxy<Atom> (Residue::*)()) &Residue::first_conformer,
         py::keep_alive<0, 1>())
    .def("sole_atom", &Residue::sole_atom)
    .def("is_water", &Residue::is_water)
    .def("trim_to_alanine", (bool (*)(Residue&)) &trim_to_alanine)
    .def("clone", [](const Residue& self) { return new Residue(self); })
    .def("__repr__", [](const Residue& self) {
        return tostr("<gemmi.Residue ", self.str(), " with ",
                     self.atoms.size(), " atoms>");
    });

  pyAtom
    .def(py::init<>())
    .def_readwrite("name", &Atom::name)
    .def_readwrite("altloc", &Atom::altloc)
    .def_readwrite("charge", &Atom::charge)
    .def_readwrite("element", &Atom::element)
    .def_readwrite("pos", &Atom::pos)
    .def_readwrite("occ", &Atom::occ)
    .def_readwrite("b_iso", &Atom::b_iso)
    .def_readwrite("aniso", &Atom::aniso)
    .def_readwrite("serial", &Atom::serial)
    .def_readwrite("flag", &Atom::flag)
    .def("is_hydrogen", &Atom::is_hydrogen)
    .def("has_altloc", &Atom::has_altloc)
    .def("b_eq", &Atom::b_eq)
    .def("clone", [](const Atom& self) { return new Atom(self); })
    .def("__repr__", [](const Atom& self) {
        std::string r = "<gemmi.Atom " + self.name;
        if (self.altloc) {
            r += '.';
            r += self.altloc;
        }
        using namespace std;  // VS2015/17 doesn't like std::snprintf
        char buf[128];
        snprintf(buf, 128, " at (%.1f, %.1f, %.1f)>",
                 self.pos.x, self.pos.y, self.pos.z);
        return r + buf;
    });

  m.def("calculate_b_est", &calculate_b_est);
  m.def("calculate_angle", &calculate_angle,
        "Input: three points. Output: angle in radians.");
  m.def("calculate_dihedral", &calculate_dihedral,
        "Input: four points. Output: dihedral angle in radians.");
  m.def("calculate_phi_psi", &calculate_phi_psi,
        py::arg("prev_residue"), py::arg("residue"), py::arg("next_residue"));
  m.def("calculate_omega", &calculate_omega,
        py::arg("residue"), py::arg("next_residue"));
  m.def("calculate_sequence_weight", &calculate_sequence_weight,
        py::arg("sequence"), py::arg("unknown")=0.);
}
