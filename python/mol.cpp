// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/entstr.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/model.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/polyheur.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/to_mmcif.hpp"
#include "gemmi/to_mmcif.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#ifdef  __INTEL_COMPILER
// warning #597: "X<T>::operator X<T>() const" will not be called for implicit
// or explicit conversions. That warning is triggered when template
// UniqProxy is expanded with const Value.
# pragma warning disable 597
#endif

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<NcsOp>)
PYBIND11_MAKE_OPAQUE(std::vector<Entity>)
using info_map_type = std::map<std::string, std::string>;
PYBIND11_MAKE_OPAQUE(info_map_type)

namespace pybind11 { namespace detail {
  template<> struct type_caster<SeqId::OptionalNum>
    : optional_caster<SeqId::OptionalNum> {};
}} // namespace pybind11::detail


void add_mol(py::module& m) {
  py::class_<ResidueInfo>(m, "ResidueInfo")
    .def_readonly("one_letter_code", &ResidueInfo::one_letter_code)
    .def_readonly("hydrogen_count", &ResidueInfo::hydrogen_count)
    .def("found", &ResidueInfo::found)
    .def("is_standard", &ResidueInfo::is_standard)
    .def("is_water", &ResidueInfo::is_water)
    .def("is_nucleic_acid", &ResidueInfo::is_nucleic_acid)
    .def("is_amino_acid", &ResidueInfo::is_amino_acid);

  m.def("find_tabulated_residue", &find_tabulated_residue, py::arg("name"),
        "Find chemical component information in the internal table.");

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

  py::class_<PolySeqItem>(m, "PolySeqItem")
    .def(py::init<int, std::string>())
    .def_readwrite("num", &PolySeqItem::num)
    .def_readwrite("mon", &PolySeqItem::mon)
    .def("__repr__", [](const PolySeqItem& self) {
        std::string r = "<gemmi.PolySeqItem " + self.mon;
        if (self.num != -1)
          r += " (" + std::to_string(self.num) + ")";
        return r + ">";
    });

  py::class_<Entity>(m, "Entity")
    .def(py::init<std::string>())
    .def_readwrite("name", &Entity::name)
    .def_readwrite("subchains", &Entity::subchains)
    .def_readwrite("entity_type", &Entity::entity_type)
    .def_readwrite("polymer_type", &Entity::polymer_type)
    .def_readwrite("poly_seq", &Entity::poly_seq)
    .def("__repr__", [](const Entity& self) {
        std::string r = "<gemmi.Entity '" + self.name + "' ";
        r += entity_type_to_string(self.entity_type);
        if (self.polymer_type != PolymerType::Unknown)
          r += " " + polymer_type_to_qstring(self.polymer_type);
        using namespace std;  // VS2015/17 doesn't like std::snprintf
        char buf[64];
        snprintf(buf, 64, " object at %p>", (void*)&self);
        return r + buf;
    });

  py::class_<NcsOp>(m, "NcsOp")
    .def_readwrite("id", &NcsOp::id)
    .def_readwrite("given", &NcsOp::given)
    .def_readonly("tr", &NcsOp::tr)
    .def("apply", &NcsOp::apply)
    .def("__repr__", [](const NcsOp& self) {
        return "<gemmi.NcsOp " + self.id +
               " |shift|=" + std::to_string(self.tr.vec.length()) +
               (self.given ? " (" : " (not ") + "given)>";
    });

  py::bind_vector<std::vector<NcsOp>>(m, "VectorNcsOp");
  py::bind_vector<std::vector<Entity>>(m, "VectorEntity");
  py::bind_map<info_map_type>(m, "InfoMap");

  py::class_<Structure>(m, "Structure")
    .def(py::init<>())
    .def_readwrite("name", &Structure::name)
    .def_readwrite("cell", &Structure::cell)
    .def_readwrite("spacegroup_hm", &Structure::spacegroup_hm)
    .def_readwrite("ncs", &Structure::ncs)
    .def_readwrite("resolution", &Structure::resolution)
    .def_readwrite("entities", &Structure::entities)
    .def_readwrite("info", &Structure::info)
    .def("get_entity",
         (Entity* (Structure::*)(const std::string&)) &Structure::get_entity,
         py::arg("subchain"), py::return_value_policy::reference_internal)
    .def("get_entity_of",
         (Entity* (Structure::*)(const ResidueSpan&)) &Structure::get_entity_of,
         py::arg("subchain"), py::return_value_policy::reference_internal)
    .def("__len__", [](const Structure& st) { return st.models.size(); })
    .def("__iter__", [](const Structure& st) {
        return py::make_iterator(st.models);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Structure& st, int index) -> Model& {
        return st.models.at(index >= 0 ? index : index + st.models.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](Structure& st, const std::string& name) -> Model& {
        return *impl::find_iter(st.models, name);
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__delitem__", &Structure::remove_model, py::arg("name"))
    .def("find_or_add_model", &Structure::find_or_add_model,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("merge_chain_parts", &Structure::merge_chain_parts,
         py::arg("min_sep")=0)
    .def("make_pdb_headers", &make_pdb_headers)
    .def("write_pdb", [](const Structure& st, const std::string& path,
                         bool seqres_records, bool ssbond_records,
                         bool link_records, bool cispep_records,
                         bool ter_records, bool numbered_ter) {
       PdbWriteOptions options;
       options.seqres_records = seqres_records;
       options.ssbond_records = ssbond_records;
       options.link_records = link_records;
       options.cispep_records = cispep_records;
       options.ter_records = ter_records;
       options.numbered_ter = numbered_ter;
       std::ofstream f(path.c_str());
       write_pdb(st, f, options);
    }, py::arg("path"),
       py::arg("seqres_records")=true, py::arg("ssbond_records")=true,
       py::arg("link_records")=true, py::arg("cispep_records")=true,
       py::arg("ter_records")=true, py::arg("numbered_ter")=true)
    .def("write_minimal_pdb",
         [](const Structure& st, const std::string& path) {
       std::ofstream f(path.c_str());
       write_minimal_pdb(st, f);
    }, py::arg("path"))
    .def("make_minimal_pdb", [](const Structure& st) -> std::string {
       std::ostringstream os;
       write_minimal_pdb(st, os);
       return os.str();
    })
    .def("make_mmcif_document", &make_mmcif_document)
    .def("add_entity_types", (void (*)(Structure&, bool)) &add_entity_types,
         py::arg("overwrite")=false)
    .def("assign_subchains", (void (*)(Structure&, bool)) &assign_subchains,
         py::arg("overwrite")=false)
    .def("ensure_entities", &ensure_entities)
    .def("deduplicate_entities", &deduplicate_entities)
    .def("setup_entities", &setup_entities)
    .def("remove_hydrogens", remove_hydrogens<Structure>)
    .def("remove_waters", remove_waters<Structure>)
    .def("remove_ligands_and_waters",
         (void (*)(Structure&)) &remove_ligands_and_waters)
    .def("remove_empty_chains", (void (*)(Structure&)) &remove_empty_chains)
    .def("__repr__", [](const Structure& self) {
        return "<gemmi.Structure " + self.name + " with " +
               std::to_string(self.models.size()) + " model(s)>";
    });

  py::class_<CRA>(m, "CRA")
    .def_readonly("chain", &CRA::chain)
    .def_readonly("residue", &CRA::residue)
    .def_readonly("atom", &CRA::atom);

  py::class_<Model>(m, "Model")
    .def(py::init<std::string>())
    .def_readwrite("name", &Model::name)
    .def("__len__", [](const Model& self) { return self.chains.size(); })
    .def("__iter__", [](const Model& self) {
        return py::make_iterator(self.chains);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Model& self, const std::string& name) -> Chain& {
        return *impl::find_iter(self.chains, name);
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](Model& self, int index) -> Chain& {
        return self.chains.at(index >= 0 ? index : index + self.chains.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("get_subchain",
         (ResidueSpan (Model::*)(const std::string&)) &Model::get_subchain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("subchains", &Model::subchains,
         py::return_value_policy::reference_internal)
    .def("find_residue_group", &Model::find_residue_group,
         py::arg("chain"), py::arg("seqid"),
         py::return_value_policy::reference_internal)
    .def("sole_residue", &Model::sole_residue,
         py::arg("chain"), py::arg("seqid"),
         py::return_value_policy::reference_internal)
    .def("find_chain", &Model::find_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("find_last_chain", &Model::find_last_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("add_chain", [](Model& self, const std::string& name) -> Chain& {
        self.chains.emplace_back(name);
        return self.chains.back();
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("remove_chain", &Model::remove_chain, py::arg("name"))
    .def("__delitem__", &Model::remove_chain, py::arg("name"))
    .def("count_atom_sites", &count_atom_sites<Model>)
    .def("count_occupancies", &count_occupancies<Model>)
    .def("__repr__", [](const Model& self) {
        return "<gemmi.Model " + self.name + " with " +
               std::to_string(self.chains.size()) + " chain(s)>";
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

  py::class_<Chain>(m, "Chain")
    .def(py::init<std::string>())
    .def_readwrite("name", &Chain::name)
    .def("__len__", [](const Chain& ch) { return ch.residues.size(); })
    .def("__iter__", [](const Chain& ch) {
        return py::make_iterator(ch.residues);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Chain& ch, const std::string& seqid) {
        return ch.find_residue_group(SeqId(seqid));
    }, py::arg("pdb_seqid"), py::keep_alive<0, 1>())
    .def("__getitem__", [](Chain& ch, int index) -> Residue& {
        return ch.residues.at(index >= 0 ? index : index + ch.residues.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__delitem__", [](Chain& ch, int index) {
        if (index < 0)
          index += ch.residues.size();
        if ((size_t) index >= ch.residues.size())
          throw py::index_error();
        ch.residues.erase(ch.residues.begin() + index);
    }, py::arg("index"))
    .def("subchains", &Chain::subchains)
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
    .def("trim_to_alanine", (void (*)(Chain&)) &trim_to_alanine)
    .def("first_conformer",
         (UniqProxy<Residue> (Chain::*)()) &Chain::first_conformer)
    .def("__repr__", [](const Chain& self) {
        return "<gemmi.Chain " + self.name +
               " with " + std::to_string(self.residues.size()) + " res>";
    });

  py::class_<ResidueSpan> residue_span(m, "ResidueSpan");
  residue_span
    .def("__len__", &ResidueSpan::size)
    .def("__iter__", [](ResidueSpan& g) { return py::make_iterator(g); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const ResidueSpan &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](ResidueSpan& g, int index) -> Residue& {
        return g.at(index >= 0 ? index : index + g.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("first_conformer", (UniqProxy<Residue, ResidueSpan> (ResidueSpan::*)())
                            &ResidueSpan::first_conformer)
    .def("length", &ResidueSpan::length)
    .def("subchain_id", &ResidueSpan::subchain_id)
    .def("check_polymer_type", &check_polymer_type)
    .def("make_one_letter_sequence", &make_one_letter_sequence)
    .def("__repr__", [](const ResidueSpan& self) {
        int N = self.size();
        std::string r = "<gemmi.ResidueSpan of " + std::to_string(N) + ": [";
        for (int i = 0; i < (N < 5 ? N : 3); ++i) {
          if (i != 0) r += ' ';
          r += self[i].str();
        }
        if (N >= 5)
          r += " ... " + self[N-1].str();
        return r + "]>";
    });

  py::class_<ResidueGroup>(m, "ResidueGroup", residue_span)
    // need to duplicate it so it is visible
    .def("__getitem__", [](ResidueGroup& g, int index) -> Residue& {
        return g.at(index >= 0 ? index : index + g.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", &ResidueGroup::by_resname,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("__delitem__", &ResidueGroup::remove_residue, py::arg("name"))
    .def("__repr__", [](const ResidueGroup& self) {
        return "<gemmi.ResidueGroup [" +
               join_str(self, ' ', [](const Residue& r) { return r.str(); }) +
               "]>";
    });

  py::class_<SeqId>(m, "SeqId")
    .def(py::init<int, char>())
    .def(py::init<const std::string&>())
    .def_readwrite("num", &SeqId::num)
    .def_readwrite("icode", &SeqId::icode)
    .def("__str__", &SeqId::str)
    .def("__repr__", [](const SeqId& self) {
        return "<gemmi.SeqId " + self.str() + ">";
    });

  py::class_<Residue>(m, "Residue")
    .def(py::init<>())
    .def_readwrite("name", &Residue::name)
    .def_readwrite("seqid", &Residue::seqid)
    .def_readwrite("segment", &Residue::segment)
    .def_readwrite("subchain", &Residue::subchain)
    .def_readwrite("entity_type", &Residue::entity_type)
    .def_readwrite("het_flag", &Residue::het_flag)
    .def_readwrite("label_seq", &Residue::label_seq)
    .def("__len__", [](const Residue& res) { return res.atoms.size(); })
    .def("__contains__", [](const Residue& res, const std::string& name) {
        return res.find_atom(name) != nullptr;
    })
    .def("__iter__", [](const Residue& res) {
        return py::make_iterator(res.atoms);
    }, py::keep_alive<0, 1>())
    // TODO: should it return a list of items with this name
    //       or should altloc be specified as part of the name ('CA.A')
    .def("__getitem__", [](Residue& res, const std::string& name) -> Atom& {
        return *impl::find_iter(res.atoms, name);
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__getitem__", [](Residue& self, int index) -> Atom& {
        return self.atoms.at(index >= 0 ? index : index + self.atoms.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__delitem__", &Residue::remove_atom, py::arg("name"))
    .def("first_conformer",
         (UniqProxy<Atom> (Residue::*)()) &Residue::first_conformer)
    .def("is_water", &Residue::is_water)
    .def("trim_to_alanine", (bool (*)(Residue&)) &trim_to_alanine)
    .def("__repr__", [](const Residue& self) {
        return "<gemmi.Residue " + self.str() +
               " with " + std::to_string(self.atoms.size()) + " atoms>";
    });

  py::class_<Atom>(m, "Atom")
    .def(py::init<>())
    .def_readwrite("name", &Atom::name)
    .def_readwrite("altloc", &Atom::altloc)
    .def_readwrite("charge", &Atom::charge)
    .def_readwrite("element", &Atom::element)
    .def_readwrite("pos", &Atom::pos)
    .def_readwrite("occ", &Atom::occ)
    .def_readwrite("b_iso", &Atom::b_iso)
    .def_readwrite("serial", &Atom::serial)
    .def("is_hydrogen", &Atom::is_hydrogen)
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

  m.def("calculate_dihedral", &calculate_dihedral,
        "Input: four points. Output: dihedral angle in radians.");
}
