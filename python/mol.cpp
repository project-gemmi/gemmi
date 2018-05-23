// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/resinfo.hpp"
#include "gemmi/model.hpp"
#include "gemmi/calculate.hpp"
#include "gemmi/modify.hpp"
#include "gemmi/polyheur.hpp"  // for extract_sequence_info
#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/to_pdb.hpp"
#include "gemmi/to_mmcif.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<NcsOp>)
using entity_map_type = std::map<std::string, Entity>;
PYBIND11_MAKE_OPAQUE(entity_map_type)
using info_map_type = std::map<std::string, std::string>;
PYBIND11_MAKE_OPAQUE(info_map_type)

namespace pybind11 { namespace detail {
  template<> struct type_caster<ResidueId::OptionalNum>
    : optional_caster<ResidueId::OptionalNum> {};
}} // namespace pybind11::detail


void add_mol(py::module& m) {
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });

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

  py::class_<Entity>(m, "Entity")
    .def(py::init<>())
    .def_readwrite("entity_type", &Entity::entity_type)
    .def_readwrite("polymer_type", &Entity::polymer_type)
    .def("__repr__", [](const Entity& self) {
        std::string r = "<gemmi.Entity " + self.type_as_string();
        if (self.polymer_type != PolymerType::Unknown)
          r += " / " + self.polymer_type_as_string();
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
  py::bind_map<entity_map_type>(m, "EntityMap");
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
    .def("get_entity_of",
         (Entity* (Structure::*)(const Chain&)) &Structure::get_entity_of,
         py::arg("chain"), py::return_value_policy::reference_internal)
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
    .def("make_pdb_headers", &make_pdb_headers)
    .def("write_pdb", [](const Structure& st, const std::string& path) {
       std::ofstream f(path.c_str());
       write_pdb(st, f);
    }, py::arg("path"))
    .def("write_minimal_pdb",
         [](const Structure& st, const std::string& path, const char* chain) {
       std::ofstream f(path.c_str());
       write_minimal_pdb(st, f, chain);
    }, py::arg("path"), py::arg("chain")=nullptr)
    .def("make_mmcif_document", &make_mmcif_document)
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
    .def("residues", &Model::residues,
         py::arg("auth_chain"), py::arg("resnum"), py::arg("icode"),
         py::return_value_policy::reference_internal)
    .def("residue", &Model::residue,
         py::arg("auth_chain"), py::arg("resnum"), py::arg("icode"),
         py::return_value_policy::reference_internal)
    .def("find_or_add_chain", &Model::find_or_add_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
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

  py::class_<Chain>(m, "Chain")
    .def(py::init<std::string>())
    .def_readwrite("name", &Chain::name)
    .def_readwrite("auth_name", &Chain::auth_name)
    .def_readwrite("entity_id", &Chain::entity_id)
    .def("__len__", [](const Chain& ch) { return ch.residues.size(); })
    .def("__iter__", [](const Chain& ch) {
        return py::make_iterator(ch.residues);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", &Chain::find_by_seqid,
         py::arg("pdb_seqid"), py::keep_alive<0, 1>())
    .def("__getitem__", [](Chain& ch, int index) -> Residue& {
        return ch.residues.at(index >= 0 ? index : index + ch.residues.size());
    }, py::arg("index"), py::keep_alive<0, 1>())
    .def("append_residues", &Chain::append_residues)
    .def("count_atom_sites", &count_atom_sites<Chain>)
    .def("count_occupancies", &count_occupancies<Chain>)
    .def("trim_to_alanine", &trim_to_alanine)
    .def("extract_sequence_info", &extract_sequence_info)
    .def("first_conformer",
         (UniqProxy<Residue> (Chain::*)()) &Chain::first_conformer)
    .def("__repr__", [](const Chain& self) {
        return "<gemmi.Chain " + self.name + " (" + self.auth_name +
               ") with " + std::to_string(self.residues.size()) + " res>";
    });

  py::class_<ResidueGroup>(m, "ResidueGroup")
    .def("__len__", &ResidueGroup::size)
    .def("__iter__", [](ResidueGroup& g) { return py::make_iterator(g); },
         py::keep_alive<0, 1>())
    .def("__bool__", [](const ResidueGroup &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](ResidueGroup& g, int index) -> Residue& {
        return g.at(index >= 0 ? index : index + g.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("__getitem__", &ResidueGroup::by_resname,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("__repr__", [](const ResidueGroup& self) {
        return "<gemmi.ResidueGroup [" +
               join_str(self, ' ', [](const Residue& r) { return r.str(); }) +
               "]>";
    });

  py::class_<UniqProxy<Atom>>(m, "FirstConformerAtoms")
    .def("__iter__", [](UniqProxy<Atom>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  py::class_<Residue>(m, "Residue")
    .def(py::init<>())
    .def_readwrite("name", &Residue::name)
    .def_readwrite("label_seq", &Residue::label_seq)
    .def_readwrite("seq_num", &Residue::seq_num)
    .def("seq_id", &Residue::seq_id)
    .def_readwrite("icode", &Residue::icode)
    .def_readwrite("segment", &Residue::segment)
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
