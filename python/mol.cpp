// Copyright 2017 Global Phasing Ltd.

#include "gemmi/model.hpp"
#include "gemmi/calculate.hpp"  // for calculate_mass, count_atom_sites
#include "gemmi/modify.hpp"     // for remove_alternative_conformations
#include "gemmi/polyheur.hpp"   // for one_letter_code, trim_to_alanine
#include "gemmi/assembly.hpp"   // for expand_ncs, HowToNameCopiedChain
#include "gemmi/select.hpp"     // for Selection
#include "tostr.hpp"

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "meta.h"

namespace py = pybind11;
using namespace gemmi;

using info_map_type = std::map<std::string, std::string>;
PYBIND11_MAKE_OPAQUE(info_map_type)

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

template<typename P, typename C> C& get_child(P& parent, int index) {
  auto& children = parent.children();
  return children[normalize_index(index, children)];
}

template<typename P, typename C> void set_child(P& parent, int index, C& child) {
  auto& children = parent.children();
  children[normalize_index(index, children)] = child;
}

template<typename P> void remove_child(P& parent, int index) {
  auto& children = parent.children();
  children.erase(children.begin() + normalize_index(index, children));
}

template<typename P> void remove_children(P& parent, py::slice slice) {
  delitem_slice(parent.children(), slice);
}

void add_mol(py::module& m) {

  // "Forward declaration" of python classes to avoid
  // C++ signatures in docstrings
  py::class_<Atom> pyAtom(m, "Atom");
  py::class_<Residue, ResidueId> pyResidue(m, "Residue");
  py::class_<Chain> pyChain(m, "Chain");
  py::class_<Model> pyModel(m, "Model");
  py::class_<ResidueSpan> pyResidueSpan(m, "ResidueSpan");
  py::class_<ResidueGroup, ResidueSpan> pyResidueGroup(m, "ResidueGroup");
  py::class_<Selection> pySelection(m, "Selection");

  py::enum_<HowToNameCopiedChain>(m, "HowToNameCopiedChain")
    .value("Short", HowToNameCopiedChain::Short)
    .value("AddNumber", HowToNameCopiedChain::AddNumber)
    .value("Dup", HowToNameCopiedChain::Dup)
    ;

  m.def("one_letter_code",
        (std::string (*)(const std::vector<std::string>&)) &one_letter_code);
  m.def("one_letter_code",
        [](const ResidueSpan& span) { return one_letter_code(span); });
  m.def("make_address", &make_address);

  py::enum_<CoorFormat>(m, "CoorFormat")
    .value("Unknown", CoorFormat::Unknown)
    .value("Detect", CoorFormat::Detect)
    .value("Pdb", CoorFormat::Pdb)
    .value("Mmcif", CoorFormat::Mmcif)
    .value("Mmjson", CoorFormat::Mmjson)
    .value("ChemComp", CoorFormat::ChemComp);

  py::bind_map<info_map_type>(m, "InfoMap");

  py::class_<CRA>(m, "CRA")
    .def_readonly("chain", &CRA::chain)
    .def_readonly("residue", &CRA::residue)
    .def_readonly("atom", &CRA::atom)
    .def("atom_matches", [](const CRA& self, const AtomAddress& addr) {
        return atom_matches(self, addr);
    })
    .def("__str__", [](const CRA& self) { return atom_str(self); })
    .def("__repr__", [](const CRA& self) {
        return tostr("<gemmi.CRA ", atom_str(self), '>');
    });

  py::class_<CraProxy>(m, "CraGenerator")
    .def("__iter__", [](CraProxy& self) { return py::make_iterator(self); },
         py::keep_alive<0, 1>());

  py::class_<Structure> structure(m, "Structure");
  structure
    .def(py::init<>())
    .def_readwrite("name", &Structure::name)
    .def_readwrite("cell", &Structure::cell)
    .def_readwrite("spacegroup_hm", &Structure::spacegroup_hm)
    .def_readwrite("ncs", &Structure::ncs)
    .def_readwrite("resolution", &Structure::resolution)
    .def_readwrite("input_format", &Structure::input_format)
    .def_readwrite("entities", &Structure::entities)
    .def_readwrite("connections", &Structure::connections)
    .def_readwrite("helices", &Structure::helices)
    .def_readwrite("sheets", &Structure::sheets)
    .def_readwrite("assemblies", &Structure::assemblies)
    .def_readwrite("meta", &Structure::meta)
    .def_readwrite("has_d_fraction", &Structure::has_d_fraction)
    .def_readwrite("info", &Structure::info)
    .def_readwrite("raw_remarks", &Structure::raw_remarks)
    .def("find_spacegroup", &Structure::find_spacegroup)
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
    .def("__setitem__", &set_child<Structure, Model>)
    .def("find_connection_by_cra", [](Structure& st, CRA cra1, CRA cra2) {
        return st.find_connection_by_cra(cra1, cra2);
    }, py::return_value_policy::reference_internal)
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
    .def("remove_empty_chains", &Structure::remove_empty_chains)
    .def("setup_cell_images", &Structure::setup_cell_images)
    .def("add_entity_types", (void (*)(Structure&, bool)) &add_entity_types,
         py::arg("overwrite")=false)
    .def("assign_subchains", &assign_subchains,
         py::arg("force")=false, py::arg("fail_if_unknown")=true)
    .def("ensure_entities", &ensure_entities)
    .def("deduplicate_entities", &deduplicate_entities)
    .def("setup_entities", &setup_entities)
    .def("assign_cis_flags", assign_cis_flags<Structure>)
    .def("remove_alternative_conformations",
         remove_alternative_conformations<Structure>)
    .def("remove_hydrogens", remove_hydrogens<Structure>)
    .def("remove_waters", remove_waters<Structure>)
    .def("remove_ligands_and_waters", remove_ligands_and_waters<Structure>)
    .def("expand_hd_mixture", &expand_hd_mixture)
    .def("collapse_hd_mixture", &collapse_hd_mixture)
    .def("assign_serial_numbers", (void (*)(Structure&)) &assign_serial_numbers)
    .def("shorten_chain_names", &shorten_chain_names)
    .def("expand_ncs", &expand_ncs, py::arg("how"))
    .def("transform_to_assembly",
         [](Structure& st, const std::string& assembly_name, HowToNameCopiedChain how) {
        return transform_to_assembly(st, assembly_name, how, nullptr);
    }, py::arg("assembly_name"), py::arg("how"))
    .def("calculate_box", &calculate_box, py::arg("margin")=0.)
    .def("calculate_fractional_box", &calculate_fractional_box, py::arg("margin")=0.)
    .def("clone", [](const Structure& self) { return new Structure(self); })
    .def("__repr__", [](const Structure& self) {
        return tostr("<gemmi.Structure ", self.name, " with ",
                     self.models.size(), " model(s)>");
    });
    add_assign_label_seq_id(structure);
    add_write(m, structure);

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
    .def("find_cra", (CRA (Model::*)(const AtomAddress&, bool)) &Model::find_cra,
         py::arg(), py::arg("ignore_segment")=false)
    .def("find_chain",
         (Chain* (Model::*)(const std::string&))&Model::find_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("find_last_chain", &Model::find_last_chain,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("add_chain",
         [](Model& self, const Chain& chain, int pos, bool unique_name) -> Chain& {
           Chain& ref = add_child<Model, Chain>(self, chain, pos);
           if (unique_name)
             ensure_unique_chain_name(self, ref);
           return ref;
    }, py::arg("chain"), py::arg("pos")=-1, py::arg("unique_name")=false,
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
    .def("remove_alternative_conformations", remove_alternative_conformations<Model>)
    .def("remove_hydrogens", remove_hydrogens<Model>)
    .def("remove_waters", remove_waters<Model>)
    .def("remove_ligands_and_waters", remove_ligands_and_waters<Model>)
    .def("has_hydrogen", &has_hydrogen<Model>)
    .def("count_atom_sites", &count_atom_sites<Model>, py::arg("sel")=nullptr)
    // deprecated, use has_hydrogen() or count_atom_sites(Selection('[H,D]'))
    .def("count_hydrogen_sites", &count_hydrogen_sites<Model>)
    .def("count_occupancies", &count_occupancies<Model>, py::arg("sel")=nullptr)
    .def("calculate_mass", &calculate_mass<Model>)
    .def("calculate_center_of_mass", [](const Model& self) {
        return calculate_center_of_mass(self).get();
    })
    .def("transform_pos_and_adp", transform_pos_and_adp<Model>, py::arg("tr"))
    .def("split_chains_by_segments", &split_chains_by_segments)
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
        return getitem_slice(ch.residues, slice);
    }, py::return_value_policy::reference_internal)
    .def("__delitem__", remove_child<Chain>, py::arg("index"))
    .def("__delitem__", remove_children<Chain>)
    .def("add_residue", add_child<Chain, Residue>,
         py::arg("residue"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("subchains", (std::vector<ResidueSpan> (Chain::*)()) &Chain::subchains)
    .def("whole", (ResidueSpan (Chain::*)()) &Chain::whole,
         py::keep_alive<0, 1>())
    .def("get_polymer", (ResidueSpan (Chain::*)()) &Chain::get_polymer,
         py::keep_alive<0, 1>())
    .def("get_ligands", (ResidueSpan (Chain::*)()) &Chain::get_ligands,
         py::keep_alive<0, 1>())
    .def("get_waters", (ResidueSpan (Chain::*)()) &Chain::get_waters,
         py::keep_alive<0, 1>())
    .def("get_subchain", (ResidueSpan (Chain::*)(const std::string&)) &Chain::get_subchain,
         py::keep_alive<0, 1>())
    .def("has_entity_types_and_subchains", &has_entity_types_and_subchains)
    .def("previous_residue", &Chain::previous_residue,
         py::return_value_policy::reference_internal)
    .def("next_residue", &Chain::next_residue,
         py::return_value_policy::reference_internal)
    .def("append_residues", &Chain::append_residues,
         py::arg("new_residues"), py::arg("min_sep")=0)
    .def("count_atom_sites", &count_atom_sites<Chain>, py::arg("sel")=nullptr)
    .def("count_occupancies", &count_occupancies<Chain>, py::arg("sel")=nullptr)
    .def("calculate_mass", &calculate_mass<Chain>)
    .def("calculate_center_of_mass", [](const Chain& self) {
          return calculate_center_of_mass(self).get();
     })
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
    .def("transform_pos_and_adp", transform_pos_and_adp<ResidueSpan>)
    .def("__repr__", [](const ResidueSpan& self) {
        int N = (int) self.size();
        std::string r = "<gemmi.ResidueSpan of " + std::to_string(N) + ": ";
        if (N > 0) {
          r += self.front().subchain;
          if (self.back().subchain != self.front().subchain)
            r += " - " + self.back().subchain;
          r += ' ';
        }
        r += '[';
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
    .def("__bool__", [](const AtomGroup &g) -> bool { return !g.empty(); })
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

  pyResidue
    .def(py::init<>())
    .def_readwrite("subchain", &Residue::subchain)
    .def_readonly("entity_id", &Residue::entity_id)
    .def_readwrite("label_seq", &Residue::label_seq)
    .def_readwrite("entity_type", &Residue::entity_type)
    .def_readwrite("het_flag", &Residue::het_flag)
    .def_readwrite("flag", &Residue::flag)
    .def_property_readonly("sifts_unp", [](const Residue& self) {
        return py::make_tuple(self.sifts_unp.res,
                              self.sifts_unp.acc_index,
                              self.sifts_unp.num);
    })
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
    .def("find_atom", [](Residue& self, const std::string& name, char altloc, Element el) {
           return self.find_atom(name, altloc, el);
         }, py::arg("name"), py::arg("altloc"), py::arg("el")=Element(El::X),
            py::return_value_policy::reference_internal)
    .def("remove_atom",
         [](Residue& self, const std::string& name, char altloc, Element el) {
           self.atoms.erase(self.find_atom_iter(name, altloc, el));
    }, py::arg("name"), py::arg("altloc"), py::arg("el")=Element(El::X))
    .def("add_atom", add_child<Residue, Atom>,
         py::arg("atom"), py::arg("pos")=-1,
         py::return_value_policy::reference_internal)
    .def("first_conformer",
         (UniqProxy<Atom> (Residue::*)()) &Residue::first_conformer,
         py::keep_alive<0, 1>())
    .def("sole_atom", &Residue::sole_atom)
    .def("get_ca", &Residue::get_ca, py::return_value_policy::reference_internal)
    .def("get_p", &Residue::get_p, py::return_value_policy::reference_internal)
    .def("is_water", &Residue::is_water)
    .def("remove_hydrogens", &remove_hydrogens<Residue>)
    .def("trim_to_alanine", (bool (*)(Residue&)) &trim_to_alanine)
    .def("clone", [](const Residue& self) { return new Residue(self); })
    .def("__repr__", [](const Residue& self) {
        return tostr("<gemmi.Residue ", self.str(), " with ",
                     self.atoms.size(), " atoms>");
    });

  py::enum_<CalcFlag>(m, "CalcFlag")
    .value("NotSet", CalcFlag::NotSet)
    .value("Determined", CalcFlag::Determined)
    .value("Calculated", CalcFlag::Calculated)
    .value("Dummy", CalcFlag::Dummy);

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
    .def_readwrite("fraction", &Atom::fraction)
    .def_readwrite("calc_flag", &Atom::calc_flag)
    .def_readwrite("flag", &Atom::flag)
    .def_readwrite("tls_group_id", &Atom::tls_group_id)
    .def("is_hydrogen", &Atom::is_hydrogen)
    .def("has_altloc", &Atom::has_altloc)
    .def("b_eq", &Atom::b_eq)
    .def("padded_name", &Atom::padded_name)
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
  m.def("find_best_plane", &find_best_plane,
        py::arg("atoms"));
  m.def("get_distance_from_plane", &get_distance_from_plane,
        py::arg("pos"), py::arg("coeff"));
  m.def("calculate_omega", &calculate_omega,
        py::arg("residue"), py::arg("next_residue"));
  m.def("calculate_sequence_weight", &calculate_sequence_weight,
        py::arg("sequence"), py::arg("unknown")=0.);
  m.def("make_assembly", [](const Assembly& assembly, const Model& model,
                            HowToNameCopiedChain how) {
        return make_assembly(assembly, model, how, nullptr);
  });

  // select.hpp
  py::class_<FilterProxy<Selection, Model>> pySelectionModelsProxy(m, "SelectionModelsProxy");
  py::class_<FilterProxy<Selection, Chain>> pySelectionChainsProxy(m, "SelectionChainsProxy");
  py::class_<FilterProxy<Selection, Residue>> pySelectionResiduesProxy(m, "SelectionResiduesProxy");
  py::class_<FilterProxy<Selection, Atom>> pySelectionAtomsProxy(m, "SelectionAtomsProxy");

  pySelection
    .def(py::init<>())
    .def(py::init<const std::string&>())
    .def("models", &Selection::models)
    .def("chains", &Selection::chains)
    .def("residues", &Selection::residues)
    .def("atoms", &Selection::atoms)
    .def("first_in_model", &Selection::first_in_model,
         py::keep_alive<1, 2>())
    .def("first", &Selection::first, py::return_value_policy::reference,
         py::keep_alive<1, 2>())
    .def("str", &Selection::str)
    .def("set_residue_flags", &Selection::set_residue_flags)
    .def("set_atom_flags", &Selection::set_atom_flags)
    .def("copy_model_selection", &Selection::copy_selection<Model>)
    .def("copy_structure_selection", &Selection::copy_selection<Structure>)
    .def("remove_selected", &Selection::remove_selected<Structure>)
    .def("remove_selected", &Selection::remove_selected<Model>)
    .def("remove_not_selected", &Selection::remove_not_selected<Structure>)
    .def("remove_not_selected", &Selection::remove_not_selected<Model>)
    .def("__repr__", [](const Selection& self) {
        return "<gemmi.Selection CID: " + self.str() + ">";
    });

  pySelectionModelsProxy
    .def("__iter__", [](FilterProxy<Selection, Model>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  pySelectionChainsProxy
    .def("__iter__", [](FilterProxy<Selection, Chain>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  pySelectionResiduesProxy
    .def("__iter__", [](FilterProxy<Selection, Residue>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

  pySelectionAtomsProxy
    .def("__iter__", [](FilterProxy<Selection, Atom>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
}
