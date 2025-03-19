// Copyright 2017 Global Phasing Ltd.

#include "gemmi/model.hpp"
#include "gemmi/calculate.hpp"  // for calculate_mass, count_atom_sites
#include "gemmi/modify.hpp"     // for remove_alternative_conformations
#include "gemmi/polyheur.hpp"   // for one_letter_code, trim_to_alanine
#include "gemmi/assembly.hpp"   // for expand_ncs, HowToNameCopiedChain
#include "gemmi/select.hpp"     // for Selection
#include "gemmi/sprintf.hpp"    // for snprintf_z

#include "common.h"
#include "serial.h"  // for getstate, setstate
#include "make_iterator.h"
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/array.h>  // for calculate_phi_psi, find_best_plane, ...
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/variant.h>
#include "meta.h"  // IWYU pragma: export

using namespace gemmi;

using info_map_type = std::map<std::string, std::string>;
NB_MAKE_OPAQUE(info_map_type)

namespace {

// cf. returns_references_to in nanobind docs
struct returns_references {
  static void precall(PyObject **, size_t, nb::detail::cleanup_list *) {}

  static void postcall(PyObject **args, size_t, nb::handle ret) {
    if (!nb::isinstance<nb::sequence>(ret))
      throw std::runtime_error("return value should be a sequence");
    for (nb::handle nurse : ret)
      nb::detail::keep_alive(nurse.ptr(), args[0]);
  }
};

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

template<typename P> void remove_children(P& parent, const nb::slice& slice) {
  delitem_slice(parent.children(), slice);
}

}  // anonymous namespace

void add_mol(nb::module_& m) {

  // "Forward declaration" of python classes to avoid
  // C++ signatures in docstrings
  nb::class_<Atom> pyAtom(m, "Atom");
  nb::class_<Residue, ResidueId> pyResidue(m, "Residue");
  nb::class_<Chain> pyChain(m, "Chain");
  nb::class_<Model> pyModel(m, "Model");
  nb::class_<ResidueSpan> pyResidueSpan(m, "ResidueSpan");
  nb::class_<ResidueGroup, ResidueSpan> pyResidueGroup(m, "ResidueGroup");
  nb::class_<Selection> pySelection(m, "Selection");

  nb::enum_<HowToNameCopiedChain>(m, "HowToNameCopiedChain")
    .value("Short", HowToNameCopiedChain::Short)
    .value("AddNumber", HowToNameCopiedChain::AddNumber)
    .value("Dup", HowToNameCopiedChain::Dup)
    ;

  m.def("make_address", &make_address);

  nb::enum_<CoorFormat>(m, "CoorFormat")
    .value("Unknown", CoorFormat::Unknown)
    .value("Detect", CoorFormat::Detect)
    .value("Pdb", CoorFormat::Pdb)
    .value("Mmcif", CoorFormat::Mmcif)
    .value("Mmjson", CoorFormat::Mmjson)
    .value("ChemComp", CoorFormat::ChemComp);

  nb::bind_map<info_map_type, rv_ri>(m, "InfoMap");

  nb::class_<CRA>(m, "CRA")
    .def_ro("chain", &CRA::chain)
    .def_ro("residue", &CRA::residue)
    .def_ro("atom", &CRA::atom)
    .def("atom_matches", [](const CRA& self, const AtomAddress& addr) {
        return atom_matches(self, addr);
    })
    .def("__str__", [](const CRA& self) { return atom_str(self); })
    .def("__repr__", [](const CRA& self) {
        return cat("<gemmi.CRA ", atom_str(self), '>');
    });

  nb::class_<CraProxy>(m, "CraGenerator")
    .def("__iter__", [](CraProxy& self) { return usual_iterator(self, self); },
         nb::keep_alive<0, 1>());

  nb::class_<Structure> structure(m, "Structure");
  structure
    .def(nb::init<>())
    .def_rw("name", &Structure::name)
    .def_rw("cell", &Structure::cell)
    .def_rw("spacegroup_hm", &Structure::spacegroup_hm)
    .def_rw("ncs", &Structure::ncs)
    .def_rw("resolution", &Structure::resolution)
    .def_rw("input_format", &Structure::input_format)
    .def_rw("entities", &Structure::entities)
    .def_rw("connections", &Structure::connections)
    .def_rw("cispeps", &Structure::cispeps)
    .def_rw("mod_residues", &Structure::mod_residues)
    .def_rw("helices", &Structure::helices)
    .def_rw("sheets", &Structure::sheets)
    .def_rw("assemblies", &Structure::assemblies)
    .def_rw("meta", &Structure::meta)
    .def_rw("has_d_fraction", &Structure::has_d_fraction)
    .def_rw("has_origx", &Structure::has_origx)
    .def_ro("origx", &Structure::origx)
    .def_rw("info", &Structure::info)
    .def_rw("raw_remarks", &Structure::raw_remarks)
    .def("find_spacegroup", &Structure::find_spacegroup, nb::rv_policy::reference)
    .def("get_entity",
         (Entity* (Structure::*)(const std::string&)) &Structure::get_entity,
         nb::arg("subchain"), nb::rv_policy::reference_internal)
    .def("get_entity_of", [](Structure& st, ResidueSpan& span) {
        return st.get_entity_of(span);
    }, nb::arg("subchain"), nb::rv_policy::reference_internal)
    .def("__len__", [](const Structure& st) { return st.models.size(); })
    .def("__iter__", [](Structure& st) { return usual_iterator(st, st.models); },
         nb::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Structure, Model>, nb::arg("index"),
         nb::rv_policy::reference_internal)
    .def("__delitem__", &remove_child<Structure>, nb::arg("index"))
    .def("__delitem__", &remove_children<Structure>)
    .def("__setitem__", &set_child<Structure, Model>)
    .def("find_connection_by_cra", [](Structure& st, CRA cra1, CRA cra2, bool ignore_segment) {
        return st.find_connection_by_cra(cra1, cra2, ignore_segment);
    }, nb::arg("cra1"), nb::arg("cra2"), nb::arg("ignore_segment")=false,
       nb::rv_policy::reference_internal)
    .def("find_connection", &Structure::find_connection,
         nb::arg("partner1"), nb::arg("partner2"),
         nb::rv_policy::reference_internal)
    .def("find_or_add_model", &Structure::find_or_add_model,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("add_model", add_child<Structure, Model>,
         nb::arg("model"), nb::arg("pos")=-1,
         nb::rv_policy::reference_internal)
    .def("renumber_models", &Structure::renumber_models)
    .def("merge_chain_parts", &Structure::merge_chain_parts,
         nb::arg("min_sep")=0)
    .def("remove_empty_chains", &Structure::remove_empty_chains)
    .def("setup_cell_images", &Structure::setup_cell_images)
    .def("add_entity_types", (void (*)(Structure&, bool)) &add_entity_types,
         nb::arg("overwrite")=false)
    .def("add_entity_ids", &add_entity_ids,
         nb::arg("overwrite")=false)
    .def("add_conect", &Structure::add_conect,
         nb::arg("serial1"), nb::arg("serial2"), nb::arg("order"))
    .def("clear_conect", [](Structure& self) { self.conect_map.clear(); })
    .def_ro("conect_map", &Structure::conect_map)
    .def("assign_subchains", &assign_subchains,
         nb::arg("force")=false, nb::arg("fail_if_unknown")=true)
    .def("ensure_entities", &ensure_entities)
    .def("deduplicate_entities", &deduplicate_entities)
    .def("setup_entities", &setup_entities)
    .def("assign_het_flags", assign_het_flags<Structure>, nb::arg("flag")='R')
    .def("remove_waters", remove_waters<Structure>)
    .def("remove_ligands_and_waters", remove_ligands_and_waters<Structure>)
    .def("shorten_ccd_codes", &shorten_ccd_codes)
    .def("restore_full_ccd_codes", &restore_full_ccd_codes)
    .def_rw("shortened_ccd_codes", &Structure::shortened_ccd_codes)

    // modify.hpp
    .def("remove_alternative_conformations",
         remove_alternative_conformations<Structure>)
    .def("remove_hydrogens", remove_hydrogens<Structure>)
    .def("store_deuterium_as_fraction", &store_deuterium_as_fraction)
    .def("standardize_crystal_frame", &standardize_crystal_frame)
    .def("assign_serial_numbers", (void (*)(Structure&, bool)) &assign_serial_numbers,
         nb::arg("numbered_ter")=false)
    .def("rename_chain", &rename_chain)
    .def("rename_residues", &rename_residues)
    // assembly.hpp
    .def("shorten_chain_names", &shorten_chain_names)
    .def("expand_ncs", &expand_ncs, nb::arg("how"), nb::arg("merge_dist")=0.2)
    .def("transform_to_assembly", &transform_to_assembly,
       nb::arg("assembly_name"), nb::arg("how"), nb::arg("logging")=nb::none(),
       nb::arg("keep_spacegroup")=false, nb::arg("merge_dist")=0.2)
    // calculate.hpp
    .def("calculate_box", &calculate_box, nb::arg("margin")=0.)
    .def("calculate_fractional_box", &calculate_fractional_box, nb::arg("margin")=0.)

    .def("clone", [](const Structure& self) { return new Structure(self); })
    .def("__getstate__", &getstate<Structure>)
    .def("__setstate__", &setstate<Structure>)
    .def("__repr__", [](const Structure& self) {
        return cat("<gemmi.Structure ", self.name, " with ",
                   self.models.size(), " model(s)>");
    });
    add_assign_label_seq_id(structure);
    add_write(m, structure);

  pyModel
    .def(nb::init<int>())
    .def_rw("num", &Model::num)
    .def("__len__", [](const Model& self) { return self.chains.size(); })
    .def("__iter__", [](Model& self) { return usual_iterator(self, self.chains); },
         nb::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Model, Chain>, nb::arg("index"),
         nb::rv_policy::reference_internal)
    .def("__getitem__", [](Model& self, const std::string& name) -> Chain& {
        return *impl::find_iter(self.chains, name);
    }, nb::arg("name"), nb::rv_policy::reference_internal)
    .def("all", (CraProxy (Model::*)()) &Model::all, nb::keep_alive<0, 1>())
    .def("get_subchain",
         (ResidueSpan (Model::*)(const std::string&)) &Model::get_subchain,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("subchains", (std::vector<ResidueSpan> (Model::*)()) &Model::subchains,
         nb::call_policy<returns_references>())
    .def("find_residue_group", &Model::find_residue_group,
         nb::arg("chain"), nb::arg("seqid"),
         nb::keep_alive<0, 1>())
    .def("sole_residue", &Model::sole_residue,
         nb::arg("chain"), nb::arg("seqid"),
         nb::rv_policy::reference_internal)
    .def("get_all_residue_names", &Model::get_all_residue_names)
    .def("find_cra", (CRA (Model::*)(const AtomAddress&, bool)) &Model::find_cra,
         nb::arg(), nb::arg("ignore_segment")=false)
    .def("find_chain",
         (Chain* (Model::*)(const std::string&))&Model::find_chain,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("find_last_chain", &Model::find_last_chain,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("add_chain",
         [](Model& self, const Chain& chain, int pos, bool unique_name) -> Chain& {
           Chain& ref = add_child<Model, Chain>(self, chain, pos);
           if (unique_name)
             ensure_unique_chain_name(self, ref);
           return ref;
    }, nb::arg("chain"), nb::arg("pos")=-1, nb::arg("unique_name")=false,
       nb::rv_policy::reference_internal)
    // for compatibility with older gemmi version
    .def("add_chain", [](Model& self, const std::string& name) -> Chain& {
        self.chains.emplace_back(name);
        return self.chains.back();
     }, nb::arg("name"), nb::rv_policy::reference_internal)
    .def("remove_chain", &Model::remove_chain, nb::arg("name"))
    .def("__delitem__", &Model::remove_chain, nb::arg("name"))
    .def("__delitem__", remove_child<Model>, nb::arg("index"))
    .def("__delitem__", remove_children<Model>)
    .def("remove_alternative_conformations", remove_alternative_conformations<Model>)
    .def("remove_hydrogens", remove_hydrogens<Model>)
    .def("remove_waters", remove_waters<Model>)
    .def("remove_ligands_and_waters", remove_ligands_and_waters<Model>)
    .def("has_hydrogen", &has_hydrogen<Model>)
    .def("count_atom_sites", &count_atom_sites<Model>, nb::arg("sel")=nb::none())
    .def("count_occupancies", &count_occupancies<Model>, nb::arg("sel")=nb::none())
    .def("calculate_mass", &calculate_mass<Model>)
    .def("calculate_center_of_mass", [](const Model& self) {
        return calculate_center_of_mass(self).get();
    })
    .def("calculate_b_iso_range", &calculate_b_iso_range<Model>)
    .def("calculate_b_aniso_range", &calculate_b_aniso_range)
    .def("transform_pos_and_adp", transform_pos_and_adp<Model>, nb::arg("tr"))
    .def("split_chains_by_segments", &split_chains_by_segments)
    .def("clone", [](const Model& self) { return new Model(self); })
    .def("__getstate__", &getstate<Model>)
    .def("__setstate__", &setstate<Model>)
    .def("__repr__", [](const Model& self) {
        return cat("<gemmi.Model ", self.num, " with ", self.chains.size(), " chain(s)>");
    });

  nb::class_<UniqProxy<Residue>>(m, "FirstConformerRes")
    .def("__iter__", [](UniqProxy<Residue>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  nb::class_<UniqProxy<Residue, ResidueSpan>>(m, "FirstConformerResSpan")
    .def("__iter__", [](UniqProxy<Residue, ResidueSpan>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  nb::class_<UniqProxy<Atom>>(m, "FirstConformerAtoms")
    .def("__iter__", [](UniqProxy<Atom>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  nb::class_<ResidueSpan::GroupingProxy>(m, "ResidueSpanGroups")
    .def("__iter__", [](ResidueSpan::GroupingProxy& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  pyChain
    .def(nb::init<std::string>())
    .def_rw("name", &Chain::name)
    .def("__len__", [](const Chain& ch) { return ch.residues.size(); })
    .def("__iter__", [](Chain& ch) { return usual_iterator(ch, ch.residues); },
         nb::keep_alive<0, 1>())
    .def("__getitem__", [](Chain& ch, const std::string& seqid) {
        return ch.find_residue_group(SeqId(seqid));
    }, nb::arg("pdb_seqid"), nb::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Chain, Residue>, nb::arg("index"),
         nb::rv_policy::reference_internal)
    .def("__getitem__", [](Chain &ch, const nb::slice& slice) -> nb::list {
        return getitem_slice(ch.residues, slice);
    }, nb::rv_policy::reference_internal)
    .def("__delitem__", remove_child<Chain>, nb::arg("index"))
    .def("__delitem__", remove_children<Chain>)
    .def("add_residue", add_child<Chain, Residue>,
         nb::arg("residue"), nb::arg("pos")=-1,
         nb::rv_policy::reference_internal)
    .def("subchains", (std::vector<ResidueSpan> (Chain::*)()) &Chain::subchains,
         nb::call_policy<returns_references>())
    .def("whole", (ResidueSpan (Chain::*)()) &Chain::whole,
         nb::keep_alive<0, 1>())
    .def("get_polymer", (ResidueSpan (Chain::*)()) &Chain::get_polymer,
         nb::keep_alive<0, 1>())
    .def("get_ligands", (ResidueSpan (Chain::*)()) &Chain::get_ligands,
         nb::keep_alive<0, 1>())
    .def("get_waters", (ResidueSpan (Chain::*)()) &Chain::get_waters,
         nb::keep_alive<0, 1>())
    .def("get_subchain", (ResidueSpan (Chain::*)(const std::string&)) &Chain::get_subchain,
         nb::keep_alive<0, 1>())
    .def("previous_residue", &Chain::previous_residue,
         nb::rv_policy::reference_internal)
    .def("next_residue", &Chain::next_residue,
         nb::rv_policy::reference_internal)
    .def("append_residues", &Chain::append_residues,
         nb::arg("new_residues"), nb::arg("min_sep")=0)
    .def("count_atom_sites", &count_atom_sites<Chain>, nb::arg("sel")=nb::none())
    .def("count_occupancies", &count_occupancies<Chain>, nb::arg("sel")=nb::none())
    .def("calculate_mass", &calculate_mass<Chain>)
    .def("calculate_center_of_mass", [](const Chain& self) {
          return calculate_center_of_mass(self).get();
     })
    .def("trim_to_alanine", (void (*)(Chain&)) &trim_to_alanine)
    .def("first_conformer",
         (UniqProxy<Residue> (Chain::*)()) &Chain::first_conformer,
         nb::keep_alive<0, 1>())
    .def("clone", [](const Chain& self) { return new Chain(self); })
    .def("__getstate__", &getstate<Chain>)
    .def("__setstate__", &setstate<Chain>)
    .def("__repr__", [](const Chain& self) {
        return cat("<gemmi.Chain ", self.name, " with ", self.residues.size(), " res>");
    });

  pyResidueSpan
    .def("__len__", &ResidueSpan::size)
    .def("__iter__", [](ResidueSpan& g) { return usual_iterator(g, g); },
         nb::keep_alive<0, 1>())
    .def("__bool__", [](const ResidueSpan &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](ResidueSpan& g, int index) -> Residue& {
        return g[normalize_index(index, g)];
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def("__getitem__", [](ResidueSpan& self, const std::string& seqid) {
        return self.find_residue_group(SeqId(seqid));
    }, nb::arg("pdb_seqid"), nb::keep_alive<0, 1>())
    .def("__delitem__", [](ResidueSpan& self, int index) {
        self.erase(self.begin() + normalize_index(index, self));
    }, nb::arg("index"))
    .def("add_residue", add_item<ResidueSpan, Residue>,
         nb::arg("residue"), nb::arg("pos")=-1,
         nb::rv_policy::reference_internal)
    .def("first_conformer",
         (UniqProxy<Residue, ResidueSpan> (ResidueSpan::*)())
                                                &ResidueSpan::first_conformer,
         nb::keep_alive<0, 1>())
    .def("residue_groups", &ResidueSpan::residue_groups, nb::keep_alive<0, 1>())
    .def("length", &ResidueSpan::length)
    .def("subchain_id", &ResidueSpan::subchain_id)
    .def("extract_sequence", &ResidueSpan::extract_sequence)
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
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def("__getitem__", &ResidueGroup::by_resname,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("__delitem__", &ResidueGroup::remove_residue, nb::arg("name"))
    .def("__repr__", [](const ResidueGroup& self) {
        return "<gemmi.ResidueGroup [" +
               join_str(self, ' ', [](const Residue& r) { return r.str(); }) +
               "]>";
    });

  nb::class_<AtomGroup>(m, "AtomGroup")
    .def("__len__", &AtomGroup::size)
    .def("__iter__", [](AtomGroup& g) { return usual_iterator(g, g); },
         nb::keep_alive<0, 1>())
    .def("__bool__", [](const AtomGroup &g) -> bool { return !g.empty(); })
    .def("__getitem__", [](AtomGroup& g, int index) -> Atom& {
        return g[normalize_index(index, g)];
    }, nb::arg("index"), nb::rv_policy::reference_internal)
    .def("__getitem__", &AtomGroup::by_altloc,
         nb::arg("altloc"), nb::rv_policy::reference_internal)
    .def("name", &AtomGroup::name)
    .def("__repr__", [](const AtomGroup& self) {
        return cat("<gemmi.AtomGroup ", self.name(), ", sites: ", self.size(), '>');
    });

  pyResidue
    .def(nb::init<>())
    .def_rw("subchain", &Residue::subchain)
    .def_rw("entity_id", &Residue::entity_id)
    .def_rw("label_seq", &Residue::label_seq, nb::arg().none())
    .def_rw("entity_type", &Residue::entity_type)
    .def_rw("het_flag", &Residue::het_flag)
    .def_rw("flag", &Residue::flag)
    .def_prop_ro("sifts_unp", [](const Residue& self) {
        return nb::make_tuple(self.sifts_unp.res,
                              self.sifts_unp.acc_index,
                              self.sifts_unp.num);
    })
    .def("__len__", [](const Residue& res) { return res.atoms.size(); })
    .def("__contains__", [](const Residue& res, const std::string& name) {
        return res.find_atom(name, '*') != nullptr;
    })
    .def("__iter__", [](Residue& res) {
        return usual_iterator(res, res.atoms);
    }, nb::keep_alive<0, 1>())
    .def("__getitem__", &get_child<Residue, Atom>, nb::arg("index"),
         nb::rv_policy::reference_internal)
    .def("__getitem__", &Residue::get,
         nb::arg("name"), nb::rv_policy::reference_internal)
    .def("__delitem__", remove_child<Residue>, nb::arg("index"))
    .def("__delitem__", remove_children<Residue>)
    .def("find_atom", [](Residue& self, const std::string& name, char altloc, Element el) {
           return self.find_atom(name, altloc, el);
         }, nb::arg("name"), nb::arg("altloc"), nb::arg("el")=Element(El::X),
            nb::rv_policy::reference_internal)
    .def("remove_atom",
         [](Residue& self, const std::string& name, char altloc, Element el) {
           self.atoms.erase(self.find_atom_iter(name, altloc, el));
    }, nb::arg("name"), nb::arg("altloc"), nb::arg("el")=Element(El::X))
    .def("add_atom", add_child<Residue, Atom>,
         nb::arg("atom"), nb::arg("pos")=-1,
         nb::rv_policy::reference_internal)
    .def("first_conformer",
         (UniqProxy<Atom> (Residue::*)()) &Residue::first_conformer,
         nb::keep_alive<0, 1>())
    .def("sole_atom", &Residue::sole_atom)
    .def("get_ca", &Residue::get_ca, nb::rv_policy::reference_internal)
    .def("get_p", &Residue::get_p, nb::rv_policy::reference_internal)
    .def("is_water", &Residue::is_water)
    .def("remove_hydrogens", &remove_hydrogens<Residue>)
    .def("trim_to_alanine", (bool (*)(Residue&)) &trim_to_alanine)
    .def("recommended_het_flag", &recommended_het_flag)
    .def("clone", [](const Residue& self) { return new Residue(self); })
    .def("__getstate__", &getstate<Residue>)
    .def("__setstate__", &setstate<Residue>)
    .def("__repr__", [](const Residue& self) {
        return cat("<gemmi.Residue ", self.str(), " with ", self.atoms.size(), " atoms>");
    });

  nb::enum_<CalcFlag>(m, "CalcFlag")
    .value("NotSet", CalcFlag::NotSet)
    .value("NoHydrogen", CalcFlag::NoHydrogen)
    .value("Determined", CalcFlag::Determined)
    .value("Calculated", CalcFlag::Calculated)
    .value("Dummy", CalcFlag::Dummy);

  pyAtom
    .def(nb::init<>())
    .def_rw("name", &Atom::name)
    .def_rw("altloc", &Atom::altloc)
    .def_rw("charge", &Atom::charge)
    .def_rw("element", &Atom::element)
    .def_prop_rw("pos",
        [](Atom& self) -> Position& { return self.pos; },
        [](Atom& self, const std::variant<Position, std::array<double,3>>& v) {
          if (const Position* p = std::get_if<0>(&v)) {
            self.pos = *p;
          } else {
            const std::array<double,3>& t = *std::get_if<1>(&v);
            self.pos.x = t[0];
            self.pos.y = t[1];
            self.pos.z = t[2];
          }
        })
    .def_rw("occ", &Atom::occ)
    .def_rw("b_iso", &Atom::b_iso)
    .def_rw("aniso", &Atom::aniso)
    .def_rw("serial", &Atom::serial)
    .def_rw("fraction", &Atom::fraction)
    .def_rw("calc_flag", &Atom::calc_flag)
    .def_rw("flag", &Atom::flag)
    .def_rw("tls_group_id", &Atom::tls_group_id)
    .def("is_hydrogen", &Atom::is_hydrogen)
    .def("has_altloc", &Atom::has_altloc)
    .def("b_eq", &Atom::b_eq)
    .def("padded_name", &Atom::padded_name)
    .def("clone", [](const Atom& self) { return new Atom(self); })
    .def("__getstate__", &getstate<Atom>)
    .def("__setstate__", &setstate<Atom>)
    .def("__repr__", [](const Atom& self) {
        std::string r = "<gemmi.Atom " + self.name;
        if (self.altloc) {
            r += '.';
            r += self.altloc;
        }
        char buf[128];
        snprintf_z(buf, 128, " at (%.1f, %.1f, %.1f)>",
                   self.pos.x, self.pos.y, self.pos.z);
        return r + buf;
    });

  m.def("calculate_b_est", &calculate_b_est);
  m.def("calculate_angle", &calculate_angle,
        "Input: three points. Output: angle in radians.");
  m.def("calculate_dihedral", &calculate_dihedral,
        "Input: four points. Output: dihedral angle in radians.");
  m.def("calculate_omega", &calculate_omega,
        nb::arg("residue"), nb::arg("next_residue"));
  m.def("calculate_phi_psi", &calculate_phi_psi,
        nb::arg("prev_residue").none(), nb::arg("residue"), nb::arg("next_residue").none());
  m.def("find_best_plane", &find_best_plane,
        nb::arg("atoms"));
  m.def("get_distance_from_plane", &get_distance_from_plane,
        nb::arg("pos"), nb::arg("coeff"));
  m.def("parse_triplet_as_ftransform", &parse_triplet_as_ftransform);
  m.def("calculate_u_from_tls", &calculate_u_from_tls);
  m.def("make_assembly", &make_assembly,
        nb::arg("assembly"), nb::arg("model"), nb::arg("how"), nb::arg("logging")=nb::none());
  m.def("expand_ncs_model", &expand_ncs_model);
  m.def("merge_atoms_in_expanded_model", &merge_atoms_in_expanded_model,
        nb::arg("model"), nb::arg("cell"), nb::arg("max_dist")=0.2,
        nb::arg("compare_serial")=true);

  // select.hpp
  nb::class_<FilterProxy<Selection, Model>> pySelectionModelsProxy(m, "SelectionModelsProxy");
  nb::class_<FilterProxy<Selection, Chain>> pySelectionChainsProxy(m, "SelectionChainsProxy");
  nb::class_<FilterProxy<Selection, Residue>> pySelectionResiduesProxy(m, "SelectionResiduesProxy");
  nb::class_<FilterProxy<Selection, Atom>> pySelectionAtomsProxy(m, "SelectionAtomsProxy");

  pySelection
    .def(nb::init<>())
    .def(nb::init<const std::string&>())
    .def("models", &Selection::models)
    .def("chains", &Selection::chains)
    .def("residues", &Selection::residues)
    .def("atoms", &Selection::atoms)
    .def("first_in_model", &Selection::first_in_model,
         nb::keep_alive<1, 2>())
    .def("first", &Selection::first, nb::rv_policy::reference,
         nb::keep_alive<1, 2>())
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
    })
    .def("str", &Selection::str);

  pySelectionModelsProxy
    .def("__iter__", [](FilterProxy<Selection, Model>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  pySelectionChainsProxy
    .def("__iter__", [](FilterProxy<Selection, Chain>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  pySelectionResiduesProxy
    .def("__iter__", [](FilterProxy<Selection, Residue>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());

  pySelectionAtomsProxy
    .def("__iter__", [](FilterProxy<Selection, Atom>& self) {
        return usual_iterator(self, self);
    }, nb::keep_alive<0, 1>());
}
