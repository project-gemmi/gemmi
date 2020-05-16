// Copyright 2018 Global Phasing Ltd.

#include "gemmi/chemcomp.hpp"
#include "gemmi/gzread.hpp" // read_cif_gz
#include "gemmi/linkhunt.hpp"
#include "gemmi/monlib.hpp"
#include "gemmi/topo.hpp"
#include "gemmi/seqalign.hpp"  // for align_string_sequences
#include "gemmi/labelseq.hpp"  // for align_sequence_to_polymer
#include "gemmi/select.hpp"

#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace gemmi;

PYBIND11_MAKE_OPAQUE(std::vector<Restraints::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<ChemComp::Atom>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::Bond>)
PYBIND11_MAKE_OPAQUE(std::vector<Topo::ExtraLink>)
using monomers_type = std::map<std::string, ChemComp>;
using links_type = std::map<std::string, ChemLink>;
using modifications_type = std::map<std::string, ChemMod>;
PYBIND11_MAKE_OPAQUE(monomers_type)
PYBIND11_MAKE_OPAQUE(links_type)
PYBIND11_MAKE_OPAQUE(modifications_type)

void add_monlib(py::module& m) {

  py::class_<ChemMod> pyChemMod(m, "ChemMod");
  py::class_<ChemLink> chemlink(m, "ChemLink");
  py::class_<ChemComp> chemcomp(m, "ChemComp");

  py::class_<ChemLink::Side> pyChemLinkSide(chemlink, "Side");

  py::class_<Restraints> restraints(m, "Restraints");
  py::class_<Restraints::Bond> pyRestraintsBond(restraints, "Bond");

  py::class_<Topo> topo(m, "Topo");
  py::class_<Topo::Bond> pyTopoBond(topo, "Bond");
  py::class_<Topo::ExtraLink> pyTopoExtraLink(topo, "ExtraLink");

  py::class_<ChemComp::Atom> pyChemCompAtom(chemcomp, "Atom");

  py::bind_vector<std::vector<Restraints::Bond>>(m, "RestraintsBonds");
  py::bind_vector<std::vector<ChemComp::Atom>>(m, "ChemCompAtoms");
  py::bind_vector<std::vector<Topo::Bond>>(m, "TopoBonds");
  py::bind_vector<std::vector<Topo::ExtraLink>>(m, "TopoExtraLinks");
  py::bind_map<monomers_type>(m, "ChemCompMap");
  py::bind_map<links_type>(m, "ChemLinkMap");
  py::bind_map<modifications_type>(m, "ChemModMap");

  py::enum_<BondType>(m, "BondType")
    .value("Unspec", BondType::Unspec)
    .value("Single", BondType::Single)
    .value("Double", BondType::Double)
    .value("Triple", BondType::Triple)
    .value("Aromatic", BondType::Aromatic)
    .value("Deloc", BondType::Deloc)
    .value("Metal", BondType::Metal);


  py::class_<Restraints::AtomId>(restraints, "AtomId")
    .def_readwrite("comp", &Restraints::AtomId::comp)
    .def_readwrite("atom", &Restraints::AtomId::atom)
    .def("get_from",
         (Atom* (Restraints::AtomId::*)(Residue&, Residue*, char) const)
         &Restraints::AtomId::get_from,
         py::arg("res1"), py::arg("res2"), py::arg("altloc"),
         py::return_value_policy::reference)
    ;
  pyRestraintsBond
    .def_readwrite("id1", &Restraints::Bond::id1)
    .def_readwrite("id2", &Restraints::Bond::id2)
    .def_readwrite("type", &Restraints::Bond::type)
    .def_readwrite("aromatic", &Restraints::Bond::aromatic)
    .def_readwrite("value", &Restraints::Bond::value)
    .def_readwrite("esd", &Restraints::Bond::esd)
    .def("lexicographic_str", &Restraints::Bond::lexicographic_str)
    .def("__repr__", [](const Restraints::Bond& self) {
        return "<gemmi.Restraints.Bond " + self.str() + ">";
    });
  restraints
    .def_readwrite("bonds", &Restraints::bonds)
    .def("empty", &Restraints::empty)
    .def("get_bond", &Restraints::get_bond,
         py::return_value_policy::reference_internal)
    .def("get_bond", [](Restraints& self,
                        const std::string& a1, const std::string& a2)
                                                      -> Restraints::Bond& {
            auto it = self.find_bond(a1, a2);
            if (it == self.bonds.end())
              fail("Bond restraint not found: " + a1 + "-" + a2);
            return *it;
         }, py::return_value_policy::reference_internal)
    ;

  pyChemCompAtom
    .def_readonly("id", &ChemComp::Atom::id)
    .def_readonly("el", &ChemComp::Atom::el)
    .def_readonly("charge", &ChemComp::Atom::charge)
    .def_readonly("chem_type", &ChemComp::Atom::chem_type)
    .def("is_hydrogen", &ChemComp::Atom::is_hydrogen)
    ;
  chemcomp
    .def_readonly("name", &ChemComp::name)
    .def_readonly("group", &ChemComp::group)
    .def_readonly("atoms", &ChemComp::atoms)
    .def_readonly("rt", &ChemComp::rt)
    .def("get_atom", &ChemComp::get_atom)
    .def("remove_hydrogens", &ChemComp::remove_hydrogens)
    ;
  m.def("make_chemcomp_from_block", &make_chemcomp_from_block);

  chemlink
    .def_readwrite("id", &ChemLink::id)
    .def_readwrite("name", &ChemLink::name)
    .def_readwrite("side1", &ChemLink::side1)
    .def_readwrite("side2", &ChemLink::side2)
    .def_readwrite("rt", &ChemLink::rt)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemLink " + self.id + ">";
    });

  py::enum_<ChemLink::Group>(pyChemLinkSide, "Group")
      .value("Peptide",  ChemLink::Group::Peptide)
      .value("PPeptide", ChemLink::Group::PPeptide)
      .value("MPeptide", ChemLink::Group::MPeptide)
      .value("Pyranose", ChemLink::Group::Pyranose)
      .value("DnaRna",   ChemLink::Group::DnaRna)
      .value("Null",     ChemLink::Group::Null);

  pyChemLinkSide
    .def_readwrite("comp", &ChemLink::Side::comp)
    .def_readwrite("mod", &ChemLink::Side::mod)
    .def_readwrite("group", &ChemLink::Side::group)
    .def("__repr__", [](const ChemLink::Side& self) {
        return "<gemmi.ChemLink.Side " + self.comp + "/" +
               ChemLink::group_str(self.group) + ">";
    });

  pyChemMod
    .def_readwrite("id", &ChemMod::id)
    .def("__repr__", [](const ChemLink& self) {
        return "<gemmi.ChemMod " + self.id + ">";
    });

  py::class_<MonLib>(m, "MonLib")
    .def(py::init<>())
    .def_readonly("monomers", &MonLib::monomers)
    .def_readonly("links", &MonLib::links)
    .def_readonly("modifications", &MonLib::modifications)
    .def("find_link", &MonLib::find_link, py::arg("link_id"),
         py::return_value_policy::reference_internal)
    .def("add_monomer_if_present", &MonLib::add_monomer_if_present)
    .def("add_monomers_if_present", &MonLib::add_monomers_if_present)
    .def("__repr__", [](const MonLib& self) {
        return "<gemmi.MonLib with " +
               std::to_string(self.monomers.size()) + " monomers, " +
               std::to_string(self.links.size()) + " links, " +
               std::to_string(self.modifications.size()) + " modifications>";
    });

  pyTopoBond
    .def_readonly("restr", &Topo::Bond::restr)
    .def_readonly("atoms", &Topo::Bond::atoms)
    .def("calculate", &Topo::Bond::calculate)
    .def("calculate_z", &Topo::Bond::calculate_z)
    ;
  pyTopoExtraLink
    .def_readonly("res1", &Topo::ExtraLink::res1)
    .def_readonly("res2", &Topo::ExtraLink::res2)
    .def_readonly("alt1", &Topo::ExtraLink::alt1)
    .def_readonly("alt2", &Topo::ExtraLink::alt2)
    .def_readonly("link_id", &Topo::ExtraLink::link_id)
    ;
  topo
    .def(py::init<>())
    .def("prepare_refmac_topology",
         [](Topo& self, Structure& st, MonLib& monlib) {
        self.initialize_refmac_topology(st, st.first_model(), monlib);
        self.finalize_refmac_topology(monlib);
    })
    .def_readonly("bonds", &Topo::bonds)
    .def_readonly("extras", &Topo::extras)
    ;

  m.def("read_monomer_lib", [](const std::string& monomer_dir,
                               const std::vector<std::string>& resnames) {
    return read_monomer_lib(monomer_dir, resnames, gemmi::read_cif_gz);
  });
  m.def("read_monomer_cif", [](const std::string& path) {
    return read_monomer_cif(path, gemmi::read_cif_gz);
  });

  py::class_<BondIndex>(m, "BondIndex")
    .def(py::init<const Model&>(), py::keep_alive<1, 2>())
    .def("add_link", &BondIndex::add_link)
    .def("add_monomer_bonds", &BondIndex::add_monomer_bonds)
    .def("are_linked", &BondIndex::are_linked)
    .def("graph_distance", &BondIndex::graph_distance,
         py::arg("a"), py::arg("b"), py::arg("same_index"),
         py::arg("max_distance")=4)
    ;

  py::class_<ContactSearch> contactsearch(m, "ContactSearch");
  py::class_<ContactSearch::Result> pyContactSearchResult(contactsearch, "Result");
  py::enum_<ContactSearch::Ignore> pyContactSearchIgnore(contactsearch, "Ignore");

  contactsearch
    .def(py::init<float>())
    .def_readwrite("search_radius", &ContactSearch::search_radius)
    .def_readwrite("ignore", &ContactSearch::ignore)
    .def_readwrite("twice", &ContactSearch::twice)
    .def_readwrite("special_pos_cutoff_sq", &ContactSearch::special_pos_cutoff_sq)
    .def_readwrite("min_occupancy", &ContactSearch::min_occupancy)
    .def("setup_atomic_radii", &ContactSearch::setup_atomic_radii)
    .def("get_radius", [](const ContactSearch& self, Element el) {
        return self.get_radius(el.elem);
    })
    .def("set_radius", [](ContactSearch& self, Element el, float r) {
        self.set_radius(el.elem, r);
    })
    .def("find_contacts", &ContactSearch::find_contacts)
    ;

  pyContactSearchIgnore
    .value("Nothing", ContactSearch::Ignore::Nothing)
    .value("SameResidue", ContactSearch::Ignore::SameResidue)
    .value("AdjacentResidues", ContactSearch::Ignore::AdjacentResidues)
    .value("SameChain", ContactSearch::Ignore::SameChain)
    .value("SameAsu", ContactSearch::Ignore::SameAsu);

  pyContactSearchResult
    .def_readonly("partner1", &ContactSearch::Result::partner1)
    .def_readonly("partner2", &ContactSearch::Result::partner2)
    .def_readonly("image_idx", &ContactSearch::Result::image_idx)
    .def_property_readonly("dist", [](ContactSearch::Result& self) {
        return std::sqrt(self.dist_sq);
    })
    ;

  py::class_<LinkHunt> linkhunt(m, "LinkHunt");
  py::class_<LinkHunt::Match> pyLinkHuntMatch(linkhunt, "Match");
  linkhunt
    .def(py::init<>())
    .def("index_chem_links", &LinkHunt::index_chem_links,
         py::arg("monlib"), py::keep_alive<1, 2>())
    .def("find_possible_links", &LinkHunt::find_possible_links,
         py::arg("st"), py::arg("bond_margin"), py::arg("radius_margin"),
         py::arg("ignore")=ContactSearch::Ignore::SameResidue)
    ;

  pyLinkHuntMatch
    .def_readonly("chem_link", &LinkHunt::Match::chem_link)
    .def_readonly("chem_link_count", &LinkHunt::Match::chem_link_count)
    .def_readonly("cra1", &LinkHunt::Match::cra1)
    .def_readonly("cra2", &LinkHunt::Match::cra2)
    .def_readonly("same_image", &LinkHunt::Match::same_image)
    .def_readonly("bond_length", &LinkHunt::Match::bond_length)
    .def_readonly("conn", &LinkHunt::Match::conn)
    ;
}

void add_select(py::module& m) {
  py::class_<Selection> pySelection(m, "Selection");
  py::class_<FilterProxy<Selection, Model>> pySelectionModelsProxy(m, "SelectionModelsProxy");
  py::class_<FilterProxy<Selection, Chain>> pySelectionChainsProxy(m, "SelectionChainsProxy");
  py::class_<FilterProxy<Selection, Residue>> pySelectionResidusProxy(m, "SelectionResidusProxy");
  py::class_<FilterProxy<Selection, Atom>> pySelectionAtomsProxy(m, "SelectionAtomsProxy");

  m.def("parse_cid", &parse_cid);
  pySelection
    .def(py::init<>())
    .def("models", &Selection::models)
    .def("chains", &Selection::chains)
    .def("residues", &Selection::residues)
    .def("atoms", &Selection::atoms)
    .def("first_in_model", &Selection::first_in_model,
         py::keep_alive<1, 2>())
    .def("first", &Selection::first, py::return_value_policy::reference,
         py::keep_alive<1, 2>())
    .def("to_cid", &Selection::to_cid)
    .def("set_residue_flags", &Selection::set_residue_flags)
    .def("set_atom_flags", &Selection::set_atom_flags)
    .def("copy_model_selection", &Selection::copy_selection<Model>)
    .def("copy_structure_selection", &Selection::copy_selection<Structure>)
    .def("__repr__", [](const Selection& self) {
        return "<gemmi.Selection CID: " + self.to_cid() + ">";
    });

    pySelectionModelsProxy
    .def("__iter__", [](FilterProxy<Selection, Model>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

    pySelectionChainsProxy
    .def("__iter__", [](FilterProxy<Selection, Chain>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

    pySelectionResidusProxy
    .def("__iter__", [](FilterProxy<Selection, Residue>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());

    pySelectionAtomsProxy
    .def("__iter__", [](FilterProxy<Selection, Atom>& self) {
        return py::make_iterator(self);
    }, py::keep_alive<0, 1>());
}

void add_alignment(py::module& m) {
  py::class_<AlignmentResult>(m, "AlignmentResult")
    .def_readonly("score", &AlignmentResult::score)
    .def_readonly("match_count", &AlignmentResult::match_count)
    .def_readonly("match_string", &AlignmentResult::match_string)
    .def("cigar_str", &AlignmentResult::cigar_str)
    .def("calculate_identity", &AlignmentResult::calculate_identity,
         py::arg("which")=0)
    .def("add_gaps", &AlignmentResult::add_gaps, py::arg("s"), py::arg("which"))
    .def("formatted", &AlignmentResult::formatted)
    ;

  py::class_<AlignmentScoring>(m, "AlignmentScoring")
    .def(py::init<>())
    .def_readwrite("match", &AlignmentScoring::match)
    .def_readwrite("mismatch", &AlignmentScoring::mismatch)
    .def_readwrite("gapo", &AlignmentScoring::gapo)
    .def_readwrite("gape", &AlignmentScoring::gape)
    ;

  m.def("prepare_blosum62_scoring", &prepare_blosum62_scoring);
  m.def("align_string_sequences", &align_string_sequences,
        py::arg("query"), py::arg("target"), py::arg("free_gapo"),
        py::arg_v("scoring", AlignmentScoring(), "gemmi.AlignmentScoring()"));
  m.def("align_sequence_to_polymer",
        [](const std::vector<std::string>& full_seq,
           const ResidueSpan& polymer, PolymerType polymer_type,
           AlignmentScoring& sco) {
      return align_sequence_to_polymer(full_seq, polymer, polymer_type, sco);
  }, py::arg("full_seq"), py::arg("polymer"), py::arg("polymer_type"),
     py::arg_v("scoring", AlignmentScoring(), "gemmi.AlignmentScoring()"));
}
