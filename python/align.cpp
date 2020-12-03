// Copyright 2018 Global Phasing Ltd.

#include "gemmi/align.hpp"     // for align_sequence_to_polymer
#include "gemmi/seqalign.hpp"  // for align_string_sequences
#include "gemmi/select.hpp"

#include "common.h"
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;


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

  // sequence alignment
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

  // structure superposition
  py::enum_<SupSelect>(m, "SupSelect")
    .value("CaP", SupSelect::CaP)
    .value("All", SupSelect::All);

  py::class_<SupResult>(m, "SupResult")
    .def_readonly("rmsd", &SupResult::rmsd)
    .def_readonly("count", &SupResult::count)
    .def_readonly("center1", &SupResult::center1)
    .def_readonly("center2", &SupResult::center2)
    .def_readonly("transform", &SupResult::transform)
    .def("apply", &apply_superposition)
    ;

  m.def("calculate_superposition",
        [](const ResidueSpan& fixed, const ResidueSpan& movable,
           PolymerType ptype, SupSelect sel, char altloc, bool current_rmsd) {
          return calculate_superposition(fixed, movable, ptype, sel, altloc, current_rmsd);
        }, py::arg("fixed"), py::arg("movable"), py::arg("ptype"), py::arg("sel"),
           py::arg("altloc")='\0', py::arg("current_rmsd")=false);

  m.def("superpose_positions",
        [](std::vector<Position> pos1, std::vector<Position> pos2,
           const std::vector<double>& weight) {
          return superpose_positions(pos1.data(), pos2.data(), pos1.size(),
                                     weight.empty() ? nullptr : weight.data());
        }, py::arg("pos1"), py::arg("pos2"), py::arg("weight")=std::vector<int>{});
}

void add_assign_label_seq_id(py::class_<Structure>& structure) {
  structure
    .def("assign_label_seq_id", &assign_label_seq_id, py::arg("force")=false);
}
