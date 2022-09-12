// Copyright 2018 Global Phasing Ltd.

#include "gemmi/align.hpp"     // for align_sequence_to_polymer
#include "gemmi/seqalign.hpp"  // for align_string_sequences

#include "common.h"
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;


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
    .value("MainChain", SupSelect::MainChain)
    .value("All", SupSelect::All);

  py::class_<SupResult>(m, "SupResult")
    .def_readonly("rmsd", &SupResult::rmsd)
    .def_readonly("count", &SupResult::count)
    .def_readonly("center1", &SupResult::center1)
    .def_readonly("center2", &SupResult::center2)
    .def_readonly("transform", &SupResult::transform)
    ;

  m.def("calculate_current_rmsd",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           SupSelect sel, char altloc) {
          return calculate_current_rmsd(fixed, movable, ptype, sel, altloc);
        }, py::arg("fixed"), py::arg("movable"), py::arg("ptype"), py::arg("sel"),
           py::arg("altloc")='\0');
  m.def("calculate_superposition",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           SupSelect sel, int trim_cycles, double trim_cutoff, char altloc) {
          return calculate_superposition(fixed, movable, ptype, sel,
                                         trim_cycles, trim_cutoff, altloc);
        }, py::arg("fixed"), py::arg("movable"), py::arg("ptype"), py::arg("sel"),
           py::arg("trim_cycles")=0, py::arg("trim_cutoff")=2.0,
           py::arg("altloc")='\0');
  m.def("calculate_superpositions_in_moving_window",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           double radius) {
          return calculate_superpositions_in_moving_window(fixed, movable, ptype, radius);
        }, py::arg("fixed"), py::arg("movable"), py::arg("ptype"), py::arg("radius")=10.0);

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
