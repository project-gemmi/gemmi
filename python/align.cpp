// Copyright 2018 Global Phasing Ltd.

#include "gemmi/align.hpp"     // for align_sequence_to_polymer
#include "gemmi/seqalign.hpp"  // for align_string_sequences

#include "common.h"
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

using namespace gemmi;

void add_alignment(nb::module_& m) {

  // sequence alignment
  nb::class_<AlignmentResult>(m, "AlignmentResult")
    .def_ro("score", &AlignmentResult::score)
    .def_ro("match_count", &AlignmentResult::match_count)
    .def_ro("match_string", &AlignmentResult::match_string)
    .def("cigar_str", &AlignmentResult::cigar_str)
    .def("calculate_identity", &AlignmentResult::calculate_identity,
         nb::arg("which")=0)
    .def("add_gaps", &AlignmentResult::add_gaps, nb::arg("s"), nb::arg("which"))
    .def("formatted", &AlignmentResult::formatted)
    ;

  nb::class_<AlignmentScoring>(m, "AlignmentScoring")
    .def("__init__", [](AlignmentScoring* t, char what) {
          const AlignmentScoring* s = AlignmentScoring::simple();
          if (what == 'p')
            s = AlignmentScoring::partial_model();
          else if (what == 'b')
            s = AlignmentScoring::blosum62();
          new(t) AlignmentScoring(*s);
    }, nb::arg("what")='s')
    .def_rw("match", &AlignmentScoring::match)
    .def_rw("mismatch", &AlignmentScoring::mismatch)
    .def_rw("gapo", &AlignmentScoring::gapo)
    .def_rw("gape", &AlignmentScoring::gape)
    .def_rw("good_gapo", &AlignmentScoring::good_gapo)
    .def_rw("bad_gapo", &AlignmentScoring::bad_gapo)
    ;

  m.def("align_string_sequences", &align_string_sequences,
        nb::arg("query"), nb::arg("target"), nb::arg("target_gapo"),
        nb::arg("scoring")=nb::none());
  m.def("align_sequence_to_polymer",
        [](const std::vector<std::string>& full_seq, const ResidueSpan& polymer,
           PolymerType polymer_type, AlignmentScoring* scoring) {
      return align_sequence_to_polymer(full_seq, polymer, polymer_type, scoring);
  }, nb::arg("full_seq"), nb::arg("polymer"),
     nb::arg("polymer_type"), nb::arg("scoring")=nb::none());

  // structure superposition
  nb::enum_<SupSelect>(m, "SupSelect")
    .value("CaP", SupSelect::CaP)
    .value("MainChain", SupSelect::MainChain)
    .value("All", SupSelect::All);

  nb::class_<SupResult>(m, "SupResult")
    .def_ro("rmsd", &SupResult::rmsd)
    .def_ro("count", &SupResult::count)
    .def_ro("center1", &SupResult::center1)
    .def_ro("center2", &SupResult::center2)
    .def_ro("transform", &SupResult::transform)
    ;

  m.def("calculate_current_rmsd",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           SupSelect sel, char altloc) {
          return calculate_current_rmsd(fixed, movable, ptype, sel, altloc);
        }, nb::arg("fixed"), nb::arg("movable"), nb::arg("ptype"), nb::arg("sel"),
           nb::arg("altloc")='\0');
  m.def("calculate_superposition",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           SupSelect sel, int trim_cycles, double trim_cutoff, char altloc) {
          return calculate_superposition(fixed, movable, ptype, sel,
                                         trim_cycles, trim_cutoff, altloc);
        }, nb::arg("fixed"), nb::arg("movable"), nb::arg("ptype"), nb::arg("sel"),
           nb::arg("trim_cycles")=0, nb::arg("trim_cutoff")=2.0,
           nb::arg("altloc")='\0');
  m.def("calculate_superpositions_in_moving_window",
        [](const ResidueSpan& fixed, const ResidueSpan& movable, PolymerType ptype,
           double radius) {
          return calculate_superpositions_in_moving_window(fixed, movable, ptype, radius);
        }, nb::arg("fixed"), nb::arg("movable"), nb::arg("ptype"), nb::arg("radius")=10.0);

  m.def("superpose_positions",
        [](std::vector<Position> pos1, std::vector<Position> pos2,
           const std::vector<double>& weight) {
          if (pos1.size() != pos2.size())
            fail("superpose_positions: pos1 and pos2 must have equal lengths");
          if (!weight.empty() && weight.size() != pos1.size())
            fail("superpose_positions: weights must be empty or of the same length as pos1/pos2");
          return superpose_positions(pos1.data(), pos2.data(), pos1.size(),
                                     weight.empty() ? nullptr : weight.data());
        }, nb::arg("pos1"), nb::arg("pos2"), nb::arg("weight")=std::vector<int>{});
}

void add_assign_label_seq_id(nb::class_<Structure>& structure) {
  structure
    .def("assign_label_seq_id", &assign_label_seq_id, nb::arg("force")=false)
    .def("clear_sequences", &clear_sequences)
    .def("assign_best_sequences", &assign_best_sequences, nb::arg("fasta_sequences"))
    ;
}
