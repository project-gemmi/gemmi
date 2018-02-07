// Copyright 2017 Global Phasing Ltd.

#include "gemmi/elem.hpp"
#include "gemmi/unitcell.hpp"
#include "gemmi/model.hpp"
#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/to_pdb.hpp"

#include <cstdio>
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace gemmi;

static std::string triple(double x, double y, double z) {
  using namespace std;  // VS2015/17 doesn't like std::snprintf
  char buf[128];
  snprintf(buf, 128, "%g, %g, %g", x, y, z);
  return std::string(buf);
}

namespace pybind11 { namespace detail {
  template<> struct type_caster<ResidueId::OptionalNum>
    : optional_caster<ResidueId::OptionalNum> {};
}} // namespace pybind11::detail


void add_mol(py::module& m) {
  py::class_<Position>(m, "Position")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Position::x)
    .def_readwrite("y", &Position::y)
    .def_readwrite("z", &Position::z)
    .def("__repr__", [](const Position& self) {
        return "<gemmi.Position(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<Fractional>(m, "Fractional")
    .def(py::init<double,double,double>())
    .def_readwrite("x", &Fractional::x)
    .def_readwrite("y", &Fractional::y)
    .def_readwrite("z", &Fractional::z)
    .def("__repr__", [](const Fractional& self) {
        return "<gemmi.Fractional(" + triple(self.x, self.y, self.z) + ")>";
    });
  py::class_<UnitCell>(m, "UnitCell")
    .def(py::init<>())
    .def(py::init([](double a, double b, double c,
                     double alpha, double beta, double gamma) {
           UnitCell cell;
           cell.set(a, b, c, alpha, beta, gamma);
           return cell;
         }), py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("alpha"), py::arg("beta"), py::arg("gamma"))
    .def_readonly("a", &UnitCell::a)
    .def_readonly("b", &UnitCell::b)
    .def_readonly("c", &UnitCell::c)
    .def_readonly("alpha", &UnitCell::alpha)
    .def_readonly("beta", &UnitCell::beta)
    .def_readonly("gamma", &UnitCell::gamma)
    .def_readonly("volume", &UnitCell::volume)
    .def("set", &UnitCell::set)
    .def("fractionalize", &UnitCell::fractionalize)
    .def("orthogonalize", &UnitCell::orthogonalize)
    .def("__repr__", [](const UnitCell& self) {
        return "<gemmi.UnitCell(" + triple(self.a, self.b, self.c)
             + ", " + triple(self.alpha, self.beta, self.gamma) + ")>";
    });

  m.def("calculate_dihedral", &calculate_dihedral,
        "Input: four points. Output: dihedral angle in radians.");
  py::class_<Element>(m, "Element")
    .def(py::init<const std::string &>())
    .def(py::init<int>())
    .def_property_readonly("name", &Element::name)
    .def_property_readonly("weight", &Element::weight)
    .def_property_readonly("atomic_number", &Element::atomic_number)
    .def("__repr__", [](const Element& self) {
        return "<gemmi.Element: " + std::string(self.name()) + ">";
    });

  py::class_<Structure>(m, "Structure")
    .def(py::init<>())
    .def_readwrite("name", &Structure::name)
    .def_readwrite("cell", &Structure::cell)
    .def_readwrite("sg_hm", &Structure::sg_hm)
    .def("__len__", [](const Structure& st) { return st.models.size(); })
    .def("__iter__", [](const Structure& st) {
        return py::make_iterator(st.models);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Structure& st, int index) -> Model& {
        return st.models.at(index >= 0 ? index : index + st.models.size());
    }, py::arg("index"), py::return_value_policy::reference_internal)
    .def("find_or_add_model", &Structure::find_or_add_model,
         py::arg("name"), py::return_value_policy::reference_internal)
    .def("write_pdb", [](const Structure& st, const std::string& path) {
       std::ofstream f(path.c_str());
       write_pdb(st, f);
    }, py::arg("path"))
    .def("write_minimal_pdb",
         [](const Structure& st, const std::string& path, const char* chain) {
       std::ofstream f(path.c_str());
       write_minimal_pdb(st, f, chain);
    }, py::arg("path"), py::arg("chain")=nullptr);

  py::class_<Model>(m, "Model")
    .def(py::init<std::string>())
    .def_readwrite("name", &Model::name)
    .def("__len__", [](const Model& mdl) { return mdl.chains.size(); })
    .def("__iter__", [](const Model& mdl) {
        return py::make_iterator(mdl.chains);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Model& mdl, const std::string& name) -> Chain& {
        Chain* ch = mdl.find_chain(name);
        if (!ch)
          throw py::key_error("chain '" + name + "' does not exist");
        return *ch;
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("find_or_add_chain", &Model::find_or_add_chain,
         py::arg("name"), py::return_value_policy::reference_internal);

  py::class_<Chain>(m, "Chain")
    .def(py::init<std::string>())
    .def_readwrite("name", &Chain::name)
    .def_readwrite("auth_name", &Chain::auth_name)
    .def("__len__", [](const Chain& ch) { return ch.residues.size(); })
    .def("__iter__", [](const Chain& ch) {
        return py::make_iterator(ch.residues);
    }, py::keep_alive<0, 1>())
    .def("append_residues", &Chain::append_residues)
    // we need __getitem__ here, but what should be the argument?
    // label_seq_id+label_comp_id or seq-num+icode+resname
    ;

  py::class_<Residue>(m, "Residue")
    .def(py::init<>())
    .def_readwrite("name", &Residue::name)
    .def_readwrite("label_seq", &Residue::label_seq)
    .def_readwrite("seq_num", &Residue::seq_num)
    .def_property("icode",
        [](const Residue& r) { return r.icode ? std::string(1, r.icode) : ""; },
        [](Residue& r, const char* c) { r.icode = c ? *c : '\0'; })
    .def_readwrite("segment", &Residue::segment)
    .def("__len__", [](const Residue& res) { return res.atoms.size(); })
    .def("__iter__", [](const Residue& res) {
        return py::make_iterator(res.atoms);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](Residue& res, const std::string& name) -> Atom& {
        Atom* atom = res.find_atom(name);
        if (!atom)
          throw py::key_error("residue has no atom '" + name + "'");
        return *atom;
    }, py::arg("name"), py::return_value_policy::reference_internal)
    .def("__repr__", [](const Residue& self) {
        std::string r = "<gemmi.Residue " + self.name + " " + self.seq_id();
        if (self.label_seq)
          r += " (" + self.label_seq.str() + ")";
        return r + " with " + std::to_string(self.atoms.size()) + " atoms>";
    });

  py::class_<Atom>(m, "Atom")
    .def(py::init<>())
    .def_readwrite("name", &Atom::name)
    .def_readwrite("altloc", &Atom::altloc)
    .def_readwrite("charge", &Atom::charge)
    .def_readwrite("element", &Atom::element)
    .def_readwrite("pos", &Atom::pos)
    .def_readwrite("occ", &Atom::occ)
    .def_readwrite("b_iso", &Atom::b_iso);
}
