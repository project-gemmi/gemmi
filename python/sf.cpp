// Copyright 2020 Global Phasing Ltd.

#include "common.h"
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "gemmi/it92.hpp"
#include "gemmi/c4322.hpp"
#include "gemmi/neutron92.hpp"
#include "gemmi/sfcalc.hpp"   // for StructureFactorCalculator
#include "gemmi/dencalc.hpp"  // for DensityCalculator
#include "gemmi/fprime.hpp"   // for add_cl_fprime_for_all_elements

namespace py = pybind11;
using gemmi::Element;

template<typename Table>
void add_sfcalc(py::module& m, const char* name, bool with_mb) {
  using SFC = gemmi::StructureFactorCalculator<Table>;
  py::class_<SFC> sfc(m, name);
  sfc
    .def(py::init<const gemmi::UnitCell&>())
    .def_readwrite("addends", &SFC::addends)
    .def("calculate_sf_from_model", &SFC::calculate_sf_from_model)
    .def("calculate_sf_from_small_structure", &SFC::calculate_sf_from_small_structure);
  if (with_mb)
    sfc
      .def("mott_bethe_factor", &SFC::mott_bethe_factor)
      .def("calculate_mb_z", &SFC::calculate_mb_z,
           py::arg("model"), py::arg("hkl"), py::arg("only_h")=false);
}

template<typename Table>
void add_dencalc(py::module& m, const char* name) {
  using DenCalc = gemmi::DensityCalculator<Table, float>;
  py::class_<DenCalc>(m, name)
    .def(py::init<>())
    .def_readonly("grid", &DenCalc::grid)
    .def_readwrite("d_min", &DenCalc::d_min)
    .def_readwrite("rate", &DenCalc::rate)
    .def_readwrite("blur", &DenCalc::blur)
    .def_readwrite("cutoff", &DenCalc::cutoff)
    .def_readwrite("addends", &DenCalc::addends)
    .def("set_refmac_compatible_blur", &DenCalc::set_refmac_compatible_blur)
    .def("put_model_density_on_grid", &DenCalc::put_model_density_on_grid)
    .def("initialize_grid", &DenCalc::initialize_grid)
    .def("add_model_density_to_grid", &DenCalc::add_model_density_to_grid)
    .def("add_atom_density_to_grid", &DenCalc::add_atom_density_to_grid)
    .def("add_c_contribution_to_grid", &DenCalc::add_c_contribution_to_grid)
    .def("set_grid_cell_and_spacegroup", &DenCalc::set_grid_cell_and_spacegroup)
    .def("reciprocal_space_multiplier", &DenCalc::reciprocal_space_multiplier)
    .def("mott_bethe_factor", &DenCalc::mott_bethe_factor)
    .def("estimate_radius", [](const DenCalc &self, const gemmi::Atom &atom){
        double b;
        if (!atom.aniso.nonzero())
            b = atom.b_iso + self.blur;
        else {
            gemmi::SMat33<double> aniso_b = atom.aniso.scaled(gemmi::u_to_b()).added_kI(self.blur);
            // rough estimate, so we don't calculate eigenvalues
            b = std::max(std::max(aniso_b.u11, aniso_b.u22), aniso_b.u33);
        }
        auto coef = Table::get(atom.element);
        auto precal = coef.precalculate_density_iso(b, self.addends.get(atom.element));
        return self.estimate_radius(precal, b);
    })
    ;
}

void add_sf(py::module& m) {
  py::class_<gemmi::Addends>(m, "Addends")
    .def("set", &gemmi::Addends::set)
    .def("get", &gemmi::Addends::get)
    .def("clear", &gemmi::Addends::clear)
    .def("add_cl_fprime", [](gemmi::Addends& self, double energy) {
        gemmi::add_cl_fprime_for_all_elements(&self.values[1], energy);
    }, py::arg("energy"))
    .def("subtract_z", &gemmi::Addends::subtract_z,
         py::arg("except_hydrogen")=false)
    ;

  using IT92 = gemmi::IT92<double>;
  using C4322 = gemmi::C4322<double>;
  using Neutron92 = gemmi::Neutron92<double>;
  add_sfcalc<IT92>(m, "StructureFactorCalculatorX", true);
  add_sfcalc<C4322>(m, "StructureFactorCalculatorE", false);
  add_sfcalc<Neutron92>(m, "StructureFactorCalculatorN", false);
  add_dencalc<IT92>(m, "DensityCalculatorX");
  add_dencalc<C4322>(m, "DensityCalculatorE");
  add_dencalc<Neutron92>(m, "DensityCalculatorN");
  m.def("IT92_normalize", &IT92::normalize);
}
