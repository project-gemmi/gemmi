
#pragma once
#include <pybind11/pybind11.h>

void add_elem(pybind11::module& m); // elem.cpp
void add_symmetry(pybind11::module& m); // sym.cpp
void add_ccp4(pybind11::module& m); // ccp4.cpp
void add_grid(pybind11::module& m); // grid.cpp
void add_recgrid(pybind11::module& m); // recgrid.cpp
void add_unitcell(pybind11::module& m); // unitcell.cpp
void add_hkl(pybind11::module& m); // hkl.cpp
void add_meta(pybind11::module& m); // meta.cpp
void add_mol(pybind11::module& m); // mol.cpp
void add_mtz(pybind11::module& m); // mtz.cpp
void add_cif(pybind11::module& cif); // cif.cpp
void add_cif_read(pybind11::module& cif); // read.cpp
void add_read_structure(pybind11::module& m); // read.cpp
void add_small(pybind11::module& m); // read.cpp
void add_chemcomp(pybind11::module& m); // chemcomp.cpp
void add_monlib(pybind11::module& m); // monlib.cpp
void add_topo(pybind11::module& m); // topo.cpp
void add_alignment(pybind11::module& m); // align.cpp
void add_scaling(pybind11::module& m); // scaling.cpp
void add_search(pybind11::module& m); // search.cpp
void add_sf(pybind11::module& m); // sf.cpp
void add_custom(pybind11::module& m); // custom.cpp

// defined in write.cpp
namespace gemmi {
  struct Structure;
  namespace cif { struct Document; }
}

void add_write(pybind11::module& m, pybind11::class_<gemmi::Structure>& structure);
// defined in align.cpp
void add_assign_label_seq_id(pybind11::class_<gemmi::Structure>& structure);

// defined in read.cpp
void cif_parse_string(gemmi::cif::Document& doc, const std::string& data);
void cif_parse_file(gemmi::cif::Document& doc, const std::string& filename);


template<typename T> int normalize_index(int index, const T& container) {
  if (index < 0)
    index += (int) container.size();
  if ((size_t) index >= container.size())
    throw pybind11::index_error();
  return index;
}

template<typename Item>
void delitem_at_index(std::vector<Item>& items, pybind11::ssize_t idx) {
  items.erase(items.begin() + idx);
}

template<typename Item>
void delitem_range(std::vector<Item>& items, pybind11::ssize_t start, pybind11::ssize_t end) {
  items.erase(items.begin() + start, items.begin() + end);
}

template<typename Items>
void delitem_slice(Items& items, const pybind11::slice& slice) {
  namespace py = pybind11;
  py::ssize_t start, stop, step, slength;
  if (!slice.compute((py::ssize_t)items.size(), &start, &stop, &step, &slength))
    throw py::error_already_set();
  if (step == 1) {
    delitem_range(items, start, start + slength);
  } else {
    for (int i = 0; i < slength; ++i)
      delitem_at_index(items, start + (step > 0 ? slength - 1 - i : i) * step);
  }
}

template<typename Items>
pybind11::list getitem_slice(Items& items, const pybind11::slice& slice) {
  namespace py = pybind11;
  py::ssize_t start, stop, step, slength;
  if (!slice.compute((py::ssize_t)items.size(), &start, &stop, &step, &slength))
    throw py::error_already_set();
  py::list l;
  for (py::ssize_t i = 0; i < slength; ++i)
    l.append(py::cast(&items[start + i * step]));
  return l;
}
