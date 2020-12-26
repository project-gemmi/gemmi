
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
void add_cif(pybind11::module& cif); // cif.cpp
void add_read_structure(pybind11::module& m); // read.cpp
void add_cif_read(pybind11::module& cif); // read.cpp
void add_monlib(pybind11::module& m); // monlib.cpp
void add_alignment(pybind11::module& m); // align.cpp
void add_select(pybind11::module& m); // align.cpp
void add_search(pybind11::module& m); // search.cpp
void add_sf(pybind11::module& m); // sf.cpp

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
