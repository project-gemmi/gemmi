
#pragma once

#if defined(__clang__)
  #pragma clang diagnostic push
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wpedantic"
#endif

#include <nanobind/nanobind.h>  // IWYU pragma: export
#include <gemmi/logger.hpp>     // for Logger

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#if NB_VERSION_MAJOR < 2 || (NB_VERSION_MAJOR == 2 && NB_VERSION_MINOR < 2)
  #error Required nanobind version >= 2.2
#endif

namespace nb = nanobind;
constexpr auto rv_ri = nb::rv_policy::reference_internal;

void add_elem(nb::module_& m); // elem.cpp
void add_xds(nb::module_& m); // elem.cpp
void add_symmetry(nb::module_& m); // sym.cpp
void add_ccp4(nb::module_& m); // ccp4.cpp
void add_grid(nb::module_& m); // grid.cpp
void add_recgrid(nb::module_& m); // recgrid.cpp
void add_unitcell(nb::module_& m); // unitcell.cpp
void add_hkl(nb::module_& m); // hkl.cpp
void add_meta(nb::module_& m); // meta.cpp
void add_mol(nb::module_& m); // mol.cpp
void add_mtz(nb::module_& m); // mtz.cpp
void add_cif(nb::module_& cif); // cif.cpp
void add_cif_read(nb::module_& cif); // read.cpp
void add_read_structure(nb::module_& m); // read.cpp
void add_small(nb::module_& m); // read.cpp
void add_chemcomp(nb::module_& m); // chemcomp.cpp
void add_monlib(nb::module_& m); // monlib.cpp
void add_topo(nb::module_& m); // topo.cpp
void add_alignment(nb::module_& m); // align.cpp
void add_scaling(nb::module_& m); // scaling.cpp
void add_search(nb::module_& m); // search.cpp
void add_sf(nb::module_& m); // sf.cpp
void add_custom(nb::module_& m); // custom.cpp

// defined in write.cpp
namespace gemmi {
  struct Structure;
  namespace cif { struct Document; }
}

void add_write(nb::module_& m, nb::class_<gemmi::Structure>& structure);
// defined in align.cpp
void add_assign_label_seq_id(nb::class_<gemmi::Structure>& structure);

// convert pythonic index to C++ index
inline int c_index(int index, size_t size) {
  if (index < 0)
    index += (int) size;
  if ((size_t) index >= size)
    throw nb::index_error();
  return index;
}

template<typename T> int normalize_index(int index, const T& container) {
  return c_index(index, container.size());
}

// specialized for cif::Table in cif.cpp
template<typename Vec>
void delitem_at_index(Vec& items, size_t idx) {
  items.erase(items.begin() + idx);
}

// specialized for cif::Table in cif.cpp
template<typename Vec>
void delitem_range(Vec& items, size_t start, size_t end) {
  items.erase(items.begin() + start, items.begin() + end);
}

template<typename Items>
void delitem_slice(Items& items, const nb::slice& slice) {
  auto [start, stop, step, length] = slice.compute(items.size());
  if (step == 1) {
    delitem_range(items, start, start + length);
  } else {
    for (size_t i = 0; i < length; ++i)
      delitem_at_index(items, start + (step > 0 ? length - 1 - i : i) * step);
  }
}

template<typename Items>
nb::list getitem_slice(Items& items, const nb::slice& slice) {
  auto [start, stop, step, length] = slice.compute(items.size());
  nb::list l;
  for (size_t i = 0; i < length; ++i)
    l.append(nb::cast(&items[start + i * step]));
  return l;
}

// for numpy __array__ method
inline nb::object handle_numpy_array_args(const nb::object& o, nb::handle dtype, nb::handle copy) {
  if (dtype.is_none() || dtype.is(o.attr("dtype"))) {
    if (copy.ptr() != Py_True)
      return o;
    dtype = o.attr("dtype");
  }
  if (copy.ptr() == Py_False)  // astype() would copy even with copy=False
    throw nb::value_error("Unable to avoid copy while creating an array as requested.");
  return o.attr("astype")(dtype);
}

namespace nanobind { namespace detail {
template <> struct type_caster<gemmi::Logger> {
  NB_TYPE_CASTER(gemmi::Logger, const_name("object"))
  bool from_python(handle src, uint8_t, cleanup_list *) noexcept {
    value = {};
    if (PyTuple_Check(src.ptr()) && PyTuple_Size(src.ptr()) == 2) {
      value.threshold = (int) PyLong_AsLong(PyTuple_GetItem(src.ptr(), 1));
      if (value.threshold == -1 && PyErr_Occurred())
        return false;
      src = PyTuple_GetItem(src.ptr(), 0);
    }
    if (src.is_none()) {
      // nothing
    } else if (nb::hasattr(src, "write") && nb::hasattr(src, "flush")) {
      value.callback = {[obj=nb::borrow(src)](const std::string& s) {
        obj.attr("write")(s + "\n");
        obj.attr("flush")();
      }};
    } else if (PyCallable_Check(src.ptr())) {
      value.callback = {[obj=nb::borrow(src)](const std::string& s) { obj(s); }};
    } else {
      return false;
    }
    return true;
  }
};
}} // namespace nanobind::detail
