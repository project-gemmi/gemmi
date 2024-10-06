#include "common.h"
#include <nanobind/make_iterator.h>

// This wrapper simplifies the call to nb::make_iterator (by assuming
// the name "iterator") and changes the default rv_policy to what it was
// in pybind11 and nanobind<2.0.
template<nb::rv_policy Policy = nb::rv_policy::reference_internal, typename S, typename T>
auto usual_iterator(const S&, T& value) {
  return nb::make_iterator<Policy>(nb::type<S>(), "iterator", value);
}
