#pragma once
#include <pybind11/numpy.h>

template<typename T>
pybind11::array_t<T> py_array_from_vector(std::vector<T>&& original_vec) {
  auto v = new std::vector<T>(std::move(original_vec));
  pybind11::capsule cap(v, [](void* p) { delete (std::vector<T>*) p; });
  return pybind11::array_t<T>(v->size(), v->data(), cap);
}
