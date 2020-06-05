
#pragma once
#include <pybind11/pybind11.h>

template<typename T> int normalize_index(int index, const T& container) {
  if (index < 0)
    index += (int) container.size();
  if ((size_t) index >= container.size())
    throw pybind11::index_error();
  return index;
}


