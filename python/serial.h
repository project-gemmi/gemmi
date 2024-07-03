// Copyright 2024 Global Phasing Ltd.

#pragma once
#include "common.h"
#include "gemmi/serialize.hpp"
#include "../third_party/serializer.h"

template<typename T>
nb::bytes getstate(const T& self) {
  std::vector<unsigned char> data;
  zpp::serializer::memory_output_archive out(data);
  // awkward zpp design: calling zpp::serializer::*_archive::operator()()
  out(self);
  return nb::bytes(data.data(), data.size());
}

template<typename T>
void setstate(T& self, const nb::bytes& state) {
  new(&self) T();  // probably not needed
  const unsigned char* ptr = (unsigned char*) state.data();
  zpp::serializer::memory_view_input_archive(ptr, state.size())(self);
}
