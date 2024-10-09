#pragma once
#include <array>
#include <vector>
#include "common.h"
#include <nanobind/ndarray.h>

template<typename T>
using cpu_array = nb::ndarray<T, nb::shape<-1>, nb::device::cpu>;
template<typename T>
using cpu_c_array = nb::ndarray<T, nb::shape<-1>, nb::device::cpu, nb::c_contig>;
using cpu_miller_array = nb::ndarray<int, nb::shape<-1,3>, nb::device::cpu>;
using cpu_c_miller_array = nb::ndarray<int, nb::shape<-1,3>, nb::device::cpu, nb::c_contig>;

template<typename T>
auto numpy_array_from_vector(std::vector<T>&& original_vec) {
  using V = std::vector<T>;
  auto v = new V(std::move(original_vec));
  nb::capsule owner(v, [](void* p) noexcept { delete static_cast<V*>(p); });
  return nb::ndarray<nb::numpy, T, nb::shape<-1>>(v->data(), {v->size()}, owner);
}

template<typename T, size_t N>
auto py_array2d_from_vector(std::vector<std::array<T,N>>&& original_vec) {
  using V = std::vector<std::array<T,N>>;
  auto v = new V(std::move(original_vec));
  nb::capsule owner(v, [](void* p) noexcept { delete static_cast<V*>(p); });
  return nb::ndarray<nb::numpy, T, nb::shape<-1, N>>(v->data(), {v->size(), N}, owner);
}

// to be used with rv_policy::reference_internal
template<typename T, typename S>
auto vector_member_array(std::vector<S>& vec, T S::*ptr) {
  constexpr int64_t stride = static_cast<int64_t>(sizeof(S) / sizeof(T));
  static_assert(stride * sizeof(T) == sizeof(S),
                "vector_member_array(): problem with stride");
  return nb::ndarray<nb::numpy, T, nb::shape<-1>>(
          &(vec.data()->*ptr),
          {vec.size()},
          nb::handle(),
          {stride});
}

template<typename T>
auto make_numpy_array(std::initializer_list<size_t> size,
                      std::initializer_list<int64_t> strides={}) {
  size_t total_size = 1;
  for (size_t i : size)
    total_size *= i;
  T* c_array = new T[total_size];
  nb::capsule owner(c_array, [](void* p) noexcept { delete [] static_cast<T*>(p); });
  return nb::ndarray<nb::numpy, T>(c_array, size, owner, strides);
}

template<typename Ret, typename Obj, typename Func>
auto miller_function(const Obj& obj, Func func, const cpu_miller_array& hkl) {
  auto h = hkl.view();
  size_t n = h.shape(0);
  auto result = make_numpy_array<Ret>({n});
  Ret* rptr = result.data();
  for (size_t i = 0; i < n; ++i)
    rptr[i] = (obj.*func)({h(i, 0), h(i, 1), h(i, 2)});
  return result;
}

