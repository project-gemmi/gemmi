
#include <pybind11/numpy.h>

template<typename Ret, typename Obj, typename Func>
pybind11::array_t<Ret>
miller_function(const Obj& obj, Func func, pybind11::array_t<int> hkl) {
  auto h = hkl.unchecked<2>();
  if (h.shape(1) != 3)
    throw std::domain_error("error: the size of the second dimension != 3");
  auto result = pybind11::array_t<Ret>(h.shape(0));
  pybind11::buffer_info rbuf = result.request();
  Ret* rptr = (Ret*) rbuf.ptr;
  for (pybind11::ssize_t i = 0; i < h.shape(0); ++i)
    rptr[i] = (obj.*func)({{h(i, 0), h(i, 1), h(i, 2)}});
  return result;
}
