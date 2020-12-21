
namespace gemmi {
  std::ostream& operator<< (std::ostream& os, const Entity& ent);
}

namespace pybind11 { namespace detail {
  template<> struct type_caster<gemmi::SeqId::OptionalNum>
    : optional_caster<gemmi::SeqId::OptionalNum> {};
}} // namespace pybind11::detail

PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly::Gen>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly::Operator>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly>)
