
namespace gemmi {
  std::ostream& operator<< (std::ostream& os, const Entity& ent);
}

namespace pybind11 { namespace detail {
  template<> struct type_caster<gemmi::SeqId::OptionalNum>
    : optional_caster<gemmi::SeqId::OptionalNum> {};
}} // namespace pybind11::detail

PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Helix>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Sheet>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Sheet::Strand>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly::Gen>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly::Operator>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Assembly>)

PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Connection>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::NcsOp>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::Entity>)

PYBIND11_MAKE_OPAQUE(std::vector<gemmi::ReflectionsInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::DiffractionInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::TlsGroup::Selection>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::BasicRefinementInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::RefinementInfo::Restr>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::TlsGroup>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::ExperimentInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::CrystalInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::RefinementInfo>)
PYBIND11_MAKE_OPAQUE(std::vector<gemmi::SoftwareItem>)
