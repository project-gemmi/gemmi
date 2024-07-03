
namespace gemmi {
  std::ostream& operator<< (std::ostream& os, const Entity& ent);
}

namespace nanobind { namespace detail {
  template<> struct type_caster<gemmi::SeqId::OptionalNum> {
    NB_TYPE_CASTER(gemmi::SeqId::OptionalNum, optional_name(const_name("int")))
    using Caster =  make_caster<int>;

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
      if (src.is_none()) {
        value = Value();
        return true;
      }
      Caster caster;
      if (!caster.from_python(src, flags_for_local_caster<int>(flags), cleanup))
        return false;
      value = caster.operator cast_t<int>();
      return true;
    }

    static handle from_cpp(Value value, rv_policy policy, cleanup_list *cleanup) noexcept {
      if (!value)
        return none().release();
      return Caster::from_cpp(*value, policy, cleanup);
    }
  };
}} // namespace nanobind::detail

NB_MAKE_OPAQUE(std::vector<gemmi::Helix>)
NB_MAKE_OPAQUE(std::vector<gemmi::Sheet>)
NB_MAKE_OPAQUE(std::vector<gemmi::Sheet::Strand>)
NB_MAKE_OPAQUE(std::vector<gemmi::Assembly::Gen>)
NB_MAKE_OPAQUE(std::vector<gemmi::Assembly::Operator>)
NB_MAKE_OPAQUE(std::vector<gemmi::Assembly>)

NB_MAKE_OPAQUE(std::vector<gemmi::Connection>)
NB_MAKE_OPAQUE(std::vector<gemmi::NcsOp>)
NB_MAKE_OPAQUE(std::vector<gemmi::Entity>)
