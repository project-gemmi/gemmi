#pragma once
#include "common.h"
#include <nanobind/stl/detail/nb_optional.h>
#include "gemmi/metadata.hpp"

namespace nanobind { namespace detail {
  template<> struct type_caster<gemmi::SeqId::OptionalNum>
    : optional_caster<gemmi::SeqId::OptionalNum> {};
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
