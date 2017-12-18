
#include <cstring>
#include "gemmi/symmetry.hpp"

extern "C" {
#include "symmetry.h"
}

using gemmi::SpaceGroup;
using gemmi::GroupOps;

static const cSpaceGroup* as_c(const SpaceGroup* sg) {
  return reinterpret_cast<const cSpaceGroup*>(sg);
}
static cGroupOps* as_c(GroupOps* sg) {
  return reinterpret_cast<cGroupOps*>(sg);
}

static const SpaceGroup* as_cpp(const cSpaceGroup* sg) {
  return reinterpret_cast<const SpaceGroup*>(sg);
}
static GroupOps* as_cpp(cGroupOps* sg) {
  return reinterpret_cast<GroupOps*>(sg);
}

extern "C" {

const cSpaceGroup* find_spacegroup_by_name(const char* name) {
  return as_c(gemmi::find_spacegroup_by_name(name));
}

const cSpaceGroup* find_spacegroup_by_number(int n) {
  return as_c(gemmi::find_spacegroup_by_number(n));
}

int SpaceGroup_number(const cSpaceGroup* sg) {
  return sg ? as_cpp(sg)->number : 0;
}

const char* SpaceGroup_hm(const cSpaceGroup* sg) {
  return sg ? as_cpp(sg)->hm : nullptr;
}

const char* SpaceGroup_hall(const cSpaceGroup* sg) {
  return sg ? as_cpp(sg)->hall : nullptr;
}

void SpaceGroup_short_name(const cSpaceGroup* sg, char* dest) {
  if (sg && dest)
    std::strcpy(dest, as_cpp(sg)->short_name().c_str());
}

cGroupOps* SpaceGroup_operations(const cSpaceGroup* sg) {
  GroupOps* ops = new GroupOps;
  *ops = as_cpp(sg)->operations();
  return as_c(ops);
}

int GroupOps_order(cGroupOps* ops) {
  return as_cpp(ops)->order();
}

int GroupOps_free(cGroupOps* ops) {
  delete as_cpp(ops);
}

}

// vim:sw=2:ts=2:et
