
#include <cstring>
#include "gemmi/symmetry.hpp"

extern "C" {
#include "symmetry.h"
}

using gemmi::SpaceGroup;
using gemmi::GroupOps;

static const geSpaceGroup* as_c(const SpaceGroup* sg) {
  return reinterpret_cast<const geSpaceGroup*>(sg);
}
static geGroupOps* as_c(GroupOps* sg) {
  return reinterpret_cast<geGroupOps*>(sg);
}

static const SpaceGroup* as_cpp(const geSpaceGroup* sg) {
  return reinterpret_cast<const SpaceGroup*>(sg);
}
static GroupOps* as_cpp(geGroupOps* sg) {
  return reinterpret_cast<GroupOps*>(sg);
}

extern "C" {

const geSpaceGroup* find_spacegroup_by_name(const char* name) {
  return as_c(gemmi::find_spacegroup_by_name(name));
}

const geSpaceGroup* find_spacegroup_by_number(int n) {
  return as_c(gemmi::find_spacegroup_by_number(n));
}

int geSpaceGroup_number(const geSpaceGroup* sg) {
  return sg ? as_cpp(sg)->number : 0;
}

const char* geSpaceGroup_hm(const geSpaceGroup* sg) {
  return sg ? as_cpp(sg)->hm : nullptr;
}

const char* geSpaceGroup_hall(const geSpaceGroup* sg) {
  return sg ? as_cpp(sg)->hall : nullptr;
}

void geSpaceGroup_short_name(const geSpaceGroup* sg, char* dest) {
  if (sg && dest)
    std::strcpy(dest, as_cpp(sg)->short_name().c_str());
}

geGroupOps* geSpaceGroup_operations(const geSpaceGroup* sg) {
  GroupOps* ops = new GroupOps;
  *ops = as_cpp(sg)->operations();
  return as_c(ops);
}

int geGroupOps_order(geGroupOps* ops) {
  return as_cpp(ops)->order();
}

void geGroupOps_free(geGroupOps* ops) {
  delete as_cpp(ops);
}

}

// vim:sw=2:ts=2:et
