
#include <cstring>
#include "gemmi/symmetry.hpp"

extern "C" {
#include "symmetry.h"
}

using gemmi::SpaceGroup;
using gemmi::GroupOps;
using gemmi::Op;

static const geSpaceGroup* as_c(const SpaceGroup* sg) {
  return reinterpret_cast<const geSpaceGroup*>(sg);
}
static geGroupOps* as_c(GroupOps* sg) {
  return reinterpret_cast<geGroupOps*>(sg);
}
static geOp* as_c(Op* op) {
  return reinterpret_cast<geOp*>(op);
}

static const SpaceGroup* as_cpp(const geSpaceGroup* sg) {
  return reinterpret_cast<const SpaceGroup*>(sg);
}
static GroupOps* as_cpp(geGroupOps* sg) {
  return reinterpret_cast<GroupOps*>(sg);
}
static Op* as_cpp(geOp* op) {
  return reinterpret_cast<Op*>(op);
}

extern "C" {

const geSpaceGroup* ge_find_spacegroup_by_name(const char* name) {
  return as_c(gemmi::find_spacegroup_by_name(name));
}

const geSpaceGroup* ge_find_spacegroup_by_number(int n) {
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

// char[16] for dest is more than enough
void geSpaceGroup_short_name(const geSpaceGroup* sg, char* dest) {
  if (sg && dest) {
    std::strncpy(dest, as_cpp(sg)->short_name().c_str(), 15);
    dest[15] = '\0';
  }
}

geGroupOps* geSpaceGroup_operations(const geSpaceGroup* sg) {
  GroupOps* ops = new GroupOps;
  *ops = as_cpp(sg)->operations();
  return as_c(ops);
}

int geGroupOps_order(geGroupOps* ops) {
  return as_cpp(ops)->order();
}

geOp* geGroupOps_get_op(geGroupOps* ops, int n) {
  Op* op = new Op;
  *op = as_cpp(ops)->get_op(n);
  return as_c(op);
}

void geGroupOps_free(geGroupOps* ops) {
  delete as_cpp(ops);
}

// char[32] should be enough for dest
void geOp_triplet(geOp* op, char* dest) {
  if (op && dest) {
    std::strncpy(dest, as_cpp(op)->triplet().c_str(), 31);
    dest[31] = '\0';
  }
}

void geOp_free(geOp* op) {
  delete as_cpp(op);
}

}
