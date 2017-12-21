
#include "symmetry.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  int i, j;
  for (i = 1; i < argc; ++i) {
    char buffer[32];
    const geSpaceGroup* sg = ge_find_spacegroup_by_name(argv[i]);
    if (sg) {
      geGroupOps* ops = geSpaceGroup_operations(sg);
      int n = geGroupOps_order(ops);
      geSpaceGroup_short_name(sg, buffer);
      printf("space group %d  %s  or  %s    order: %d\n",
             geSpaceGroup_number(sg), geSpaceGroup_hm(sg), buffer, n);
      for (j = 0; j < n; ++j) {
          geOp* op = geGroupOps_get_op(ops, j);
          geOp_triplet(op, buffer);
          geOp_free(op);
          printf("   %s\n", buffer);
      }
      geGroupOps_free(ops);
    } else {
      printf("%s: not found\n", argv[i]);
    }
  }
}
