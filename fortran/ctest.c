
#include "symmetry.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  int i;
  for (i = 1; i < argc; ++i) {
    char short_hm[16];
    const cSpaceGroup* sg = find_spacegroup_by_name(argv[i]);
    if (!sg) {
      printf("n/a\n");
      continue;
    }
    SpaceGroup_short_name(sg, short_hm);
    printf("space group %d  %s  or  %s",
            SpaceGroup_number(sg), SpaceGroup_hm(sg), short_hm);
    cGroupOps* ops = SpaceGroup_operations(sg);
    printf("   order: %d\n", GroupOps_order(ops));
    GroupOps_free(ops);
  }
}
