
#include "symmetry.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  int i;
  for (i = 1; i < argc; ++i) {
    char short_hm[16];
    const geSpaceGroup* sg = find_spacegroup_by_name(argv[i]);
    if (!sg) {
      printf("n/a\n");
      continue;
    }
    geSpaceGroup_short_name(sg, short_hm);
    printf("space group %d  %s  or  %s",
           geSpaceGroup_number(sg), geSpaceGroup_hm(sg), short_hm);
    geGroupOps* ops = geSpaceGroup_operations(sg);
    printf("   order: %d\n", geGroupOps_order(ops));
    geGroupOps_free(ops);
  }
}
