

struct cSpaceGroup;
struct cGroupOps;
typedef struct cSpaceGroup cSpaceGroup;
typedef struct cGroupOps cGroupOps;

const cSpaceGroup* find_spacegroup_by_name(const char* name);
const cSpaceGroup* find_spacegroup_by_number(int n);
int SpaceGroup_number(const cSpaceGroup* sg);
const char* SpaceGroup_hm(const cSpaceGroup* sg);
const char* SpaceGroup_hall(const cSpaceGroup* sg);
void SpaceGroup_short_name(const cSpaceGroup* sg, char* dest);
cGroupOps* SpaceGroup_operations(const cSpaceGroup* sg);
int GroupOps_order(cGroupOps* ops);
void GroupOps_free(cGroupOps* ops);
