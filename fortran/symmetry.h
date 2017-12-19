

struct geSpaceGroup;
struct geGroupOps;
typedef struct geSpaceGroup geSpaceGroup;
typedef struct geGroupOps geGroupOps;

const geSpaceGroup* find_spacegroup_by_name(const char* name);
const geSpaceGroup* find_spacegroup_by_number(int n);
int geSpaceGroup_number(const geSpaceGroup* sg);
const char* geSpaceGroup_hm(const geSpaceGroup* sg);
const char* geSpaceGroup_hall(const geSpaceGroup* sg);
void geSpaceGroup_short_name(const geSpaceGroup* sg, char* dest);
geGroupOps* geSpaceGroup_operations(const geSpaceGroup* sg);
int geGroupOps_order(geGroupOps* ops);
void geGroupOps_free(geGroupOps* ops);
