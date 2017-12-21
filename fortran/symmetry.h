

struct geSpaceGroup;
struct geGroupOps;
struct geOp;
typedef struct geSpaceGroup geSpaceGroup;
typedef struct geGroupOps geGroupOps;
typedef struct geOp geOp;

const geSpaceGroup* ge_find_spacegroup_by_name(const char* name);
const geSpaceGroup* ge_find_spacegroup_by_number(int n);
int geSpaceGroup_number(const geSpaceGroup* sg);
const char* geSpaceGroup_hm(const geSpaceGroup* sg);
const char* geSpaceGroup_hall(const geSpaceGroup* sg);
void geSpaceGroup_short_name(const geSpaceGroup* sg, char* dest);
geGroupOps* geSpaceGroup_operations(const geSpaceGroup* sg);
int geGroupOps_order(geGroupOps* ops);
geOp* geGroupOps_get_op(geGroupOps* ops, int n);
void geGroupOps_free(geGroupOps* ops);
void geOp_triplet(geOp* op, char* dest);
void geOp_free(geOp* op);
