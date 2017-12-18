
//struct cGridF;
struct cGridC;
//typedef struct cGridF cGridF;
typedef struct cGridC cGridC;

cGridC* GridC_init(int nx, int ny, int nz);
void GridC_set_unit_cell(cGridC* grid, double a, double b, double c,
                         double alpha, double beta, double gamma);
void GridC_mask_atom(cGridC* grid, double x, double y, double z, double radius);
void GridC_apply_space_group(cGridC* grid, int ccp4_num);
signed char* GridC_data(cGridC* grid);
void GridC_free(cGridC* grid);
