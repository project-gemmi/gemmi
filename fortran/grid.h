
//struct geGridF;
struct geGridC;
//typedef struct geGridF geGridF;
typedef struct geGridC geGridC;

geGridC* GridC_init(int nx, int ny, int nz);
void GridC_set_unit_cell(geGridC* grid, double a, double b, double c,
                         double alpha, double beta, double gamma);
void GridC_mask_atom(geGridC* grid, double x, double y, double z,
                     double radius);
void GridC_apply_space_group(geGridC* grid, int ccp4_num);
signed char* GridC_data(geGridC* grid);
void GridC_free(geGridC* grid);
