
#include <stdint.h>

//struct geMap;  // gemmi::Ccp4<float> -- corresponds to ccp4 mode 2
struct geMask;  // gemmi::Ccp4<int8_t> -- corresponds to ccp4 mode 0
//typedef struct geGrid2 geGrid2;
typedef struct geMask geMask;

geMask* geMask_init(int nx, int ny, int nz);
void geMask_set_unit_cell(geMask* grid, double a, double b, double c,
                          double alpha, double beta, double gamma);
void geMask_mask_atom(geMask* grid, double x, double y, double z,
                      double radius);
void geMask_apply_space_group(geMask* grid, int ccp4_num);
int8_t* geMask_data(geMask* grid);
int8_t geMask_get_value(geMask* grid, int u, int v, int w);
void geMask_free(geMask* grid);
void geMask_update_ccp4_header(geMask* grid, int mode, bool update_stats);
void geMask_write_ccp4_map(geMask* grid, const char* path);
