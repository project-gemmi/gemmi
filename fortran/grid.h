
#include <stdint.h>

//struct geGrid2;  // gemmi::Grid<float> -- corresponds to ccp4 mode 2
struct geGrid0;  // gemmi::Grid<int8_t> -- corresponds to ccp4 mode 0
//typedef struct geGrid2 geGrid2;
typedef struct geGrid0 geGrid0;

geGrid0* geGrid0_init(int nx, int ny, int nz);
void geGrid0_set_unit_cell(geGrid0* grid, double a, double b, double c,
                           double alpha, double beta, double gamma);
void geGrid0_mask_atom(geGrid0* grid, double x, double y, double z,
                       double radius);
void geGrid0_apply_space_group(geGrid0* grid, int ccp4_num);
int8_t* geGrid0_data(geGrid0* grid);
int8_t geGrid0_get_value(geGrid0* grid, int u, int v, int w);
void geGrid0_free(geGrid0* grid);
void geGrid0_update_ccp4_header(geGrid0* grid, int mode, bool update_stats);
void geGrid0_write_ccp4_map(geGrid0* grid, const char* path);
