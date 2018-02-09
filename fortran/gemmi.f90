
module gemmi
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: spacegroup, groupops, op, grid0, find_spacegroup_by_name, &
            find_spacegroup_by_number

  interface
    ! functions from symmetry.h

    ! const geSpaceGroup* find_spacegroup_by_name(const char* name);
    type(c_ptr) function c_find_spacegroup_by_name(name) &
    bind(C, name="ge_find_spacegroup_by_name")
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*)
    end function

    ! const geSpaceGroup* find_spacegroup_by_number(int n);
    type(c_ptr) function c_find_spacegroup_by_number(n) &
    bind(C, name="ge_find_spacegroup_by_number")
      use iso_c_binding
      integer(c_int), value :: n
    end function

    ! int geSpaceGroup_number(const geSpaceGroup* sg);
    integer(c_int) function c_spacegroup_number(sg) &
    bind(C, name="geSpaceGroup_number")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
    end function

    ! const char* geSpaceGroup_hm(const geSpaceGroup* sg);
    type(c_ptr) function c_spacegroup_hm(sg) &
    bind(C, name="geSpaceGroup_hm")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
    end function

    ! const char* geSpaceGroup_hall(const geSpaceGroup* sg);
    type(c_ptr) function c_spacegroup_hall(sg) &
    bind(C, name="geSpaceGroup_hall")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
    end function

    ! void geSpaceGroup_short_name(const geSpaceGroup* sg, char* dest);
    subroutine c_spacegroup_short_name(sg, dest) &
    bind(C, name="geSpaceGroup_short_name")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
      character(kind=c_char), intent(out) :: dest(*)
    end subroutine

    ! geGroupOps* geSpaceGroup_operations(const geSpaceGroup* sg);
    type(c_ptr) function c_spacegroup_operations(sg) &
    bind(C, name="geSpaceGroup_operations")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
    end function

    ! int geGroupOps_order(geGroupOps* ops);
    integer(c_int) function c_groupops_order(ops) &
    bind(C, name="geGroupOps_order")
      use iso_c_binding
      type(c_ptr), intent(in), value :: ops
    end function

    ! geOp* geGroupOps_get_op(geGroupOps* ops, int n);
    type(c_ptr) function c_groupops_get_op(ops, n) &
    bind(C, name="geGroupOps_get_op")
      use iso_c_binding
      type(c_ptr), intent(in), value :: ops
      integer(c_int), intent(in), value :: n
    end function

    ! void geGroupOps_free(geGroupOps* ops);
    subroutine c_groupops_free(ops) bind(C, name="geGroupOps_free")
      use iso_c_binding
      type(c_ptr), intent(in), value :: ops
    end subroutine

    ! void geOp_triplet(geOp* op, char* dest);
    subroutine c_op_triplet(op, dest) bind(C, name="geOp_triplet")
      use iso_c_binding
      type(c_ptr), intent(in), value :: op
      character(kind=c_char), intent(out) :: dest(*)
    end subroutine

    ! void geOp_free(geOp* op);
    subroutine c_op_free(op) bind(C, name="geOp_free")
      use iso_c_binding
      type(c_ptr), intent(in), value :: op
    end subroutine

    ! functions from grid.h

    !geGrid0* geGrid0_init(int nx, int ny, int nz);
    type(c_ptr) function c_grid0_init(nx, ny, nz) &
    bind(C, name="geGrid0_init")
      use iso_c_binding
      integer(c_int), intent(in), value :: nx, ny, nz
    end function

    !void geGrid0_set_unit_cell(geGrid0* grid, double a, double b, double c,
    !                           double alpha, double beta, double gamma);
    subroutine c_grid0_set_unit_cell(grid, a, b, c, alpha, beta, gamma) &
    bind(C, name="geGrid0_set_unit_cell")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      real(c_double), intent(in), value :: a, b, c, alpha, beta, gamma
    end subroutine

    !void geGrid0_mask_atom(geGrid0* grid, double x, double y, double z,
    !                       double radius);
    subroutine c_grid0_mask_atom(grid, x, y, z, radius) &
    bind(C, name="geGrid0_mask_atom")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      real(c_double), intent(in), value :: x, y, z, radius
    end subroutine

    !void geGrid0_apply_space_group(geGrid0* grid, int ccp4_num);
    subroutine c_grid0_apply_space_group(grid, ccp4_num) &
    bind(C, name="geGrid0_apply_space_group")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      integer(c_int), value :: ccp4_num
    end subroutine

    !int8_t* geGrid0_data(geGrid0* grid);
    type(c_ptr) function c_grid0_data(grid) bind(C, name="geGrid0_data")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
    end function

    !int8_t geGrid0_get_value(geGrid0* grid, int u, int v, int w);
    integer(c_int8_t) function c_grid0_get_value(grid, u, v, w) &
    bind(C, name="geGrid0_get_value")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      integer(c_int), value :: u, v, w
    end function

    !void geGrid0_free(geGrid0* grid);
    subroutine c_grid0_free(grid) bind(C, name="geGrid0_free")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
    end subroutine

    !void geGrid0_update_ccp4_header(geGrid0* grid, int mode, int update_stats);
    subroutine c_grid0_update_ccp4_header(grid, mode, update_stats) &
    bind(C, name="geGrid0_update_ccp4_header")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      integer(c_int), value :: mode
      logical(c_bool), value :: update_stats
    end subroutine

    !void geGrid0_write_ccp4_map(geGrid0* grid, const char* path);
    subroutine c_grid0_write_ccp4_map(grid, path) &
    bind(C, name="geGrid0_write_ccp4_map")
      use iso_c_binding
      type(c_ptr), intent(in), value :: grid
      character(kind=c_char), intent(in) :: path(*)
    end subroutine

    ! C runtime functions

    ! size_t strlen(const char *s);
    integer(c_size_t) function strlen(s) bind(C, name='strlen')
      use iso_c_binding
      type(c_ptr), intent(in), value :: s
    end function strlen
  end interface

  type spacegroup
    private
    type(c_ptr) :: ptr = c_null_ptr
  contains
    procedure :: num => spacegroup_number
    procedure :: hall => spacegroup_hall
    procedure :: hm => spacegroup_hm
    procedure :: short_name => spacegroup_short_name
    procedure :: operations => spacegroup_operations
  end type

  type groupops
    private
    type(c_ptr) :: ptr = c_null_ptr
  contains
    procedure :: order => groupops_order
    procedure :: get_op => groupops_get_op
    procedure :: free => groupops_free
  end type

  type op
    private
    type(c_ptr) :: ptr = c_null_ptr
  contains
    procedure :: triplet => op_triplet
    procedure :: free => op_free
  end type

  type grid0
    private
    type(c_ptr) :: ptr = c_null_ptr
  contains
    procedure :: init => grid0_init
    procedure :: set_unit_cell => grid0_set_unit_cell
    procedure :: mask_atom => grid0_mask_atom
    procedure :: apply_space_group => grid0_apply_space_group
    procedure :: data => grid0_data
    procedure :: get_value => grid0_get_value
    procedure :: free => grid0_free
    procedure :: update_ccp4_header => grid0_update_ccp4_header
    procedure :: write_ccp4_map => grid0_write_ccp4_map
  end type

contains

  type(spacegroup) function find_spacegroup_by_name(str)
    character(len=*), intent(in) :: str
    find_spacegroup_by_name%ptr = c_find_spacegroup_by_name(str//c_null_char)
  end function

  type(spacegroup) function find_spacegroup_by_number(num)
    integer, intent(in) :: num
    find_spacegroup_by_number%ptr = c_find_spacegroup_by_number(num)
  end function

  integer function spacegroup_number(this)
    class(spacegroup), intent(in) :: this
    spacegroup_number = c_spacegroup_number(this%ptr)
  end function

  function spacegroup_hall(this)
    character, pointer, dimension(:) :: spacegroup_hall
    class(spacegroup), intent(in) :: this
    type(c_ptr) :: c_string
    c_string = c_spacegroup_hall(this%ptr);
    call c_f_pointer(c_string, spacegroup_hall, [strlen(c_string)])
  end function

  function spacegroup_hm(this)
    character, pointer, dimension(:) :: spacegroup_hm
    class(spacegroup), intent(in) :: this
    type(c_ptr) :: c_string
    c_string = c_spacegroup_hm(this%ptr)
    call c_f_pointer(c_string, spacegroup_hm, [strlen(c_string)])
  end function

  character(16) function spacegroup_short_name(this)
    class(spacegroup), intent(in) :: this
    character(len=16, kind=c_char) :: string
    call c_spacegroup_short_name(this%ptr, string)
    spacegroup_short_name = string
  end function

  type(groupops) function spacegroup_operations(this)
    class(spacegroup), intent(in) :: this
    spacegroup_operations%ptr = c_spacegroup_operations(this%ptr)
  end function

  integer function groupops_order(this)
    class(groupops), intent(in) :: this
    groupops_order = c_groupops_order(this%ptr)
  end function

  type(op) function groupops_get_op(this, n)
    class(groupops), intent(in) :: this
    integer, intent(in), value :: n
    groupops_get_op%ptr = c_groupops_get_op(this%ptr, n-1)
  end function

  subroutine groupops_free(this)
    class(groupops) :: this
    call c_groupops_free(this%ptr)
    this%ptr = c_null_ptr
  end subroutine

  character(32) function op_triplet(this)
    class(op), intent(in) :: this
    character(len=32, kind=c_char) :: string
    call c_op_triplet(this%ptr, string)
    op_triplet = string
  end function

  subroutine op_free(this)
    class(op) :: this
    call c_op_free(this%ptr)
    this%ptr = c_null_ptr
  end subroutine


  subroutine grid0_init(this, nx, ny, nz)
    class(grid0) :: this
    integer, intent(in) :: nx, ny, nz
    this%ptr = c_grid0_init(nx, ny, nz)
  end subroutine

  subroutine grid0_set_unit_cell(this, a, b, c, alpha, beta, gamma)
    class(grid0) :: this
    double precision, intent(in) :: a, b, c, alpha, beta, gamma
    call c_grid0_set_unit_cell(this%ptr, a, b, c, alpha, beta, gamma)
  end subroutine

  subroutine grid0_mask_atom(this, x, y, z, radius)
    class(grid0) :: this
    double precision, intent(in) :: x, y, z, radius
    call c_grid0_mask_atom(this%ptr, x, y, z, radius)
  end subroutine

  subroutine grid0_apply_space_group(this, ccp4_num)
    class(grid0) :: this
    integer, intent(in) :: ccp4_num
    call c_grid0_apply_space_group(this%ptr, ccp4_num)
  end subroutine

  type(c_ptr) function grid0_data(this)
    class(grid0), intent(in) :: this
    grid0_data = c_grid0_data(this%ptr)
  end function

  integer function grid0_get_value(this, u, v, w)
    class(grid0), intent(in) :: this
    integer, intent(in) :: u, v, w
    grid0_get_value = c_grid0_get_value(this%ptr, u, v, w)
  end function

  subroutine grid0_free(this)
    class(grid0) :: this
    call c_grid0_free(this%ptr)
    this%ptr = c_null_ptr
  end subroutine

  subroutine grid0_update_ccp4_header(this, mode, update_stats)
    class(grid0) :: this
    integer, intent(in) :: mode
    logical, intent(in) :: update_stats
    call c_grid0_update_ccp4_header(this%ptr, mode, &
      logical(update_stats .eqv. .true., kind=c_bool))
  end subroutine

  subroutine grid0_write_ccp4_map(this, path)
    class(grid0), intent(in) :: this
    character(len=*), intent(in) :: path
    call c_grid0_write_ccp4_map(this%ptr, path//c_null_char)
  end subroutine

end module gemmi
! vim:sw=2:ts=2:et
