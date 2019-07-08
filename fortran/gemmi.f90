
module gemmi
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: spacegroup, groupops, op, mask, find_spacegroup_by_name, &
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

    !geMask* geMask_init(int nx, int ny, int nz);
    type(c_ptr) function c_mask_init(nx, ny, nz) &
    bind(C, name="geMask_init")
      use iso_c_binding
      integer(c_int), intent(in), value :: nx, ny, nz
    end function

    !void geMask_set_unit_cell(geMask* mask, double a, double b, double c,
    !                           double alpha, double beta, double gamma);
    subroutine c_mask_set_unit_cell(mask, a, b, c, alpha, beta, gamma) &
    bind(C, name="geMask_set_unit_cell")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
      real(c_double), intent(in), value :: a, b, c, alpha, beta, gamma
    end subroutine

    !void geMask_mask_atom(geMask* mask, double x, double y, double z,
    !                       double radius);
    subroutine c_mask_mask_atom(mask, x, y, z, radius) &
    bind(C, name="geMask_mask_atom")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
      real(c_double), intent(in), value :: x, y, z, radius
    end subroutine

    !void geMask_apply_space_group(geMask* mask, int ccp4_num);
    subroutine c_mask_apply_space_group(mask, ccp4_num) &
    bind(C, name="geMask_apply_space_group")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
      integer(c_int), value :: ccp4_num
    end subroutine

    !int8_t* geMask_data(geMask* mask);
    type(c_ptr) function c_mask_data(mask) bind(C, name="geMask_data")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
    end function

    !int8_t geMask_get_value(geMask* mask, int u, int v, int w);
    integer(c_int8_t) function c_mask_get_value(mask, u, v, w) &
    bind(C, name="geMask_get_value")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
      integer(c_int), value :: u, v, w
    end function

    !void geMask_free(geMask* mask);
    subroutine c_mask_free(mask) bind(C, name="geMask_free")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
    end subroutine

    !void geMask_update_ccp4_header(geMask* mask, int mode, int update_stats);
    subroutine c_mask_update_ccp4_header(mask, mode, update_stats) &
    bind(C, name="geMask_update_ccp4_header")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
      integer(c_int), value :: mode
      logical(c_bool), value :: update_stats
    end subroutine

    !void geMask_write_ccp4_map(geMask* mask, const char* path);
    subroutine c_mask_write_ccp4_map(mask, path) &
    bind(C, name="geMask_write_ccp4_map")
      use iso_c_binding
      type(c_ptr), intent(in), value :: mask
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

  type mask
    private
    type(c_ptr) :: ptr = c_null_ptr
  contains
    procedure :: init => mask_init
    procedure :: set_unit_cell => mask_set_unit_cell
    procedure :: mask_atom => mask_mask_atom
    procedure :: apply_space_group => mask_apply_space_group
    procedure :: data => mask_data
    procedure :: get_value => mask_get_value
    procedure :: free => mask_free
    procedure :: update_ccp4_header => mask_update_ccp4_header
    procedure :: write_ccp4_map => mask_write_ccp4_map
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


  subroutine mask_init(this, nx, ny, nz)
    class(mask) :: this
    integer, intent(in) :: nx, ny, nz
    this%ptr = c_mask_init(nx, ny, nz)
  end subroutine

  subroutine mask_set_unit_cell(this, a, b, c, alpha, beta, gamma)
    class(mask) :: this
    double precision, intent(in) :: a, b, c, alpha, beta, gamma
    call c_mask_set_unit_cell(this%ptr, a, b, c, alpha, beta, gamma)
  end subroutine

  subroutine mask_mask_atom(this, x, y, z, radius)
    class(mask) :: this
    double precision, intent(in) :: x, y, z, radius
    call c_mask_mask_atom(this%ptr, x, y, z, radius)
  end subroutine

  subroutine mask_apply_space_group(this, ccp4_num)
    class(mask) :: this
    integer, intent(in) :: ccp4_num
    call c_mask_apply_space_group(this%ptr, ccp4_num)
  end subroutine

  type(c_ptr) function mask_data(this)
    class(mask), intent(in) :: this
    mask_data = c_mask_data(this%ptr)
  end function

  integer function mask_get_value(this, u, v, w)
    class(mask), intent(in) :: this
    integer, intent(in) :: u, v, w
    mask_get_value = c_mask_get_value(this%ptr, u, v, w)
  end function

  subroutine mask_free(this)
    class(mask) :: this
    call c_mask_free(this%ptr)
    this%ptr = c_null_ptr
  end subroutine

  subroutine mask_update_ccp4_header(this, mode, update_stats)
    class(mask) :: this
    integer, intent(in) :: mode
    logical, intent(in) :: update_stats
    call c_mask_update_ccp4_header(this%ptr, mode, &
                                    logical(update_stats, kind=c_bool))
  end subroutine

  subroutine mask_write_ccp4_map(this, path)
    class(mask), intent(in) :: this
    character(len=*), intent(in) :: path
    call c_mask_write_ccp4_map(this%ptr, path//c_null_char)
  end subroutine

end module gemmi
