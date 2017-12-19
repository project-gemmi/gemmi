
module gemmi
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: spacegroup, groupops, find_spacegroup_by_name, &
            find_spacegroup_by_number

  interface
    ! const geSpaceGroup* find_spacegroup_by_name(const char* name);
    type(c_ptr) function c_find_spacegroup_by_name(name) &
    bind(C, name="find_spacegroup_by_name")
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*)
    end function

    ! const geSpaceGroup* find_spacegroup_by_number(int n);
    type(c_ptr) function c_find_spacegroup_by_number(n) &
    bind(C, name="find_spacegroup_by_number")
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

    ! void geGroupOps_free(geGroupOps* ops);
    subroutine c_groupops_free(ops) bind(C, name="geGroupOps_free")
      use iso_c_binding
      type(c_ptr), intent(in), value :: ops
    end subroutine

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
    procedure :: free => groupops_free
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
    ! TODO fix it
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

  subroutine groupops_free(this)
    implicit none
    class(groupops) :: this
    call c_groupops_free(this%ptr)
    this%ptr = c_null_ptr
  end subroutine

end module gemmi
! vim:sw=2:ts=2:et
