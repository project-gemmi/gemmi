
module gemmi
  use iso_c_binding
  implicit none
  interface
    ! const cSpaceGroup* find_spacegroup_by_name(const char* name);
    type(c_ptr) function C_find_spacegroup_by_name(name) &
    bind(C, name="find_spacegroup_by_name")
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*)
    end function

    ! const cSpaceGroup* find_spacegroup_by_number(int n);
    type(c_ptr) function C_find_spacegroup_by_number(n) &
    bind(C, name="find_spacegroup_by_number")
      use iso_c_binding
      integer(c_int), value :: n
    end function

    ! int SpaceGroup_number(const cSpaceGroup* sg);
    integer(c_int) function C_SpaceGroup_number(sg) &
    bind(C, name="SpaceGroup_number")
      use iso_c_binding
      type(c_ptr), intent(in), value :: sg
    end function

    ! const char* SpaceGroup_hm(const cSpaceGroup* sg);

    ! const char* SpaceGroup_hall(const cSpaceGroup* sg);

    ! void SpaceGroup_short_name(const cSpaceGroup* sg, char* dest);

    ! cGroupOps* SpaceGroup_operations(const cSpaceGroup* sg);

    ! int GroupOps_order(cGroupOps* ops);

    ! void GroupOps_free(cGroupOps* ops);
    subroutine C_GroupOps_free(ops) bind(C, name="GroupOps_free")
      use iso_c_binding
      type(c_ptr), intent(in) :: ops
    end subroutine
  end interface


end module
