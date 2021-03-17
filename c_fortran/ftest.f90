
program main
  use gemmi_mod
  implicit none
  type(unitcell) :: cell
  type(mtz) :: mtz1
  character(len=256) :: arg
  integer :: i
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    mtz1 = read_mtz_file(trim(arg))
    print *, "mtz title: ", mtz1%get_title()
    print *, "reflections: ", mtz1%get_nreflections()
    print *, "space group: ", mtz1%get_spacegroup_name()
    cell = mtz1%get_cell(-1)
    print *, "unit cell:", cell%get_a(), cell%get_b(), cell%get_c(), &
             cell%get_alpha(), cell%get_beta(), cell%get_gamma()
    !call mtz1%set_title('new title')
    !call mtz1%write_to_file("out.mtz")
  end do
end program main

