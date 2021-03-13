
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
    cell = mtz1%get_cell(-1)
    print *, "mtz unit cell:", cell%get_a(), cell%get_b(), cell%get_c(), &
             cell%get_alpha(), cell%get_beta(), cell%get_gamma()
  end do
end program main

