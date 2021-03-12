
program main
  use gemmi_mod
  implicit none
  type(unitcell) :: cell
  type(mtz) :: mtz1
  character(len=256) :: arg
  integer :: i, j, num
  !sg = find_spacegroup_by_number(5)
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    mtz1 = read_mtz_file(trim(arg))
    cell = mtz1%get_cell(-1)
    print *, "mtz unit cell:", cell%get_a(), cell%get_b(), cell%get_c(), &
             cell%get_alpha(), cell%get_beta(), cell%get_gamma()
    !write(*, *) "mtz unit cell: ", get_cell(-1)mtz1%spacegroup%xhm()
    !sg = find_spacegroup_by_name(trim(arg))
  end do
end program main

!program main
!  use gemmi
!  implicit none
!  type(spacegroup) :: sg
!  type(groupops) :: ops
!  type(op) :: op_
!  character(len=256) :: arg
!  integer :: i, j, num
!  do i = 1, command_argument_count()
!    call get_command_argument(i, arg)
!    sg = find_spacegroup_by_name(trim(arg))
!    num = sg%num()
!    if (num .eq. 0) then
!      write(*, *) "n/a"
!      cycle
!    end if
!    ops = sg%operations()
!    write(*, *) "space group", sg%num(), "  ", sg%hm(), &
!                " short: ", sg%short_name()
!    do j = 1, ops%order()
!      op_ = ops%get_op(j)
!      write(*, *) j, ": ", op_%triplet()
!      call op_%free()
!    end do
!    call ops%free()
!  end do
!end program main
