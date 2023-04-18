
program main
  use gemmi_mod
  implicit none
  type(spacegroup) :: sg
  type(groupops) :: ops
  type(op) :: op_
  character(len=256) :: arg
  integer :: i, j
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    sg = find_spacegroup_by_name(trim(arg))
    if (.not. sg%associated()) then
      write(*, *) "n/a"
      cycle
    end if
    ops = sg%operations()
    write(*, *) "space group", sg%get_number(), "  ", sg%xhm(), &
                " short: ", sg%short_name()
    do j = 1, ops%order()
      op_ = ops%get_op(j - 1)
      write(*, *) j, ": ", op_%triplet()
      !call op_%free()
    end do
    !call ops%free()
  end do
end program main
