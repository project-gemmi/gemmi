
program main
  use gemmi
  implicit none
  type(spacegroup) :: sg
  type(groupops) :: ops
  character(len=256) :: arg
  integer :: i, num
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    sg = find_spacegroup_by_name(trim(arg))
    num = sg%num()
    if (num .eq. 0) then
      write(*, *) "n/a"
      cycle
    end if
    ops = sg%operations()
    write(*, *) sg%num(), "  ", sg%hm(), "  order:", ops%order()
    call ops%free()
  end do
end program main

! vim:sw=2:ts=2:et
