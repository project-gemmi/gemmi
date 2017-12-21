program main
  use gemmi
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  type(grid0) :: grid
  type(c_ptr) :: values

  call grid%init(180, 180, 480)
  call grid%set_unit_cell(80.518d0, 80.518d0, 200.952d0, 90d0, 90d0, 120d0)

  ! do i=1,natoms
  !   grid%mask_atom(atom_x, atom_y, atom_z, 3.0)
  ! end do
  call grid%mask_atom(30d0, 40d0, 55d0, 3.0d0)

  call grid%apply_space_group(180)

  values = grid%data()

  ! for testing
  call grid%prepare_ccp4_header(0) ! arg is mode: 0 or 2
  call grid%write_ccp4_map("test.map")

  call grid%free()
end program main
! vim:sw=2:ts=2:et
