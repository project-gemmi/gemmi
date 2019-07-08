program main
  use gemmi, only : mask
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  type(mask) :: grid
  type(c_ptr) :: values

  call grid%init(180, 180, 480)
  call grid%set_unit_cell(80.518d0, 80.518d0, 200.952d0, 90d0, 90d0, 120d0)

  ! do i=1,natoms
  !   grid%mask_atom(atom_x, atom_y, atom_z, 3.0)
  ! end do
  call grid%mask_atom(30d0, 40d0, 55d0, 0.5d0)

  call grid%apply_space_group(180)

  values = grid%data()

  print *, grid%get_value(0, 0, 0), grid%get_value(119, 103, 131)

  ! for testing
  call grid%update_ccp4_header(0, .true.) ! args: mode (0 or 2), update_stats
  call grid%write_ccp4_map("test.map")

  call grid%free()
end program main
