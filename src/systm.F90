module system
use specie
use atom
use mesh
use states

implicit none

type system_type
  integer :: natoms
  type(atom_type), pointer :: atom(:)

  integer :: nspecies
  type(specie_type), pointer :: specie(:)

  type(mesh_type) :: m
  type(states_type) :: st
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  integer :: i
  real(r8) :: val_charge

  sub_name = 'system_init'; call push_sub()

  s%nspecies = specie_init(s%specie)
  s%natoms = atom_init(s%atom, s%nspecies, s%specie)
  call mesh_init(s%m, s%natoms, s%atom)

  !  find total charge of the system
  val_charge = 0
  do i = 1, s%natoms
    val_charge = val_charge - s%atom(i)%spec%Z_val
  enddo

  call states_init(s%st, val_charge)

  call pop_sub()
end subroutine system_init

subroutine system_end(s)
  type(system_type), intent(inout) :: s

  sub_name = 'system_end'; call push_sub()

  call states_end(s%st)
  call mesh_end(s%m)
  call atom_end(s%natoms, s%atom)
  call specie_end(s%nspecies, s%specie)

  call pop_sub()
end subroutine system_end

end module system
