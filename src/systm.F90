module system
use specie
use atom

implicit none

type system_type
  integer :: natoms
  type(atom_type), pointer :: atom(:)

  integer :: nspecies
  type(specie_type), pointer :: specie(:)
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  sub_name = 'system_init'; call push_sub()

  s%nspecies = specie_init(s%specie)
  s%natoms = atom_init(s%atom, s%nspecies, s%specie)

  call atom_end(s%natoms, s%atom)
  call specie_end(s%nspecies, s%specie)

  call pop_sub()
end subroutine system_init

end module system
