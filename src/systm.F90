module system
use io
use specie
use atom
use mesh
use states

implicit none

type system_type
  character(len=20) :: sysname ! the name of the system we are running

  integer :: natoms
  type(atom_type), pointer :: atom(:)
  integer :: ncatoms  ! For QM+MM calculations
  type(atom_classical_type), pointer :: catom(:)

  real(r8) :: eii, kinetic_energy ! the ion-ion energy

  integer :: nspecies
  type(specie_type), pointer :: specie(:)

  type(mesh_type) :: m
  type(states_type), pointer :: st
end type system_type

contains

subroutine system_init(s)
  type(system_type), intent(out) :: s

  integer :: i
  real(r8) :: val_charge
  logical :: l

  sub_name = 'system_init'; call push_sub()

  call oct_parse_str('SystemName', 'system', s%sysname)
  s%nspecies = specie_init(s%specie)
  call atom_init(s%natoms, s%atom, s%ncatoms, s%catom, s%nspecies, s%specie)
  call oct_parse_logical("OutputCoordinates", .false., l)
  if(l) then
    call geom_write_xyz(s%sysname, s%natoms, s%atom, s%ncatoms, s%catom)
  end if

  call mesh_init(s%m, s%natoms, s%atom)

  !  find total charge of the system
  val_charge = 0
  do i = 1, s%natoms
    val_charge = val_charge - s%atom(i)%spec%Z_val
  enddo

  allocate(s%st)
  call states_init(s%st, s%m, val_charge)

  call pop_sub()
end subroutine system_init

subroutine system_end(s)
  type(system_type), intent(inout) :: s

  sub_name = 'system_end'; call push_sub()

  if(associated(s%st)) then
    call states_end(s%st)
    deallocate(s%st); nullify(s%st)
  end if
  call mesh_end(s%m)
  call atom_end(s%natoms, s%atom, s%ncatoms, s%catom)
  call specie_end(s%nspecies, s%specie)

  call pop_sub()
end subroutine system_end

end module system
