module system
use specie
use atom
use mesh
use states

implicit none

type system_type
  character(len=20) :: sysname ! the name of the system we are running

  integer :: natoms
  type(atom_type), pointer :: atom(:)
  real(r8) :: eii ! the ion-ion energy

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

  s%sysname = fdf_string('SystemName', 'system')
  s%nspecies = specie_init(s%specie)
  s%natoms = atom_init(s%atom, s%nspecies, s%specie)
  call mesh_init(s%m, s%natoms, s%atom)

  !  find total charge of the system
  val_charge = 0
  do i = 1, s%natoms
    val_charge = val_charge - s%atom(i)%spec%Z_val
  enddo

  call states_init(s%st, s%m, val_charge)

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

subroutine ion_ion_energy(sys)
  type(system_type), intent(inout) :: sys
  
  real(r8) :: r
  integer :: i, j

  sub_name = 'ion_ion_energy'; call push_sub()

  ! calculate the ion-ion energy
  sys%eii = 0.0_r8

  do i = 2, sys%natoms
    do j = 1, i - 1
      r = sqrt(sum((sys%atom(i)%x - sys%atom(j)%x)**2))
      sys%eii = sys%eii + &
           sys%atom(i)%spec%Z_val*sys%atom(j)%spec%Z_val/r
    end do
  end do

  call pop_sub()
end subroutine ion_ion_energy

end module system
