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
    call geom_write_xyz(s)
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

subroutine kinetic_energy(sys)
  type(system_type), intent(inout) :: sys

  integer :: i

  sub_name = 'kinetic_energy'; call push_sub()  

  sys%kinetic_energy = 0.0_r8
  do i = 1, sys%natoms
     sys%kinetic_energy = sys%kinetic_energy + 0.5_r8*sys%atom(i)%spec%weight* &
         sum(sys%atom(i)%v(:)**2)
  enddo

  call pop_sub()
end subroutine kinetic_energy

subroutine geom_write_xyz(sys)
  type(system_type), intent(IN) :: sys

  integer i, iunit
  
  ! xyz format, for easy plot in rasmol
#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    call io_assign(iunit)
    open(iunit, file=trim(sys%sysname)//'.xyz', status='unknown')
    write(iunit, '(i4)') sys%natoms
    write(iunit, '(1x)')
    do i = 1, sys%natoms
      write(iunit, '(6x,a,2x,3f12.6)') &
           sys%atom(i)%spec%label, sys%atom(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)

    if(sys%ncatoms > 0) then
      call io_assign(iunit)
      open(iunit, file=trim(sys%sysname)//'_classical.xyz', status='unknown')
      write(iunit, '(i4)') sys%ncatoms
      write(iunit, '(1x)')
      do i = 1, sys%ncatoms
        write(iunit, '(6x,a1,2x,3f12.6,a,f12.6)') &
             sys%catom(i)%label(1:1), sys%catom(i)%x(:)/units_out%length%factor, &
             " # ", sys%catom(i)%charge
      end do
      call io_close(iunit)
    end if

#ifdef HAVE_MPI
  end if
#endif
  
  return
end subroutine geom_write_xyz


end module system
