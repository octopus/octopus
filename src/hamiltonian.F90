#include "config.h"

module hamiltonian
use global
use spline
use fft
use units
use system
use hartree
use xc

implicit none

type hamiltonian_type
  integer :: ispin ! How to handle spin (duplicated in states_type)
  integer :: np  ! number of points (duplicated in mesh)

  integer :: soc ! spin-orbit coupling?

  integer :: vpsl_space           ! How should the local potential be calculated
  real(r8), pointer :: Vpsl(:)    ! the external potential
  real(r8), pointer :: Vhartree(:)! the external potential
  real(r8), pointer :: Vxc(:,:)   ! xc potential
  real(r8), pointer :: rho_core(:)! core charge for nl core corrections

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec, epot

  ! System under the independent particle approximation, or not.
  logical :: ip_app

  ! should we include the classical point charges
  logical :: classic_pot
  real(r8), pointer :: Vclassic(:)! potential created by classic point charges

  ! hartree potential structure
  type(hartree_type) :: hart
  type(xc_type) :: xc

end type hamiltonian_type

contains

subroutine hamiltonian_init(h, sys)
  type(hamiltonian_type), intent(out) :: h
  type(system_type), intent(inout) :: sys

  sub_name = 'hamiltonian_init'; call push_sub()

  ! Duplicate this two variables
  h%ispin = sys%st%ispin
  h%np  = sys%m%np

  ! allocate potentials and density of the cores
  allocate(h%Vpsl(h%np), h%Vhartree(h%np), h%Vxc(h%np, h%ispin), h%rho_core(h%np))
  h%Vhartree = 0._r8
  h%Vxc = 0._r8

  if(sys%ncatoms > 0) then
    h%classic_pot = fdf_boolean("ClassicPotential", .false.)
    if(h%classic_pot) allocate(h%Vclassic(h%np))
  end if

  ! should we calculate the local pseudopotentials in Fourier space?
  h%vpsl_space = fdf_integer('LocalPotentialSpace', 0)
  if(h%vpsl_space < 0 .or. h%vpsl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vpsl_space, &
         "' is not a valid LocalPotentialSpace"
    message(2) = '(0 <= LocalPotentialSpace <=1)'
    call write_fatal(2)
  end if

  if(h%vpsl_space == 1) then
    call specie_fourier_init(sys%nspecies, sys%specie, sys%m)
  end if

  ! Should we treat the particles as independent?
  h%ip_app = fdf_boolean("NonInteractingElectrons", .false.)
  if(h%ip_app) then
    message(1) = 'Info: Treating the electrons as non-interacting'
    call write_info(1)
  else
    ! initilize hartree and xc modules
    call hartree_init(h%hart, sys%m)
    call xc_init(h%xc, sys%m, h%ispin)
  end if

  call pop_sub()
end subroutine hamiltonian_init

subroutine hamiltonian_end(h)
  type(hamiltonian_type) :: h

  sub_name = 'hamiltonian_end'; call push_sub()

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhartree, h%Vxc)
    nullify(h%Vpsl, h%Vhartree, h%Vxc)
  end if

  if(h%classic_pot .and. associated(h%Vclassic)) then
    deallocate(h%Vclassic); nullify(h%Vclassic)
  end if

  if(.not.h%ip_app) then
    call hartree_end(h%hart)
    call xc_end(h%xc)
  end if

  call pop_sub()
end subroutine hamiltonian_end

subroutine hamiltonian_setup(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  real(r8), allocatable :: v_aux(:,:)
  integer :: i

  h%epot = 0._r8 ! The energy coming from the potentials

  if(.not.h%ip_app) then
    call hartree_solve(h%hart, sys%m, h%Vhartree, sys%st%rho)
    do i = 1, sys%st%nspin
      h%epot = h%epot - 0.5_r8*dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vhartree)
    end do
    
    allocate(v_aux(h%np, sys%st%nspin))
    call R_FUNC(xc_pot)(h%xc, sys%m, sys%st, h%hart, h%rho_core, &
         h%Vxc, v_aux, h%ex, h%ec)
    h%Vxc = h%Vxc + v_aux
    do i = 1, sys%st%nspin
      h%epot = h%epot - dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vxc(:, i))
    end do
    deallocate(v_aux)
    
  end if


end subroutine hamiltonian_setup

! This subroutine calculates the total energy of the system. Basically, it
! adds up the KS eigenvalues, and then it substracts the whatever double
! counts exist (see TDDFT theory for details).
subroutine hamiltonian_energy(h, sys, iunit, reduce)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys
  integer, intent(in) :: iunit
  logical, intent(in), optional :: reduce

  integer :: ik
  real(r8) :: s, e
#ifdef HAVE_MPI
  integer :: ierr
#endif 

  sub_name = 'hamiltonian_energy'; call push_sub()

  e = sum(sys%st%occ(sys%st%st_start:sys%st%st_end, :)*&
       sys%st%eigenval(sys%st%st_start:sys%st%st_end, :))
#ifdef HAVE_MPI
  if(present(reduce) .and. reduce) then
    call MPI_ALLREDUCE(e, s, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    e = s
  end if
#endif
  
  h%etot = e + sys%eii + h%epot + h%ex + h%ec
  
  if(iunit>0) then
    write(message(1), '(6x,a, f15.8)')'Ion-ion     = ', sys%eii / units_out%energy%factor
    write(message(2), '(6x,a, f15.8)')'Eigenvalues = ', e       / units_out%energy%factor
    write(message(3), '(6x,a, f15.8)')'Potentials  = ', h%epot  / units_out%energy%factor
    write(message(4), '(6x,a, f15.8)')'Exchange    = ', h%ex    / units_out%energy%factor
    write(message(5), '(6x,a, f15.8)')'Correlation = ', h%ec    / units_out%energy%factor
    write(message(6), '(6x,a, f15.8)')'Total       = ', h%etot  / units_out%energy%factor
    call write_info(6, iunit)
  end if

  call pop_sub()
  return
end subroutine hamiltonian_energy

#include "h_external_pot.F90"

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_inc.F90"

end module hamiltonian
