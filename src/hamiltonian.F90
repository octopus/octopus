#include "config.h"

module hamiltonian
use global
use spline
use fft
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
  real(r8), pointer :: Vhxc(:,:)  ! Hartre+xc potential
  real(r8), pointer :: rho_core(:)! core charge for nl core corrections

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec, epot

  ! System under the independent particle approximation, or not.
  logical :: ip_app

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
  allocate(h%Vpsl(h%np), h%Vhxc(h%np, h%ispin), h%rho_core(h%np))

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
    deallocate(h%Vpsl, h%Vhxc)
    nullify(h%Vpsl, h%Vhxc)
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

  real(r8), allocatable :: v_aux1(:,:), v_aux2(:,:)
  integer :: i, is

  allocate(v_aux1(h%np, h%ispin), v_aux2(h%np, h%ispin))

  is = min(H%ispin, 2)
  h%epot = 0._r8 ! The energy coming from the potentials

  do i = 1, is
    h%Vhxc(:, i) = H%vpsl(:)
  end do
  if(is > 2) h%Vhxc(:, 2:) = 0._r8

  if(.not.h%ip_app) then
    call hartree_solve(h%hart, sys%m, v_aux1(:,1), sys%st%rho)
    do i = 1, is
      h%Vhxc(:, i) = h%Vhxc(:, i) + v_aux1(:, 1)
      h%epot = h%epot - 0.5_r8*dmesh_dp(sys%m, sys%st%rho(:, is), v_aux1(:, 1))
    end do
    
    call R_FUNC(xc_pot)(h%xc, sys%m, sys%st, h%hart, h%rho_core, &
         v_aux1, v_aux2, h%ex, h%ec)
    h%Vhxc = h%Vhxc + v_aux1 + v_aux2
    do i = 1, is
      h%epot = h%epot - dmesh_dp(sys%m, sys%st%rho(:, is), v_aux1(:, is))
      h%epot = h%epot - dmesh_dp(sys%m, sys%st%rho(:, is), v_aux2(:, is))
    end do
    
  end if

  deallocate(v_aux1, v_aux2)

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
#if defined(HAVE_MPI) && defined(MPI_TD)
  integer :: ierr
#endif 

  sub_name = 'tenergy'; call push_sub()

  e = sum(sys%st%occ(sys%st%st_start:sys%st%st_end, :)*&
       sys%st%eigenval(sys%st%st_start:sys%st%st_end, :))
#if defined(HAVE_MPI) && defined(MPI_TD)
  if(present(reduce) .and. reduce) then
    call MPI_ALLREDUCE(e, s, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    e = s
  end if
#endif
  
  h%etot = h%etot + sys%eii + h%epot + h%ex + h%ec
  
  if(iunit>0) then
    write(message(1), '(6x,a, f15.8)')'Ion-ion     = ', sys%eii
    write(message(2), '(6x,a, f15.8)')'Eigenvalues = ', e
    write(message(3), '(6x,a, f15.8)')'Potentials  = ', h%epot
    write(message(4), '(6x,a, f15.8)')'Exchange    = ', h%ex
    write(message(5), '(6x,a, f15.8)')'Correlation = ', h%ec
    write(message(6), '(6x,a, f15.8)')'Total       = ', h%etot
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
