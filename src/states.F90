#include "config.h"

module states
use global
use fdf
use math
use mesh

implicit none

type states_type
  integer :: dim ! Dimension of the state
  integer :: nst ! Number of states in each irreducible subspace
  integer :: nik ! Number of irreducible subspaces
  integer :: ispin ! spin mode

  ! pointers to the wavefunctions
  real(r8), pointer :: dpsi(:,:,:,:)
  complex(r8), pointer :: zpsi(:,:,:,:)

  ! the density (after all we are doing DFT :)
  real(r8), pointer :: rho(:,:)

  real(r8), pointer :: eigenval(:,:)
  real(r8), pointer :: occ(:,:)

  integer :: st_start, st_end ! needed for some parallel parts
end type states_type

contains

subroutine states_init(st, m, val_charge)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  real(r8), intent(in) :: val_charge

  real(r8) :: excess_charge
  integer :: nempty

  sub_name = 'states_init'; call push_sub()

  st%ispin = fdf_integer('SpinComponents', 1)
  if (st%ispin /= 1 .and. st%ispin /= 2 .and. st%ispin /= 4) then
    write(message(1),'(a,i4,a)') "Input: '", st%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 4)'
    call write_fatal(2)
  end if

#ifdef PERIODIC_1D
  st%nik = fdf_integer('NumberKPoints', 1)
  if (st%nik < 1) then
    write(message(1),'(a,i4,a)') "Input: '", st%nik,"' is not a valid NumberKPoints"
    message(2) = '(NumberKPoints >= 1)'
    call write_fatal(2)
  end if
#else
  st%nik = 1
#endif
  
  excess_charge = fdf_double('ExcessCharge', 0.0_r8)

  nempty = fdf_integer('EmptyStates', 0)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid EmptyStates"
    message(2) = '(0 <= EmptyStates)'
    call write_fatal(2)
  end if
  
  st%nst = -int((val_charge + excess_charge)/2._r8)
  if(st%nst*2._r8 < -(val_charge + excess_charge)) &
       st%nst = st%nst + 1
  st%nst = st%nst + nempty

  select case(st%ispin)
  case(1)
    st%dim = 1
  case(2)
    st%dim = 1
    st%nik = st%nik*2
  case(4)
    st%dim = 2
    st%nik = st%nik*2
  end select

  ! we now allocate some arrays
  allocate(st%rho(m%np, st%ispin), &
       st%occ(st%nst, st%nik), st%eigenval(st%nst, st%nik))

  st%st_start = 1; st%st_end = st%nst

  call pop_sub()
end subroutine states_init

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  sub_name = 'states_end'; call push_sub()

  if(associated(st%rho)) then
    deallocate(st%rho); nullify(st%rho)
    deallocate(st%occ); nullify(st%occ)
    deallocate(st%eigenval); nullify(st%eigenval)
  end if

  if(associated(st%dpsi)) then
    deallocate(st%dpsi); nullify(st%dpsi)
  end if

  if(associated(st%zpsi)) then
    deallocate(st%zpsi); nullify(st%zpsi)
  end if

  call pop_sub()
end subroutine states_end

! generate a hydrogen s-wavefunction around a random point
subroutine states_generate_random(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m

  integer, save :: iseed = 123
  integer :: ist, ik, id, i
  real(r8) :: a(3), rnd, f, r

  sub_name = 'states_generate_random'; call push_sub()

  st%R_FUNC(psi)(0, :, :, :) = 0.0_r8
  do ik = 1, st%nik
    do ist = 1, st%nst
      do id = 1, st%dim
        call quickrnd(iseed, rnd)
        a(1) = 2.0_r8*(2*rnd - 1)
        call quickrnd(iseed, rnd)
        a(2) = 2.0_r8*(2*rnd - 1)
        call quickrnd(iseed, rnd)
        a(3) = 2.0_r8*(2*rnd - 1)
#ifdef COMPLEX_WFNS
        call quickrnd(iseed, rnd)
        f = M_PI*(2*rnd - 1)
#endif
        do i = 1, m%np
          call mesh_r(m, i, r, a=a)
#ifdef COMPLEX_WFNS
          st%R_FUNC(psi)(i, id, ist, ik) = exp(-0.5_r8*r*r + M_zI*f)
#else
          st%R_FUNC(psi)(i, id, ist, ik) = exp(-0.5_r8*r*r)
#endif
        end do
      end do
    end do

    call R_FUNC(states_gram_schmidt)(st%nst, m, st%dim, st%R_FUNC(psi)(:,:,:,ik))
  end do
  st%eigenval = 0._r8

  call pop_sub()
end subroutine states_generate_random

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
