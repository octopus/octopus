!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "config_F90.h"

module mix
use global
use mesh
use states

implicit none

private
public :: mix_type, mix_init, mix_end, mix_dens, dcalcdens, zcalcdens

integer, parameter :: LINEAR   = 0, &
                      ANDERSON = 1, &
                      BROYDEN  = 2

! Anderson Mixing
!!$ integer, parameter :: A_IN = 1, A_OUT = 2

type mix_type
  integer  :: type_of_mixing                
  integer  :: norbitals               !  Number of orbitals
  real(r8) :: weight                  !  Weight of each component to perform 
                                      !  the integrals.

  real(r8) :: alpha                   !  rho_new = (1-a)*rho_in + a*rho_out

  ! Anderson Mixing
!!$  real(r8) :: anderson_mixing_alpha
!!$  real(r8), pointer :: densities(:,:,:,:)
!!$  integer :: anderson_number

  ! Broyden Mixing
  integer :: broyden_number
  real(r8), pointer :: df(:,:,:), dv(:,:,:)
end type mix_type

contains

! Initialization...
subroutine mix_init(smix, m, st)
  type(mix_type), intent(out) :: smix
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st

  sub_name = 'mix_init'; call push_sub()
  
  ! check input parameters
  call oct_parse_int(C_string("TypeOfMixing"), 0, smix%type_of_mixing)
  if(smix%type_of_mixing == 1) then
    message(1) = 'Anderson mixing (TypeOfMixing = 1) currently not implemented)'
    call write_fatal(1)
  end if
  if(smix%type_of_mixing < 0 .or. smix%type_of_mixing > 2) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  call oct_parse_double(C_string("Mixing"), 0.3_r8, smix%alpha);
  if(smix%alpha <= 0.0_r8 .or. smix%alpha > 1.0_r8) then
    write(message(1), '(a, f14.6,a)') "Input: '", smix%alpha, &
         "' is not a valid Mixing"
    message(2) = '(0 < Mixing <= 1)'
    call write_fatal(2)
  end if

  if(smix%type_of_mixing == 2) then
    call oct_parse_int(C_string("BroydenNumber"), 3, smix%broyden_number)
    if(smix%broyden_number <= 1 .or. smix%broyden_number > 5) then
      write(message(1), '(a, i4,a)') "Input: '", smix%broyden_number, &
           "' is not a valid BroydenNumber"
      message(2) = '(1 < BroydenNumber <= 5)'
      call write_fatal(2)
    end if
  end if

  if(smix%type_of_mixing == BROYDEN) then
    allocate(smix%df(m%np, st%nspin, smix%broyden_number + 1))
    allocate(smix%dv(m%np, st%nspin, smix%broyden_number + 1))
    smix%df = 0._r8
    smix%dv = 0._r8
  end if

  call pop_sub()
  return
end subroutine mix_init

subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix

  sub_name = 'mix_end'; call push_sub()

  if(smix%type_of_mixing == BROYDEN) then
    if(associated(smix%df)) then
      deallocate(smix%df); nullify(smix%df);
      deallocate(smix%dv); nullify(smix%dv);
    end if
  end if

  call pop_sub()
end subroutine mix_end

subroutine mix_dens(smix, iter, st, m, dist)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)      :: iter
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  real(r8), intent(out)    :: dist
  
  real(r8), allocatable :: rhoout(:,:), dummy(:)
  integer :: is, errorflag
    
  sub_name = 'mix_dens'; call push_sub()

  allocate(rhoout(m%np, st%nspin))
  call R_FUNC(calcdens)(st, m%np, rhoout)

  dist = 0._r8
  allocate(dummy(m%np))
  do is = 1, st%nspin
    dummy = (st%rho(:,is) - rhoout(:,is))**2
    dist = dist + dmesh_integrate(m, dummy)
  end do
  dist = sqrt(dist)
  deallocate(dummy)
    
  select case(smix%type_of_mixing)
  case(LINEAR)
    call mix_linear(smix, m%np, st%nspin, st%rho, rhoout)
!!$  case(ANDERSON)
!!$      call anderson_mix(rho,rhoout, iter, errorflag)
  case(BROYDEN)
    call mix_broyden(st%nspin, m%np, smix, st%rho, rhoout, iter, errorflag)
  end select
    
  deallocate(rhoout)

  call pop_sub()
end subroutine mix_dens

! Performs the linear mixing...
subroutine mix_linear(smix, np, nspin, rho, rhoout)
  type(mix_type), intent(IN) :: smix
  integer, intent(in) :: np, nspin
  real(r8), intent(inout) :: rho(np, nspin)
  real(r8), intent(IN) :: rhoout(np, nspin)

  call dscal(np*nspin, 1.0_r8 - smix%alpha, rho, 1)
  call daxpy(np*nspin, smix%alpha, rhoout, 1, rho, 1)

  return
end subroutine mix_linear

! Broyden mixing...
subroutine mix_broyden(nspin, np, smix, vin, vout, iter, errorflag)
  integer, intent(in) :: nspin
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: iter, np
  integer, intent(out)    :: errorflag
  real(r8), intent(inout) :: vout(np, nspin), & ! the new density; 
       vin(np, nspin)     ! the old H+xc density

  integer, parameter :: maxter = 5

  integer  :: is, i, j, iwork(maxter), info, iter_used, ipos
  real(r8) :: beta(maxter, maxter), gamma, work(maxter), w(maxter), w0  
  real(r8), external :: ddot, dnrm2 ! BLAS routines

  if (iter.lt.1) then
    errorflag = -1
    return
  endif
  if (smix%broyden_number.gt.maxter) then
    errorflag = -2
    return
  endif
  
  w0 = 0.01_r8
  w  = maxter*1.0_r8

  vout = vout - vin
  iter_used = min(iter - 1, smix%broyden_number)
  
  ! ipos is the position in which results from the present iteraction 
  ! are stored. ipos = iter until ipos = smix%broyden_number, then back to 1, 2,...
  ipos = mod(iter - 2, smix%broyden_number) + 1
  
  if(iter.gt.1) then
    smix%df(:,:, ipos) = vout - smix%df(:,:, ipos)
    smix%dv(:,:, ipos) = vin  - smix%dv(:,:, ipos)

    do is = 1, nspin
      gamma = 1._r8/DNRM2(np, smix%df(1, is, ipos), 1)
      call dscal (np, gamma, smix%df(1, is, ipos), 1)
      call dscal (np, gamma, smix%dv(1, is, ipos), 1)
    end do
  endif
    
  ! save values for next iteration
  i = mod(iter - 1, smix%broyden_number) + 1
  smix%df(:,:, smix%broyden_number + 1) = vout
  smix%dv(:,:, smix%broyden_number + 1) = vin

  is_loop: do is = 1, nspin
    beta = 0._r8
    do i = 1, iter_used
      do j = i + 1, iter_used
        beta(i, j) = w(i)*w(j)*DDOT(np, smix%df(1, is, j), 1, smix%df(1, is, i), 1)
        beta(j, i) = beta(i, j)
      end do
      beta(i, i) = w0**2 + w(i)**2
    end do

    ! invert matrix beta
    call dsytrf('u', iter_used, beta, maxter, iwork, work, maxter, info)
    if(info .ne. 0) then
      errorflag = 1
      return
    end if

    call dsytri('u', iter_used, beta, maxter, iwork, work, info)
    if(info .ne. 0) then
      errorflag = 2
      return
    end if
  
    ! complete the matrix
    do i = 1, iter_used
      do j = i + 1, iter_used
        beta(j, i) = beta(i, j)
      end do
    end do
  
    do i = 1, iter_used
      work(i) = ddot(np, smix%df(1, is, i), 1, vout(1, is), 1)
    end do
  
    call daxpy(np, smix%alpha, vout(1, is), 1, vin(1, is), 1)
  
    do i = 1, iter_used
      gamma = 0.d0
      do j = 1, iter_used
        gamma = gamma + beta(j, i)*w(j)*work(j)
      end do
      vin(:, is) = vin(:, is) - w(i)*gamma*(smix%alpha*smix%df(:, is, i) + smix%dv(:, is, i))
    end do
  
  end do is_loop

  ! put values in right positions
  i = mod(iter - 1, smix%broyden_number) + 1
  smix%df(:,:, i) = smix%df(:,:, smix%broyden_number + 1)
  smix%dv(:,:, i) = smix%dv(:,:, smix%broyden_number + 1)

  return
end subroutine mix_broyden
  
#include "undef.F90"
#include "real.F90"
#include "mix_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mix_inc.F90"
#include "undef.F90"

end module mix
