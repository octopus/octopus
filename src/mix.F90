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
use oct_parser

implicit none

private
public :: mix_type, mix_init, mix_end, mixing

integer, parameter :: LINEAR    = 0, &
                      ANDERSON  = 1, &
                      BROYDEN   = 2

! Anderson Mixing
!!$ integer, parameter :: A_IN = 1, A_OUT = 2

type mix_type
  private
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
subroutine mix_init(smix, np, nspin)
  type(mix_type), intent(out) :: smix
  integer, intent(in) :: np
  integer, intent(in) :: nspin

  call push_sub('mix_init')
  
  ! check input parameters
  call oct_parse_int("TypeOfMixing", 2, smix%type_of_mixing)
  if(smix%type_of_mixing == 1) then
    message(1) = 'Anderson mixing (TypeOfMixing = 1) currently not implemented)'
    call write_fatal(1)
  end if
  if(smix%type_of_mixing < 0 .or. smix%type_of_mixing > 2) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  call oct_parse_double("Mixing", 0.3_r8, smix%alpha);
  if(smix%alpha <= 0.0_r8 .or. smix%alpha > 1.0_r8) then
    write(message(1), '(a, f14.6,a)') "Input: '", smix%alpha, &
         "' is not a valid Mixing"
    message(2) = '(0 < Mixing <= 1)'
    call write_fatal(2)
  end if

  if(smix%type_of_mixing == BROYDEN) then
    call oct_parse_int("BroydenNumber", 3, smix%broyden_number)
    if(smix%broyden_number <= 1 .or. smix%broyden_number > 5) then
      write(message(1), '(a, i4,a)') "Input: '", smix%broyden_number, &
           "' is not a valid BroydenNumber"
      message(2) = '(1 < BroydenNumber <= 5)'
      call write_fatal(2)
    end if
    allocate(smix%df(np, nspin, smix%broyden_number + 1))
    allocate(smix%dv(np, nspin, smix%broyden_number + 1))
    smix%df = 0._r8
    smix%dv = 0._r8
    write(message(1), '(a)') 'Info: Broyden mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)
  end if

  call pop_sub()
  return
end subroutine mix_init

subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix

  call push_sub('mix_end')

  if(smix%type_of_mixing == BROYDEN) then
    if(associated(smix%df)) then
      deallocate(smix%df); nullify(smix%df);
      deallocate(smix%dv); nullify(smix%dv);
    end if
  end if

  call pop_sub()
end subroutine mix_end

subroutine mixing(smix, iter, np, nspin, vin, vout, vnew)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)      :: iter, np, nspin
  real(r8), dimension(np, nspin), intent(in) :: vin, vout
  real(r8), dimension(np, nspin), intent(out) :: vnew

  integer :: is, errorflag
    
  call push_sub('mixing')

  select case(smix%type_of_mixing)
  case(LINEAR)
    call mix_linear(smix, np, nspin, vin, vout, vnew)
!!$  case(ANDERSON)
!!$      call anderson_mix(rho,rhoout, iter, errorflag)
  case(BROYDEN)
    call mix_broyden(smix, np, nspin, vin, vout, vnew, iter, errorflag)
    if(errorflag .ne. 0) then
      write(message(1), '(a,i3)') 'mix_broyden returned error ', errorflag
      call write_fatal(1)
    endif
  end select

  call pop_sub()
end subroutine mixing

! Performs the linear mixing...
subroutine mix_linear(smix, np, nspin, vin, vout, vnew)
  type(mix_type), intent(IN) :: smix
  integer, intent(in) :: np, nspin
  real(r8), dimension(np, nspin),  intent(in) :: vin, vout
  real(r8), dimension(np, nspin), intent(out) :: vnew

  vnew = vin
  call dscal(np*nspin, 1.0_r8 - smix%alpha, vnew, 1)
  call daxpy(np*nspin, smix%alpha, vnew, 1, vnew, 1)

  return
end subroutine mix_linear

! Broyden mixing...
subroutine mix_broyden(smix, np, nspin, vin, vout, vnew, iter, errorflag)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: nspin, np
  integer, intent(in)     :: iter
  integer, intent(out)    :: errorflag
  real(r8), dimension(np, nspin), intent(in)  :: vin, vout 
  real(r8), dimension(np, nspin), intent(out) ::vnew

  integer, parameter :: maxter = 5

  integer  :: is, i, j, iwork(maxter), info, iter_used, ipos
  real(r8) :: beta(maxter, maxter), gamma, work(maxter), w(maxter), w0, r(np, nspin)
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

  r = vout - vin
  iter_used = min(iter - 1, smix%broyden_number)
  
  ! ipos is the position in which results from the present iteraction 
  ! are stored. ipos = iter until ipos = smix%broyden_number, then back to 1, 2,...
  ipos = mod(iter - 2, smix%broyden_number) + 1
  
  if(iter.gt.1) then
    smix%df(:,:, ipos) = r - smix%df(:,:, ipos)
    smix%dv(:,:, ipos) = vin  - smix%dv(:,:, ipos)

    do is = 1, nspin
      gamma = dnrm2(np, smix%df(1, is, ipos), 1)
      if(gamma > 1e-8_r8) then
         gamma = 1.0_r8/gamma
      else
         gamma = 1.0_r8
      endif
      !gamma = 1._r8/DNRM2(np, smix%df(1, is, ipos), 1)
      call dscal (np, gamma, smix%df(1, is, ipos), 1)
      call dscal (np, gamma, smix%dv(1, is, ipos), 1)
    end do
  endif
    
  ! save values for next iteration
  i = mod(iter - 1, smix%broyden_number) + 1
  smix%df(:,:, smix%broyden_number + 1) = r
  smix%dv(:,:, smix%broyden_number + 1) = vin

  is_loop: do is = 1, nspin
    beta = 0._r8
    do i = 1, iter_used
      do j = i + 1, iter_used
        beta(i, j) = w(i)*w(j)*DDOT(np, smix%df(1:np, is, j), 1, smix%df(1:np, is, i), 1)
        beta(j, i) = beta(i, j)
      end do
      beta(i, i) = w0**2 + w(i)**2
    end do

    ! invert matrix beta
    call dsytrf('u', iter_used, beta, maxter, iwork, work, maxter, info)
    if(info .ne. 0) then
      write(message(1), '(a, i3)') 'mix_broyden: dsytrf returned info = ', info
      call write_fatal(1)
    end if

    call dsytri('u', iter_used, beta, maxter, iwork, work, info)
    if(info .ne. 0) then
      write(message(1), '(a, i3)') 'mix_broyden: dsytri returned info = ', info
      call write_fatal(1)
    end if
  
    ! complete the matrix
    do i = 1, iter_used
      do j = i + 1, iter_used
        beta(j, i) = beta(i, j)
      end do
    end do
  
    do i = 1, iter_used
      work(i) = ddot(np, smix%df(1, is, i), 1, r(1, is), 1)
    end do

    call daxpy(np, smix%alpha, r(1, is), 1, vnew(1, is), 1)
  
    do i = 1, iter_used
      gamma = M_ZERO
      do j = 1, iter_used
        gamma = gamma + beta(j, i)*w(j)*work(j)
      end do
      vnew(:, is) = vnew(:, is) - w(i)*gamma*(smix%alpha*smix%df(:, is, i) + smix%dv(:, is, i))
    end do
  
  end do is_loop

  ! put values in right positions
  i = mod(iter - 1, smix%broyden_number) + 1
  smix%df(:,:, i) = smix%df(:,:, smix%broyden_number + 1)
  smix%dv(:,:, i) = smix%dv(:,:, smix%broyden_number + 1)

  errorflag = 0; return
end subroutine mix_broyden

end module mix
