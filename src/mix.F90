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
use linalg

implicit none

private
public :: mix_type, mix_init, mix_end, mixing

integer, parameter :: LINEAR  = 0, &
                      GRPULAY = 1, &
                      BROYDEN = 2

type mix_type
  private
  integer  :: type_of_mixing

  real(r8) :: alpha               !  vnew = (1-alpha)*vin + alpha*vout

  integer :: ns    ! number of steps used to extrapolate the new vector

  real(r8), pointer :: df(:,:,:)    => NULL(), &
                       dv(:,:,:)    => NULL(), &
                       f_old(:,:)   => NULL(), &
                       vin_old(:,:) => NULL()

  integer :: last_ipos = 0
end type mix_type

contains

! Initialization...
subroutine mix_init(smix, np, nc)
  type(mix_type), intent(out) :: smix
  integer, intent(in) :: np
  integer, intent(in) :: nc

  call push_sub('mix_init')

  ! check input parameters
  call oct_parse_int("TypeOfMixing", 2, smix%type_of_mixing)
  if(smix%type_of_mixing < 0 .or. smix%type_of_mixing > 2) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  if (smix%type_of_mixing == LINEAR .or. smix%type_of_mixing == BROYDEN) then
    call oct_parse_double("Mixing", 0.3_r8, smix%alpha)
    if(smix%alpha <= 0.0_r8 .or. smix%alpha > 1.0_r8) then
      write(message(1), '(a, f14.6,a)') "Input: '",smix%alpha,"' is not a valid Mixing"
      message(2) = '(0 < Mixing <= 1)'
      call write_fatal(2)
    end if
  end if

  if (smix%type_of_mixing == GRPULAY .or. smix%type_of_mixing == BROYDEN) then
    call oct_parse_int("MixNumberSteps", 3, smix%ns)
    if(smix%ns <= 1) then
      write(message(1), '(a, i4,a)') "Input: '", smix%ns, &
                     "' is not a valid MixNumberSteps"
      message(2) = '(1 < MixNumberSteps)'
      call write_fatal(2)
    end if
  end if

  select case (smix%type_of_mixing)
  case (GRPULAY)
    write(message(1), '(a)') 'Info: GR-Pulay mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(np, smix%ns + 1, nc),&
             smix%dv(np, smix%ns + 1, nc),&
             smix%vin_old(np, nc),    &
             smix%f_old(np, nc)        )
    smix%df = 0._r8; smix%dv = 0._r8; smix%vin_old = 0._r8; smix%f_old = 0._r8

  case (BROYDEN)
    write(message(1), '(a)') 'Info: Broyden mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(np, smix%ns, nc),&
             smix%dv(np, smix%ns, nc),&
             smix%vin_old(np, nc),    &
             smix%f_old(np, nc)        )
    smix%df = 0._r8; smix%dv = 0._r8; smix%vin_old = 0._r8; smix%f_old = 0._r8

  end select

  call pop_sub()

  return
end subroutine mix_init

subroutine mix_end(smix)
  type(mix_type), intent(in) :: smix

  call push_sub('mix_end')

  if (associated(smix%df))      deallocate(smix%df)
  if (associated(smix%dv))      deallocate(smix%dv)
  if (associated(smix%vin_old)) deallocate(smix%vin_old)
  if (associated(smix%f_old))   deallocate(smix%f_old)

  call pop_sub()
end subroutine mix_end

subroutine mixing(smix, iter, np, nc, vin, vout, vnew)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)      :: iter, np, nc
  real(r8), dimension(np, nc), intent(in) :: vin, vout
  real(r8), dimension(np, nc), intent(out) :: vnew

  call push_sub('mixing')

  if (iter.lt.1) then
    message(1) = 'Wrong number of iterations in suboutine mixing.'
    call write_fatal(1)
  end if

  select case (smix%type_of_mixing)
  case (LINEAR)
    call mix_linear(smix%alpha, np, nc, vin, vout, vnew)

  case (BROYDEN)
    call mix_broyden(smix, np, nc, vin, vout, vnew, iter)

  case (GRPULAY)
    call mix_grpulay(smix, np, nc, vin, vout, vnew, iter)

  end select

  call pop_sub()
end subroutine mixing

! Performs the linear mixing...
subroutine mix_linear(alpha, np, nc, vin, vout, vnew)
  real(r8), intent(in) :: alpha
  integer, intent(in) :: np, nc
  real(r8), dimension(np, nc), intent(in) :: vin, vout
  real(r8), dimension(np, nc), intent(out) :: vnew

  call dcopy(np*nc, vin(1, 1), 1, vnew(1, 1), 1)
  call dscal(np*nc, 1.0_r8 - alpha, vnew(1, 1), 1)
  call daxpy(np*nc, alpha, vout(1, 1), 1, vnew(1, 1), 1)

  return
end subroutine mix_linear

! Broyden mixing...
subroutine mix_broyden(smix, np, nc, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: np, nc
  integer, intent(in)     :: iter
  real(r8), dimension(np, nc), intent(in)  :: vin, vout 
  real(r8), dimension(np, nc), intent(out) :: vnew

  integer :: i, ipos, iter_used
  real(r8) :: gamma, f(np, nc)

  f = vout - vin
  if(iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1

    smix%df(:, ipos, :) = f - smix%f_old
    smix%dv(:, ipos, :) = vin - smix%vin_old

    do i = 1,nc
      gamma = dnrm2(np, smix%df(1, ipos, i), 1)
      if(gamma > 1e-8_r8) then
        gamma = 1.0_r8/gamma
      else
        gamma = 1.0_r8
      endif
      call dscal (np, gamma, smix%df(1, ipos, i), 1)
      call dscal (np, gamma, smix%dv(1, ipos, i), 1)
    end do

    smix%last_ipos = ipos
  end if
  ! Store residual and vin for next iteration
  smix%vin_old = vin
  smix%f_old = f

  ! extrapotate new vector
  iter_used = min(iter - 1, smix%ns)
  do i = 1,nc
    call broyden_extrapolation(smix%alpha, np, vin(:, i), vout(:, i), vnew(:, i), iter_used, &
                               f(:, i), smix%df(:, 1:iter_used, i), smix%dv(:, 1:iter_used, i))
  end do

end subroutine mix_broyden

subroutine broyden_extrapolation(alpha, np, vin, vout, vnew, iter_used, f, df, dv)
  real(r8), intent(in)  :: alpha
  integer, intent(in)   :: np, iter_used
  real(r8), intent(in)  :: vin(np), vout(np), f(np), df(np, iter_used), dv(np, iter_used)
  real(r8), intent(out) :: vnew(np)

  real(r8), parameter :: w0 = 0.01_r8

  integer  :: i, j
  real(r8) :: beta(iter_used, iter_used), gamma, work(iter_used), w(iter_used)

  if (iter_used == 0) then
    ! linear mixing...
    call dcopy(np, vin(1), 1, vnew(1), 1)
    call daxpy(np, alpha, f(1), 1, vnew(1), 1) 
    return
  end if

  w  = 5.0_r8

  ! compute matrix beta
  beta = 0._r8
  do i = 1, iter_used
    do j = i + 1, iter_used
      beta(i, j) = w(i)*w(j)*DDOT(np, df(1, j), 1, df(1, i), 1)
      beta(j, i) = beta(i, j)
    end do
    beta(i, i) = w0**2 + w(i)**2
  end do

  ! invert matrix beta
  call dsyinvert(iter_used, iter_used, beta)

  do i = 1, iter_used
    work(i) = ddot(np, df(1, i), 1, f(1), 1)
  end do

  ! linear mixing term
  call dcopy(np, vin(1), 1, vnew(1), 1)
  call daxpy(np, alpha, f(1), 1, vnew(1), 1)

  ! other terms
  do i = 1, iter_used
    gamma = M_ZERO
    do j = 1, iter_used
      gamma = gamma + beta(j, i)*w(j)*work(j)
    end do
    vnew = vnew - w(i)*gamma*(alpha*df(:, i) + dv(:, i))
  end do

end subroutine broyden_extrapolation


! Guaranteed-reduction Pulay
subroutine mix_grpulay(smix, np, nc, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: np, nc
  integer, intent(in)     :: iter
  real(r8), dimension(np, nc), intent(in)  :: vin, vout 
  real(r8), dimension(np, nc), intent(out) :: vnew

  integer :: ipos, iter_used, i
  real(r8) :: f(np, nc)

  f = vout - vin

  ! we only extrapolate a new vector every two iterations
  select case (mod(iter, 2_i4))
  case (1)
    ! Store df and dv from current iteration
    if (iter > 1) then
      ipos = smix%last_ipos
      smix%df(:, ipos, :) = f - smix%f_old
      smix%dv(:, ipos, :) = vin - smix%vin_old
    end if

    ! Store residual and vin for next extrapolation
    smix%vin_old = vin
    smix%f_old = f

    ! we need the output vector for vout. So lets do vnew = vout to get that information
    vnew = vout
  case (0)
    ! Store df and dv from current iteration in arrays df and dv so that we can use them
    ! for the extrapolation. Next iterations they will be lost.
    ipos = mod(smix%last_ipos, smix%ns + 1) + 1
    smix%df(:, ipos, :) = f - smix%f_old
    smix%dv(:, ipos, :) = vin - smix%vin_old
    smix%last_ipos = ipos

    ! extrapotate new vector
    iter_used = min(iter/2, smix%ns + 1)
    do i = 1,nc
      call pulay_extrapolation(np, vin(:, i), vout(:, i), vnew(:, i), iter_used, &
                               f(:, i), smix%df(:, 1:iter_used, i), smix%dv(:, 1:iter_used, i))
    end do
  end select

end subroutine mix_grpulay

subroutine pulay_extrapolation(np, vin, vout, vnew, iter_used, f, df, dv)
  integer, intent(in)   :: np, iter_used
  real(r8), intent(in)  :: vin(np), vout(np), f(np), df(np, iter_used), dv(np, iter_used)
  real(r8), intent(out) :: vnew(np)

  integer :: i, j
  real(r8) :: a(iter_used, iter_used), alpha

  ! set matrix A
  a = 0._r8
  do i = 1, iter_used
    do j = i + 1, iter_used
      a(i, j) = ddot(np, df(1, j), 1, df(1, i), 1)
      a(j, i) = a(i, j)
    end do
    a(i, i) = ddot(np, df(1, i), 1, df(1, i), 1)
  end do
  if (all(a < 1.0E-8)) then
    ! residuals are too small. Do not mix.
    vnew = vout
    return
  end if

  ! invert matrix A
  call dsyinvert(iter_used, iter_used, a)

  ! compute new density
  vnew = vin
  do i = 1,iter_used
    alpha = 0._r8
    do j = 1,iter_used
      alpha = alpha - a(j, i)*DDOT(np, df(1, j), 1, f(1), 1)
    end do
    vnew = vnew + alpha * dv(:, i)
  end do

end subroutine pulay_extrapolation

end module mix
