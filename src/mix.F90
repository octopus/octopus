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

#include "global.h"

module mix
use global
use lib_oct_parser
use lib_adv_alg

implicit none

private
public :: mix_type, mix_init, mix_end, mixing

integer, parameter :: LINEAR  = 0, &
                      GRPULAY = 1, &
                      BROYDEN = 2

type mix_type
  private
  integer  :: type_of_mixing

  FLOAT :: alpha               !  vnew = (1-alpha)*vin + alpha*vout

  integer :: ns    ! number of steps used to extrapolate the new vector

  FLOAT, pointer :: df(:,:,:),  dv(:,:,:), f_old(:,:), vin_old(:,:)

  integer :: last_ipos

  logical :: icomp
  integer :: nc
end type mix_type

contains

! Initialization...
subroutine mix_init(smix, np, nv)
  type(mix_type), intent(out) :: smix
  integer, intent(in) :: np
  integer, intent(in) :: nv

  integer :: i

  call push_sub('mix_init')

  ! check input parameters
  call loct_parse_int("TypeOfMixing", 2, smix%type_of_mixing)
  if(smix%type_of_mixing < 0 .or. smix%type_of_mixing > 2) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  if (smix%type_of_mixing == LINEAR .or. smix%type_of_mixing == BROYDEN) then
    call loct_parse_float("Mixing", CNST(0.3), smix%alpha)
    if(smix%alpha <= M_ZERO .or. smix%alpha > M_ONE) then
      write(message(1), '(a, f14.6,a)') "Input: '",smix%alpha,"' is not a valid Mixing"
      message(2) = '(0 < Mixing <= 1)'
      call write_fatal(2)
    end if
  end if

  if (smix%type_of_mixing == GRPULAY .or. smix%type_of_mixing == BROYDEN) then
    call loct_parse_int("MixNumberSteps", 3, smix%ns)
    if(smix%ns <= 1) then
      write(message(1), '(a, i4,a)') "Input: '", smix%ns, &
                     "' is not a valid MixNumberSteps"
      message(2) = '(1 < MixNumberSteps)'
      call write_fatal(2)
    end if

    if (nv /= 1) then
      call loct_parse_logical("MixComponentsIndep", .false., smix%icomp)
      if (smix%icomp) then
        smix%nc = 1
      else
        smix%nc = nv
      end if
    else
      smix%icomp = .true.
    end if
  end if


  select case (smix%type_of_mixing)
  case (GRPULAY)
    write(message(1), '(a)') 'Info: GR-Pulay mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    if (smix%icomp) then
      allocate(smix%df(np, smix%ns + 1, nv), smix%vin_old(np, nv), &
               smix%dv(np, smix%ns + 1, nv), smix%f_old  (np, nv)   )
    else
      allocate(smix%df(np*nv, smix%ns + 1, 1), smix%vin_old(np*nv, 1), &
               smix%dv(np*nv, smix%ns + 1, 1), smix%f_old  (np*nv, 1)   )
    end if
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  case (BROYDEN)
    write(message(1), '(a)') 'Info: Broyden mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    if (smix%icomp) then
      allocate(smix%df(np, smix%ns, nv), smix%vin_old(np, nv), &
               smix%dv(np, smix%ns, nv), smix%f_old  (np, nv)   )
    else
      allocate(smix%df(np*nv, smix%ns, 1), smix%vin_old(np*nv, 1), &
               smix%dv(np*nv, smix%ns, 1), smix%f_old  (np*nv, 1)   )
    end if
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  end select

  smix%last_ipos = 0

  call pop_sub()
end subroutine mix_init

subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix

  call push_sub('mix_end')

  if (associated(smix%df))      deallocate(smix%df)
  if (associated(smix%dv))      deallocate(smix%dv)
  if (associated(smix%vin_old)) deallocate(smix%vin_old)
  if (associated(smix%f_old))   deallocate(smix%f_old)

  call pop_sub()
end subroutine mix_end

subroutine mixing(smix, iter, np, nv, vin, vout, vnew)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)      :: iter, np, nv
  FLOAT, dimension(np, nv), intent(in) :: vin, vout
  FLOAT, dimension(np, nv), intent(out) :: vnew

  integer :: i
  FLOAT, allocatable :: vint(:,:), voutt(:,:), vnewt(:,:)

  call push_sub('mixing')

  if (iter.lt.1) then
    message(1) = 'Wrong number of iterations in suboutine mixing.'
    call write_fatal(1)
  end if

  select case (smix%type_of_mixing)
  case (LINEAR)
    call mix_linear(smix%alpha, np, nv, vin, vout, vnew)

  case (BROYDEN, GRPULAY)
    if (smix%icomp) then
      select case(smix%type_of_mixing)
      case (BROYDEN)
        call mix_broyden(smix, np, nv, vin, vout, vnew, iter)
      case (GRPULAY)
        call mix_grpulay(smix, np, nv, vin, vout, vnew, iter)
      end select

    else
      ! Build total vectors
      allocate(vint(np*nv, 1), voutt(np*nv, 1), vnewt(np*nv, 1))
      do i = 1, nv
        vint((i-1)*np+1:i*np, 1) = vin(:, i)
        voutt((i-1)*np+1:i*np, 1) = vout(:, i)
      end do

      select case(smix%type_of_mixing)
      case (BROYDEN)
        call mix_broyden(smix, np*nv, 1, vint, voutt, vnewt, iter)
      case (GRPULAY)
        call mix_grpulay(smix, np*nv, 1, vint, voutt, vnewt, iter)
      end select

      ! recover vnew from vnewt
      do i = 1, nv
        vnew(:, i) = vnewt((i-1)*np+1:i*np, 1)
      end do
      deallocate(vint, voutt, vnewt)

    end if
  end select

  call pop_sub()
end subroutine mixing

! Performs the linear mixing...
subroutine mix_linear(alpha, np, nv, vin, vout, vnew)
  FLOAT, intent(in) :: alpha
  integer, intent(in) :: np, nv
  FLOAT, dimension(np, nv), intent(in) :: vin, vout
  FLOAT, dimension(np, nv), intent(out) :: vnew

  vnew = vin*(M_ONE - alpha) + alpha*vout

end subroutine mix_linear

! Broyden mixing...
subroutine mix_broyden(smix, np, nv, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: np, nv
  integer, intent(in)     :: iter
  FLOAT, dimension(np, nv), intent(in)  :: vin, vout 
  FLOAT, dimension(np, nv), intent(out) :: vnew

  integer :: i, ipos, iter_used
  FLOAT :: gamma, f(np, nv)

  f = vout - vin
  if(iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1

    smix%df(:, ipos, :) = f - smix%f_old
    smix%dv(:, ipos, :) = vin - smix%vin_old

    do i = 1,nv
      gamma = sqrt(dot_product(smix%df(:, ipos, i), smix%df(:, ipos, i)))
      if(gamma > CNST(1e-8)) then
        gamma = M_ONE/gamma
      else
        gamma = M_ONE
      endif
      smix%df(:, ipos, i) = smix%df(:, ipos, i)*gamma
      smix%dv(:, ipos, i) = smix%dv(:, ipos, i)*gamma
    end do

    smix%last_ipos = ipos
  end if
  ! Store residual and vin for next iteration
  smix%vin_old = vin
  smix%f_old = f

  ! extrapotate new vector
  iter_used = min(iter - 1, smix%ns)
  do i = 1,nv
    call broyden_extrapolation(smix%alpha, np, vin(:, i), vout(:, i), vnew(:, i), iter_used, &
                               f(:, i), smix%df(:, 1:iter_used, i), smix%dv(:, 1:iter_used, i))
  end do

end subroutine mix_broyden

subroutine broyden_extrapolation(alpha, np, vin, vout, vnew, iter_used, f, df, dv)
  FLOAT, intent(in)  :: alpha
  integer, intent(in)   :: np, iter_used
  FLOAT, intent(in)  :: vin(np), vout(np), f(np), df(np, iter_used), dv(np, iter_used)
  FLOAT, intent(out) :: vnew(np)

  FLOAT, parameter :: w0 = CNST(0.01)

  integer  :: i, j
  FLOAT :: beta(iter_used, iter_used), gamma, work(iter_used), w(iter_used)

  if (iter_used == 0) then
    ! linear mixing...
    vnew = vin + alpha*f
    return
  end if

  w  = M_FIVE

  ! compute matrix beta
  beta = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      beta(i, j) = w(i)*w(j)*dot_product(df(:, j), df(:, i))
      beta(j, i) = beta(i, j)
    end do
    beta(i, i) = w0**2 + w(i)**2
  end do

  ! invert matrix beta
  call lalg_invert(iter_used, beta)

  do i = 1, iter_used
    work(i) = dot_product(df(:, i), f)
  end do

  ! linear mixing term
  vnew = vin + alpha*f

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
subroutine mix_grpulay(smix, np, nv, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)     :: np, nv
  integer, intent(in)     :: iter
  FLOAT, dimension(np, nv), intent(in)  :: vin, vout 
  FLOAT, dimension(np, nv), intent(out) :: vnew

  integer :: ipos, iter_used, i
  FLOAT :: f(np, nv)

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
    do i = 1, nv
      call pulay_extrapolation(np, vin(:, i), vout(:, i), vnew(:, i), iter_used, &
                               f(:, i), smix%df(:, 1:iter_used, i), smix%dv(:, 1:iter_used, i))
    end do
  end select

end subroutine mix_grpulay

subroutine pulay_extrapolation(np, vin, vout, vnew, iter_used, f, df, dv)
  integer, intent(in)   :: np, iter_used
  FLOAT, intent(in)  :: vin(np), vout(np), f(np), df(np, iter_used), dv(np, iter_used)
  FLOAT, intent(out) :: vnew(np)

  integer :: i, j
  FLOAT :: a(iter_used, iter_used), alpha

  ! set matrix A
  a = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      a(i, j) = dot_product(df(:, j), df(:, i))
      a(j, i) = a(i, j)
    end do
    a(i, i) = dot_product(df(:, i), df(:, i))
  end do
  if (all(a < 1.0E-8)) then
    ! residuals are too small. Do not mix.
    vnew = vout
    return
  end if

  ! invert matrix A
  call lalg_invert(iter_used, a)

  ! compute new vector
  vnew = vin
  do i = 1,iter_used
    alpha = M_ZERO
    do j = 1,iter_used
      alpha = alpha - a(j, i)*dot_product(df(:, j), f)
    end do
    vnew = vnew + alpha * dv(:, i)
  end do

end subroutine pulay_extrapolation

end module mix
