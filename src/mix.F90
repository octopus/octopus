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
use lib_basic_alg
use lib_adv_alg

implicit none

private
public :: mix_type, mix_init, mix_end, mixing

integer, parameter, public :: &
   MIX_LINEAR  = 0, &
   MIX_GRPULAY = 1, &
   MIX_BROYDEN = 2

type mix_type
  private
  integer  :: type_of_mixing

  FLOAT :: alpha               !  vnew = (1-alpha)*vin + alpha*vout

  integer :: ns    ! number of steps used to extrapolate the new vector

  FLOAT, pointer :: df(:, :, :, :),  dv(:, :, :, :), f_old(:, :, :), vin_old(:, :, :)

  integer :: last_ipos
end type mix_type

contains

! ---------------------------------------------------------
subroutine mix_init(smix, dim, np, nspin, def_)
  type(mix_type),    intent(out) :: smix
  integer,           intent(in)  :: dim, np, nspin
  integer, optional, intent(in)  :: def_

  integer :: def

  call push_sub('mix_init')

  def = MIX_BROYDEN
  if(present(def_)) def = def_

  ! check input parameters
  call loct_parse_int("TypeOfMixing", def, smix%type_of_mixing)
  if(smix%type_of_mixing < MIX_LINEAR .or. smix%type_of_mixing > MIX_BROYDEN) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  if (smix%type_of_mixing == MIX_LINEAR .or. smix%type_of_mixing == MIX_BROYDEN) then
    call loct_parse_float("Mixing", CNST(0.3), smix%alpha)
    if(smix%alpha <= M_ZERO .or. smix%alpha > M_ONE) then
      write(message(1), '(a, f14.6,a)') "Input: '",smix%alpha,"' is not a valid Mixing"
      message(2) = '(0 < Mixing <= 1)'
      call write_fatal(2)
    end if
  end if

  if (smix%type_of_mixing == MIX_GRPULAY .or. smix%type_of_mixing == MIX_BROYDEN) then
    call loct_parse_int("MixNumberSteps", 3, smix%ns)
    if(smix%ns <= 1) then
      write(message(1), '(a, i4,a)') "Input: '", smix%ns, &
                     "' is not a valid MixNumberSteps"
      message(2) = '(1 < MixNumberSteps)'
      call write_fatal(2)
    end if
  end if


  select case (smix%type_of_mixing)
  case (MIX_GRPULAY)
    write(message(1), '(a)') 'Info: GR-Pulay mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(dim, np, nspin, smix%ns + 1), smix%vin_old(dim, np, nspin), &
             smix%dv(dim, np, nspin, smix%ns + 1), smix%f_old  (dim, np, nspin)   )
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  case (MIX_BROYDEN)
    write(message(1), '(a)') 'Info: Broyden mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(dim, np, nspin, smix%ns), smix%vin_old(dim, np, nspin), &
             smix%dv(dim, np, nspin, smix%ns), smix%f_old  (dim, np, nspin)   )
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  end select

  smix%last_ipos = 0

  call pop_sub()
end subroutine mix_init


! ---------------------------------------------------------
subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix
  call push_sub('mix_end')

  if (associated(smix%df))      deallocate(smix%df)
  if (associated(smix%dv))      deallocate(smix%dv)
  if (associated(smix%vin_old)) deallocate(smix%vin_old)
  if (associated(smix%f_old))   deallocate(smix%f_old)

  call pop_sub()
end subroutine mix_end


! ---------------------------------------------------------
subroutine mixing(smix, iter, dim, np, nspin, vin, vout, vnew)
  type(mix_type), intent(inout) :: smix
  integer,        intent(in)    :: iter, dim, np, nspin
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)

  call push_sub('mixing')

  if (iter.lt.1) then
    message(1) = 'Wrong number of iterations in suboutine mixing.'
    call write_fatal(1)
  end if

  select case (smix%type_of_mixing)
  case (MIX_LINEAR)
    call mixing_linear(smix%alpha, np, vin, vout, vnew)

  case (MIX_BROYDEN)
    call mixing_broyden(smix, dim, np, nspin, vin, vout, vnew, iter)

  case (MIX_GRPULAY)
    call mixing_grpulay(smix, dim, np, nspin, vin, vout, vnew, iter)

  end select

  call pop_sub()
end subroutine mixing


! ---------------------------------------------------------
subroutine mixing_linear(alpha, np, vin, vout, vnew)
  FLOAT,   intent(in)  :: alpha
  integer, intent(in)  :: np
  FLOAT,   intent(in)  :: vin(:, :, :), vout(:, :, :)
  FLOAT,   intent(out) :: vnew(:, :, :)

  vnew = vin*(M_ONE - alpha) + alpha*vout

end subroutine mixing_linear


! ---------------------------------------------------------
subroutine mixing_broyden(smix, dim, np, nspin, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer,        intent(in)    :: dim, np, nspin, iter
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)

  integer :: i, ipos, iter_used
  FLOAT :: gamma
  FLOAT, allocatable :: f(:, :, :)

  allocate(f(dim, np, nspin))

  f(1:dim, 1:np, 1:nspin) = vout(1:dim, 1:np, 1:nspin) - vin(1:dim, 1:np, 1:nspin)
  if(iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1

    call lalg_copy(dim, np, nspin, f(:, :, :), smix%df(:, :, :, ipos))
    call lalg_copy(dim, np, nspin, vin(:, :, :), smix%dv(:, :, :, ipos))
    call lalg_axpy(dim, np, nspin, -M_ONE, smix%f_old(:, :, :),   smix%df(:, :, :, ipos))
    call lalg_axpy(dim, np, nspin, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))

    gamma = sqrt(sum(smix%df(1:dim, 1:np, 1:nspin, ipos)*smix%df(1:dim, 1:np, 1:nspin, ipos)))

    if(gamma > CNST(1e-8)) then
      gamma = M_ONE/gamma
    else
      gamma = M_ONE
    endif
    call lalg_scal(dim, np, nspin, gamma, smix%df(:, :, :, ipos))
    call lalg_scal(dim, np, nspin, gamma, smix%dv(:, :, :, ipos))

    smix%last_ipos = ipos
  end if

  ! Store residual and vin for next iteration
  smix%vin_old = vin
  smix%f_old   = f

  ! extrapotate new vector
  iter_used = min(iter - 1, smix%ns)
  call broyden_extrapolation(smix%alpha, dim, np, nspin, vin, vout, vnew, iter_used, f, &
     smix%df(1:dim, 1:np, 1:nspin, 1:iter_used), &
     smix%dv(1:dim, 1:np, 1:nspin, 1:iter_used))

  deallocate(f)

end subroutine mixing_broyden


! ---------------------------------------------------------
subroutine broyden_extrapolation(alpha, dim, np, nspin, vin, vout, vnew, iter_used, f, df, dv)
  FLOAT,   intent(in)  :: alpha
  integer, intent(in)  :: dim, np, nspin, iter_used
  FLOAT,   intent(in)  :: vin(:, :, :), vout(:, :, :), f(:, :, :)
  FLOAT,   intent(in)  :: df(:, :, :, :), dv(:, :, :, :)
  FLOAT,   intent(out) :: vnew(:, :, :)


  FLOAT, parameter :: w0 = CNST(0.01)

  integer  :: i, j, ispin, idim
  FLOAT :: beta(iter_used, iter_used), gamma, work(iter_used), w(iter_used)

  if (iter_used == 0) then
    ! linear mixing...
    vnew(1:dim, 1:np, 1:nspin) = vin(1:dim, 1:np, 1:nspin) + alpha*f(1:dim, 1:np, 1:nspin)
    return
  end if

  w  = M_FIVE

  ! compute matrix beta
  beta = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      beta(i, j) = w(i)*w(j)*sum(df(1:dim, 1:np, 1:nspin, j)*df(1:dim, 1:np, 1:nspin, i))
      beta(j, i) = beta(i, j)
    end do
    beta(i, i) = w0**2 + w(i)**2
  end do

  ! invert matrix beta
  call lalg_invert(iter_used, beta)

  do i = 1, iter_used
    work(i) = sum(df(1:dim, 1:np, 1:nspin, i)*f(1:dim, 1:np, 1:nspin)) ! dot_product(df(:, :, :, i), f)
  end do


  ! linear mixing term
  vnew(1:dim, 1:np, 1:nspin) = vin(1:dim, 1:np, 1:nspin) + alpha*f(1:dim, 1:np, 1:nspin)

  ! other terms
  do i = 1, iter_used
    gamma = M_ZERO
    do j = 1, iter_used
      gamma = gamma + beta(j, i)*w(j)*work(j)
    end do
    vnew(1:dim, 1:np, 1:nspin) = vnew(1:dim, 1:np, 1:nspin) - w(i)*gamma*(alpha*df(1:dim, 1:np, 1:nspin, i) + &
                                                                          dv(1:dim, 1:np, 1:nspin, i))
  end do

end subroutine broyden_extrapolation


! ---------------------------------------------------------
! Guaranteed-reduction Pulay
! ---------------------------------------------------------
subroutine mixing_grpulay(smix, dim, np, nspin, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  integer,        intent(in)    :: dim, np, nspin
  integer,        intent(in)    :: iter
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)

  integer :: ipos, iter_used, i
  FLOAT, allocatable :: f(:, :, :)

  allocate(f(dim, np, nspin))
  f = vout - vin

  ! we only extrapolate a new vector every two iterations
  select case (mod(iter, 2_i4))
  case (1)
    ! Store df and dv from current iteration
    if (iter > 1) then
      ipos = smix%last_ipos
      call lalg_copy(dim, np, nspin, f(:, :, :), smix%df(:, :, :, ipos))
      call lalg_copy(dim, np, nspin, vin(:, :, :), smix%dv(:, :, :, ipos))
      call lalg_axpy(dim, np, nspin, -M_ONE, smix%f_old(:, :, :), smix%df(:, :, :, ipos))
      call lalg_axpy(dim, np, nspin, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))
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
    call lalg_copy(dim, np, nspin, f(:, :, :), smix%df(:, :, :, ipos))
    call lalg_copy(dim, np, nspin, vin(:, :, :), smix%dv(:, :, :, ipos))
    call lalg_axpy(dim, np, nspin, -M_ONE, smix%f_old(:, :, :), smix%df(:, :, :, ipos))
    call lalg_axpy(dim, np, nspin, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))

    smix%last_ipos = ipos

    ! extrapotate new vector
    iter_used = min(iter/2, smix%ns + 1)
    call pulay_extrapolation(dim, np, nspin, vin, vout, vnew, iter_used, f, &
                                 smix%df(1:dim, 1:np, 1:nspin, 1:iter_used), &
                                 smix%dv(1:dim, 1:np, 1:nspin, 1:iter_used))
  end select

  deallocate(f)
end subroutine mixing_grpulay


! ---------------------------------------------------------
subroutine pulay_extrapolation(dim, np, nspin, vin, vout, vnew, iter_used, f, df, dv)
  integer, intent(in)   :: dim, np, nspin, iter_used 
  FLOAT, intent(in)  :: vin(:, :, :), vout(:, :, :), f(:, :, :), df(:, :, :, :), dv(:, :, :, :)
  FLOAT, intent(out) :: vnew(:, :, :)

  integer :: i, j
  FLOAT :: a(iter_used, iter_used), alpha

  ! set matrix A
  a = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      a(i, j) = sum(df(:, :, :, j)*df(:, :, :, i))
      a(j, i) = a(i, j)
    end do
    a(i, i) = sum(df(:, :, :, i)*df(:, :, :, i))
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
       alpha = alpha - a(j, i)*sum(df(:, :, :, j)*f(:, :, :))
    end do
    vnew(:, :, :) = vnew(:, :, :) + alpha * dv(:, :, :, i)
  end do

end subroutine pulay_extrapolation

end module mix
