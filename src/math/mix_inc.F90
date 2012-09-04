!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$

! ---------------------------------------------------------
subroutine X(mixing)(smix, iter, vin, vout, vnew, dotp)
  type(mix_t),  intent(inout) :: smix
  integer,      intent(in)    :: iter
  R_TYPE,       intent(in)    :: vin(:, :, :), vout(:, :, :)
  R_TYPE,       intent(out)   :: vnew(:, :, :)
  interface
    R_TYPE function dotp(x, y)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  
  PUSH_SUB(X(mixing))
  
  ASSERT(iter >= 1)
  
  select case (smix%type_of_mixing)
  case (MIX_LINEAR)
    call X(mixing_linear)(smix%alpha, smix%d1, smix%d2, smix%d3, vin, vout, vnew)
    
  case (MIX_BROYDEN)
    call X(mixing_broyden)(smix, vin, vout, vnew, iter, dotp)
    
  case (MIX_GRPULAY)
    call X(mixing_grpulay)(smix, vin, vout, vnew, iter, dotp)
    
  end select
  
  POP_SUB(X(mixing))
end subroutine X(mixing)


! ---------------------------------------------------------
subroutine X(mixing_linear)(alpha, d1, d2, d3, vin, vout, vnew)
  FLOAT,        intent(in) :: alpha
  integer,      intent(in) :: d1, d2, d3
  R_TYPE,       intent(in) :: vin(:, :, :), vout(:, :, :)
  R_TYPE,       intent(out):: vnew(:, :, :)

  PUSH_SUB(X(mixing_linear))
  
  vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3)*(M_ONE - alpha) + alpha*vout(1:d1, 1:d2, 1:d3)
  
  POP_SUB(X(mixing_linear))
end subroutine X(mixing_linear)


! ---------------------------------------------------------
subroutine X(mixing_broyden)(smix, vin, vout, vnew, iter, dotp)
  type(mix_t), intent(inout) :: smix
  R_TYPE,         intent(in) :: vin(:, :, :), vout(:, :, :)
  R_TYPE,         intent(out):: vnew(:, :, :)
  integer,        intent(in) :: iter
  interface
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface

  integer :: ipos, iter_used, i, j, d1, d2, d3
  R_TYPE :: gamma
  R_TYPE, allocatable :: f(:, :, :)

  PUSH_SUB(X(mixing_broyden))

  d1 = smix%d1
  d2 = smix%d2
  d3 = smix%d3
  
  SAFE_ALLOCATE(f(1:d1, 1:d2, 1:d3))
  
  f(1:d1, 1:d2, 1:d3) = vout(1:d1, 1:d2, 1:d3) - vin(1:d1, 1:d2, 1:d3)
  if(iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1
    
    call lalg_copy(d1, d2, d3, f(:, :, :), smix%X(df)(:, :, :, ipos))
    call lalg_copy(d1, d2, d3, vin(:, :, :), smix%X(dv)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(f_old)(:, :, :),   smix%X(df)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(vin_old)(:, :, :), smix%X(dv)(:, :, :, ipos))

    
    gamma = M_ZERO
    do i = 1, d2
      do j = 1, d3
        gamma = gamma + dotp(smix%X(df)(:, i, j, ipos), smix%X(df)(:, i, j, ipos))
      end do
    end do
    gamma = sqrt(gamma)
    
    if(abs(gamma) > CNST(1e-8)) then
      gamma = M_ONE/gamma
    else
      gamma = M_ONE
    end if
    call lalg_scal(d1, d2, d3, gamma, smix%X(df)(:, :, :, ipos))
    call lalg_scal(d1, d2, d3, gamma, smix%X(dv)(:, :, :, ipos))
    
    smix%last_ipos = ipos
  end if
  
  ! Store residual and vin for next iteration
  smix%X(vin_old)(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3)
  smix%X(f_old)  (1:d1, 1:d2, 1:d3) = f  (1:d1, 1:d2, 1:d3)
  
  ! extrapolate new vector
  iter_used = min(iter - 1, smix%ns)
  call X(broyden_extrapolation)(smix%alpha, d1, d2, d3, vin, vnew, iter_used, f, &
       smix%X(df), smix%X(dv), dotp)
  
  SAFE_DEALLOCATE_A(f)
  
  POP_SUB(X(mixing_broyden))
end subroutine X(mixing_broyden)


! ---------------------------------------------------------
subroutine X(broyden_extrapolation)(alpha, d1, d2, d3, vin, vnew, iter_used, f, df, dv, dotp)
  FLOAT,    intent(in)       :: alpha
  integer,  intent(in)       :: d1, d2, d3, iter_used
  R_TYPE,   intent(in)       :: vin(:, :, :), f(:, :, :)
  R_TYPE,   intent(in)       :: df(:, :, :, :), dv(:, :, :, :)
  R_TYPE,   intent(out)      :: vnew(:, :, :)
  interface
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  
  FLOAT, parameter :: w0 = CNST(0.01)
  integer  :: i, j, k, l
  R_TYPE    :: gamma, x
  R_TYPE, allocatable :: beta(:, :), work(:), w(:)

  PUSH_SUB(X(broyden_extrapolation))
  
  if (iter_used == 0) then
    ! linear mixing...
    vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3) + alpha*f(1:d1, 1:d2, 1:d3)
    POP_SUB(X(broyden_extrapolation))
    return
  end if
  
  SAFE_ALLOCATE(beta(1:iter_used, 1:iter_used))
  SAFE_ALLOCATE(work(1:iter_used))
  SAFE_ALLOCATE(   w(1:iter_used))

  w  = M_FIVE
  
  ! compute matrix beta
  beta = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      beta(i, j) = M_ZERO
      do k = 1, d2
        do l = 1, d3
          beta(i, j) = beta(i, j) + w(i)*w(j)*dotp(df(:, k, l, j), df(:, k, l, i))
        end do
      end do
      beta(j, i) = beta(i, j)
    end do
    beta(i, i) = w0**2 + w(i)**2
  end do
  
  ! invert matrix beta
  x = lalg_inverter(iter_used, beta)
  
  do i = 1, iter_used
    work(i) = M_ZERO
    do k = 1, d2
      do l = 1, d3
        work(i) = work(i) + dotp(df(:, k, l, i), f(:, k, l))
      end do
    end do
  end do
  
  ! linear mixing term
  vnew(1:d1, 1:d2, 1:d3) = vin(1:d1, 1:d2, 1:d3) + alpha*f(1:d1, 1:d2, 1:d3)
  
  ! other terms
  do i = 1, iter_used
    gamma = M_ZERO
    do j = 1, iter_used
      gamma = gamma + beta(j, i)*w(j)*work(j)
    end do
    vnew(1:d1, 1:d2, 1:d3) = vnew(1:d1, 1:d2, 1:d3) - w(i)*gamma*(alpha*df(1:d1, 1:d2, 1:d3, i) + &
        dv(1:d1, 1:d2, 1:d3, i))
  end do
  
  SAFE_DEALLOCATE_A(beta)
  SAFE_DEALLOCATE_A(work)
  SAFE_DEALLOCATE_A(w)
  
  POP_SUB(X(broyden_extrapolation))
end subroutine X(broyden_extrapolation)


! ---------------------------------------------------------
! Guaranteed-reduction Pulay
! ---------------------------------------------------------
subroutine X(mixing_grpulay)(smix, vin, vout, vnew, iter, dotp)
  type(mix_t), intent(inout) :: smix
  integer,      intent(in)   :: iter
  R_TYPE,        intent(in)   :: vin(:, :, :), vout(:, :, :)
  R_TYPE,        intent(out)  :: vnew(:, :, :)
  interface
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  
  integer :: d1, d2, d3
  integer :: ipos, iter_used
  R_TYPE, allocatable :: f(:, :, :)
    
  PUSH_SUB(X(mixing_grpulay))

  d1 = smix%d1
  d2 = smix%d2
  d3 = smix%d3
  
  SAFE_ALLOCATE(f(1:d1, 1:d2, 1:d3))
  f(1:d1, 1:d2, 1:d3) = vout(1:d1, 1:d2, 1:d3) - vin(1:d1, 1:d2, 1:d3)
  
  ! we only extrapolate a new vector every two iterations
  select case (mod(iter, 2))
  case (1)
    ! Store df and dv from current iteration
    if (iter > 1) then
      ipos = smix%last_ipos
      call lalg_copy(d1, d2, d3, f(:, :, :), smix%X(df)(:, :, :, ipos))
      call lalg_copy(d1, d2, d3, vin(:, :, :), smix%X(dv)(:, :, :, ipos))
      call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(f_old)(:, :, :), smix%X(df)(:, :, :, ipos))
      call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(vin_old)(:, :, :), smix%X(dv)(:, :, :, ipos))
    end if
    
    ! Store residual and vin for next extrapolation
    smix%X(vin_old) = vin
    smix%X(f_old) = f
    
    ! we need the output vector for vout. So let`s do vnew = vout to get that information
    vnew = vout
  case (0)
    ! Store df and dv from current iteration in arrays df and dv so that we can use them
    ! for the extrapolation. Next iterations they will be lost.
    ipos = mod(smix%last_ipos, smix%ns + 1) + 1
    call lalg_copy(d1, d2, d3, f(:, :, :), smix%X(df)(:, :, :, ipos))
    call lalg_copy(d1, d2, d3, vin(:, :, :), smix%X(dv)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(f_old)(:, :, :), smix%X(df)(:, :, :, ipos))
    call lalg_axpy(d1, d2, d3, R_TOTYPE(-M_ONE), smix%X(vin_old)(:, :, :), smix%X(dv)(:, :, :, ipos))
    
    smix%last_ipos = ipos
    
    ! extrapolate new vector
    iter_used = min(iter/2, smix%ns + 1)
    call X(pulay_extrapolation)(d2, d3, vin, vout, vnew, iter_used, f, &
         smix%X(df)(1:d1, 1:d2, 1:d3, 1:iter_used), &
         smix%X(dv)(1:d1, 1:d2, 1:d3, 1:iter_used), dotp)

  end select
  
  SAFE_DEALLOCATE_A(f)
  
  POP_SUB(X(mixing_grpulay))
end subroutine X(mixing_grpulay)


! ---------------------------------------------------------
subroutine X(pulay_extrapolation)(d2, d3, vin, vout, vnew, iter_used, f, df, dv, dotp)
  integer, intent(in) :: d2, d3
  integer, intent(in)   :: iter_used
  R_TYPE,  intent(in)  :: vin(:, :, :), vout(:, :, :), f(:, :, :), df(:, :, :, :), dv(:, :, :, :)
  R_TYPE,  intent(out) :: vnew(:, :, :)
  interface
    R_TYPE function dotp(x, y) result(res)
      R_TYPE, intent(in) :: x(:)
      R_TYPE, intent(in) :: y(:)
    end function dotp
  end interface
  
  integer :: i, j, k, l
  R_TYPE :: alpha
  R_TYPE, allocatable :: a(:, :)
  
  PUSH_SUB(X(pulay_extrapolation))
  
  SAFE_ALLOCATE(a(1:iter_used, 1:iter_used))
  
  ! set matrix A
  a = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      a(i, j) = M_ZERO
      do k = 1, d2
        do l = 1, d3
          a(i, j) = a(i, j) + dotp(df(:, k, l, j), df(:, k, l, i))
        end do
      end do
      a(j, i) = a(i, j)
    end do
    a(i, i) = M_ZERO
    do k = 1, d2
      do l = 1, d3
        a(i, i) = a(i, i) + dotp(df(:, k, l, i), df(:, k, l, i))
      end do
    end do
  end do
  if (all(abs(a) < 1.0E-8)) then
      ! residuals are too small. Do not mix.
    vnew = vout
    POP_SUB(X(pulay_extrapolation))
    return
  end if
  
  alpha = lalg_inverter(iter_used, a)
  
  ! compute new vector
  vnew = vin
  do i = 1,iter_used
    alpha = M_ZERO
    do j = 1,iter_used
      do k = 1, d2
        do l = 1, d3
          alpha = alpha - a(i, j) * dotp(df(:, k, l, j), f(:, k, l))
        end do
      end do
    end do
    vnew(:, :, :) = vnew(:, :, :) + alpha * dv(:, :, :, i)
  end do
  
  SAFE_DEALLOCATE_A(a)
  
  POP_SUB(X(pulay_extrapolation))
end subroutine X(pulay_extrapolation)



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
