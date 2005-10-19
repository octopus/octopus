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
!!
!! $Id$

#include "global.h"

module mix
use global
use messages
use syslabels
use mesh
use mesh_function
use lib_oct_parser
use lib_basic_alg
use lib_adv_alg
use varinfo

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
subroutine mix_init(smix, m, d2, d3, def_)
  type(mix_type),    intent(out) :: smix
  type(mesh_type),   intent(in)  :: m
  integer,           intent(in)  :: d2, d3
  integer, optional, intent(in)  :: def_

  integer :: def

  call push_sub('mix.mix_init')

  def = MIX_BROYDEN
  if(present(def_)) def = def_


  !%Variable TypeOfMixing
  !%Type integer
  !%Section 8 SCF
  !%Description
  !% The scheme scheme used to produce, at each iteration in the self consistent cycle
  !% that attempts to solve the Kohn-Sham equations, the input density from the value
  !% of the input and output densities of previous iterations.
  !%Option linear 0
  !% Simple linear mixing.
  !%Option gr_pulay 1
  !% "Guaranteed-reduction" Pulay scheme (FIXME: Add reference).
  !%Option broyden 2
  !% Broyden scheme (FIXME: Add reference).
  !%End
  call loct_parse_int(check_inp('TypeOfMixing'), def, smix%type_of_mixing)
  if(.not.varinfo_valid_option('TypeOfMixing', smix%type_of_mixing)) call input_error('TypeOfMixing')
  
  !%Variable Mixing
  !%Type float
  !%Section 8 SCF
  !%Description
  !% Both the linear and the Broyden scheme depend on a "mixing parameter", set by this variable.
  !%End
  if (smix%type_of_mixing == MIX_LINEAR .or. smix%type_of_mixing == MIX_BROYDEN) then
    call loct_parse_float(check_inp('Mixing'), CNST(0.3), smix%alpha)
    if(smix%alpha <= M_ZERO .or. smix%alpha > M_ONE) call input_error('Mixing')
  end if

  !%Variable MixNumberSteps
  !%Type integer
  !%Section 8 SCF
  !%Description
  !% In the Broyden and in the GR-Pulay scheme, the new input density or potential is constructed 
  !% from the values of densities/potentials of previous a given number of previous iterations.
  !% This number is set by this variable.
  !%End
  if (smix%type_of_mixing == MIX_GRPULAY .or. smix%type_of_mixing == MIX_BROYDEN) then
    call loct_parse_int(check_inp('MixNumberSteps'), 3, smix%ns)
    if(smix%ns <= 1) call input_error('MixNumberSteps')
  end if


  select case (smix%type_of_mixing)
  case (MIX_GRPULAY)
    write(message(1), '(a)') 'Info: GR-Pulay mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(m%np, d2, d3, smix%ns + 1), smix%vin_old(m%np, d2, d3), &
             smix%dv(m%np, d2, d3, smix%ns + 1), smix%f_old  (m%np, d2, d3)   )
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  case (MIX_BROYDEN)
    write(message(1), '(a)') 'Info: Broyden mixing used. It can (i) boost your convergence, '
    write(message(2), '(a)') '      (ii) do nothing special, or (iii) totally screw up the run.'
    write(message(3), '(a)') '      Good luck!'
    call write_info(3)

    allocate(smix%df(m%np, d2, d3, smix%ns), smix%vin_old(m%np, d2, d3), &
             smix%dv(m%np, d2, d3, smix%ns), smix%f_old  (m%np, d2, d3)   )
    smix%df = M_ZERO; smix%dv = M_ZERO; smix%vin_old = M_ZERO; smix%f_old = M_ZERO

  end select

  smix%last_ipos = 0

  call pop_sub()
end subroutine mix_init


! ---------------------------------------------------------
subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix
  call push_sub('mix.mix_end')

  if (associated(smix%df))      deallocate(smix%df)
  if (associated(smix%dv))      deallocate(smix%dv)
  if (associated(smix%vin_old)) deallocate(smix%vin_old)
  if (associated(smix%f_old))   deallocate(smix%f_old)

  call pop_sub()
end subroutine mix_end


! ---------------------------------------------------------
subroutine mixing(smix, m, iter, d2, d3, vin, vout, vnew)
  type(mix_type), intent(inout) :: smix
  type(mesh_type), intent(in)   :: m
  integer,        intent(in)    :: iter, d2, d3
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)

  call push_sub('mix.mixing')

  if (iter.lt.1) then
    message(1) = 'Wrong number of iterations in suboutine mixing.'
    call write_fatal(1)
  end if

  select case (smix%type_of_mixing)
  case (MIX_LINEAR)
    call mixing_linear(smix%alpha, vin, vout, vnew)

  case (MIX_BROYDEN)
    call mixing_broyden(smix, m, d2, d3, vin, vout, vnew, iter)

  case (MIX_GRPULAY)
    call mixing_grpulay(smix, m, d2, d3, vin, vout, vnew, iter)

  end select

  call pop_sub()
end subroutine mixing


! ---------------------------------------------------------
subroutine mixing_linear(alpha, vin, vout, vnew)
  FLOAT,   intent(in)  :: alpha
  FLOAT,   intent(in)  :: vin(:, :, :), vout(:, :, :)
  FLOAT,   intent(out) :: vnew(:, :, :)

  vnew = vin*(M_ONE - alpha) + alpha*vout

end subroutine mixing_linear


! ---------------------------------------------------------
subroutine mixing_broyden(smix, m, d2, d3, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  type(mesh_type), intent(in)   :: m
  integer,        intent(in)    :: d2, d3, iter
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)


  integer :: ipos, iter_used, i, j
  FLOAT :: gamma
  FLOAT, allocatable :: f(:, :, :)

  allocate(f(m%np, d2, d3))

  f(1:m%np, 1:d2, 1:d3) = vout(1:m%np, 1:d2, 1:d3) - vin(1:m%np, 1:d2, 1:d3)
  if(iter > 1) then
    ! Store df and dv from current iteration
    ipos = mod(smix%last_ipos, smix%ns) + 1

    call lalg_copy(m%np, d2, d3, f(:, :, :), smix%df(:, :, :, ipos))
    call lalg_copy(m%np, d2, d3, vin(:, :, :), smix%dv(:, :, :, ipos))
    call lalg_axpy(m%np, d2, d3, -M_ONE, smix%f_old(:, :, :),   smix%df(:, :, :, ipos))
    call lalg_axpy(m%np, d2, d3, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))

    gamma = M_ZERO
    do i = 1, d2
       do j = 1, d3
          gamma = gamma + dmf_integrate(m, smix%df(:, i, j, ipos)*smix%df(:, i, j, ipos))
       enddo
    enddo
    gamma = sqrt(gamma)

    if(gamma > CNST(1e-8)) then
      gamma = M_ONE/gamma
    else
      gamma = M_ONE
    endif
    call lalg_scal(m%np, d2, d3, gamma, smix%df(:, :, :, ipos))
    call lalg_scal(m%np, d2, d3, gamma, smix%dv(:, :, :, ipos))

    smix%last_ipos = ipos
  end if

  ! Store residual and vin for next iteration
  smix%vin_old = vin
  smix%f_old   = f

  ! extrapolate new vector
  iter_used = min(iter - 1, smix%ns)
  call broyden_extrapolation(smix%alpha, m, d2, d3, vin, vnew, iter_used, f, &
     smix%df(1:m%np, 1:d2, 1:d3, 1:iter_used), &
     smix%dv(1:m%np, 1:d2, 1:d3, 1:iter_used))

  deallocate(f)

end subroutine mixing_broyden


! ---------------------------------------------------------
subroutine broyden_extrapolation(alpha, m, d2, d3, vin, vnew, iter_used, f, df, dv)
  FLOAT,   intent(in)         :: alpha
  type(mesh_type), intent(in) :: m
  integer, intent(in)         :: d2, d3, iter_used
  FLOAT,   intent(in)         :: vin(:, :, :), f(:, :, :)
  FLOAT,   intent(in)         :: df(:, :, :, :), dv(:, :, :, :)
  FLOAT,   intent(out)        :: vnew(:, :, :)


  FLOAT, parameter :: w0 = CNST(0.01)

  integer  :: i, j, k, l
  FLOAT :: beta(iter_used, iter_used), gamma, work(iter_used), w(iter_used), x

  if (iter_used == 0) then
    ! linear mixing...
    vnew(1:m%np, 1:d2, 1:d3) = vin(1:m%np, 1:d2, 1:d3) + alpha*f(1:m%np, 1:d2, 1:d3)
    return
  end if

  w  = M_FIVE

  ! compute matrix beta
  beta = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      beta(i, j) = M_ZERO
      do k = 1, d2
         do l = 1, d3
            beta(i, j) = beta(i, j) + w(i)*w(j)*dmf_integrate(m, df(:, k, l, j)*df(:, k, l, i))
         enddo
      enddo
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
           work(i) = work(i) + dmf_integrate(m, df(:, k, l, i)*f(:, k, l))
        enddo
     enddo
  end do


  ! linear mixing term
  vnew(1:m%np, 1:d2, 1:d3) = vin(1:m%np, 1:d2, 1:d3) + alpha*f(1:m%np, 1:d2, 1:d3)

  ! other terms
  do i = 1, iter_used
    gamma = M_ZERO
    do j = 1, iter_used
      gamma = gamma + beta(j, i)*w(j)*work(j)
    end do
    vnew(1:m%np, 1:d2, 1:d3) = vnew(1:m%np, 1:d2, 1:d3) - w(i)*gamma*(alpha*df(1:m%np, 1:d2, 1:d3, i) + &
                                                                          dv(1:m%np, 1:d2, 1:d3, i))
  end do

end subroutine broyden_extrapolation


! ---------------------------------------------------------
! Guaranteed-reduction Pulay
! ---------------------------------------------------------
subroutine mixing_grpulay(smix, m, d2, d3, vin, vout, vnew, iter)
  type(mix_type), intent(inout) :: smix
  type(mesh_type), intent(in)   :: m
  integer,        intent(in)    :: d2, d3
  integer,        intent(in)    :: iter
  FLOAT,          intent(in)    :: vin(:, :, :), vout(:, :, :)
  FLOAT,          intent(out)   :: vnew(:, :, :)

  integer :: ipos, iter_used
  FLOAT, allocatable :: f(:, :, :)

  allocate(f(m%np, d2, d3))
  f = vout - vin

  ! we only extrapolate a new vector every two iterations
  select case (mod(iter, 2_i4))
  case (1)
    ! Store df and dv from current iteration
    if (iter > 1) then
      ipos = smix%last_ipos
      call lalg_copy(m%np, d2, d3, f(:, :, :), smix%df(:, :, :, ipos))
      call lalg_copy(m%np, d2, d3, vin(:, :, :), smix%dv(:, :, :, ipos))
      call lalg_axpy(m%np, d2, d3, -M_ONE, smix%f_old(:, :, :), smix%df(:, :, :, ipos))
      call lalg_axpy(m%np, d2, d3, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))
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
    call lalg_copy(m%np, d2, d3, f(:, :, :), smix%df(:, :, :, ipos))
    call lalg_copy(m%np, d2, d3, vin(:, :, :), smix%dv(:, :, :, ipos))
    call lalg_axpy(m%np, d2, d3, -M_ONE, smix%f_old(:, :, :), smix%df(:, :, :, ipos))
    call lalg_axpy(m%np, d2, d3, -M_ONE, smix%vin_old(:, :, :), smix%dv(:, :, :, ipos))

    smix%last_ipos = ipos

    ! extrapotate new vector
    iter_used = min(iter/2, smix%ns + 1)
    call pulay_extrapolation(m, d2, d3, vin, vout, vnew, iter_used, f, &
                             smix%df(1:m%np, 1:d2, 1:d3, 1:iter_used), &
                             smix%dv(1:m%np, 1:d2, 1:d3, 1:iter_used))
  end select

  deallocate(f)
end subroutine mixing_grpulay


! ---------------------------------------------------------
subroutine pulay_extrapolation(m, d2, d3, vin, vout, vnew, iter_used, f, df, dv)
  type(mesh_type), intent(in) :: m
  integer, intent(in) :: d2, d3
  integer, intent(in)   :: iter_used
  FLOAT, intent(in)  :: vin(:, :, :), vout(:, :, :), f(:, :, :), df(:, :, :, :), dv(:, :, :, :)
  FLOAT, intent(out) :: vnew(:, :, :)

  integer :: i, j, k, l
  FLOAT :: a(iter_used, iter_used), alpha

  ! set matrix A
  a = M_ZERO
  do i = 1, iter_used
    do j = i + 1, iter_used
      a(i, j) = M_ZERO
      do k = 1, d2
         do l = 1, d3
            a(i, j) = a(i, j) + dmf_integrate(m, df(:, k, l, j)*df(:, k, l, i))
         enddo
      enddo
      a(j, i) = a(i, j)
    end do
    a(i, i) = M_ZERO
    do k = 1, d2
       do l = 1, d3
          a(i, i) = a(i, i) + dmf_integrate(m, df(:, k, l, i)*df(:, k, l, i))
       enddo
    enddo
  end do
  if (all(a < 1.0E-8)) then
    ! residuals are too small. Do not mix.
    vnew = vout
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
             alpha = alpha - a(i, j)*dmf_integrate(m, df(:, k, l, j)*f(:, k, l))
          enddo
       enddo
    end do
    vnew(:, :, :) = vnew(:, :, :) + alpha * dv(:, :, :, i)
  end do

end subroutine pulay_extrapolation

end module mix
