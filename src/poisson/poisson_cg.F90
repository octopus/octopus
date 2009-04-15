!! Copyright (C) 2004-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module poisson_cg_m
  use derivatives_m
  use global_m
  use math_m
  use mesh_m
  use messages_m
  use poisson_corrections_m
  use profiling_m
  use solvers_m

  implicit none

  private
  public ::           &
    poisson_cg_init,  &
    poisson_cg_end,   &
    poisson_cg1,      &
    poisson_cg2

  FLOAT, public :: threshold

contains


  ! ---------------------------------------------------------
  subroutine poisson_cg_init(m, ml, thr)
    type(mesh_t), intent(in) :: m
    integer, intent(in) :: ml
    FLOAT,   intent(in) :: thr

    call push_sub('poisson_cg.poisson_cg1_init')
    threshold = thr
    call pop_sub()
  end subroutine poisson_cg_init


  ! ---------------------------------------------------------
  subroutine poisson_cg_end

  end subroutine poisson_cg_end


  ! ---------------------------------------------------------
  subroutine poisson_cg1(m, corrector, der, pot, rho)
    type(mesh_t),      target, intent(in)    :: m
    type(poisson_corr_t), intent(inout) :: corrector
    type(derivatives_t), target, intent(in)    :: der
    FLOAT,                     intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                     intent(in)    :: rho(:) ! rho(m%np)


    integer :: iter
    FLOAT :: res
    FLOAT, allocatable :: wk(:), lwk(:), zk(:), pk(:)

    call push_sub('poisson_cg.poisson_cg1')

    ALLOCATE( wk(m%np_part), m%np_part)
    ALLOCATE(lwk(m%np_part), m%np_part)
    ALLOCATE( zk(m%np_part), m%np_part)
    ALLOCATE( pk(m%np_part), m%np_part)

    ! build initial guess for the potential
    wk(1:m%np) = pot(1:m%np)
    call boundary_conditions(corrector, m, rho, wk)
    call dderivatives_lapl(der, wk, lwk, .true.)

    zk(1:m%np) = -M_FOUR*M_PI*rho(1:m%np) - lwk(1:m%np)
    SAFE_DEALLOCATE_A(wk)
    SAFE_DEALLOCATE_A(lwk) ! they are no longer needed

    der_pointer  => der
    mesh_pointer => m
    pk = zk
    iter = 400
    call dconjugate_gradients(m%np_part, pk, zk, &
      internal_laplacian_op, internal_dotp, iter, res, threshold)
    if(res >= threshold) then
      message(1) = 'Conjugate gradients Poisson solver did not converge.'
      write(message(2), '(a,i8)')    '  Iter = ',iter
      write(message(3), '(a,e14.6)') '  Res = ', res
      call write_warning(3)
    end if
    nullify(der_pointer, mesh_pointer)
    pot(1:m%np) = pot(1:m%np) + pk(1:m%np)

    SAFE_DEALLOCATE_A(zk)
    SAFE_DEALLOCATE_A(pk)
    call pop_sub()
  end subroutine poisson_cg1


  ! ---------------------------------------------------------
  subroutine poisson_cg2(m, der, pot, rho)
    implicit none
    type(mesh_t), target,      intent(in) :: m
    type(derivatives_t), target, intent(in) :: der
    FLOAT,                  intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                     intent(in) :: rho(:) ! rho(m%np)

    integer :: iter
    FLOAT, allocatable :: rhs(:), x(:)
    FLOAT :: res

    call push_sub('poisson_cg.poisson_cg2')

    iter = 400
    der_pointer  => der
    mesh_pointer => m

    ALLOCATE(rhs(m%np_part), m%np_part)
    ALLOCATE(x(m%np_part), m%np_part)

    rhs(1:m%np)         = - M_FOUR*M_PI*rho(1:m%np)
    rhs(m%np:m%np_part) = M_ZERO

    x(1:m%np)           = pot(1:m%np)
    x(m%np:m%np_part)   = M_ZERO
    
    call dconjugate_gradients(m%np_part, x, rhs, &
      internal_laplacian_op, internal_dotp, iter, res, threshold)

    if(res >= threshold) then
      message(1) = 'Conjugate gradients Poisson solver did not converge.'
      write(message(2), '(a,i8)')    '  Iter = ',iter
      write(message(3), '(a,e14.6)') '  Res = ', res
      call write_warning(3)
    end if
    pot(1:m%np) = x(1:m%np)

    nullify(der_pointer, mesh_pointer)
    SAFE_DEALLOCATE_A(rhs)
    SAFE_DEALLOCATE_A(x)
    call pop_sub()
  end subroutine poisson_cg2



end module poisson_cg_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
