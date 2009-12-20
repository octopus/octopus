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

#include "global.h"

module poisson_corrections_m
  use derivatives_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::                   &
    poisson_corrections_init, &
    poisson_corrections_end,  &
    correct_rho,              &
    boundary_conditions,      &
    poisson_corr_t

  public ::                   &
    internal_laplacian_op,    &
    internal_dotp,            &
    der_pointer,              &
    mesh_pointer

  FLOAT, parameter :: alpha_ = M_FIVE

  type poisson_corr_t 
     integer :: maxl
     FLOAT, pointer :: phi(:, :)
     FLOAT, pointer :: aux(:, :)
     FLOAT, pointer :: gaussian(:)
  end type poisson_corr_t

  type(derivatives_t), pointer :: der_pointer
  type(mesh_t),        pointer :: mesh_pointer
     
contains

  ! ---------------------------------------------------------
  subroutine poisson_corrections_init(this, ml, m)
    type(poisson_corr_t), intent(inout) :: this
    integer, intent(in)      :: ml
    type(mesh_t), intent(in) :: m

    FLOAT :: alpha, gamma, ylm, gylm(1:MAX_DIM), r, x(MAX_DIM)
    integer :: i, l, add_lm, lldfac, j, mm

    this%maxl = ml

    add_lm = 0
    do l = 0, this%maxl
      do mm = -l, l
        add_lm = add_lm + 1
      end do
    end do
    
    SAFE_ALLOCATE(this%phi(1:m%np, 1:add_lm))
    SAFE_ALLOCATE(this%aux(1:m%np, 1:add_lm))
    SAFE_ALLOCATE(this%gaussian(1:m%np))

    alpha = alpha_*m%spacing(1)
    do i = 1, m%np
      call mesh_r(m, i, r, x = x)
      this%gaussian(i) = exp(-(r/alpha)**2)
      add_lm = 1
      do l = 0, this%maxl
        lldfac = 1; do j = 1, 2*l+1, 2; lldfac = lldfac * j; end do
        gamma = ( sqrt(M_PI)*2**(l+3) ) / lldfac
        do mm = -l, l
          call grylmr(x(1), x(2), x(3), l, mm, ylm, gylm)
          if(r > M_EPSILON) then
            this%phi(i, add_lm) = gamma*isubl(l, r/alpha)*ylm/r**(l+1)
            this%aux(i, add_lm) = r**l*ylm
          else
            this%phi(i, add_lm) = gamma*ylm / alpha
            if(l == 0) then
              this%aux(i, add_lm) = ylm
            else
              this%aux(i, add_lm) = M_ZERO
            end if
          end if
          add_lm = add_lm + 1
        end do
      end do
    end do

  contains

    ! ---------------------------------------------------------
    FLOAT function isubl( l, x)
      integer, intent(in) :: l
      FLOAT,   intent(in) :: x

      isubl = M_HALF*loct_gamma(l + M_HALF)*(M_ONE - loct_incomplete_gamma(l+M_HALF, x**2) )
    end function isubl

  end subroutine poisson_corrections_init


  ! ---------------------------------------------------------
  subroutine poisson_corrections_end(this)
    type(poisson_corr_t), intent(inout) :: this

    SAFE_DEALLOCATE_P(this%phi)
    SAFE_DEALLOCATE_P(this%aux)
    SAFE_DEALLOCATE_P(this%gaussian)
  end subroutine poisson_corrections_end


  ! ---------------------------------------------------------
  subroutine correct_rho(this, m, rho, rho_corrected, vh_correction)
    type(poisson_corr_t), intent(in) :: this
    type(mesh_t), intent(in)  :: m
    FLOAT,        intent(in)  :: rho(:)
    FLOAT,        intent(out) :: rho_corrected(:)
    FLOAT,        intent(out) :: vh_correction(:)

    integer :: i, add_lm, l, mm, lldfac, j
    FLOAT   :: alpha
    FLOAT, allocatable :: mult(:)
    FLOAT, allocatable :: betal(:)
    type(profile_t), save :: prof

    call push_sub('poisson_corrections.correct_rho')
    call profiling_in(prof, "POISSON_CORRECT")

    SAFE_ALLOCATE(mult(1:(this%maxl+1)**2))
    call get_multipoles(this, m, rho, this%maxl, mult)

    alpha = alpha_*m%spacing(1)

    SAFE_ALLOCATE(betal(1:(this%maxl+1)**2))
    add_lm = 1
    do l = 0, this%maxl
      do mm = -l, l
        lldfac = 1
        do j = 1, 2*l+1, 2
          lldfac = lldfac*j
        end do
        betal(add_lm) = (2**(l+2))/( alpha**(2*l+3) * sqrt(M_PI) * lldfac )
        add_lm = add_lm + 1
      end do
    end do

    call lalg_copy(m%np, rho, rho_corrected)
    vh_correction = M_ZERO
    add_lm = 1
    do l = 0, this%maxl
      do mm = -l, l
        forall(i = 1:m%np)
          rho_corrected(i) = rho_corrected(i) - mult(add_lm)*betal(add_lm)*this%aux(i, add_lm)*this%gaussian(i)
        end forall
        call lalg_axpy(m%np, mult(add_lm), this%phi(:, add_lm), vh_correction)
        add_lm = add_lm + 1
      end do
    end do

    SAFE_DEALLOCATE_A(mult)
    SAFE_DEALLOCATE_A(betal)

    call profiling_out(prof)
    call pop_sub()
  end subroutine correct_rho


  ! ---------------------------------------------------------
  subroutine get_multipoles(this, m, rho, ml, mult)
    type(poisson_corr_t), intent(in) :: this
    type(mesh_t), intent(in)  :: m
    FLOAT,           intent(in)  :: rho(:)  ! rho(m%np)
    integer,         intent(in)  :: ml
    FLOAT,           intent(out) :: mult((ml+1)**2)

    integer :: add_lm, l, mm

    call push_sub('poisson_corrections.get_multipoles')

    mult(:) = M_ZERO
    add_lm = 1
    do l = 0, ml
      do mm = -l, l
        mult(add_lm) = dmf_dotp(m, rho, this%aux(:, add_lm))
        add_lm = add_lm + 1
      end do
    end do

    call pop_sub()
  end subroutine get_multipoles

  ! ---------------------------------------------------------
  subroutine internal_laplacian_op(x, y)
    FLOAT, intent(inout) :: x(:)
    FLOAT, intent(out)   :: y(:)

    call push_sub('poisson_corrections.internal_laplacian_op')
    call dderivatives_lapl(der_pointer, x, y)
    call pop_sub()

  end subroutine internal_laplacian_op


  ! ---------------------------------------------------------
  FLOAT function internal_dotp(x, y) result(res)
    FLOAT, intent(inout) :: x(:)
    FLOAT, intent(in)    :: y(:)

    res = dmf_dotp(mesh_pointer, x, y)
  end function internal_dotp


  ! ---------------------------------------------------------
  subroutine boundary_conditions(this, m, rho, pot)
    type(poisson_corr_t), intent(in) :: this    
    type(mesh_t), intent(in)  :: m
    FLOAT,        intent(in)  :: rho(:)  ! rho(m%np)
    FLOAT,        intent(inout) :: pot(:)  ! pot(m%np_part)

    integer :: i, add_lm, l, mm, bp_lower
    FLOAT   :: x(MAX_DIM), r, s1, sa, gylm(1:MAX_DIM)
    FLOAT, allocatable :: mult(:)

    SAFE_ALLOCATE(mult(1:(this%maxl+1)**2))

    call get_multipoles(this, m, rho, this%maxl, mult)

    bp_lower = m%np + 1
    if(m%parallel_in_domains) bp_lower = bp_lower + m%vp%np_ghost(m%vp%partno)

    pot(bp_lower:m%np_part) = M_ZERO
    do i = bp_lower, m%np_part ! boundary conditions
      call mesh_r(m, i, r, x=x)
      add_lm = 1
      do l = 0, this%maxl
        s1 = M_FOUR*M_PI/((M_TWO*l + M_ONE)*r**(l + 1))
        do mm = -l, l
          call grylmr(x(1), x(2), x(3), l, mm, sa, gylm)
          pot(i) = pot(i) + sa * mult(add_lm) * s1
          add_lm = add_lm+1
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(mult)
  end subroutine boundary_conditions


end module poisson_corrections_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
