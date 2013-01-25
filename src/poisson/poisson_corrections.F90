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
  use datasets_m
  use derivatives_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use nl_operator_m
  use parser_m
  use profiling_m
  use simul_box_m

  implicit none

  private
  public ::                   &
    poisson_corrections_init, &
    poisson_corrections_end,  &
    correct_rho,              &
    poisson_boundary_conditions, &
    poisson_corr_t,           &
    internal_laplacian_op,    &
    internal_dotp,            &
    der_pointer,              &
    mesh_pointer

  FLOAT, parameter :: alpha_ = M_FIVE

  type poisson_corr_t
    integer :: method
    integer :: maxl
    FLOAT, pointer :: phi(:, :)
    FLOAT, pointer :: aux(:, :)
    FLOAT, pointer :: gaussian(:)
  end type poisson_corr_t

  type(derivatives_t), pointer :: der_pointer
  type(mesh_t),        pointer :: mesh_pointer
  
  integer, parameter  ::     &
    CORR_MULTIPOLE = 1,     &
    CORR_EXACT     = 3

contains

  ! ---------------------------------------------------------
  subroutine poisson_corrections_init(this, ml, mesh)
    type(poisson_corr_t), intent(inout) :: this
    integer, intent(in)      :: ml
    type(mesh_t), intent(in) :: mesh

    FLOAT :: alpha, gamma, ylm, gylm(1:MAX_DIM), rr, xx(MAX_DIM)
    integer :: ip, ll, add_lm, lldfac, jj, mm

    PUSH_SUB(poisson_corrections_init)

    !%Variable PoissonSolverBoundaries
    !%Type integer
    !%Section Hamiltonian::Poisson
    !%Default multipole
    !%Description
    !% For finite systems, some Poisson solvers (<tt>multigrid</tt>,
    !% <tt>cg_corrected</tt>, and <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>)
    !% require the calculation of the
    !% boundary conditions with an auxiliary method. This variable selects that method.
    !%Option multipole 1
    !% A multipole expansion of the density is used to approximate the potential on the boundaries.
    !%Option exact 3
    !% An exact integration of the Poisson equation is done over the boundaries. This option is
    !% experimental, and not implemented for domain parallelization.
    !%End
    call parse_integer(datasets_check('PoissonSolverBoundaries'), CORR_MULTIPOLE, this%method)

    select case(this%method)
    case(CORR_MULTIPOLE)
      this%maxl = ml

      add_lm = 0
      do ll = 0, this%maxl
        do mm = -ll, ll
          add_lm = add_lm + 1
        end do
      end do

      SAFE_ALLOCATE(this%phi(1:mesh%np, 1:add_lm))
      SAFE_ALLOCATE(this%aux(1:mesh%np, 1:add_lm))
      SAFE_ALLOCATE(this%gaussian(1:mesh%np))

      alpha = alpha_ * mesh%spacing(1)
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, coords = xx)
        this%gaussian(ip) = exp(-(rr/alpha)**2)
        add_lm = 1
        do ll = 0, this%maxl
          lldfac = 1
          do jj = 1, 2*ll+1, 2
            lldfac = lldfac * jj
          end do
          gamma = ( sqrt(M_PI)*2**(ll+3) ) / lldfac
          do mm = -ll, ll
            call grylmr(xx(1), xx(2), xx(3), ll, mm, ylm, gylm)
            if(rr > M_EPSILON) then
              this%phi(ip, add_lm) = gamma*isubl(ll, rr/alpha)*ylm/rr**(ll+1)
              this%aux(ip, add_lm) = rr**ll*ylm
            else
              this%phi(ip, add_lm) = gamma*ylm / alpha
              if(ll == 0) then
                this%aux(ip, add_lm) = ylm
              else
                this%aux(ip, add_lm) = M_ZERO
              end if
            end if
            add_lm = add_lm + 1
          end do
        end do
      end do

    case(CORR_EXACT)
      call messages_experimental('Exact Poisson solver boundaries')
      if(mesh%parallel_in_domains) call messages_not_implemented('Exact Poisson solver boundaries with domain parallelization')

    end select

    POP_SUB(poisson_corrections_init)

  contains

    ! ---------------------------------------------------------
    FLOAT function isubl( ll, xx)
      integer, intent(in) :: ll
      FLOAT,   intent(in) :: xx

      ! no push_sub, called too frequently
      isubl = M_HALF*loct_gamma(ll + M_HALF)*(M_ONE - loct_incomplete_gamma(ll+M_HALF, xx**2) )

    end function isubl

  end subroutine poisson_corrections_init


  ! ---------------------------------------------------------
  subroutine poisson_corrections_end(this)
    type(poisson_corr_t), intent(inout) :: this

    PUSH_SUB(poisson_corrections_end)

    select case(this%method)
    case(CORR_MULTIPOLE)
      SAFE_DEALLOCATE_P(this%phi)
      SAFE_DEALLOCATE_P(this%aux)
      SAFE_DEALLOCATE_P(this%gaussian)
    case(CORR_EXACT)
    end select
    
    POP_SUB(poisson_corrections_end)
  end subroutine poisson_corrections_end


  ! ---------------------------------------------------------
  subroutine correct_rho(this, der, rho, rho_corrected, vh_correction)
    type(poisson_corr_t), intent(in)  :: this
    type(derivatives_t),  intent(in)  :: der
    FLOAT,                intent(in)  :: rho(:)
    FLOAT,                intent(out) :: rho_corrected(:)
    FLOAT,                intent(out) :: vh_correction(:)

    integer :: ip, add_lm, ll, mm, lldfac, jj, ip2
    FLOAT   :: alpha, vv, rr
    FLOAT, allocatable :: mult(:)
    FLOAT, allocatable :: betal(:)
    type(profile_t), save :: prof

    PUSH_SUB(correct_rho)
    call profiling_in(prof, "POISSON_CORRECT")

    ASSERT(ubound(vh_correction, dim = 1) == der%mesh%np_part)

    select case(this%method)
    case(CORR_MULTIPOLE)

      SAFE_ALLOCATE(mult(1:(this%maxl+1)**2))
      call get_multipoles(this, der%mesh, rho, this%maxl, mult)

      alpha = alpha_ * der%mesh%spacing(1)

      SAFE_ALLOCATE(betal(1:(this%maxl+1)**2))
      add_lm = 1
      do ll = 0, this%maxl
        do mm = -ll, ll
          lldfac = 1
          do jj = 1, 2*ll+1, 2
            lldfac = lldfac*jj
          end do
          betal(add_lm) = (2**(ll+2))/( alpha**(2*ll+3) * sqrt(M_PI) * lldfac )
          add_lm = add_lm + 1
        end do
      end do

      call lalg_copy(der%mesh%np, rho, rho_corrected)
      vh_correction = M_ZERO
      add_lm = 1
      do ll = 0, this%maxl
        do mm = -ll, ll
          forall(ip = 1:der%mesh%np)
            rho_corrected(ip) = rho_corrected(ip) - mult(add_lm)*betal(add_lm)*this%aux(ip, add_lm)*this%gaussian(ip)
          end forall
          call lalg_axpy(der%mesh%np, mult(add_lm), this%phi(:, add_lm), vh_correction)
          add_lm = add_lm + 1
        end do
      end do

      SAFE_DEALLOCATE_A(mult)
      SAFE_DEALLOCATE_A(betal)

    case(CORR_EXACT)

      forall(ip = 1:der%mesh%np) vh_correction(ip) = M_ZERO
      
      do ip = der%mesh%np + 1, der%mesh%np_part
        vv = M_ZERO
        do ip2 = 1, der%mesh%np
          rr = sqrt(sum((der%mesh%x(ip, 1:MAX_DIM) - der%mesh%x(ip2, 1:MAX_DIM))**2))
          vv = vv + rho(ip2)/rr
        end do
        vh_correction(ip) = der%mesh%volume_element*vv
      end do

      ASSERT(.not. nl_operator_compact_boundaries(der%lapl))

      call dderivatives_lapl(der, vh_correction, rho_corrected, set_bc = .false.)
 
      forall(ip = 1:der%mesh%np) rho_corrected(ip) = rho(ip) + CNST(1.0)/(CNST(4.0)*M_PI)*rho_corrected(ip)

    end select

    call profiling_out(prof)
    POP_SUB(correct_rho)
  end subroutine correct_rho


  ! ---------------------------------------------------------
  subroutine get_multipoles(this, mesh, rho, ml, mult)
    type(poisson_corr_t), intent(in)  :: this
    type(mesh_t),         intent(in)  :: mesh
    FLOAT,                intent(in)  :: rho(:)  ! rho(mesh%np)
    integer,              intent(in)  :: ml
    FLOAT,                intent(out) :: mult((ml+1)**2)

    integer :: add_lm, ll, mm

    PUSH_SUB(get_multipoles)

    mult(:) = M_ZERO
    add_lm = 1
    do ll = 0, ml
      do mm = -ll, ll
        mult(add_lm) = dmf_dotp(mesh, rho, this%aux(:, add_lm))
        add_lm = add_lm + 1
      end do
    end do

    POP_SUB(get_multipoles)
  end subroutine get_multipoles

  ! ---------------------------------------------------------
  subroutine internal_laplacian_op(xx, yy)
    FLOAT, intent(inout) :: xx(:)
    FLOAT, intent(out)   :: yy(:)

    PUSH_SUB(internal_laplacian_op)
    call dderivatives_lapl(der_pointer, xx, yy)
    POP_SUB(internal_laplacian_op)

  end subroutine internal_laplacian_op


  ! ---------------------------------------------------------
  FLOAT function internal_dotp(xx, yy) result(res)
    FLOAT, intent(inout) :: xx(:)
    FLOAT, intent(in)    :: yy(:)

    PUSH_SUB(internal_dotp)

    res = dmf_dotp(mesh_pointer, xx, yy)
    POP_SUB(internal_dotp)
  end function internal_dotp


  ! ---------------------------------------------------------
  subroutine poisson_boundary_conditions(this, mesh, rho, pot)
    type(poisson_corr_t), intent(in)    :: this    
    type(mesh_t),         intent(in)    :: mesh
    FLOAT,                intent(in)    :: rho(:)  !< rho(mesh%np)
    FLOAT,                intent(inout) :: pot(:)  !< pot(mesh%np_part)

    integer :: ip, add_lm, ll, mm, bp_lower
    FLOAT   :: xx(MAX_DIM), rr, s1, sa, gylm(1:MAX_DIM)
    FLOAT, allocatable :: mult(:)

    PUSH_SUB(poisson_boundary_conditions)

    SAFE_ALLOCATE(mult(1:(this%maxl+1)**2))

    call get_multipoles(this, mesh, rho, this%maxl, mult)

    bp_lower = mesh%np + 1
    if(mesh%parallel_in_domains) then
      bp_lower = bp_lower + mesh%vp%np_ghost(mesh%vp%partno)
    endif

    pot(bp_lower:mesh%np_part) = M_ZERO
    do ip = bp_lower, mesh%np_part ! boundary conditions
      call mesh_r(mesh, ip, rr, coords = xx)
      add_lm = 1
      do ll = 0, this%maxl
        s1 = M_FOUR*M_PI/((M_TWO*ll + M_ONE)*rr**(ll + 1))
        do mm = -ll, ll
          call grylmr(xx(1), xx(2), xx(3), ll, mm, sa, gylm)
          pot(ip) = pot(ip) + sa * mult(add_lm) * s1
          add_lm = add_lm+1
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(mult)
    POP_SUB(poisson_boundary_conditions)
  end subroutine poisson_boundary_conditions

end module poisson_corrections_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
