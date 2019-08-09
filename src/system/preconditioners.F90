!! Copyright (C) 2006 X. Andrade, M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module preconditioners_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use messages_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use stencil_star_oct_m
  use simul_box_oct_m
  use varinfo_oct_m

  implicit none
  private

  integer, public, parameter ::     &
    PRE_NONE      = 0,              &
    PRE_FILTER    = 1,              &
    PRE_JACOBI    = 2,              &
    PRE_POISSON   = 3,              &
    PRE_MULTIGRID = 7

  public ::                            &
    preconditioner_t,                  &
    preconditioner_init,               &
    preconditioner_null,               &
    preconditioner_end,                &
    preconditioner_is_multigrid,       &
    dpreconditioner_apply,             &
    zpreconditioner_apply,             &
    dpreconditioner_apply_batch,       &
    zpreconditioner_apply_batch,       &
    preconditioner_obsolete_variables

  type preconditioner_t
    private
    integer :: which

    type(nl_operator_t) :: op
    FLOAT, pointer      :: diag_lapl(:) !< diagonal of the laplacian
    integer             :: npre, npost, nmiddle
  end type preconditioner_t

contains

  ! ---------------------------------------------------------
  subroutine preconditioner_init(this, namespace, gr)
    type(preconditioner_t), intent(out)    :: this
    type(namespace_t),      intent(in)     :: namespace
    type(grid_t),           intent(in)     :: gr

    FLOAT :: alpha, default_alpha
    FLOAT :: vol
    integer :: default
    integer :: maxp, is, ns, ip, ip2

    PUSH_SUB(preconditioner_init)

    !%Variable Preconditioner
    !%Type integer
    !%Section SCF::Eigensolver
    !%Description
    !% Which preconditioner to use in order to solve the Kohn-Sham
    !% equations or the linear-response equations. The default is
    !% pre_filter, except for curvilinear coordinates, where no
    !% preconditioner is applied by default.
    !%Option no 0
    !% Do not apply preconditioner.
    !%Option pre_filter 1
    !% Filter preconditioner.
    !%Option pre_jacobi 2
    !% Jacobi preconditioner. Only the local part of the pseudopotential is used.
    !% Not very helpful.
    !%Option pre_poisson 3
    !% Uses the full Laplacian as preconditioner. The inverse is calculated through
    !% the solution of the Poisson equation. This is, of course, very slow.
    !%Option pre_multigrid 7
    !% Multigrid preconditioner.
    !%End

    if(gr%mesh%use_curvilinear) then
      default = PRE_NONE
    else
      default = PRE_FILTER
    end if

    call parse_variable(namespace, 'Preconditioner', default, this%which)
    if(.not.varinfo_valid_option('Preconditioner', this%which)) call messages_input_error('Preconditioner')
    call messages_print_var_option(stdout, 'Preconditioner', this%which)

    select case(this%which)
    case(PRE_FILTER)
      ! the smoothing has a star stencil like the laplacian
      call nl_operator_init(this%op, "Preconditioner")
      call stencil_star_get_lapl(this%op%stencil, gr%mesh%sb%dim, 1)
      call nl_operator_build(gr%mesh, this%op, gr%mesh%np, const_w = .not. gr%mesh%use_curvilinear)

      !%Variable PreconditionerFilterFactor
      !%Type float
      !%Section SCF::Eigensolver
      !%Description
      !% This variable controls how much filter preconditioner is
      !% applied. A value of 1.0 means no preconditioning, 0.5 is the
      !% standard.
      !%
      !% The default is 0.5, except for periodic systems where the
      !% default is 0.6.
      !%
      !% If you observe that the first eigenvectors are not converging
      !% properly, especially for periodic systems, you should
      !% increment this value.
      !%
      !% The allowed range for this parameter is between 0.5 and 1.0.
      !% For other values, the SCF may converge to wrong results.
      !%End
      default_alpha = CNST(0.5)
      if(simul_box_is_periodic(gr%sb)) default_alpha = CNST(0.6)

      call parse_variable(namespace, 'PreconditionerFilterFactor', default_alpha, alpha)

      call messages_print_var_value(stdout, 'PreconditionerFilterFactor', alpha)

      ! check for correct interval of alpha
      if (alpha < CNST(0.5) .or. alpha > CNST(1.0)) then
        call messages_input_error('PreconditionerFilterFactor')
      end if

      ns = this%op%stencil%size

      if (this%op%const_w) then
        maxp = 1
      else
        maxp = gr%mesh%np
      end if

      do ip = 1,maxp

        if(gr%mesh%use_curvilinear) vol = sum(gr%mesh%vol_pp(ip + this%op%ri(1:ns, this%op%rimap(ip))))

        do is = 1, ns
          if(is /= this%op%stencil%center) then
            this%op%w(is, ip) = M_HALF*(M_ONE - alpha)/gr%mesh%sb%dim
          else
            this%op%w(is, ip) = alpha
          end if
          this%op%w(is, ip) = this%op%w(is, ip) * M_TWO * M_TWO*(M_ONE - alpha)/gr%mesh%sb%dim
          ip2 = ip + this%op%ri(is, this%op%rimap(ip))
          if(gr%mesh%use_curvilinear) this%op%w(is, ip) = this%op%w(is, ip)*(ns*gr%mesh%vol_pp(ip2)/vol)
        end do
      end do

      call nl_operator_update_weights(this%op)

    case(PRE_JACOBI, PRE_MULTIGRID)
      SAFE_ALLOCATE(this%diag_lapl(1:gr%mesh%np))
      call derivatives_lapl_diag(gr%der, this%diag_lapl)
      call lalg_scal(gr%mesh%np, -M_HALF, this%diag_lapl(:))
    end select

    if(this%which == PRE_MULTIGRID) then
      !%Variable PreconditionerIterationsPre
      !%Type integer
      !%Section SCF::Eigensolver
      !%Description
      !% This variable is the number of pre-smoothing iterations for the multigrid
      !% preconditioner. The default is 1.
      !%End
      call parse_variable(namespace, 'PreconditionerIterationsPre', 1, this%npre)

      !%Variable PreconditionerIterationsMiddle
      !%Type integer
      !%Section SCF::Eigensolver
      !%Description
      !% This variable is the number of smoothing iterations on the coarsest grid for the multigrid
      !% preconditioner. The default is 1.
      !%End
      call parse_variable(namespace, 'PreconditionerIterationsMiddle', 1, this%nmiddle)

      !%Variable PreconditionerIterationsPost
      !%Type integer
      !%Section SCF::Eigensolver
      !%Description
      !% This variable is the number of post-smoothing iterations for the multigrid
      !% preconditioner. The default is 2.
      !%End
      call parse_variable(namespace, 'PreconditionerIterationsPost', 2, this%npost)
    end if

    POP_SUB(preconditioner_init)
  end subroutine preconditioner_init


  ! ---------------------------------------------------------
  subroutine preconditioner_null(this)
    type(preconditioner_t), intent(inout) :: this

    PUSH_SUB(preconditioner_null)
    this%which = PRE_NONE

    POP_SUB(preconditioner_null)
  end subroutine preconditioner_null


  ! ---------------------------------------------------------
  subroutine preconditioner_end(this)
    type(preconditioner_t), intent(inout) :: this

    PUSH_SUB(preconditioner_end)

    select case(this%which)
    case(PRE_FILTER)
      call nl_operator_end(this%op)

    case(PRE_JACOBI, PRE_MULTIGRID)
      SAFE_DEALLOCATE_P(this%diag_lapl)
    end select

    call preconditioner_null(this)
    POP_SUB(preconditioner_end)
  end subroutine preconditioner_end


  ! ---------------------------------------------------------
  logical pure function preconditioner_is_multigrid(this) result(req)
    type(preconditioner_t), intent(in) :: this

    req = (this%which == PRE_MULTIGRID)
  end function preconditioner_is_multigrid

  ! ---------------------------------------------------------
  subroutine preconditioner_obsolete_variables(namespace, old_prefix, new_prefix)
    type(namespace_t),   intent(in)    :: namespace
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix

    call messages_obsolete_variable(namespace, trim(old_prefix)//'Preconditioner', trim(new_prefix)//'Preconditioner')
  end subroutine preconditioner_obsolete_variables

#include "undef.F90"
#include "complex.F90"
#include "preconditioners_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "preconditioners_inc.F90"

end module preconditioners_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
