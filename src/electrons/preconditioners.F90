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
  use boundaries_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use stencil_star_oct_m
  use space_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m

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
    dpreconditioner_apply,             &
    zpreconditioner_apply,             &
    dpreconditioner_apply_batch,       &
    zpreconditioner_apply_batch,       &
    preconditioner_obsolete_variables

  type preconditioner_t
    private
    integer :: which

    type(nl_operator_t), allocatable :: op_array(:)  !< this array is necessary for derivatives_get_lapl() to work
    type(nl_operator_t), pointer, public :: op   !< pointer to access op_array(1) simply as op

    FLOAT, allocatable  :: diag_lapl(:) !< diagonal of the laplacian
    integer             :: npre, npost, nmiddle

    type(multigrid_t) :: mgrid  ! multigrid object

    type(derivatives_t), pointer :: der => null()
  end type preconditioner_t

contains

  ! ---------------------------------------------------------
  subroutine preconditioner_init(this, namespace, gr, mc, space)
    type(preconditioner_t), target, intent(out)    :: this
    type(namespace_t),              intent(in)     :: namespace
    type(grid_t), target,           intent(in)     :: gr
    type(multicomm_t),              intent(in)     :: mc
    type(space_t),                  intent(in)     :: space

    FLOAT :: alpha, default_alpha, omega
    FLOAT :: vol
    integer :: default
    integer :: maxp, is, ns, ip, ip2
    character(len=32) :: name

    PUSH_SUB(preconditioner_init)

    SAFE_ALLOCATE(this%op_array(1))
    this%op => this%op_array(1)

    this%der => gr%der

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
    if(.not.varinfo_valid_option('Preconditioner', this%which)) call messages_input_error(namespace, 'Preconditioner')
    call messages_print_var_option(stdout, 'Preconditioner', this%which)

    select case(this%which)
    case(PRE_FILTER)
      ! the smoothing is performed uing the same stencil as the Laplacian
      name = "Preconditioner"
      call derivatives_get_lapl(gr%der, this%op_array, name, 1)

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
      if(space%is_periodic()) default_alpha = CNST(0.6)

      call parse_variable(namespace, 'PreconditionerFilterFactor', default_alpha, alpha)

      call messages_print_var_value(stdout, 'PreconditionerFilterFactor', alpha)

      ! check for correct interval of alpha
      if (alpha < CNST(0.5) .or. alpha > CNST(1.0)) then
        call messages_input_error(namespace, 'PreconditionerFilterFactor')
      end if

      ns = this%op%stencil%size

      if (this%op%const_w) then
        maxp = 1
      else
        maxp = gr%mesh%np
      end if

      !We change the weights to be the one of the kinetic energy operator
      this%op%w = -M_HALF * this%op%w
      
      SAFE_ALLOCATE(this%diag_lapl(1:this%op%np))
      call dnl_operator_operate_diag(this%op, this%diag_lapl)

      do ip = 1, maxp

        if(gr%mesh%use_curvilinear) vol = sum(gr%mesh%vol_pp(ip + this%op%ri(1:ns, this%op%rimap(ip))))

        !The filter preconditioner is given by two iterations of the Relaxation Jacobi method
        !This leads to \tilde{\psi} =  \omega D^{-1}(2\psi - \omega D^{-1} (-0.5\Laplacian) \psi),
        ! where \omega = 2(1-\alpha), and D is the diagonal part of (-0.5\Laplacian)
        !In order to have this to work in all cases, such as different spacings, nonorthogonal cells, ...
        !We directly apply this formula to renormalize the weights of the Laplacian 
        !and get the correct preconditioner

        omega = M_TWO * (M_ONE-alpha)

        do is = 1, ns
          this%op%w(is, ip) = - omega / this%diag_lapl(ip) * this%op%w(is, ip)
          if(is == this%op%stencil%center) then
            this%op%w(is, ip) = this%op%w(is, ip) + M_TWO
          end if
          this%op%w(is, ip) = this%op%w(is, ip) * omega / this%diag_lapl(ip)
         
          if(gr%mesh%use_curvilinear) then
            ip2 = ip + this%op%ri(is, this%op%rimap(ip))
            this%op%w(is, ip) = this%op%w(is, ip)*(ns*gr%mesh%vol_pp(ip2)/vol)
          end if
        end do
      end do

      SAFE_DEALLOCATE_A(this%diag_lapl)

      call nl_operator_output_weights(this%op)

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

      call multigrid_init(this%mgrid, namespace, space, gr%cv, gr%mesh, gr%der, gr%stencil, mc, used_for_preconditioner = .true.)
    end if

    POP_SUB(preconditioner_init)
  end subroutine preconditioner_init


  ! ---------------------------------------------------------
  subroutine preconditioner_null(this)
    type(preconditioner_t), intent(inout) :: this

    PUSH_SUB(preconditioner_null)
    this%which = PRE_NONE

    nullify(this%der)

    POP_SUB(preconditioner_null)
  end subroutine preconditioner_null


  ! ---------------------------------------------------------
  subroutine preconditioner_end(this)
    type(preconditioner_t), intent(inout) :: this

    PUSH_SUB(preconditioner_end)

    select case (this%which)
    case (PRE_FILTER)
      call nl_operator_end(this%op)

    case (PRE_JACOBI)
      SAFE_DEALLOCATE_A(this%diag_lapl)

    case (PRE_MULTIGRID)
      SAFE_DEALLOCATE_A(this%diag_lapl)
      call multigrid_end(this%mgrid)

    end select

    nullify(this%op)
    SAFE_DEALLOCATE_A(this%op_array)

    call preconditioner_null(this)
    POP_SUB(preconditioner_end)
  end subroutine preconditioner_end

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
