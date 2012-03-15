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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module preconditioners_m
  use batch_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_basic_m
  use parser_m
  use mesh_m
  use messages_m
  use multigrid_m
  use nl_operator_m
  use poisson_m
  use profiling_m
  use stencil_star_m
  use simul_box_m
  use varinfo_m

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
    integer :: which

    type(nl_operator_t) :: op
    FLOAT, pointer      :: diag_lapl(:) !< diagonal of the laplacian
  end type preconditioner_t
  
contains

  ! --------------------------------------------------------- 
  subroutine preconditioner_init(this, gr, prefix)
    type(preconditioner_t), intent(out)    :: this 
    type(grid_t),           intent(in)     :: gr
    character(len=*), optional, intent(in) :: prefix

    FLOAT, parameter :: alpha = M_HALF
    FLOAT :: vol
    character(len=256) :: prefix_
    integer :: default
    integer :: maxp, is, ns, ip, ip2
    
    PUSH_SUB(preconditioner_init)

    !%Variable Preconditioner
    !%Type integer
    !%Default filter
    !%Section SCF::Eigensolver
    !%Description
    !% Which preconditioner to use in order to solve the Kohn-Sham equations or
    !% the linear-response equations. May apply prefix of linear-response (<i>e.g.</i>
    !% <tt>EM</tt>, <tt>KdotP</tt>, <tt>VM</tt>) to differentiate from choice for ground state.
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
    prefix_ = ""
    if(present(prefix)) prefix_ = prefix

    if(gr%mesh%use_curvilinear) then
      default = PRE_NONE
    else
      default = PRE_FILTER
    end if

    if (parse_isdef(datasets_check(trim(prefix_)//'Preconditioner')) /= 0 ) then 
      call parse_integer(datasets_check(trim(prefix_)//'Preconditioner'), default, this%which)
      if(.not.varinfo_valid_option('Preconditioner', this%which)) &
        call input_error('Preconditioner')
      call messages_print_var_option(stdout, 'Preconditioner', this%which, prefix_)
    else
      call parse_integer(datasets_check('Preconditioner'), default, this%which)
      if(.not.varinfo_valid_option('Preconditioner', this%which)) &
        call input_error('Preconditioner')
      call messages_print_var_option(stdout, 'Preconditioner', this%which)
    endif  

    select case(this%which)
    case(PRE_FILTER)
      ! the smoothing has a star stencil like the laplacian
      call nl_operator_init(this%op, "Preconditioner")
      call stencil_star_get_lapl(this%op%stencil, gr%mesh%sb%dim, 1)
      call nl_operator_build(gr%mesh, this%op, gr%mesh%np, const_w = .not. gr%mesh%use_curvilinear)
      
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
            this%op%w_re(is, ip) = M_HALF*(M_ONE - alpha)/gr%mesh%sb%dim
          else
            this%op%w_re(is, ip) = alpha
          end if
          ip2 = ip + this%op%ri(is, this%op%rimap(ip))
          if(gr%mesh%use_curvilinear) this%op%w_re(is, ip) = this%op%w_re(is, ip)*(ns*gr%mesh%vol_pp(ip2)/vol)
        end do
      end do
      
      call nl_operator_update_weights(this%op)

    case(PRE_JACOBI, PRE_MULTIGRID)
      SAFE_ALLOCATE(this%diag_lapl(1:gr%mesh%np))
      call derivatives_lapl_diag(gr%der, this%diag_lapl)
      call lalg_scal(gr%mesh%np, -M_HALF, this%diag_lapl(:))
    end select

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
  subroutine preconditioner_obsolete_variables(old_prefix, new_prefix)
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix

    call messages_obsolete_variable(trim(old_prefix)//'Preconditioner', trim(new_prefix)//'Preconditioner')
  end subroutine preconditioner_obsolete_variables

#include "undef.F90"
#include "complex.F90"
#include "preconditioners_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "preconditioners_inc.F90"

end module preconditioners_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
