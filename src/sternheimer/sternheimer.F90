!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: em_resp.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module sternheimer_m
  use datasets_m
  use density_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use linear_response_m
  use linear_solver_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use multigrid_m
  use pert_m
  use poisson_m
  use preconditioners_m
  use profiling_m
  use restart_m
  use scf_tol_m
  use simul_box_m
  use smear_m
  use states_m
  use states_calc_m
  use states_dim_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use types_m
  use v_ks_m
  use xc_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::                       &
       sternheimer_t,             &
       sternheimer_init,          &
       sternheimer_end,           &
       dsternheimer_solve,        & 
       zsternheimer_solve,        &
       dsternheimer_solve_order2, & 
       zsternheimer_solve_order2, &
       sternheimer_add_fxc,       &
       sternheimer_add_hartree,   &
       dsternheimer_calc_hvar,    &
       zsternheimer_calc_hvar,    &
       dsternheimer_set_rhs,      &
       zsternheimer_set_rhs,      &
       sternheimer_have_rhs,      &
       sternheimer_unset_rhs,     &
       dsternheimer_set_inhomog,  &
       zsternheimer_set_inhomog,  &
       sternheimer_have_inhomog,  &
       sternheimer_unset_inhomog, &
       sternheimer_has_converged, &
       swap_sigma,                &
       wfs_tag_sigma,                     &
       sternheimer_obsolete_variables
  type sternheimer_t
     private
     type(linear_solver_t) :: solver
     type(mix_t)           :: mixer
     type(scf_tol_t)       :: scf_tol
     FLOAT, pointer        :: fxc(:,:,:)    ! linear change of the XC potential (fxc)
     FLOAT, pointer        :: drhs(:, :, :, :) ! precomputed bare perturbation on RHS
     CMPLX, pointer        :: zrhs(:, :, :, :)
     FLOAT, pointer        :: dinhomog(:, :, :, :, :) ! fixed inhomogeneous term on RHS
     CMPLX, pointer        :: zinhomog(:, :, :, :, :)
     logical               :: add_fxc
     logical               :: add_hartree
     logical               :: ok
     logical               :: occ_response
     logical               :: last_occ_response
     logical               :: preorthogonalization
  end type sternheimer_t
  
  type(profile_t), save :: prof, prof_hvar

contains
  
  !-----------------------------------------------------------
  subroutine sternheimer_init(this, sys, hm, prefix, wfs_are_cplx, &
    set_ham_var, set_occ_response, set_last_occ_response, set_default_solver)
    type(sternheimer_t), intent(out)   :: this
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    character(len=*),    intent(in)    :: prefix
    logical,             intent(in)    :: wfs_are_cplx
    integer, optional,   intent(in)    :: set_ham_var
    logical, optional,   intent(in)    :: set_occ_response
    logical, optional,   intent(in)    :: set_last_occ_response
    integer, optional,   intent(in)    :: set_default_solver

    integer :: ham_var, default_solver
    logical :: default_preorthog

    PUSH_SUB(sternheimer_init)

    if(simul_box_is_periodic(sys%gr%mesh%sb)) call messages_experimental("Sternheimer equation for periodic systems")
    if(sys%st%smear%method .eq. SMEAR_FIXED_OCC) then
      call messages_experimental("Sternheimer equation for arbitrary occupations")
    endif

    if(wfs_are_cplx) then
      call mix_init(this%mixer, sys%gr%mesh%np, sys%st%d%nspin, 1, func_type = TYPE_CMPLX)
    else
      call mix_init(this%mixer, sys%gr%mesh%np, sys%st%d%nspin, 1, func_type = TYPE_FLOAT)
    endif

    if(present(set_occ_response)) then
       this%occ_response = set_occ_response
    else
       this%occ_response = .false.
    endif

    !%Variable Preorthogonalization
    !%Type logical 
    !%Section Linear Response::Sternheimer 
    !%Description 
    !% Whether initial linear-response wavefunctions should be orthogonalized 
    !% or not against the occupied states, at the start of each SCF cycle.
    !% Default is true only if <tt>SmearingFunction = semiconducting</tt>,
    !% or if the <tt>Occupations</tt> block specifies all full or empty states,
    !% and we are not solving for linear response in the unoccupied subspace only.
    !%End 
    default_preorthog = (sys%st%smear%method == SMEAR_SEMICONDUCTOR .or. sys%st%smear%integral_occs) &
                        .and. .not. this%occ_response
    if (parse_isdef(datasets_check(trim(prefix)//'Preorthogonalization')) /= 0) then 
      call parse_logical(datasets_check(trim(prefix)//'Preorthogonalization'), default_preorthog, this%preorthogonalization) 
    else 
      call parse_logical(datasets_check('Preorthogonalization'), default_preorthog, this%preorthogonalization) 
    end if

    !%Variable HamiltonianVariation
    !%Type integer
    !%Default hartree+fxc
    !%Section Linear Response::Sternheimer
    !%Description
    !% The terms to be considered in the variation of the
    !% Hamiltonian. The external potential (V_ext) is always considered. The default is to include
    !% also the exchange-correlation and Hartree terms, which fully
    !% takes into account local fields.
    !% Just <tt>hartree</tt> gives you the random-phase approximation (RPA).
    !% If you want to choose the exchange-correlation kernel, use the variable
    !% <tt>XCKernel</tt>. For <tt>kdotp</tt> and magnetic <tt>em_resp</tt> modes,
    !% or if <tt>TheoryLevel = independent_particles</tt>, 
    !% the value <tt>V_ext_only</tt> is used and this variable is ignored.
    !%Option V_ext_only 0
    !% Neither Hartree nor XC potentials included.
    !%Option hartree 1
    !% The variation of the Hartree potential only.
    !%Option fxc 2
    !% The exchange-correlation kernel (the variation of the
    !% exchange-correlation potential) only.
    !%End

    if(present(set_ham_var)) then
      ham_var = set_ham_var
    else if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      if (parse_isdef(datasets_check(trim(prefix)//'HamiltonianVariation')) /= 0) then
        call parse_integer(datasets_check(trim(prefix)//'HamiltonianVariation'), 3, ham_var)
      else
        call parse_integer(datasets_check('HamiltonianVariation'), 3, ham_var)
      end if
    end if

    if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      this%add_fxc = ((ham_var / 2) == 1)
      this%add_hartree = (mod(ham_var, 2) == 1)
    else
      this%add_fxc = .false. 
      this%add_hartree = .false.
    end if
    
    if(present(set_last_occ_response)) then
       this%last_occ_response = set_last_occ_response
    else
       this%last_occ_response = .false.
    endif

    message(1) = "Variation of the Hamiltonian in Sternheimer equation: V_ext"
    if(this%add_hartree) write(message(1), '(2a)') trim(message(1)), ' + hartree'
    if(this%add_fxc)     write(message(1), '(2a)') trim(message(1)), ' + fxc'

    message(2) = "Solving Sternheimer equation for"
    if (this%occ_response) then
       write(message(2), '(2a)') trim(message(2)), ' full linear response.'
    else
       write(message(2), '(2a)') trim(message(2)), ' linear response in unoccupied subspace only.'
    endif

    message(3) = "Sternheimer preorthogonalization:"
    if (this%preorthogonalization) then
       write(message(3), '(2a)') trim(message(3)), ' yes'
    else
       write(message(3), '(2a)') trim(message(3)), ' no'
    endif
    call messages_info(3) 

    if(present(set_default_solver)) then
      default_solver = set_default_solver
    else
      if(conf%devel_version) then
        default_solver = LS_QMR_DOTP
      else
        if(states_are_real(sys%st)) then
          default_solver = LS_QMR_SYMMETRIC
          ! in this case, it is equivalent to LS_QMR_DOTP
        else
          default_solver = LS_QMR_SYMMETRIZED
        endif
      endif
    endif
    call linear_solver_init(this%solver, sys%gr, prefix, default_solver)

    if(this%solver%solver == LS_MULTIGRID .or. preconditioner_is_multigrid(this%solver%pre)) then
      if(.not. associated(sys%gr%mgrid)) then
        SAFE_ALLOCATE(sys%gr%mgrid)
        call multigrid_init(sys%gr%mgrid, sys%geo, sys%gr%cv, sys%gr%mesh, sys%gr%der, sys%gr%stencil)
      end if
    end if

    ! will not converge for non-self-consistent calculation unless LRTolScheme = fixed
    if (ham_var == 0) then
      call scf_tol_init(this%scf_tol, prefix, tol_scheme = 0) ! fixed
    else
      call scf_tol_init(this%scf_tol, prefix)
    end if

    if(this%add_fxc) call sternheimer_build_fxc(this, sys%gr%mesh, sys%st, sys%ks)

    nullify(this%drhs)
    nullify(this%zrhs)
    nullify(this%dinhomog)
    nullify(this%zinhomog)
    POP_SUB(sternheimer_init)
  end subroutine sternheimer_init


  !-----------------------------------------------------------
  subroutine sternheimer_end(this)
    type(sternheimer_t), intent(inout) :: this

    PUSH_SUB(sternheimer_end)

    call linear_solver_end(this%solver)
    call scf_tol_end(this%scf_tol)
    call mix_end(this%mixer)

    if (this%add_fxc) then
      SAFE_DEALLOCATE_P(this%fxc)
    end if

    POP_SUB(sternheimer_end)
  end subroutine sternheimer_end


  !-----------------------------------------------------------
  subroutine sternheimer_build_fxc(this, mesh, st, ks)
    type(sternheimer_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    type(v_ks_t),        intent(in)    :: ks

    FLOAT, allocatable :: rho(:, :)

    PUSH_SUB(sternheimer_build_fxc)

    SAFE_ALLOCATE(this%fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
    this%fxc = M_ZERO

    SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
    call states_total_density(st, mesh, rho)
    call xc_get_fxc(ks%xc, mesh, rho, st%d%ispin, this%fxc)
    SAFE_DEALLOCATE_A(rho)

    POP_SUB(sternheimer_build_fxc)

  end subroutine sternheimer_build_fxc


  !-----------------------------------------------------------
  logical function sternheimer_add_fxc(this) result(rr)
    type(sternheimer_t), intent(in) :: this
    rr = this%add_fxc
  end function sternheimer_add_fxc


  !-----------------------------------------------------------
  logical function sternheimer_add_hartree(this) result(rr)
    type(sternheimer_t), intent(in) :: this
    rr = this%add_hartree
  end function sternheimer_add_hartree


  !-----------------------------------------------------------
  logical function sternheimer_has_converged(this) result(rr)
    type(sternheimer_t), intent(in) :: this
    rr = this%ok
  end function sternheimer_has_converged

  !-----------------------------------------------------------
  logical pure function sternheimer_have_rhs(this) result(have)
    type(sternheimer_t), intent(in) :: this
    have = associated(this%drhs) .or. associated(this%zrhs)
  end function sternheimer_have_rhs

  !-----------------------------------------------------------
  subroutine sternheimer_unset_rhs(this)
    type(sternheimer_t), intent(inout) :: this
    
    nullify(this%drhs)
    nullify(this%zrhs)
  end subroutine sternheimer_unset_rhs

  !-----------------------------------------------------------
  logical pure function sternheimer_have_inhomog(this) result(have)
    type(sternheimer_t), intent(in) :: this
    have = associated(this%dinhomog) .or. associated(this%zinhomog)
  end function sternheimer_have_inhomog

  !-----------------------------------------------------------
  subroutine sternheimer_unset_inhomog(this)
    type(sternheimer_t), intent(inout) :: this
    
    nullify(this%dinhomog)
    nullify(this%zinhomog)
  end subroutine sternheimer_unset_inhomog

  !-----------------------------------------------------------
  integer pure function swap_sigma(sigma)
    integer, intent(in) :: sigma
    
    if(sigma == 1) then
      swap_sigma = 2
    else
      swap_sigma = 1
    endif

  end function swap_sigma

! ---------------------------------------------------------
  character(len=100) function wfs_tag_sigma(base_name, isigma) result(str)
    character(len=*), intent(in) :: base_name
    integer,          intent(in) :: isigma

    character :: sigma_char

    PUSH_SUB(wfs_tag_sigma)

    select case(isigma)
    case(1)
      sigma_char = '+'
    case(2)
      sigma_char = '-'
    case default 
      write(message(1),'(a,i2)') "Illegal integer isigma passed to wfs_tag_sigma: ", isigma
      call messages_fatal(1)
    end select

    str = trim(base_name) // sigma_char

    POP_SUB(wfs_tag_sigma)

  end function wfs_tag_sigma

  ! --------------------------------------------------------

  subroutine sternheimer_obsolete_variables(old_prefix, new_prefix)
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix
    
    call messages_obsolete_variable(trim(old_prefix)//'Preorthogonalization', trim(new_prefix)//'Preorthogonalization')
    call messages_obsolete_variable(trim(old_prefix)//'HamiltonianVariation', trim(new_prefix)//'HamiltonianVariation')

    call linear_solver_obsolete_variables(old_prefix, new_prefix)
    call scf_tol_obsolete_variables(old_prefix, new_prefix)

  end subroutine sternheimer_obsolete_variables
  
#include "complex.F90"
#include "sternheimer_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "sternheimer_inc.F90"

end module sternheimer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
