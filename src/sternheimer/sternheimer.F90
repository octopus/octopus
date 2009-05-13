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
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use loct_parser_m
  use XC_F90(lib_m)
  use linear_solver_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use multigrid_m
  use mix_m
  use h_sys_output_m
  use preconditioners_m
  use poisson_m
  use profiling_m
  use pert_m
  use restart_m
  use simul_box_m
  use scf_tol_m
  use smear_m
  use states_m
  use states_dim_m
  use states_calc_m
  use string_m
  use system_m
  use units_m
  use v_ks_m
  use xc_m
  use xc_OEP_kernel_m

  implicit none

  private
  public ::                       &
       sternheimer_t,             &
       sternheimer_init,          &
       sternheimer_end,           &
       dsternheimer_solve,        & 
       zsternheimer_solve,        &
       sternheimer_add_fxc,       &
       sternheimer_add_hartree,   &
       dsternheimer_calc_hvar,    &
       zsternheimer_calc_hvar,    &
       dsternheimer_seT_rhs,      &
       zsternheimer_seT_rhs,      &
       sternheimer_have_rhs,      &
       sternheimer_unset_rhs,     &
       sternheimer_has_converged
  
  type sternheimer_t
     private
     type(linear_solver_t) :: solver
     type(mix_t)           :: mixer
     type(scf_tol_t)       :: scf_tol
     FLOAT, pointer        :: fxc(:,:,:)    ! linear change of the xc potential (fxc)
     FLOAT, pointer        :: drhs(:, :, :, :)
     CMPLX, pointer        :: zrhs(:, :, :, :)
     logical               :: add_fxc
     logical               :: add_hartree
     logical               :: ok
     logical               :: hmermitian 
     logical               :: occ_response
     logical               :: preorthogonalization
     logical               :: oep_kernel
  end type sternheimer_t
  
  type(profile_t), save :: prof, prof_hvar

contains
  
  !-----------------------------------------------------------
  subroutine sternheimer_init(this, sys, hm, prefix, hermitian, set_ham_var, set_occ_response)
    type(sternheimer_t), intent(out)   :: this
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    character(len=*),    intent(in)    :: prefix
    logical, optional,   intent(in)    :: hermitian
    integer, optional,   intent(in)    :: set_ham_var
    logical, optional,   intent(in)    :: set_occ_response

    integer :: ham_var
    logical :: default

    if(simul_box_is_periodic(sys%gr%mesh%sb)) call messages_devel_version("Sternheimer equation for periodic systems")

    !%Variable Preorthogonalization
    !%Type logical 
    !%Default true
    !%Section Linear Response::Sternheimer 
    !%Description 
    !% Whether initial linear-response wavefunctions should be orthogonalized 
    !% or not against the occupied states, at the start of each SCF cycle.
    !%End 
    default = sys%st%smear%method == SMEAR_SEMICONDUCTOR
    if (loct_parse_isdef(datasets_check(trim(prefix)//'Preorthogonalization')) /= 0) then 
      call loct_parse_logical(datasets_check(trim(prefix)//'Preorthogonalization'), default, this%preorthogonalization) 
    else 
      call loct_parse_logical(datasets_check('Preorthogonalization'), default, this%preorthogonalization) 
    end if

    !%Variable HamiltonianVariation
    !%Type integer
    !%Default hartree+fxc
    !%Section Linear Response::Sternheimer
    !%Description
    !% The terms to be considered in the variation of the
    !% Hamiltonian. V_ext is always considered. The default is to include
    !% also the exchange, correlation and Hartree terms. If you want
    !% to choose the exchange and correlation kernel use the variable
    !% XCKernel. For kdotp and magnetic em_resp modes, the value V_ext_only
    !% is used and this variable is ignored.
    !%Option V_ext_only 0
    !% Neither Hartree nor xc potentials included.
    !%Option hartree 1
    !% The variation of the Hartree potential only.
    !%Option fxc 2
    !% The exchange and correlation kernel, the variation of the
    !% exchange and correlation potential only.
    !%End

    if(present(set_ham_var)) then
      ham_var = set_ham_var
    else if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      if (loct_parse_isdef(datasets_check(trim(prefix)//'HamiltonianVariation')) /= 0) then
        call loct_parse_int(datasets_check(trim(prefix)//'HamiltonianVariation'), 3, ham_var)
      else
        call loct_parse_int(datasets_check('HamiltonianVariation'), 3, ham_var)
      end if
    end if

    if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
      this%add_fxc = ((ham_var / 2) == 1)
      this%add_hartree = (mod(ham_var, 2) == 1)
    else
      this%add_fxc = .false. 
      this%add_hartree = .false.
    end if
    
    if(present(set_occ_response)) then
       this%occ_response = set_occ_response
    else
       this%occ_response = .false.
    endif

    message(1) = "Variation of the Hamiltonian in Sternheimer equation: V_ext"
    if(this%add_hartree) write(message(1), '(2a)') trim(message(1)), ' + hartree'
    if(this%add_fxc)     write(message(1), '(2a)') trim(message(1)), ' + fxc'

    message(2) = "Solving Sternheimer equation for"
    if (this%occ_response) then
       write(message(2), '(2a)') trim(message(2)), ' full linear response'
    else
       write(message(2), '(2a)') trim(message(2)), ' linear response in unoccupied subspace only'
    endif

    message(3) = "Sternheimer preorthogonalization:"
    if (this%preorthogonalization) then
       write(message(3), '(2a)') trim(message(3)), ' yes'
    else
       write(message(3), '(2a)') trim(message(3)), ' no'
    endif
    call write_info(3) 

    call linear_solver_init(this%solver, sys%gr, prefix)

    if(this%solver%solver == LS_MULTIGRID .or. preconditioner_is_multigrid(this%solver%pre)) then
      if(.not. associated(sys%gr%mgrid)) then
        SAFE_ALLOCATE(sys%gr%mgrid)
        call multigrid_init(sys%gr%mgrid, sys%geo, sys%gr%cv, sys%gr%mesh, sys%gr%der, sys%gr%stencil)
      end if
      call hamiltonian_mg_init(hm, sys%gr)
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
  end subroutine sternheimer_init


  !-----------------------------------------------------------
  subroutine sternheimer_end(this)
    type(sternheimer_t), intent(inout) :: this

    call linear_solver_end(this%solver)
    call scf_tol_end(this%scf_tol)

    if (this%add_fxc) then
      SAFE_DEALLOCATE_P(this%fxc)
    end if

  end subroutine sternheimer_end


  !-----------------------------------------------------------
  subroutine sternheimer_build_fxc(this, mesh, st, ks)
    type(sternheimer_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(in)    :: st
    type(v_ks_t),        intent(in)    :: ks

    FLOAT, allocatable :: rho(:, :)

    call push_sub('sternheimer.sternheimer_build_fxc')

    SAFE_ALLOCATE(this%fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
    this%fxc = M_ZERO

    if( iand(ks%xc%kernel_family, XC_FAMILY_OEP) == 0 ) then 
      this%oep_kernel = .false.

      SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
      call states_total_density(st, mesh, rho)
      call xc_get_fxc(ks%xc, mesh, rho, st%d%ispin, this%fxc)
      SAFE_DEALLOCATE_A(rho)
    else
      this%oep_kernel = .true.

      call xc_oep_kernel_init(ks%oep)
    end if

    call pop_sub()

  end subroutine sternheimer_build_fxc


  !-----------------------------------------------------------
  logical function sternheimer_add_fxc(this) result(r)
    type(sternheimer_t), intent(in) :: this
    r = this%add_fxc
  end function sternheimer_add_fxc


  !-----------------------------------------------------------
  logical function sternheimer_add_hartree(this) result(r)
    type(sternheimer_t), intent(in) :: this
    r = this%add_hartree
  end function sternheimer_add_hartree


  !-----------------------------------------------------------
  logical function sternheimer_has_converged(this) result(r)
    type(sternheimer_t), intent(in) :: this
    r = this%ok
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
