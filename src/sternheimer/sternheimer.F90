!! Copyright (C) 2004-2012 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca), David Strubbe
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

module sternheimer_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use density_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use linear_response_oct_m
  use linear_solver_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mix_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pert_oct_m
  use photon_mode_oct_m
  use poisson_oct_m
  use preconditioners_oct_m
  use profiling_oct_m
  use restart_oct_m
  use scf_tol_oct_m
  use smear_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use types_oct_m
  use wfs_elec_oct_m
  use xc_oct_m
  use xc_f03_lib_m

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
       sternheimer_build_kxc,     &  
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
       sternheimer_unset_kxc,     & 
       sternheimer_has_converged, &
       swap_sigma,                &
       wfs_tag_sigma,             &
       sternheimer_obsolete_variables, &
       dcalc_hvar,                &
       zcalc_hvar,                &
       dcalc_kvar,                &
       zcalc_kvar

  character(len=*), public, parameter :: EM_RESP_PHOTONS_DIR = "em_resp_photons/"

  type sternheimer_t
     private
     type(linear_solver_t) :: solver
     type(mix_t)           :: mixer
     type(scf_tol_t)       :: scf_tol
     FLOAT, allocatable    :: fxc(:,:,:)    !< linear change of the XC potential (fxc)
     FLOAT, allocatable    :: kxc(:,:,:,:)  !< quadratic change of the XC potential (kxc)
     FLOAT, pointer        :: drhs(:, :, :, :) => NULL() !< precomputed bare perturbation on RHS
     CMPLX, pointer        :: zrhs(:, :, :, :) => NULL()
     FLOAT, pointer        :: dinhomog(:, :, :, :, :) => NULL() !< fixed inhomogeneous term on RHS
     CMPLX, pointer        :: zinhomog(:, :, :, :, :) => NULL()
     logical               :: add_fxc
     logical               :: add_hartree
     logical               :: ok
     logical               :: occ_response
     logical               :: last_occ_response
     logical               :: occ_response_by_sternheimer
     logical               :: preorthogonalization
     logical, public       :: enable_el_pt_coupling  !< switch on photoncoupling
     FLOAT                 :: domega          !< current frequency for which we solve the freq.-dep. equation
     CMPLX                 :: zomega          !< current frequency for which we solve the freq.-dep. equation
     FLOAT, allocatable, public  :: dphoton_coord_q(:) !< canonical photon coordinate
     CMPLX, allocatable, public  :: zphoton_coord_q(:) !< canonical photon coordinate
     FLOAT                 :: el_pt_eta      !< broadening for photonic subsystem
     FLOAT, allocatable :: omg2_lmda_r(:)
     FLOAT, allocatable :: lambda_dot_r(:)
     type(photon_mode_t)   :: pt_modes
  end type sternheimer_t
  
  type(profile_t), save :: prof, prof_hvar, prof_hvar_photons

contains
  
  !-----------------------------------------------------------
  subroutine sternheimer_init(this, namespace, space, gr, st, hm, xc, mc, wfs_are_cplx, set_ham_var, set_occ_response, &
    set_last_occ_response, occ_response_by_sternheimer)
    type(sternheimer_t),      intent(out)   :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(inout) :: gr
    type(states_elec_t),      intent(in)    :: st
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(xc_t),               intent(in)    :: xc
    type(multicomm_t),        intent(in)    :: mc
    logical,                  intent(in)    :: wfs_are_cplx
    integer,        optional, intent(in)    :: set_ham_var
    logical,        optional, intent(in)    :: set_occ_response
    logical,        optional, intent(in)    :: set_last_occ_response
    logical,        optional, intent(in)    :: occ_response_by_sternheimer

    integer :: ip, nm, ii
    integer :: ham_var, iunit
    logical :: default_preorthog

    PUSH_SUB(sternheimer_init)

    if(st%smear%method  ==  SMEAR_FIXED_OCC) then
      call messages_experimental("Sternheimer equation for arbitrary occupations")
    end if
    if(st%smear%method  ==  SMEAR_SEMICONDUCTOR .and. &
      (abs(st%smear%ef_occ) > M_EPSILON) .and. abs(st%smear%ef_occ - M_ONE) > M_EPSILON) then
      write(message(1),'(a,f12.6)') 'Partial occupation at the Fermi level: ', st%smear%ef_occ
      message(2) = 'Semiconducting smearing cannot be used for Sternheimer in this situation.'
      call messages_fatal(2, namespace=namespace)
    end if

    if(wfs_are_cplx) then
      call mix_init(this%mixer, namespace, space, gr%der, gr%mesh%np, st%d%nspin, 1, func_type_= TYPE_CMPLX)
    else
      call mix_init(this%mixer, namespace, space, gr%der, gr%mesh%np, st%d%nspin, 1, func_type_= TYPE_FLOAT)
    end if

    if(present(set_occ_response)) then
       this%occ_response = set_occ_response
    else
       this%occ_response = .false.
    end if

    this%occ_response_by_sternheimer = optional_default(occ_response_by_sternheimer, .false.)

    !%Variable Preorthogonalization
    !%Type logical 
    !%Section Linear Response::Sternheimer 
    !%Description 
    !% Whether initial linear-response wavefunctions should be orthogonalized 
    !% or not against the occupied states, at the start of each SCF cycle.
    !% Default is true only if <tt>SmearingFunction = semiconducting</tt>,
    !% or if the <tt>Occupations</tt> block specifies all full or empty states,
    !% and we are not solving for linear response in the occupied subspace too.
    !%End 
    default_preorthog = (st%smear%method == SMEAR_SEMICONDUCTOR .or. &
      (st%smear%method == SMEAR_FIXED_OCC .and. st%smear%integral_occs)) &
      .and. .not. this%occ_response
    call parse_variable(namespace, 'Preorthogonalization', default_preorthog, this%preorthogonalization) 

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
    else if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call parse_variable(namespace, 'HamiltonianVariation', 3, ham_var)
    else
      ham_var = 0
    end if

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
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
    end if

    message(1) = "Variation of the Hamiltonian in Sternheimer equation: V_ext"
    if(this%add_hartree) write(message(1), '(2a)') trim(message(1)), ' + hartree'
    if(this%add_fxc)     write(message(1), '(2a)') trim(message(1)), ' + fxc'

    message(2) = "Solving Sternheimer equation for"
    if (this%occ_response) then
       write(message(2), '(2a)') trim(message(2)), ' full linear response.'
    else
       write(message(2), '(2a)') trim(message(2)), ' linear response in unoccupied subspace only.'
    end if

    message(3) = "Sternheimer preorthogonalization:"
    if (this%preorthogonalization) then
       write(message(3), '(2a)') trim(message(3)), ' yes'
    else
       write(message(3), '(2a)') trim(message(3)), ' no'
    end if
    call messages_info(3) 

    call linear_solver_init(this%solver, namespace, gr, states_are_real(st), mc, space)

    ! will not converge for non-self-consistent calculation unless LRTolScheme = fixed
    if (ham_var == 0) then
      call scf_tol_init(this%scf_tol, namespace, st%qtot, tol_scheme = 0) ! fixed
    else
      call scf_tol_init(this%scf_tol, namespace, st%qtot)
    end if

    if(this%add_fxc) call sternheimer_build_fxc(this, namespace, gr%mesh, st, xc)


    ! This variable is documented in xc_oep_init.
    call parse_variable(namespace, 'EnablePhotons', .false., this%enable_el_pt_coupling)
    call messages_print_var_value(stdout, 'EnablePhotons', this%enable_el_pt_coupling)

    if(this%enable_el_pt_coupling) then
      call photon_mode_init(this%pt_modes, namespace, gr%mesh, space%dim, M_ZERO)
      call io_mkdir(EM_RESP_PHOTONS_DIR, namespace)
      iunit = io_open(EM_RESP_PHOTONS_DIR // 'photon_modes', namespace, action='write')
      call photon_mode_write_info(this%pt_modes, iunit)
      SAFE_ALLOCATE(this%zphoton_coord_q(1:this%pt_modes%nmodes))

      nm = this%pt_modes%nmodes
      SAFE_ALLOCATE(this%omg2_lmda_r(1:gr%mesh%np))
      SAFE_ALLOCATE(this%lambda_dot_r(1:gr%mesh%np))
      do ii = 1, nm
        do ip = 1, gr%mesh%np
          this%omg2_lmda_r(ip) = - (this%pt_modes%omega(ii))**2*this%pt_modes%lambda(ii) * &
            sum(this%pt_modes%pol(ii, 1:space%dim)*gr%mesh%x(ip, 1:space%dim))
          this%lambda_dot_r(ip) = this%pt_modes%lambda(ii) * &
            sum(this%pt_modes%pol(ii, 1:space%dim)*gr%mesh%x(ip, 1:space%dim))
        end do
      end do
    end if

    !%Variable PhotonEta
    !%Type float
    !%Default 0.0000367
    !%Section Linear Response::Sternheimer
    !%Description
    !% This variable provides the value for the broadening of the photonic spectra
    !% when the coupling of electrons to photons is enabled in the frequency-dependent Sternheimer equation
    !%End
    call parse_variable(namespace, 'PhotonEta', CNST(0.0000367), this%el_pt_eta, units_inp%energy)
    call messages_print_var_value(stdout, 'PhotonEta', this%el_pt_eta, units_inp%energy)

    POP_SUB(sternheimer_init)
  end subroutine sternheimer_init

  !-----------------------------------------------------------
  subroutine sternheimer_end(this)
    type(sternheimer_t), intent(inout) :: this

    PUSH_SUB(sternheimer_end)

    SAFE_DEALLOCATE_A(this%zphoton_coord_q)
    SAFE_DEALLOCATE_A(this%omg2_lmda_r)
    SAFE_DEALLOCATE_A(this%lambda_dot_r)

    call linear_solver_end(this%solver)
    call scf_tol_end(this%scf_tol)
    call mix_end(this%mixer)

    SAFE_DEALLOCATE_A(this%fxc)

    POP_SUB(sternheimer_end)
  end subroutine sternheimer_end


  !-----------------------------------------------------------
  subroutine sternheimer_build_fxc(this, namespace, mesh, st, xc)
    type(sternheimer_t), intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    type(xc_t),          intent(in)    :: xc

    FLOAT, allocatable :: rho(:, :)

    PUSH_SUB(sternheimer_build_fxc)

    SAFE_ALLOCATE(this%fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
    this%fxc = M_ZERO

    SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
    call states_elec_total_density(st, mesh, rho)
    call xc_get_fxc(xc, mesh, namespace, rho, st%d%ispin, this%fxc)
    SAFE_DEALLOCATE_A(rho)

    POP_SUB(sternheimer_build_fxc)

  end subroutine sternheimer_build_fxc


  !-----------------------------------------------------------
  subroutine sternheimer_build_kxc(this, namespace, mesh, st, xc)
    type(sternheimer_t), intent(inout) :: this
    type(namespace_t),   intent(in)    :: namespace
    type(mesh_t),        intent(in)    :: mesh
    type(states_elec_t), intent(in)    :: st
    type(xc_t),          intent(in)    :: xc

    FLOAT, allocatable :: rho(:, :)

    PUSH_SUB(sternheimer_build_kxc)

    if(this%add_fxc) then
      SAFE_ALLOCATE(this%kxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
      this%kxc = M_ZERO

      SAFE_ALLOCATE(rho(1:mesh%np, 1:st%d%nspin))
      call states_elec_total_density(st, mesh, rho)
      call xc_get_kxc(xc, mesh, namespace, rho, st%d%ispin, this%kxc)
      SAFE_DEALLOCATE_A(rho)
    end if

    POP_SUB(sternheimer_build_kxc)

  end subroutine sternheimer_build_kxc
 
 !---------------------------------------------
  subroutine sternheimer_unset_kxc(this)
    type(sternheimer_t), intent(inout) :: this
    
    PUSH_SUB(sternheimer_unset_kxc)

    SAFE_DEALLOCATE_A(this%kxc)

    POP_SUB(sternheimer_unset_kxc)
  end subroutine sternheimer_unset_kxc

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
    
    PUSH_SUB(sternheimer_unset_rhs)

    nullify(this%drhs)
    nullify(this%zrhs)

    POP_SUB(sternheimer_unset_rhs)
  end subroutine sternheimer_unset_rhs

  !-----------------------------------------------------------
  logical pure function sternheimer_have_inhomog(this) result(have)
    type(sternheimer_t), intent(in) :: this
    have = associated(this%dinhomog) .or. associated(this%zinhomog)
  end function sternheimer_have_inhomog

  !-----------------------------------------------------------
  subroutine sternheimer_unset_inhomog(this)
    type(sternheimer_t), intent(inout) :: this
    
    PUSH_SUB(sternheimer_unset_inhomog)

    nullify(this%dinhomog)
    nullify(this%zinhomog)

    POP_SUB(sternheimer_unset_inhomog)
  end subroutine sternheimer_unset_inhomog

  !-----------------------------------------------------------
  integer pure function swap_sigma(sigma)
    integer, intent(in) :: sigma
    
    if(sigma == 1) then
      swap_sigma = 2
    else
      swap_sigma = 1
    end if

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

  subroutine sternheimer_obsolete_variables(namespace, old_prefix, new_prefix)
    type(namespace_t),   intent(in)    :: namespace
    character(len=*),    intent(in)    :: old_prefix
    character(len=*),    intent(in)    :: new_prefix
    
    PUSH_SUB(sternheimer_obsolete_variables)

    call messages_obsolete_variable(namespace, trim(old_prefix)//'Preorthogonalization', trim(new_prefix)//'Preorthogonalization')
    call messages_obsolete_variable(namespace, trim(old_prefix)//'HamiltonianVariation', trim(new_prefix)//'HamiltonianVariation')

    call linear_solver_obsolete_variables(namespace, old_prefix, new_prefix)
    call scf_tol_obsolete_variables(namespace, old_prefix, new_prefix)

    POP_SUB(sternheimer_obsolete_variables)
  end subroutine sternheimer_obsolete_variables
  
  !--------------------------------------------------------------
  subroutine calc_hvar_photons(this, mesh, st, lr_rho, nsigma, hvar)
    type(sternheimer_t),    intent(inout) :: this
    type(mesh_t),           intent(in)    :: mesh
    type(states_elec_t),    intent(in)    :: st
    integer,                intent(in)    :: nsigma
    CMPLX,                  intent(in)    :: lr_rho(:,:)
    CMPLX,                  intent(inout) :: hvar(:,:,:) !< (1:mesh%np, 1:st%d%nspin, 1:nsigma)

    CMPLX, allocatable :: s_lr_rho(:), vp_dip_self_ener(:), vp_bilinear_el_pt(:)
    CMPLX, allocatable :: first_moments(:)
    integer :: nm, is, ii
    CMPLX :: integral_result

    PUSH_SUB(calc_hvar_photons)
    call profiling_in(prof_hvar_photons, 'CALC_HVAR_PHOTONS')

    nm = this%pt_modes%nmodes

    ! photonic terms
    SAFE_ALLOCATE(s_lr_rho(1:mesh%np))
    SAFE_ALLOCATE(first_moments(1:nm))
    SAFE_ALLOCATE(vp_dip_self_ener(1:mesh%np))
    SAFE_ALLOCATE(vp_bilinear_el_pt(1:mesh%np))

    ! spin summed density
    s_lr_rho = M_ZERO
    do is = 1, st%d%nspin
      s_lr_rho = s_lr_rho + lr_rho(:, is)
    end do

    ! Compute photon q_{\alpha}s and potential for bilinear el-pt coupling
    vp_bilinear_el_pt = M_ZERO
    do ii = 1, nm
      first_moments(ii) = zmf_integrate(mesh, this%omg2_lmda_r(1:mesh%np)*s_lr_rho(1:mesh%np))

      this%zphoton_coord_q(ii) = (M_ONE/(M_TWO*(this%pt_modes%omega(ii))**2)) * &
        ((M_ONE/(this%zomega - this%pt_modes%omega(ii) + M_zI*this%el_pt_eta)) -  &
        (M_ONE/(this%zomega + this%pt_modes%omega(ii) + M_zI*this%el_pt_eta))) * &
        first_moments(ii)

      vp_bilinear_el_pt = vp_bilinear_el_pt - &
        this%pt_modes%omega(ii)*this%lambda_dot_r(1:mesh%np)*this%zphoton_coord_q(ii)
    end do

    ! Compute potential with dipole-self energy contribution
    vp_dip_self_ener = M_ZERO
    do ii = 1, nm
      integral_result = zmf_integrate(mesh, this%lambda_dot_r(1:mesh%np)*s_lr_rho(1:mesh%np))
      vp_dip_self_ener = vp_dip_self_ener + integral_result*this%lambda_dot_r(1:mesh%np)
    end do

    hvar(1:mesh%np, 1, 1) = hvar(1:mesh%np, 1, 1) + vp_dip_self_ener(1:mesh%np) + vp_bilinear_el_pt(1:mesh%np)

    SAFE_DEALLOCATE_A(s_lr_rho)
    SAFE_DEALLOCATE_A(first_moments)
    SAFE_DEALLOCATE_A(vp_dip_self_ener)
    SAFE_DEALLOCATE_A(vp_bilinear_el_pt)

    if (nsigma == 2) hvar(1:mesh%np, 1:st%d%nspin, 2) = conjg(hvar(1:mesh%np, 1:st%d%nspin, 1))

    call profiling_out(prof_hvar_photons)
    POP_SUB(calc_hvar_photons)
  end subroutine calc_hvar_photons


#include "complex.F90"
#include "sternheimer_inc.F90"

#include "undef.F90"

#include "real.F90"
#include "sternheimer_inc.F90"

end module sternheimer_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
