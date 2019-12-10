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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module propagator_etrs_oct_m
  use accel_oct_m
  use batch_oct_m
  use density_oct_m
  use exponential_oct_m
  use grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_mxll_oct_m
  use ion_dynamics_oct_m
  use lalg_basic_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use math_oct_m
  use messages_oct_m
  use mesh_function_oct_m
  use namespace_oct_m
  use parser_oct_m
  use potential_interpolation_oct_m
  use profiling_oct_m
  use propagator_base_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use states_mxll_oct_m
  use types_oct_m
  use v_ks_oct_m
  use propagation_ops_elec_oct_m
  use xc_oct_m

  implicit none

  private

  public ::                       &
    td_etrs,                      &
    td_etrs_sc,                   &
    td_aetrs,                     &
    td_caetrs

contains

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry
  subroutine td_etrs(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
    integer :: ik, ib
    type(batch_t) :: zpsib_dt
    type(density_calc_t) :: dens_calc

    PUSH_SUB(td_etrs)
    
    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      
      SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

      call propagation_ops_elec_fuse_density_exp_apply(tr%te, st, gr, hm, CNST(0.5)*dt, dt)

      call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)
      call hamiltonian_elec_update(hm, gr%mesh, namespace, time = time - dt)

    else

      call propagation_ops_elec_exp_apply(tr%te, st, gr%mesh, hm, CNST(0.5)*dt)

    end if

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, geo, &
               time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, hm, dt, time)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if
    
    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)
    
    ! propagate dt/2 with H(time - dt)
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, st, gr, hm, CNST(0.5)*dt)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_DEALLOCATE_A(vhxc_t1)
      SAFE_DEALLOCATE_A(vhxc_t2)
    end if

    POP_SUB(td_etrs)
  end subroutine td_etrs

  ! ---------------------------------------------------------
  !> Propagator with enforced time-reversal symmetry and self-consistency
  subroutine td_etrs_sc(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions, sctol, scsteps,&
        hm_mxll, gr_mxll, st_mxll, tr_mxll, mx_sctol)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions
    FLOAT,                            intent(in)    :: sctol
    integer,                optional, intent(out)   :: scsteps
    type(hamiltonian_mxll_t),  optional, intent(inout) :: hm_mxll
    type(grid_t),              optional, intent(inout) :: gr_mxll
    type(states_mxll_t),       optional, intent(inout) :: st_mxll
    type(propagator_mxll_t),   optional, intent(inout) :: tr_mxll
    FLOAT,                     optional, intent(in)    :: mx_sctol

    FLOAT :: diff, maxwell_diff, mx_dt, mx_time
    FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:), rs_state_diff(:), tmp_ma_gr(:,:), tmp_mx_gr(:,:)
    integer :: ik, ib, iter, ip
    type(batch_t), allocatable :: psi2(:, :)
    ! these are hardcoded for the moment
    integer, parameter :: niter = 10
    logical :: prop_maxwell = .false., pml_check = .false.
    CMPLX, allocatable :: rs_current_density_t1(:,:), rs_current_density_t2(:,:)
    CMPLX, allocatable :: rs_current_density_ext(:,:)
    CMPLX, allocatable :: rs_charge_density_t1(:), rs_charge_density_t2(:)
    CMPLX, allocatable :: rs_state_t1(:,:), rs_state_t2(:,:), rs_state_plane_waves(:,:)

    PUSH_SUB(td_etrs_sc)

    ASSERT(hm%theory_level /= INDEPENDENT_PARTICLES)
    
    prop_maxwell = (present(hm_mxll) .and. present(st_mxll) .and. present(gr_mxll) .and. present(tr_mxll))
    
    ! Maxwell field propagation: auxiliary arrays and RS rs_density at time t1 (1 of 3 "if (maxwell)" insertions)
    if (prop_maxwell) then

      do idim = 1, 3
        if (hm_mxll%bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then
          pml_check = .true.
        end if
      end do

      ! auxiliary arrays initialization
      SAFE_ALLOCATE(rs_state_t1(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
      SAFE_ALLOCATE(rs_state_t2(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
      SAFE_ALLOCATE(rs_state_diff(1:gr_mxll%mesh%np))
      if (hm_mxll%ma_mx_coupling_apply .or. hm_mxll%current_density_ext_flag) then
        SAFE_ALLOCATE(rs_current_density_t1(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
        rs_current_density_t1  = M_z0
        SAFE_ALLOCATE(rs_current_density_t2(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
        rs_current_density_t2  = M_z0
        if (hm_mxll%current_density_ext_flag) then
          SAFE_ALLOCATE(rs_current_density_ext(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))
          rs_current_density_ext = M_z0
        end if
        if (hm_mxll%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
          SAFE_ALLOCATE(rs_charge_density_t1(1:gr_mxll%mesh%np_part))
          rs_charge_density_t1 = M_z0
          SAFE_ALLOCATE(rs_charge_density_t2(1:gr_mxll%mesh%np_part))
          rs_charge_density_t2 = M_z0
        end if
      end if
      SAFE_ALLOCATE(rs_state_plane_waves(1:gr_mxll%mesh%np_part,1:st_mxll%d%dim))

      ! storage of the RS state of the old step time t-dt
      rs_state_t1(:,:) = st_mxll%rs_state(:,:)

      ! storage or the RS incident waves state without any coupling
      if (tr_mxll%bc_plane_waves) then
        rs_state_plane_waves(:,:) = st%rs_state_plane_waves(:,:)
      end if

      if (hm_mxll%propagation_apply) then
        if (hm_mxll%ma_mx_coupling_apply) then
          ! calculation of rs_density at time (time-dt)
          call get_rs_density(st_mxll, gr_mxll, hm_mxll, st, gr, hm, geo, hm%mx_ma_coupling_points, &
            hm%mx_ma_coupling_points_number, rs_current_density_t1, &
            rs_charge_density_t1, time-dt, tr%current_prop_test)
        end if
        ! calculation of the external RS density at time (time-dt) 
        if (hm_mxll%current_density_ext_flag) then
          call get_rs_density_ext(st_mxll, gr_mxll%mesh, time-dt, rs_current_density_ext)
          if (.not. hm_mxll%ma_mx_coupling_apply) rs_current_density_t1 = M_z0
          rs_current_density_t1 = rs_current_density_t1 + rs_current_density_ext
        end if
        ! store old convolution function of cpml calculation
        if (pml_check) then
          hm_mxll%bc%pml_conv_plus_old  = hm_mxll%bc%pml_conv_plus
          hm_mxll%bc%pml_conv_minus_old = hm_mxll%bc%pml_conv_minus
        end if
        ! calculation of the Maxwell to matter coupling points at time (time)
        call get_mx_ma_coupling_points(geo, hm, hm%mx_ma_coupling_points(:,:))
        ! get the map of each grid point to the next Maxwell to matter coupling point
        call grid_points_coupling_points_mapping(gr, hm%mx_ma_coupling_points, hm%mx_ma_coupling_points_number, &
          hm%mx_ma_coupling_points_map)
      end if

      call hamiltonian_update(hm, gr%mesh, time = time - dt)

    end if


    SAFE_ALLOCATE(vhxc_t1(1:gr%mesh%np, 1:st%d%nspin))
    SAFE_ALLOCATE(vhxc_t2(1:gr%mesh%np, 1:st%d%nspin))
    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t1)

    call messages_new_line()
    call messages_write('        Self-consistency iteration:')
    call messages_info()

    !Propagate the states to t+dt/2 and compute the density at t+dt
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, st, gr, hm, M_HALF*dt, dt)

    call v_ks_calc(ks, namespace, hm, st, geo, calc_current = .false.)

    call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)
    call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t1, hm%vhxc)

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time - dt)

    ! propagate dt/2 with H(t)

    ! first move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, geo, &
                time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, hm, dt, time)

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      call lalg_copy(gr%mesh%np, st%d%nspin, vhxc_t2, hm%vhxc)
    end if

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    if (prop_maxwell) then
          
      ! reset RS state to the old correct step time t-dt
      st_mxll%rs_state(:,:) = rs_state_t1(:,:)

      ! reset RS incident waves state without any coupling
      if (tr_mxll%bc_plane_waves) then
        st_mxll%rs_state_plane_waves(:,:) = rs_state_plane_waves(:,:)
      end if

      if (hm_mxll%propagation_apply) then
        if (hm_mxll%ma_mx_coupling_apply) then
          ! calculation of rs_density at time (time)
          call get_rs_density(st_mxll, gr_mxll, hm_mxll, st, gr, hm, geo, hm%mx_ma_coupling_points, &
            hm%mx_ma_coupling_points_number, rs_current_density_t2, rs_charge_density_t2, time, &
            tr%current_prop_test)
        end if
        ! calculation of external RS density at time (time)
        if (hm_mxll%current_density_ext_flag) then
          call get_rs_density_ext(st, gr%mesh, time, rs_current_density_ext)
          if (.not. hm%ma_mx_coupling_apply) rs_current_density_t2 = M_z0
          rs_current_density_t2 = rs_current_density_t2 + rs_current_density_ext
        end if
        ! reset old convolution function of cpml calculation
        if (pml_check) then
          hm_mxll%bc%pml_conv_plus  = hm_mxll%bc%pml_conv_plus_old
          hm_mxll%bc%pml_conv_minus = hm_mxll%bc%pml_conv_minus_old
        end if
        ! Setting electron density for Maxwell grid
        if (hm_mxll%diamag_current) then
          tmp_ma_gr = M_ZERO
          tmp_ma_gr(1:gr%mesh%np,1:st%d%nspin) = st%rho(1:gr%mesh%np,1:st%d%nspin)
          tmp_mx_gr = M_ZERO
          call dma_mesh_to_mx_mesh(st_mxll, gr_mxll, st, gr, tmp_ma_gr, tmp_mx_gr, st%d%nspin)
          st_mxll%grid_rho(1:gr%mesh%np,1:st%d%nspin) = tmp_mx_gr(1:gr%mesh%np,1:st%d%nspin)
          SAFE_DEALLOCATE_A(tmp_ma_gr)
          SAFE_DEALLOCATE_A(tmp_mx_gr)
        end if
        ! Maxwell propagation step from time (time-dt) to (time) as a first prediction
        call propagation_mxll_etrs(hm_mxll, gr_mxll, st_mxll, tr_mxll, st_mxll%rs_state, rs_current_density_t1, &
          rs_current_density_t2, rs_charge_density_t1, rs_charge_density_t2, time-dt, dt, rs_state_t1)
      end if

      ! calculation of the Maxwell to matter coupling points at time (time)
      call get_mx_ma_coupling_points(geo, hm, hm%mx_ma_coupling_points(:,:))
      ! get the map of each grid point to the next Maxwell to matter coupling point
      call gr_mxllid_points_coupling_points_mapping(gr, hm%mx_ma_coupling_points, hm%mx_ma_coupling_points_number, &
        hm%mx_ma_coupling_points_map)

      ! coupling of Maxwell fields to the matter
      if (hm%mx_ma_coupling_apply .and. (tr%current_prop_test == 0)) then
        ! transverse field calculation
        call get_vector_pot_and_transverse_field(hm%mx_ma_trans_field_calc_method, gr_mxll, hm_mxll, st_mxll, &
          tr_mxll, gr, hm, st, tr, hm_mxll%poisson_solver, time, st_mxll%rs_state, st_mxll%rs_state_trans, &
          hm_mxll%vector_potential)
        ! electric potential calculation
        call epot_generate_maxwell_coupling(hm%ep, gr, st_mxll, gr_mxll)
      else
        st_mxll%rs_state_trans = st_mxll%rs_state
      end if

    end if

    SAFE_ALLOCATE(psi2(st%group%block_start:st%group%block_end, st%d%kpt%start:st%d%kpt%end))

    ! store the state at half iteration
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call batch_copy(st%group%psib(ib, ik), psi2(ib, ik))
        if(batch_is_packed(st%group%psib(ib, ik))) call batch_pack(psi2(ib, ik), copy = .false.)
        call batch_copy_data(gr%mesh%np, st%group%psib(ib, ik), psi2(ib, ik))
      end do
    end do

    do iter = 1, niter

      call lalg_copy(gr%mesh%np, st%d%nspin, hm%vhxc, vhxc_t2)

      call propagation_ops_elec_fuse_density_exp_apply(tr%te, st, gr, hm, M_HALF*dt)

      call v_ks_calc(ks, namespace, hm, st, geo, time = time, calc_current = .false.)
      call lda_u_update_occ_matrices(hm%lda_u, gr%mesh, st, hm%hm_base, hm%energy )
      
      ! Maxwell field propagation: self cons. step of Maxwell field propagation (3 of 3 "if (maxwell)" insertions)
      if (prop_maxwell) then

        ! reset RS state to the old correct step time t-dt
        st_mxll%rs_state(:,:) = rs_state_t1(:,:)

        ! reset RS incident waves state without any coupling
        if (tr_mxll%bc_plane_waves) then
          st_mxll%rs_state_plane_waves(:,:) = rs_state_plane_waves(:,:)
        end if

        if (hm_mxll%propagation_apply) then
          if (hm_mxll%ma_mx_coupling_apply) then
            ! calculation of rs_density at time (time)
            call get_rs_density(st_mxll, gr_mxll, hm_mxll, st, gr, hm, geo, hm%mx_ma_coupling_points, &
              hm%mx_ma_coupling_points_number, rs_current_density_t2, rs_charge_density_t2, time, &
              tr%current_prop_test)
          end if
          ! calculation of external RS density at time (time)
          if (hm_mxll%current_density_ext_flag) then
            call get_rs_density_ext(st_mxll, gr_mxll%mesh, time, rs_current_density_ext)
            if (.not. hm_mxll%ma_mx_coupling_apply) rs_current_density_t2 = M_z0
            rs_current_density_t2 = rs_current_density_t2 + rs_current_density_ext
          end if
          ! reset old convolution function of cpml calculation
          if (pml_check) then
            hm_mxll%bc%pml_conv_plus  = hm_mxll%bc%pml_conv_plus_old
            hm_mxll%bc%pml_conv_minus = hm_mxll%bc%pml_conv_minus_old
          end if
          ! Setting electron density for Maxwell grid
          if (hm_mxll%diamag_current) then
            tmp_ma_gr = M_ZERO
            tmp_ma_gr(1:gr%mesh%np,1:st%d%nspin) = st%rho(1:gr%mesh%np,1:st%d%nspin)
            tmp_mx_gr = M_ZERO
            call dma_mesh_to_mx_mesh(st_mxll, gr_mxll, st, gr, tmp_ma_gr, tmp_mx_gr, st%d%nspin)
            st_mxll%grid_rho(1:gr%mesh%np,1:st%d%nspin) = tmp_mx_gr(1:gr%mesh%np,1:st%d%nspin)
            SAFE_DEALLOCATE_A(tmp_ma_gr)
            SAFE_DEALLOCATE_A(tmp_mx_gr)
          end if
          ! Maxwell propagation step from time (time-dt) to (time) as a first prediction
          call propagation_mxll_etrs(hm_mxll, gr_mxll, st_mxll, tr_mxll, st_mxll%rs_state, & 
            rs_current_density_t1, rs_current_density_t2, rs_charge_density_t1, rs_charge_density_t2, time-dt, dt,&
            rs_state_t1)
        end if

        ! calculation of the Maxwell to matter coupling points at time (time)
        call get_mx_ma_coupling_points(geo, hm, hm%mx_ma_coupling_points(:,:))
        ! get the map of each grid point to the next Maxwell to matter coupling point
        call maxwell_grid_points_coupling_points_mapping(gr, hm%mx_ma_coupling_points, &
          hm%mx_ma_coupling_points_number, hm%mx_ma_coupling_points_map)

        ! coupling of Maxwell fields to the matter
        if (hm%mx_ma_coupling_apply .and. (tr%current_prop_test == 0)) then
          ! transverse field calculation
          call get_vector_pot_and_transverse_field(hm%mx_ma_trans_field_calc_method, gr_mxll, hm_mxll, st_mxll, &
            tr_mxll, gr, hm, st, tr, hm_mxll%poisson_solver, time,st_mxll%rs_state, st_mxll%rs_state_trans, &
            hm_mxll%vector_potential)
          ! electric potential calculation
          call epot_generate_maxwell_coupling(hm%ep, gr, st_mxll, gr_mxll)
        else
          st_mxll%rs_state_trans = st_mxll%rs_state
        end if

      end if

      ! now check how much the potential changed
      do ip = 1, gr%mesh%np
        vhxc_t2(ip, 1) = sum(abs(vhxc_t2(ip, 1:st%d%nspin) - hm%vhxc(ip, 1:st%d%nspin)))
      end do
      diff = dmf_integrate(gr%mesh, vhxc_t2(:, 1))

      ! now check how much the electromagnetic field changed
      if (prop_maxwell) then
        do ip = 1, gr_mxll%mesh%np
          rs_state_diff(ip) = sum(abs(rs_state_t2(ip,1:st_mxll%d%dim) & 
            -st_mxll%rs_state(ip,1:st_mxll%d%dim)))
        end do
        maxwell_diff = dmf_integrate(gr_mxll%mesh, rs_state_diff(:))
      end if

      if (prop_maxwell) then
        call messages_write('          step ')
        call messages_write(iter)
        call messages_write(', KS-residue = ')
        call messages_write(abs(diff), fmt = '(1x,es9.2)')
        call messages_write(', Maxwell-residue = ')
        call messages_write(abs(maxwell_diff), fmt = '(1x,es9.2)')
        call messages_info()
      else
        call messages_write('          step ')
        call messages_write(iter)
        call messages_write(', residue = ')
        call messages_write(abs(diff), fmt = '(1x,es9.2)')
        call messages_info()
      end if

      if(maxwell) then
        if((diff <= sctol) .and. (maxwell_diff <= mx_sctol)) exit
      else
        if(diff <= sctol) exit
      end if

      if(iter /= niter) then
        ! we are not converged, restore the states
        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call batch_copy_data(gr%mesh%np, psi2(ib, ik), st%group%psib(ib, ik))
          end do
        end do
        ! copy current RS state for next sc step correction
        if (prop_maxwell) rs_state_t2(:,:) = st_mxll%rs_state(:,:)
      end if

    end do

    if(hm%lda_u_level /= DFT_U_NONE) then 
      call lda_u_write_U(hm%lda_u, stdout) 
      call lda_u_write_V(hm%lda_u, stdout)
    end if

    ! print an empty line
    call messages_info()

    if(present(scsteps)) scsteps = iter

    SAFE_DEALLOCATE_A(vhxc_t1)
    SAFE_DEALLOCATE_A(vhxc_t2)

    if (prop_maxwell) then
      ! auxiliary arrays deallocation
      SAFE_DEALLOCATE_A(rs_state_t1)
      SAFE_DEALLOCATE_A(rs_state_t2)
      SAFE_DEALLOCATE_A(rs_state_diff)
      if (hm_mxll%ma_mx_coupling_apply .or. hm_mxll%current_density_ext_flag) then
        SAFE_DEALLOCATE_A(rs_current_density_t1)
        SAFE_DEALLOCATE_A(rs_current_density_t2)
        if (hm_mxll%current_density_ext_flag) then
          SAFE_DEALLOCATE_A(rs_current_density_ext)
        end if
        if (hm_mxll%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
          SAFE_DEALLOCATE_A(rs_charge_density_t1)
          SAFE_DEALLOCATE_A(rs_charge_density_t2)
        end if
      end if
      SAFE_DEALLOCATE_A(rs_state_plane_waves)
    end if


    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call batch_end(psi2(ib, ik))
      end do
    end do

    SAFE_DEALLOCATE_A(psi2)

    POP_SUB(td_etrs_sc)
  end subroutine td_etrs_sc

  ! ---------------------------------------------------------
  !> Propagator with approximate enforced time-reversal symmetry
  subroutine td_aetrs(namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    PUSH_SUB(td_aetrs)

    ! propagate half of the time step with H(time - dt)
    call propagation_ops_elec_exp_apply(tr%te, st, gr%mesh, hm, M_HALF*dt)

    !Get the potentials from the interpolator
    call propagation_ops_elec_interpolate_get(gr%mesh, hm, tr%vksold)

    ! move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, &
              geo, time, ionic_scale*dt, move_ions = move_ions)

    !Propagate gauge field
    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, hm, dt, time)

    !Update Hamiltonian
    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    !Do the time propagation for the second half of the time step
    call propagation_ops_elec_fuse_density_exp_apply(tr%te, st, gr, hm, M_HALF*dt)

    POP_SUB(td_aetrs)
  end subroutine td_aetrs

  ! ---------------------------------------------------------
  !> Propagator with approximate enforced time-reversal symmetry
  subroutine td_caetrs(ks, namespace, hm, gr, st, tr, time, dt, ionic_scale, ions, geo, move_ions)
    type(v_ks_t),             target, intent(inout) :: ks
    type(namespace_t),                intent(in)    :: namespace
    type(hamiltonian_elec_t), target, intent(inout) :: hm
    type(grid_t),             target, intent(inout) :: gr
    type(states_elec_t),      target, intent(inout) :: st
    type(propagator_t),       target, intent(inout) :: tr
    FLOAT,                            intent(in)    :: time
    FLOAT,                            intent(in)    :: dt
    FLOAT,                            intent(in)    :: ionic_scale
    type(ion_dynamics_t),             intent(inout) :: ions
    type(geometry_t),                 intent(inout) :: geo
    logical,                          intent(in)    :: move_ions

    integer :: ik, ispin, ip, ist, ib
    FLOAT :: vv
    CMPLX :: phase
    type(density_calc_t)  :: dens_calc
    type(profile_t), save :: phase_prof
    integer               :: pnp, iprange
    FLOAT, allocatable    :: vold(:, :), vtauold(:, :)
    type(accel_mem_t)    :: phase_buff

    PUSH_SUB(td_caetrs)

    SAFE_ALLOCATE(vold(1:gr%mesh%np, 1:st%d%nspin))
    if (family_is_mgga_with_exc(hm%xc)) then 
      SAFE_ALLOCATE(vtauold(1:gr%mesh%np, 1:st%d%nspin))
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold, vtau = vtauold)
      call hamiltonian_elec_set_vhxc(hm, gr%mesh, vold, vtauold)
    else
      call potential_interpolation_get(tr%vksold, gr%mesh%np, st%d%nspin, 2, vold)
      call hamiltonian_elec_set_vhxc(hm, gr%mesh, vold)
    endif

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time - dt) 

    call v_ks_calc_start(ks, namespace, hm, st, geo, time = time - dt, calc_energy = .false., &
           calc_current = .false.)

    ! propagate half of the time step with H(time - dt)
    call propagation_ops_elec_exp_apply(tr%te, st, gr%mesh, hm, M_HALF*dt)

    call v_ks_calc_finish(ks, hm, namespace)

    if (family_is_mgga_with_exc(hm%xc)) then 
      call potential_interpolation_set(tr%vksold, gr%mesh%np, st%d%nspin, 1, hm%vhxc, vtau = hm%vtau)
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%v_old(:, :, 1:3), time, tr%vksold%v_old(:, :, 0))
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%vtau_old(:, :, 1:3), time, tr%vksold%vtau_old(:, :, 0))
      forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np)
        vold(ip, ispin) =  CNST(0.5)*dt*(hm%vhxc(ip, ispin) - vold(ip, ispin))
        vtauold(ip, ispin) =  CNST(0.5)*dt*(hm%vtau(ip, ispin) - vtauold(ip, ispin))
      end forall      
    else
      call potential_interpolation_set(tr%vksold, gr%mesh%np, st%d%nspin, 1, hm%vhxc)
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
         tr%vksold%v_old(:, :, 1:3), time, tr%vksold%v_old(:, :, 0))

      forall(ispin = 1:st%d%nspin, ip = 1:gr%mesh%np)
        vold(ip, ispin) =  CNST(0.5)*dt*(hm%vhxc(ip, ispin) - vold(ip, ispin))
      end forall
    end if

    ! copy vold to a cl buffer
    if(accel_is_enabled() .and. hamiltonian_elec_apply_packed(hm, gr%mesh)) then
      if(family_is_mgga_with_exc(hm%xc)) then
        call messages_not_implemented('CAETRS propagator with accel and MGGA with energy functionals')
      end if
      pnp = accel_padded_size(gr%mesh%np)
      call accel_create_buffer(phase_buff, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, pnp*st%d%nspin)
      ASSERT(ubound(vold, dim = 1) == gr%mesh%np)
      do ispin = 1, st%d%nspin
        call accel_write_buffer(phase_buff, gr%mesh%np, vold(:, ispin), offset = (ispin - 1)*pnp)
      end do
    end if

    !Get the potentials from the interpolator
    call propagation_ops_elec_interpolate_get(gr%mesh, hm, tr%vksold)

    ! move the ions to time t
    call propagation_ops_elec_move_ions(tr%propagation_ops_elec, gr, hm, st, namespace, ions, &
              geo, time, ionic_scale*dt, move_ions = move_ions)

    call propagation_ops_elec_propagate_gauge_field(tr%propagation_ops_elec, hm, dt, time)

    call propagation_ops_elec_update_hamiltonian(namespace, st, gr%mesh, hm, time)

    call density_calc_init(dens_calc, st, gr, st%rho)

    ! propagate the other half with H(t)
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_elec_dim_get_spin_index(st%d, ik)

      do ib = st%group%block_start, st%group%block_end
        if (hamiltonian_elec_apply_packed(hm, gr%mesh)) then
          call batch_pack(st%group%psib(ib, ik))
          if (hamiltonian_elec_inh_term(hm)) call batch_pack(hm%inh_st%group%psib(ib, ik))
        end if

        call profiling_in(phase_prof, "CAETRS_PHASE")
        select case(batch_status(st%group%psib(ib, ik)))
        case(BATCH_NOT_PACKED)
          do ip = 1, gr%mesh%np
            vv = vold(ip, ispin)
            phase = TOCMPLX(cos(vv), -sin(vv))
            forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
              st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip) = st%group%psib(ib, ik)%states_linear(ist)%zpsi(ip)*phase
            end forall
          end do
        case(BATCH_PACKED)
          do ip = 1, gr%mesh%np
            vv = vold(ip, ispin)
            phase = TOCMPLX(cos(vv), -sin(vv))
            forall(ist = 1:st%group%psib(ib, ik)%nst_linear)
              st%group%psib(ib, ik)%pack%zpsi(ist, ip) = st%group%psib(ib, ik)%pack%zpsi(ist, ip)*phase
            end forall
          end do
        case(BATCH_DEVICE_PACKED)
          call accel_set_kernel_arg(kernel_phase, 0, pnp*(ispin - 1))
          call accel_set_kernel_arg(kernel_phase, 1, phase_buff)
          call accel_set_kernel_arg(kernel_phase, 2, st%group%psib(ib, ik)%pack%buffer)
          call accel_set_kernel_arg(kernel_phase, 3, log2(st%group%psib(ib, ik)%pack%size(1)))

          iprange = accel_max_workgroup_size()/st%group%psib(ib, ik)%pack%size(1)

          call accel_kernel_run(kernel_phase, (/st%group%psib(ib, ik)%pack%size(1), pnp/), &
            (/st%group%psib(ib, ik)%pack%size(1), iprange/))
        end select
        call profiling_out(phase_prof)

        if (hamiltonian_elec_inh_term(hm)) then
          call exponential_apply_batch(tr%te, gr%mesh, hm, st%group%psib(ib, ik), ik, CNST(0.5)*dt, &
            inh_psib = hm%inh_st%group%psib(ib, ik))
        else
          call exponential_apply_batch(tr%te, gr%mesh, hm, st%group%psib(ib, ik), ik, CNST(0.5)*dt)
        end if

        call density_calc_accumulate(dens_calc, ik, st%group%psib(ib, ik))

        if (hamiltonian_elec_apply_packed(hm, gr%mesh)) then
          call batch_unpack(st%group%psib(ib, ik))
          if (hamiltonian_elec_inh_term(hm)) call batch_unpack(hm%inh_st%group%psib(ib, ik))
        end if
      end do
    end do

    call density_calc_end(dens_calc)

    if(accel_is_enabled() .and. hamiltonian_elec_apply_packed(hm, gr%mesh)) then
      call accel_release_buffer(phase_buff)
    end if

    SAFE_DEALLOCATE_A(vold)
    SAFE_DEALLOCATE_A(vtauold)

    POP_SUB(td_caetrs)
  end subroutine td_caetrs

end module propagator_etrs_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
