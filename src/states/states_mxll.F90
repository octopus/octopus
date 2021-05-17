!! Copyright (C) 2019 R. Jestaedt, F. Bonafe, H. Appel, A. Rubio
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

module states_mxll_oct_m
  use accel_oct_m
  use blacs_proc_grid_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use global_oct_m
  use grid_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use states_elec_dim_oct_m
  use states_elec_group_oct_m
  use states_elec_oct_m
  use tdfunction_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    states_mxll_t,                    &
    states_mxll_init,                 &
    states_mxll_allocate,             &
    states_mxll_end,                  &
    build_rs_element,                 &
    build_rs_vector,                  &
    build_rs_state,                   &
    build_rs_current_element,         &
    build_rs_current_vector,          &
    build_rs_current_state,           &
    get_electric_field_vector,        &
    get_magnetic_field_vector,        &
    get_electric_field_state,         &
    get_magnetic_field_state,         &
    get_current_element,              &
    get_current_vector,               &
    get_current_state,                &
    get_rs_state_at_point,            &
    get_divergence_field,             &
    get_poynting_vector,              &
    get_poynting_vector_plane_waves,  &
    state_diff

  type :: states_mxll_t
    ! Components are public by default
    integer                      :: dim         !< Space dimension
    integer                      :: rs_sign
    logical                      :: pack_states
    logical                      :: parallel_in_states = .false. !< Am I parallel in states?
    type(type_t), public         :: wfs_type     !< should always be complex (TYPE_CMPLX)
    integer, public              :: nst          !< Number of RS states, currently set to 1, we keep it for future uses
    logical, public              :: packed

    CMPLX, allocatable           :: rs_state_plane_waves(:,:)
    CMPLX, allocatable           :: rs_state(:,:)
    CMPLX, allocatable           :: rs_state_trans(:,:)
    CMPLX, allocatable           :: rs_state_long(:,:)

    logical                      :: rs_current_density_restart = .false.
    CMPLX, allocatable           :: rs_current_density_restart_t1(:,:)
    CMPLX, allocatable           :: rs_current_density_restart_t2(:,:)

    FLOAT, allocatable           :: ep(:)
    FLOAT, allocatable           :: mu(:)

    integer, allocatable         :: rs_state_fft_map(:,:,:)
    integer, allocatable         :: rs_state_fft_map_inv(:,:)

    FLOAT                        :: energy_rate
    FLOAT                        :: delta_energy
    FLOAT                        :: energy_via_flux_calc

    FLOAT                        :: trans_energy_rate
    FLOAT                        :: trans_delta_energy
    FLOAT                        :: trans_energy_via_flux_calc

    FLOAT                        :: plane_waves_energy_rate
    FLOAT                        :: plane_waves_delta_energy
    FLOAT                        :: plane_waves_energy_via_flux_calc

    FLOAT                        :: poynting_vector_box_surface(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO
    FLOAT                        :: poynting_vector_box_surface_plane_waves(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO
    FLOAT                        :: electric_field_box_surface(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO
    FLOAT                        :: electric_field_box_surface_plane_waves(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO
    FLOAT                        :: magnetic_field_box_surface(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO
    FLOAT                        :: magnetic_field_box_surface_plane_waves(1:2,1:MAX_DIM,1:MAX_DIM) = M_ZERO

    logical                      :: rs_state_const_external = .false.
    CMPLX, allocatable           :: rs_state_const(:)
    CMPLX, allocatable           :: rs_state_const_amp(:,:)
    type(tdf_t), allocatable     :: rs_state_const_td_function(:)

    FLOAT                        :: poynting_mean(MAX_DIM)
    FLOAT                        :: poynting_mean_plane_waves(MAX_DIM)

    integer                      :: inner_points_number
    integer, allocatable         :: inner_points_map(:)
    logical, allocatable         :: inner_points_mask(:)
    integer                      :: boundary_points_number
    integer, allocatable         :: boundary_points_map(:)
    logical, allocatable         :: boundary_points_mask(:)

    integer                      :: surface_points_number(MAX_DIM)
    integer, allocatable         :: surface_points_map(:,:,:)
    FLOAT                        :: surface_element(MAX_DIM)

    integer                      :: surface_grid_rows_number(MAX_DIM)
    integer, allocatable         :: surface_grid_points_number(:,:,:)
    integer, allocatable         :: surface_grid_points_map(:,:,:,:,:)
    integer, allocatable         :: surface_grid_center(:,:,:,:)
    FLOAT                        :: surface_grid_element(MAX_DIM)

    type(mesh_plane_t)           :: surface(2,MAX_DIM)

    integer                      :: selected_points_number
    FLOAT, allocatable           :: selected_points_coordinate(:,:)
    CMPLX, allocatable           :: selected_points_rs_state(:,:)
    CMPLX, allocatable           :: selected_points_rs_state_trans(:,:)
    FLOAT                        :: rs_state_trans_var

    FLOAT, allocatable           :: grid_rho(:,:)
    CMPLX, allocatable           :: kappa_psi(:,:)

    character(len=1024), allocatable :: user_def_e_field(:)
    character(len=1024), allocatable :: user_def_b_field(:)

    integer                      :: energy_incident_waves_calc_iter
    logical                      :: energy_incident_waves_calc

    ! external current variables
    integer                      :: external_current_number
    integer,             allocatable :: external_current_modus(:)
    character(len=1024), allocatable :: external_current_string(:,:)
    FLOAT,               allocatable :: external_current_amplitude(:,:,:)
    type(tdf_t),         allocatable :: external_current_td_function(:)
    type(tdf_t),         allocatable :: external_current_td_phase(:)
    FLOAT,               allocatable :: external_current_omega(:)
    FLOAT,               allocatable :: external_current_phase(:)

    !> used for the user-defined wavefunctions (they are stored as formula strings)
    character(len=1024), allocatable :: user_def_states(:,:,:)
    logical                     :: fromScratch
    type(mpi_grp_t)             :: mpi_grp
    type(mpi_grp_t)             :: dom_st_mpi_grp

#ifdef HAVE_SCALAPACK
    type(blacs_proc_grid_t)     :: dom_st_proc_grid
#endif
    type(distributed_t)         :: dist
    logical                     :: scalapack_compatible
    integer                     :: lnst
    integer                     :: st_start, st_end
    integer, allocatable        :: node(:)

  end type states_mxll_t

contains

  ! ---------------------------------------------------------
  subroutine states_mxll_init(st, namespace, gr)
    type(states_mxll_t), target, intent(inout) :: st
    type(namespace_t),           intent(in)    :: namespace
    type(grid_t),                intent(in)    :: gr
    type(block_t)        :: blk

    integer :: idim, nlines, ncols, il
    logical :: defaultl
    FLOAT, allocatable   :: pos(:)
    integer :: ix_max, iy_max, iz_max
    type(profile_t), save :: prof

    PUSH_SUB(states_mxll_init)

    call profiling_in(prof, 'STATES_MXLL_INIT')

    st%wfs_type = TYPE_CMPLX
    st%fromScratch = .true. ! this will be reset if restart_read is called

    ASSERT(MAX_DIM >= gr%sb%dim)
    ASSERT(gr%sb%dim == 3)
    st%dim = gr%sb%dim
    st%nst = 1

    SAFE_ALLOCATE(st%user_def_e_field(1:st%dim))
    SAFE_ALLOCATE(st%user_def_b_field(1:st%dim))

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst

    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.
    st%packed = .false.

    defaultl = .true.
    if(accel_is_enabled()) then
      defaultl = .false.
    end if
    call parse_variable(namespace, 'StatesPack', defaultl, st%pack_states)

    call messages_print_var_value(stdout, 'StatesPack', st%pack_states)

    !%Variable RiemannSilbersteinSign
    !%Type integer
    !%Default plus
    !%Section Hamiltonian
    !%Description
    !% Sign for the imaginary part of the Riemann Silberstein vector which represents the magnetic field
    !%Option plus 1
    !% Riemann Silberstein sign is plus
    !%Option minus -1
    !% Riemann Silberstein sign is minus
    !%End
    call parse_variable(namespace, 'RiemannSilbersteinSign', OPTION__RIEMANNSILBERSTEINSIGN__PLUS, st%rs_sign)

    !%Variable MaxwellFieldsCoordinate
    !%Type block
    !%Section System::Coordinates
    !%Description
    !%
    !% <tt>%Coordinates
    !% <br>&nbsp;&nbsp;    -1.0 | 2.0 |  4.0
    !% <br>&npsp;&nbsp;     0.0 | 1.0 | -2.0
    !% <br>%</tt>
    !%
    !%End

    SAFE_ALLOCATE(pos(1:st%dim))
    st%selected_points_number = 1
    if(parse_block(namespace, 'MaxwellFieldsCoordinate', blk) == 0) then
      nlines = parse_block_n(blk)
      st%selected_points_number = nlines
      SAFE_ALLOCATE(st%selected_points_coordinate(1:st%dim,1:nlines))
      SAFE_ALLOCATE(st%selected_points_rs_state(1:st%dim,1:nlines))
      SAFE_ALLOCATE(st%selected_points_rs_state_trans(1:st%dim,1:nlines))
      do il = 1, nlines
        ncols = parse_block_cols(blk,0)
        if (ncols < 3 .or. ncols > 3) then
            message(1) = 'MaxwellFieldCoordinate must have 3 columns.'
            call messages_fatal(1, namespace=namespace)
        end if
        do idim = 1, st%dim
          call parse_block_float(blk, il-1, idim-1, pos(idim), units_inp%length)
        end do
        st%selected_points_coordinate(:,il) = pos
        st%selected_points_rs_state(:,il)  = M_z0
        st%selected_points_rs_state_trans(:,il) = M_z0
      end do
    call parse_block_end(blk)
    else
      SAFE_ALLOCATE(st%selected_points_coordinate(1:st%dim, 1))
      SAFE_ALLOCATE(st%selected_points_rs_state(1:st%dim, 1))
      SAFE_ALLOCATE(st%selected_points_rs_state_trans(1:st%dim, 1))
      st%selected_points_coordinate(:,:) = M_ZERO
      st%selected_points_rs_state(:,:) = M_z0
      st%selected_points_rs_state_trans(:,:) = M_z0
    end if

    SAFE_DEALLOCATE_A(pos)

    st%surface_grid_rows_number(1) = 3
    ix_max  = st%surface_grid_rows_number(1)
    st%surface_grid_rows_number(2) = 3
    iy_max  = st%surface_grid_rows_number(2)
    st%surface_grid_rows_number(3) = 3
    iz_max  = st%surface_grid_rows_number(3)

    SAFE_ALLOCATE(st%surface_grid_center(1:2, 1:st%dim, 1:ix_max, 1:iy_max))
    SAFE_ALLOCATE(st%surface_grid_points_number(1:st%dim, 1:ix_max, 1:iy_max))

    call profiling_out(prof)

    POP_SUB(states_mxll_init)

  end subroutine states_mxll_init

  ! ---------------------------------------------------------
  !> Allocates the Maxwell states defined within a states_mxll_t structure.
  subroutine states_mxll_allocate(st, mesh)
    type(states_mxll_t),    intent(inout)   :: st
    type(mesh_t),           intent(in)      :: mesh

    type(profile_t), save :: prof

    PUSH_SUB(states_mxll_allocate)

    call profiling_in(prof, 'STATES_MXLL_ALLOCATE')

    SAFE_ALLOCATE(st%rs_state(1:mesh%np_part, 1:st%dim))
    st%rs_state(:,:) = M_z0

    SAFE_ALLOCATE(st%rs_state_trans(1:mesh%np_part, 1:st%dim))
    st%rs_state_trans(:,:) = M_z0

    SAFE_ALLOCATE(st%rs_state_long(1:mesh%np_part, 1:st%dim))
    st%rs_state_long(:,:) = M_z0

    SAFE_ALLOCATE(st%rs_state_plane_waves(1:mesh%np_part, 1:st%dim))
    st%rs_state_plane_waves(:,:) = M_z0

    SAFE_ALLOCATE(st%rs_current_density_restart_t1(1:mesh%np_part, 1:st%dim))
    st%rs_current_density_restart_t1 = M_z0

    SAFE_ALLOCATE(st%rs_current_density_restart_t2(1:mesh%np_part, 1:st%dim))
    st%rs_current_density_restart_t2 = M_z0

    SAFE_ALLOCATE(st%ep(1:mesh%np_part))
    SAFE_ALLOCATE(st%mu(1:mesh%np_part))
    st%ep = P_ep
    st%mu = P_mu

    call profiling_out(prof)

    POP_SUB(states_mxll_allocate)
  end subroutine states_mxll_allocate

  ! ---------------------------------------------------------
  subroutine states_mxll_end(st)
    type(states_mxll_t), intent(inout) :: st

    type(profile_t), save :: prof

    PUSH_SUB(states_mxll_end)

    call profiling_in(prof, 'STATES_MXLL_END')

    SAFE_DEALLOCATE_A(st%rs_state)
    SAFE_DEALLOCATE_A(st%rs_state_trans)
    SAFE_DEALLOCATE_A(st%selected_points_coordinate)
    SAFE_DEALLOCATE_A(st%selected_points_rs_state)
    SAFE_DEALLOCATE_A(st%selected_points_rs_state_trans)
    SAFE_DEALLOCATE_A(st%rs_state_long)
    SAFE_DEALLOCATE_A(st%rs_current_density_restart_t1)
    SAFE_DEALLOCATE_A(st%rs_current_density_restart_t2)
    SAFE_DEALLOCATE_A(st%user_def_e_field)
    SAFE_DEALLOCATE_A(st%user_def_b_field)

    SAFE_DEALLOCATE_A(st%rs_state_const)
    SAFE_DEALLOCATE_A(st%rs_state_const_td_function)
    SAFE_DEALLOCATE_A(st%rs_state_const_amp)
    SAFE_DEALLOCATE_A(st%rs_state_plane_waves)

    SAFE_DEALLOCATE_A(st%surface_grid_center)
    SAFE_DEALLOCATE_A(st%surface_grid_points_number)
    SAFE_DEALLOCATE_A(st%surface_grid_points_map)
    SAFE_DEALLOCATE_A(st%inner_points_map)
    SAFE_DEALLOCATE_A(st%boundary_points_map)
    SAFE_DEALLOCATE_A(st%inner_points_mask)
    SAFE_DEALLOCATE_A(st%boundary_points_mask)
    SAFE_DEALLOCATE_A(st%ep)
    SAFE_DEALLOCATE_A(st%mu)
#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_end(st%dom_st_proc_grid)
#endif
    SAFE_DEALLOCATE_A(st%external_current_modus)
    SAFE_DEALLOCATE_A(st%external_current_string)
    SAFE_DEALLOCATE_A(st%external_current_amplitude)
    SAFE_DEALLOCATE_A(st%external_current_td_function)
    SAFE_DEALLOCATE_A(st%external_current_omega)
    SAFE_DEALLOCATE_A(st%external_current_td_phase)

    call distributed_end(st%dist)
    SAFE_DEALLOCATE_A(st%node)

    call profiling_out(prof)

    POP_SUB(states_mxll_end)
  end subroutine states_mxll_end


  !----------------------------------------------------------
  subroutine build_rs_element(e_element, b_element, rs_sign, rs_element, ep_element, mu_element)
    FLOAT,             intent(in)    :: e_element, b_element
    CMPLX,             intent(inout) :: rs_element
    integer,           intent(in)    :: rs_sign
    FLOAT,   optional, intent(in)    :: ep_element
    FLOAT,   optional, intent(in)    :: mu_element

    ! no PUSH_SUB, called too often

    if (present(ep_element) .and. present(mu_element)) then
      rs_element = sqrt(ep_element/M_TWO) * e_element + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*mu_element)) * b_element
    else
      rs_element = sqrt(P_ep/M_TWO) * e_element + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*P_mu)) * b_element
    end if

  end subroutine build_rs_element


  !----------------------------------------------------------
  subroutine build_rs_vector(e_vector, b_vector, rs_sign, rs_vector, ep_element, mu_element)
    FLOAT,             intent(in)    :: e_vector(:), b_vector(:)
    CMPLX,             intent(inout) :: rs_vector(:)
    integer,           intent(in)    :: rs_sign
    FLOAT,   optional, intent(in)    :: ep_element
    FLOAT,   optional, intent(in)    :: mu_element

    ! no PUSH_SUB, called too often

    if (present(ep_element) .and. present(mu_element)) then
      rs_vector = sqrt(ep_element/M_TWO) * e_vector + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*mu_element)) * b_vector
    else
      rs_vector = sqrt(P_ep/M_TWO) * e_vector + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*P_mu)) * b_vector
    end if

  end subroutine build_rs_vector


  !----------------------------------------------------------
  subroutine build_rs_state(e_field, b_field, rs_sign, rs_state, mesh, ep_field, mu_field, np)
    FLOAT,             intent(in)    :: e_field(:,:), b_field(:,:)
    CMPLX,             intent(inout) :: rs_state(:,:)
    integer,           intent(in)    :: rs_sign
    type(mesh_t),      intent(in)    :: mesh
    FLOAT,   optional, intent(in)    :: ep_field(:)
    FLOAT,   optional, intent(in)    :: mu_field(:)
    integer, optional, intent(in)    :: np

    integer :: ip, np_
    type(profile_t), save :: prof

    PUSH_SUB(build_rs_state)

    call profiling_in(prof, 'BUILD_RS_STATE')

    np_ = optional_default(np, mesh%np)

    do ip = 1, np_
      if (present(ep_field) .and. present(mu_field)) then
        rs_state(ip, :) = sqrt(ep_field(ip)/M_TWO) * e_field(ip, :) &
                       + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*mu_field(ip))) * b_field(ip, :)
      else
        rs_state(ip, :) = sqrt(P_ep/M_TWO) * e_field(ip, :) &
                       + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*P_mu)) * b_field(ip, :)
      end if
    end do

   call profiling_out(prof)

    POP_SUB(build_rs_state)

  end subroutine build_rs_state


  !----------------------------------------------------------
  subroutine build_rs_current_element(current_element, rs_current_element, ep_element)
    FLOAT,           intent(in)    :: current_element
    CMPLX,           intent(inout) :: rs_current_element
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      rs_current_element = M_ONE/sqrt(M_TWO*ep_element) * current_element
    else
      rs_current_element = M_ONE/sqrt(M_TWO*P_ep) * current_element
    end if

  end subroutine build_rs_current_element


  !----------------------------------------------------------
  subroutine build_rs_current_vector(current_vector, rs_current_vector, ep_element)
    FLOAT,           intent(in)    :: current_vector(:)
    CMPLX,           intent(inout) :: rs_current_vector(:)
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often
    if (present(ep_element)) then
      rs_current_vector = M_ONE/sqrt(M_TWO*ep_element) * current_vector
    else
      rs_current_vector = M_ONE/sqrt(M_TWO*P_ep) * current_vector
    end if

  end subroutine build_rs_current_vector


  !----------------------------------------------------------
  subroutine build_rs_current_state(current_state, mesh, rs_current_state, ep_field, np)
    FLOAT,             intent(in)    :: current_state(:,:)
    type(mesh_t),      intent(in)    :: mesh
    CMPLX,             intent(inout) :: rs_current_state(:,:)
    FLOAT,   optional, intent(in)    :: ep_field(:)
    integer, optional, intent(in)    :: np

    integer :: ip, idim, np_, ff_dim
    type(profile_t), save :: prof

    ! no PUSH_SUB, called too often

    call profiling_in(prof, "BUILD_RS_CURRENT_STATE")

    np_ = optional_default(np, mesh%np)
    ff_dim = size(current_state, dim=2)

    if(present(ep_field)) then
      do idim = 1, ff_dim
        do ip = 1, np_
          rs_current_state(ip, idim) = M_ONE/sqrt(M_TWO*ep_field(ip)) * current_state(ip, idim)
        end do
      end do
    else
      do idim = 1, ff_dim
        do ip = 1, np_
          rs_current_state(ip, idim) = M_ONE/sqrt(M_TWO*P_ep) * current_state(ip, idim)
        end do
      end do
    end if

    call profiling_out(prof)

  end subroutine build_rs_current_state


  !----------------------------------------------------------
  subroutine get_electric_field_vector(rs_state_vector, electric_field_vector, ep_element)
    CMPLX,             intent(in)    :: rs_state_vector(:)
    FLOAT,             intent(inout) :: electric_field_vector(:)
    FLOAT,   optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      electric_field_vector(:) = sqrt(M_TWO/ep_element) * TOFLOAT(rs_state_vector(:))
    else
      electric_field_vector(:) = sqrt(M_TWO/P_ep) * TOFLOAT(rs_state_vector(:))
    end if

  end subroutine get_electric_field_vector


  !----------------------------------------------------------
  subroutine get_magnetic_field_vector(rs_state_vector, rs_sign, magnetic_field_vector, mu_element)
    CMPLX,             intent(in)    :: rs_state_vector(:)
    integer,           intent(in)    :: rs_sign
    FLOAT,             intent(inout) :: magnetic_field_vector(:)
    FLOAT,   optional, intent(in)    :: mu_element

    ! no PUSH_SUB, called too often

    if (present(mu_element)) then
      magnetic_field_vector(:) = sqrt(M_TWO*mu_element) * rs_sign * aimag(rs_state_vector(:))
    else
      magnetic_field_vector(:) = sqrt(M_TWO*P_mu) * rs_sign * aimag(rs_state_vector(:))
    end if

  end subroutine get_magnetic_field_vector


  !----------------------------------------------------------
  subroutine get_electric_field_state(rs_state, mesh, electric_field, ep_field, np)
    CMPLX,             intent(in)    :: rs_state(:,:)
    type(mesh_t),      intent(in)    :: mesh
    FLOAT,             intent(inout) :: electric_field(:,:)
    FLOAT,   optional, intent(in)    :: ep_field(:)
    integer, optional, intent(in)    :: np

    integer :: ip, np_
    type(profile_t), save :: prof

    PUSH_SUB(get_electric_field_state)

    call profiling_in(prof, 'GET_ELECTRIC_FIELD_STATE')

    np_ = optional_default(np, mesh%np)

    do ip = 1, np_
      if (present(ep_field)) then
        electric_field(ip, :) = sqrt(M_TWO/ep_field(ip)) * TOFLOAT(rs_state(ip, :))
      else
        electric_field(ip,:) = sqrt(M_TWO/P_ep) * TOFLOAT(rs_state(ip, :))
      end if
    end do

    call profiling_out(prof)

    POP_SUB(get_electric_field_state)

  end subroutine get_electric_field_state


  !----------------------------------------------------------
  subroutine get_magnetic_field_state(rs_state, mesh, rs_sign, magnetic_field, mu_field, np)
    CMPLX,             intent(in)    :: rs_state(:,:)
    type(mesh_t),      intent(in)    :: mesh
    integer,           intent(in)    :: rs_sign
    FLOAT,             intent(inout) :: magnetic_field(:,:)
    FLOAT,   optional, intent(in)    :: mu_field(:)
    integer, optional, intent(in)    :: np

    integer :: ip, np_
    type(profile_t), save :: prof

    PUSH_SUB(get_magnetic_field_state)

    call profiling_in(prof, 'GET_MAGNETIC_FIELD_STATE')

    np_ = optional_default(np, mesh%np)

    do ip = 1, np_
      if (present(mu_field)) then
        magnetic_field(ip, :) = sqrt(M_TWO*mu_field(ip)) * rs_sign * aimag(rs_state(ip, :))
      else
        magnetic_field(ip, :) = sqrt(M_TWO*P_mu) * rs_sign * aimag(rs_state(ip, :))
      end if
   end do

   call profiling_out(prof)

   POP_SUB(get_magnetic_field_state)

  end subroutine get_magnetic_field_state


  !----------------------------------------------------------
  subroutine get_current_element(rs_current_element, current_element, ep_element)
    CMPLX,           intent(in)    :: rs_current_element
    FLOAT,           intent(inout) :: current_element
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      current_element = sqrt(M_TWO*ep_element) * TOFLOAT(rs_current_element)
    else
      current_element = sqrt(M_TWO*P_ep) * TOFLOAT(rs_current_element)
    end if

  end subroutine get_current_element


  !----------------------------------------------------------
  subroutine get_current_vector(rs_current_vector, current_vector, ep_element)
    CMPLX,           intent(in)    :: rs_current_vector(:)
    FLOAT,           intent(inout) :: current_vector(:)
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      current_vector(:) = sqrt(M_TWO*ep_element) * TOFLOAT(rs_current_vector(:))
    else
      current_vector(:) = sqrt(M_TWO*P_ep) * TOFLOAT(rs_current_vector(:))
    end if

  end subroutine get_current_vector


  !----------------------------------------------------------
  subroutine get_current_state(rs_current_field, current_field, mesh, ep_field, np)
    CMPLX,             intent(in)    :: rs_current_field(:,:)
    FLOAT,             intent(inout) :: current_field(:,:)
    FLOAT,   optional, intent(in)    :: ep_field(:)
    type(mesh_t),      intent(in)    :: mesh
    integer, optional, intent(in)    :: np

    integer :: ip, np_

    PUSH_SUB(get_current_state)

    np_ = optional_default(np, mesh%np)

    do ip = 1, np_
      if (present(ep_field)) then
        current_field(ip, :) = sqrt(M_TWO*ep_field(ip)) * TOFLOAT(rs_current_field(ip, :))
      else
        current_field(ip, :) = sqrt(M_TWO*P_ep) * TOFLOAT(rs_current_field(ip, :))
      end if
    end do

    POP_SUB(get_current_state)

  end subroutine get_current_state


  !----------------------------------------------------------
  subroutine get_rs_state_at_point(rs_state_point, rs_state, pos, st, mesh)

    CMPLX,               intent(inout)   :: rs_state_point(:,:)
    CMPLX,               intent(in)      :: rs_state(:,:)
    FLOAT,               intent(in)      :: pos(:,:)
    type(states_mxll_t), intent(in)      :: st
    type(mesh_t),        intent(in)      :: mesh

    integer :: ip, pos_index, rankmin
    FLOAT   :: dmin
    CMPLX   :: ztmp(mesh%sb%dim)
    CMPLX, allocatable :: ztmp_global(:)

    PUSH_SUB(get_rs_state_at_point)

    SAFE_ALLOCATE(ztmp_global(mesh%np_global))

    do ip = 1, st%selected_points_number
      pos_index = mesh_nearest_point(mesh, pos(:,ip), dmin, rankmin)
      ztmp(:) = rs_state(pos_index, :)
      if (mesh%parallel_in_domains) then
#ifdef HAVE_MPI
        call MPI_Bcast(ztmp, st%dim, MPI_CMPLX, rankmin, mesh%mpi_grp%comm, mpi_err)
#endif
      end if
      rs_state_point(:, ip) = ztmp(:)
    end do

    SAFE_DEALLOCATE_A(ztmp_global)

    POP_SUB(get_rs_state_at_point)
  end subroutine get_rs_state_at_point


  !----------------------------------------------------------
  subroutine get_divergence_field(gr, field, field_div, charge_density)
    type(grid_t),    intent(in)    :: gr
    FLOAT,           intent(inout) :: field(:,:)
    FLOAT,           intent(inout) :: field_div(:)
    logical,         intent(in)    :: charge_density

    PUSH_SUB(get_divergence_field)

    call dderivatives_div(gr%der, field, field_div)

    if (optional_default(charge_density,.false.)) then
      field_div = P_ep * field_div
    end if

    POP_SUB(get_divergence_field)
  end subroutine get_divergence_field


  ! ---------------------------------------------------------
  subroutine get_poynting_vector(gr, st, rs_state, rs_sign, poynting_vector, ep_field, mu_field, mean_value)
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    CMPLX,                    intent(in)    :: rs_state(:,:)
    integer,                  intent(in)    :: rs_sign
    FLOAT,                    intent(inout) :: poynting_vector(:,:)
    FLOAT,          optional, intent(in)    :: ep_field(:)
    FLOAT,          optional, intent(in)    :: mu_field(:)
    FLOAT,          optional, intent(inout) :: mean_value(:)

    integer            :: ip, ip_in

    PUSH_SUB(get_poynting_vector)

    if (present(ep_field) .and. present(mu_field)) then
      do ip = 1, gr%mesh%np
        poynting_vector(ip, 1:3) = M_ONE/mu_field(ip) * sqrt(M_TWO/ep_field(ip)) &
                                * sqrt(M_TWO*mu_field(ip)) &
                                * dcross_product(TOFLOAT(rs_state(ip, 1:3)), &
                                rs_sign*aimag(rs_state(ip,1:3)))
      end do
    else
      do ip = 1, gr%mesh%np
        poynting_vector(ip, 1:3) = M_ONE/st%mu(ip) * sqrt(M_TWO/st%ep(ip)) &
                                * sqrt(M_TWO*st%mu(ip)) &
                                * dcross_product(TOFLOAT(rs_state(ip, 1:3)), &
                                rs_sign*aimag(rs_state(ip, 1:3)))
      end do
    end if

    if (present (mean_value)) then
      ASSERT(.not. gr%mesh%use_curvilinear)
      mean_value = M_ZERO
      do ip_in = 1, st%inner_points_number
        ip = st%inner_points_map(ip_in)
        mean_value(1:3)   = mean_value(1:3) + poynting_vector(ip, 1:3)
      end do
      mean_value(1:3) = mean_value(1:3) * gr%mesh%volume_element
      if(gr%mesh%parallel_in_domains) then
        call gr%mesh%allreduce(mean_value(1:3))
      end if
    end if

    POP_SUB(get_poynting_vector)
  end subroutine get_poynting_vector


  ! ---------------------------------------------------------
  subroutine get_poynting_vector_plane_waves(gr, st, rs_sign, poynting_vector, mean_value)
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,                  intent(in)    :: rs_sign
    FLOAT,                    intent(inout) :: poynting_vector(:,:)
    FLOAT,          optional, intent(inout) :: mean_value(:)

    integer            :: ip, ip_in

    PUSH_SUB(get_poynting_vector_plane_waves)

    do ip = 1, gr%mesh%np
      poynting_vector(ip, :) = M_ONE/P_mu * sqrt(M_TWO/P_ep) * sqrt(M_TWO*P_mu) &
               & * dcross_product(TOFLOAT(st%rs_state_plane_waves(ip,:)), &
               & rs_sign*aimag(st%rs_state_plane_waves(ip,:)))
    end do

    if (present(mean_value)) then
      mean_value = M_ZERO
      do ip_in = 1, st%inner_points_number
        ip = st%inner_points_map(ip_in)
        mean_value(:)   = mean_value(:) + poynting_vector(ip,:)
      end do
      mean_value(:) = mean_value(:) * gr%mesh%volume_element
      if(gr%mesh%parallel_in_domains) then
        call gr%mesh%allreduce(mean_value(:))
      end if
    end if

    POP_SUB(get_poynting_vector_plane_waves)
  end subroutine get_poynting_vector_plane_waves


  !----------------------------------------------------------
  subroutine state_diff(rs_state_old, rs_state, gr, diff)
    CMPLX,        intent(in) :: rs_state(:,:), rs_state_old(:,:)
    type(grid_t), intent(in) :: gr

    CMPLX, allocatable :: ztmp(:), ztmp_rs(:,:)
    FLOAT              :: d, diff
    integer            :: idim

    PUSH_SUB(state_diff)

    SAFE_ALLOCATE(ztmp_rs(1:gr%mesh%np,1:gr%sb%dim))
    SAFE_ALLOCATE(ztmp(1:gr%mesh%np))

    ztmp_rs(1:gr%mesh%np,:) =  rs_state(1:gr%mesh%np,:)

    ztmp(:) = TOCMPLX(M_ZERO,M_ZERO)
    do idim = 1, gr%sb%dim
      ztmp(1:gr%mesh%np) = ztmp(1:gr%mesh%np) + &
                                   abs(rs_state_old(1:gr%mesh%np, idim)-ztmp_rs(1:gr%mesh%np, idim))
      d = zmf_nrm2(gr%mesh, ztmp)
    end do
    if(d > diff) diff = d

    SAFE_DEALLOCATE_A(ztmp_rs)
    SAFE_DEALLOCATE_A(ztmp)

    POP_SUB(state_diff)
  end subroutine state_diff


end module states_mxll_oct_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
