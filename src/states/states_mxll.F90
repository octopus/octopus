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
  use blacs_proc_grid_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use geometry_oct_m
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
    states_mxll_null,                 &
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

  type states_mxll_priv_t
    private
    type(type_t) :: wfs_type             !< always complex (TYPE_CMPLX) states
  end type states_mxll_priv_t

  type states_mxll_t
    ! Components are public by default
    type(states_elec_dim_t)      :: d
    type(states_mxll_priv_t)     :: priv          !< the private components
    integer                      :: nst           !< Number of states in each irreducible subspace
    integer                      :: rs_sign
!    type(states_elec_group_t)    :: group
    logical                      :: parallel_in_states !< Am I parallel in states?

    type(batch_t), pointer       :: rsb
    type(batch_t), pointer       :: rs_transb
    type(batch_t), pointer       :: rs_longb
    type(batch_t), pointer       :: rs_curr_dens_rest1b
    type(batch_t), pointer       :: rs_curr_dens_rest2b
    
    CMPLX, pointer               :: rs_state_plane_waves(:,:)
!   CMPLX, pointer               :: rs_state(:,:)
!    CMPLX, pointer              :: rs_state_trans(:,:)
!    CMPLX, pointer              :: rs_state_long(:,:)
    
    logical                      :: rs_current_density_restart = .false.
!    CMPLX, pointer              :: rs_current_density_restart_t1(:,:)
!    CMPLX, pointer              :: rs_current_density_restart_t2(:,:)

    FLOAT, pointer               :: ep(:)
    FLOAT, pointer               :: mu(:)

    integer, pointer             :: rs_state_fft_map(:,:,:)
    integer, pointer             :: rs_state_fft_map_inv(:,:)

    FLOAT, pointer               :: energy_rate(:)
    FLOAT, pointer               :: delta_energy(:)
    FLOAT, pointer               :: energy_via_flux_calc(:)

    FLOAT, pointer               :: trans_energy_rate(:)
    FLOAT, pointer               :: trans_delta_energy(:)
    FLOAT, pointer               :: trans_energy_via_flux_calc(:)

    FLOAT, pointer               :: plane_waves_energy_rate(:)
    FLOAT, pointer               :: plane_waves_delta_energy(:)
    FLOAT, pointer               :: plane_waves_energy_via_flux_calc(:)

    FLOAT                        :: poynting_vector_box_surface(1:2,1:3,1:3) = M_ZERO
    FLOAT                        :: poynting_vector_box_surface_plane_waves(1:2,1:3,1:3) = M_ZERO
    FLOAT                        :: electric_field_box_surface(1:2,1:3,1:3) = M_ZERO
    FLOAT                        :: electric_field_box_surface_plane_waves(1:2,1:3,1:3) = M_ZERO
    FLOAT                        :: magnetic_field_box_surface(1:2,1:3,1:3) = M_ZERO
    FLOAT                        :: magnetic_field_box_surface_plane_waves(1:2,1:3,1:3) = M_ZERO

    logical                      :: rs_state_const_external = .false.
    CMPLX, pointer               :: rs_state_const(:)
    CMPLX, pointer               :: rs_state_const_amp(:,:)
    type(tdf_t), pointer         :: rs_state_const_td_function(:)

    FLOAT                        :: poynting_mean(3)
    FLOAT                        :: poynting_mean_plane_waves(3)

    integer                      :: inner_points_number
    integer, pointer             :: inner_points_map(:)
    integer                      :: boundary_points_number
    integer, pointer             :: boundary_points_map(:)

    integer                      :: surface_points_number(MAX_DIM)
    integer, pointer             :: surface_points_map(:,:,:)
    FLOAT                        :: surface_element(MAX_DIM)

    integer                      :: surface_grid_rows_number(MAX_DIM)
    integer, pointer             :: surface_grid_points_number(:,:,:)
    integer, pointer             :: surface_grid_points_map(:,:,:,:,:)
    integer, pointer             :: surface_grid_center(:,:,:,:)
    FLOAT                        :: surface_grid_element(MAX_DIM)

    type(mesh_plane_t)           :: surface(2,3)

    integer                      :: selected_points_number
    FLOAT, pointer               :: selected_points_coordinate(:,:)
    CMPLX, pointer               :: selected_points_rs_state(:,:)
    CMPLX, pointer               :: selected_points_rs_state_trans(:,:)
    FLOAT                        :: rs_state_trans_var

    FLOAT, pointer               :: grid_rho(:,:)
    CMPLX, pointer               :: kappa_psi(:,:)

    character(len=1024), pointer :: user_def_e_field(:)
    character(len=1024), pointer :: user_def_b_field(:)

    integer                      :: energy_incident_waves_calc_iter
    logical                      :: energy_incident_waves_calc

    ! external current variables
    integer                      :: external_current_number
    integer,             pointer :: external_current_modus(:)
    character(len=1024), pointer :: external_current_string(:,:)
    FLOAT,               pointer :: external_current_amplitude(:,:,:)
    type(tdf_t),         pointer :: external_current_td_function(:)
    type(tdf_t),         pointer :: external_current_td_phase(:)
    FLOAT,               pointer :: external_current_omega(:)
    FLOAT,               pointer :: external_current_phase(:)

    !> used for the user-defined wavefunctions (they are stored as formula strings)
    !! (st%d%dim, st%nst, st%d%nik)
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
    integer, pointer            :: node(:)
    logical, private            :: packed
  end type states_mxll_t

contains

  ! ---------------------------------------------------------
  subroutine states_mxll_null(st)
    type(states_mxll_t), intent(inout) :: st

    PUSH_SUB(states_mxll_null)

    call states_elec_dim_null(st%d)
!    call states_elec_group_null(st%group)
    call distributed_nullify(st%dist) 
    st%priv%wfs_type = TYPE_CMPLX
    st%d%orth_method = 0
    st%parallel_in_states = .false.
#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_nullify(st%dom_st_proc_grid)
#endif
    nullify(st%node)

    POP_SUB(states_mxll_null)
  end subroutine states_mxll_null

  ! ---------------------------------------------------------
  subroutine states_mxll_init(st, namespace, gr, geo)
    type(states_mxll_t), target, intent(inout) :: st
    type(namespace_t),           intent(in)    :: namespace
    type(grid_t),                intent(in)    :: gr
    type(geometry_t),            intent(in)    :: geo
    type(block_t)        :: blk
    integer :: idim, nlines, ncols, il
    FLOAT   :: pos(MAX_DIM)

    PUSH_SUB(states_mxll_init)

    st%fromScratch = .true. ! this will be reset if restart_read is called
    call states_mxll_null(st)
    
    st%d%dim = 3
    st%nst   = 1
    st%d%ispin = UNPOLARIZED
    st%d%nspin = 1
    st%d%spin_channels = 1
    call states_elec_choose_kpoints(st%d, gr%sb, namespace)

    SAFE_ALLOCATE(st%user_def_e_field(1:st%d%dim))
    SAFE_ALLOCATE(st%user_def_b_field(1:st%d%dim))

    st%st_start = 1
    st%st_end = st%nst
    st%lnst = st%nst

    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(1:st%nst) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.
    st%packed = .false.
    
    st%d%block_size = 1    
    call distributed_nullify(st%d%kpt, st%d%nik)

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

    st%selected_points_number = 1
    if(parse_block(namespace, 'MaxwellFieldsCoordinate', blk) == 0) then
      nlines = parse_block_n(blk)
      st%selected_points_number = nlines
      SAFE_ALLOCATE(st%selected_points_coordinate(1:st%d%dim,1:nlines))
      SAFE_ALLOCATE(st%selected_points_rs_state(1:st%d%dim,1:nlines))
      SAFE_ALLOCATE(st%selected_points_rs_state_trans(1:st%d%dim,1:nlines))
      do il=1, nlines
        ncols = parse_block_cols(blk,0)
        if (ncols < 3 .or. ncols > 3) then
            message(1) = 'MaxwellFieldCoordinate must have 3 columns.'
            call messages_fatal(1, namespace=namespace)
        end if
        do idim=1, st%d%dim
          call parse_block_float(blk, il-1, idim-1, pos(idim), units_inp%length)
        end do
        st%selected_points_coordinate(:,il) = pos
        st%selected_points_rs_state(:,il)  = M_z0
        st%selected_points_rs_state_trans(:,il) = M_z0
      end do
    else
      SAFE_ALLOCATE(st%selected_points_coordinate(1:st%d%dim,1))
      SAFE_ALLOCATE(st%selected_points_rs_state(1:st%d%dim,1))
      SAFE_ALLOCATE(st%selected_points_rs_state_trans(1:st%d%dim,1))
      st%selected_points_coordinate(:,:) = M_ZERO
      st%selected_points_rs_state(:,:) = M_z0
      st%selected_points_rs_state_trans(:,:) = M_z0
    end if

    POP_SUB(states_mxll_init)
      
  end subroutine states_mxll_init
  
  ! ---------------------------------------------------------
  !> Allocates the Maxwell states defined within a states_mxll_t structure.
  subroutine states_mxll_allocate(st, mesh)
    type(states_mxll_t),    intent(inout)   :: st
    type(mesh_t),           intent(in)      :: mesh

    PUSH_SUB(states_mxll_allocate)

    call zbatch_init(st%rsb, st%d%dim, 1, 1, mesh%np_part)
    call batch_set_zero(st%rsb)

    call zbatch_init(st%rs_transb, st%d%dim, 1, 1, mesh%np_part)
    call batch_set_zero(st%rs_transb)
 
    call zbatch_init(st%rs_longb, st%d%dim, 1, 1, mesh%np_part)
    call batch_set_zero(st%rs_longb)

    call zbatch_init(st%rs_curr_dens_rest1b, st%d%dim, 1, 1, mesh%np_part)
    call batch_set_zero(st%rs_curr_dens_rest1b)
    
    call zbatch_init(st%rs_curr_dens_rest2b, st%d%dim, 1, 1, mesh%np_part)
    call batch_set_zero(st%rs_curr_dens_rest2b)
   
!    Another alternative
!    call batch_init(st%rs_state_transb, hm%d%dim, 1, 1, st%rs_state_trans)
!    call st%rs_state_transb%end()

    POP_SUB(states_mxll_allocate)
  end subroutine states_mxll_allocate

  ! ---------------------------------------------------------
  subroutine states_mxll_end(st)
    type(states_mxll_t), intent(inout) :: st

    PUSH_SUB(states_mxll_end)

    call states_elec_dim_end(st%d)
    call st%rsb%end()
    call st%rs_transb%end()
    call st%rs_longb%end()
    call st%rs_curr_dens_rest1b%end()
    call st%rs_curr_dens_rest2b%end()

#ifdef HAVE_SCALAPACK
    call blacs_proc_grid_end(st%dom_st_proc_grid)
#endif

    call distributed_end(st%dist)
    SAFE_DEALLOCATE_P(st%node)

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

    PUSH_SUB(build_rs_state)

    do ip = 1, np_
      if (present(ep_field) .and. present(mu_field)) then
        rs_state(ip, :) = sqrt(ep_field(ip)/M_TWO) * e_field(ip, :) & 
                       + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*mu_field(ip))) * b_field(ip, :)
      else
        rs_state(ip, :) = sqrt(P_ep/M_TWO) * e_field(ip, :) &
                       + M_zI * rs_sign * sqrt(M_ONE/(M_TWO*P_mu)) * b_field(ip, :)
      end if 
    end do

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

    integer :: ip, np_

    ! no PUSH_SUB, called too often
    np_ = optional_default(np, mesh%np)
 
    do ip = 1, np_
      if (present(ep_field)) then
        rs_current_state(ip, :) = M_ONE/sqrt(M_TWO*ep_field(ip)) * current_state(ip, :)
      else
        rs_current_state(ip, :) = M_ONE/sqrt(M_TWO*P_ep) * current_state(ip, :)
      end if
    end do

  end subroutine build_rs_current_state


  !----------------------------------------------------------
  subroutine get_electric_field_vector(rs_state_vector, electric_field_vector, ep_element)
    CMPLX,             intent(in)    :: rs_state_vector(:)
    FLOAT,             intent(inout) :: electric_field_vector(:)
    FLOAT,   optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      electric_field_vector(:) = sqrt(M_TWO/ep_element) * real(rs_state_vector(:))
    else
      electric_field_vector(:) = sqrt(M_TWO/P_ep) * real(rs_state_vector(:))
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
  subroutine get_electric_field_state(rsb, mesh, electric_field, ep_field, np)
    type(batch_t),     intent(in)    :: rsb
    type(mesh_t),      intent(in)    :: mesh
    FLOAT,             intent(inout) :: electric_field(:,:)
    FLOAT,   optional, intent(in)    :: ep_field(:)
    integer, optional, intent(in)    :: np

    CMPLX, allocatable :: rs_aux(:,:)
    integer :: ip, ii, np_

    PUSH_SUB(get_electric_field_state)

    np_ = optional_default(np, mesh%np)
    SAFE_ALLOCATE(rs_aux(1:np_, 1:3))
    
    do ii = 1, 3
       call batch_get_state(rsb, np_, ii, rs_aux(:, ii))
    end do
     
    do ip = 1, np_
      if (present(ep_field)) then
        electric_field(ip, :) = sqrt(M_TWO/ep_field(ip)) * real(rs_aux(ip, :), REAL_PRECISION)
      else 
        electric_field(ip,:) = sqrt(M_TWO/P_ep) * real(rs_aux(ip, :), REAL_PRECISION)
      end if
    end do

    SAFE_DEALLOCATE_A(rs_aux)

    POP_SUB(get_electric_field_state)

  end subroutine get_electric_field_state


  !----------------------------------------------------------
  subroutine get_magnetic_field_state(rsb, mesh, rs_sign, magnetic_field, mu_field, np)
    type(batch_t),     intent(in)    :: rsb
    type(mesh_t),      intent(in)    :: mesh
    integer,           intent(in)    :: rs_sign
    FLOAT,             intent(inout) :: magnetic_field(:,:)
    FLOAT,   optional, intent(in)    :: mu_field(:)
    integer, optional, intent(in)    :: np

    CMPLX, allocatable :: rs_aux(:,:)
    integer :: ip, ii, np_

    PUSH_SUB(get_magnetic_field_state)

    np_ = optional_default(np, mesh%np)
    SAFE_ALLOCATE(rs_aux(1:np_, 1:3))

    do ii = 1, 3
      call batch_get_state(rsb, np, ii, rs_aux(:, ii))
    end do

    
    do ip = 1, np_
      if (present(mu_field)) then
        magnetic_field(ip, :) = sqrt(M_TWO*mu_field(ip)) * rs_sign * aimag(rs_aux(ip, :))
      else
        magnetic_field(ip, :) = sqrt(M_TWO*P_mu) * rs_sign * aimag(rs_aux(ip, :))
      end if
   end do

   SAFE_DEALLOCATE_A(rs_aux)

   POP_SUB(get_magnetic_field_state)

  end subroutine get_magnetic_field_state


  !----------------------------------------------------------
  subroutine get_current_element(rs_current_element, current_element, ep_element)
    CMPLX,           intent(in)    :: rs_current_element
    FLOAT,           intent(inout) :: current_element
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      current_element = sqrt(M_TWO*ep_element) * real(rs_current_element, REAL_PRECISION)
    else
      current_element = sqrt(M_TWO*P_ep) * real(rs_current_element, REAL_PRECISION)
    end if

  end subroutine get_current_element


  !----------------------------------------------------------
  subroutine get_current_vector(rs_current_vector, current_vector, ep_element)
    CMPLX,           intent(in)    :: rs_current_vector(:)
    FLOAT,           intent(inout) :: current_vector(:)
    FLOAT, optional, intent(in)    :: ep_element

    ! no PUSH_SUB, called too often

    if (present(ep_element)) then
      current_vector(:) = sqrt(M_TWO*ep_element) * real(rs_current_vector(:), REAL_PRECISION)
    else
      current_vector(:) = sqrt(M_TWO*P_ep) * real(rs_current_vector(:), REAL_PRECISION)
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
        current_field(ip, :) = sqrt(M_TWO*ep_field(ip)) * real(rs_current_field(ip, :), REAL_PRECISION)
      else
        current_field(ip, :) = sqrt(M_TWO*P_ep) * real(rs_current_field(ip, :), REAL_PRECISION)
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

    integer :: ip, pos_index_local, pos_index_global, rankmin
    FLOAT   :: dmin
    CMPLX   :: ztmp(MAX_DIM)
    CMPLX, allocatable :: ztmp_global(:)

    PUSH_SUB(get_rs_state_at_point)

    SAFE_ALLOCATE(ztmp_global(mesh%np_global))

    do ip = 1, st%selected_points_number
      call mesh_nearest_point_infos(mesh, pos(:,ip), dmin, rankmin, pos_index_local, pos_index_global)
!      pos_index = mesh_nearest_point(mesh, pos(:,ip), dmin, rankmin)
      if (mesh%parallel_in_domains) then
        ztmp(:) = rs_state(pos_index_local,:)
#ifdef HAVE_MPI
        call MPI_Bcast(ztmp, st%d%dim, MPI_CMPLX, rankmin, mesh%mpi_grp%comm, mpi_err)
        call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif
      else
        ztmp(:) = rs_state(pos_index_global, :)
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
  subroutine get_poynting_vector(gr, st, rsb, rs_sign, poynting_vector, ep_field, mu_field)
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    type(batch_t),            intent(in)    :: rsb
    integer,                  intent(in)    :: rs_sign
    FLOAT,                    intent(inout) :: poynting_vector(:,:)
    FLOAT,          optional, intent(in)    :: ep_field(:)
    FLOAT,          optional, intent(in)    :: mu_field(:)

    integer            :: ip, ii
    CMPLX, allocatable :: rs_aux(:,:)

    PUSH_SUB(get_poynting_vector)

    SAFE_ALLOCATE(rs_aux(1:gr%mesh%np, 1:3))
    do ii = 1, 3
       call batch_get_state(rsb, gr%mesh%np, ii, rs_aux(:, ii))
    end do

    if (present(ep_field) .and. present(mu_field)) then
      do ip = 1, gr%mesh%np
        poynting_vector(ip, :) = M_ONE/mu_field(ip) * sqrt(M_TWO/ep_field(ip)) &
                              * sqrt(M_TWO*mu_field(ip)) &
                              * dcross_product(real(rs_aux(ip, :), REAL_PRECISION), rs_sign*aimag(rs_aux(ip, :)))
      end do
    else
      do ip = 1, gr%mesh%np
        poynting_vector(ip,:) = M_ONE/st%mu(ip) * sqrt(M_TWO/st%ep(ip)) &
                              * sqrt(M_TWO*st%mu(ip)) &
                              * dcross_product(real(rs_aux(ip, :), REAL_PRECISION), rs_sign*aimag(rs_aux(ip, :)))
      end do
    end if

    POP_SUB(get_poynting_vector)
  end subroutine get_poynting_vector


  ! ---------------------------------------------------------
  subroutine get_poynting_vector_plane_waves(gr, st, rs_sign, poynting_vector)
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(in)    :: st
    integer,                  intent(in)    :: rs_sign
    FLOAT,                    intent(inout) :: poynting_vector(:,:)

    integer            :: ip

    PUSH_SUB(get_poynting_vector_plane_waves)

    do ip = 1, gr%mesh%np
      poynting_vector(ip, :) = M_ONE/P_mu * sqrt(M_TWO/P_ep) * sqrt(M_TWO*P_mu) &
               & * dcross_product(real(st%rs_state_plane_waves(ip,:), REAL_PRECISION), &
               & rs_sign*aimag(st%rs_state_plane_waves(ip,:)))
    end do

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
    do idim=1, gr%sb%dim
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
