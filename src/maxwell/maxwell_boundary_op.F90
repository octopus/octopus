!! Copyright (C) 2019 R. Jestaedt, H. Appel, F. Bonafe, M. Oliveira, N. Tancogne-Dejean
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

#include "global.h"

module maxwell_boundary_op_oct_m
  use box_sphere_oct_m
  use box_parallelepiped_oct_m
  use cube_function_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use maxwell_function_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use string_oct_m
  use system_linear_medium_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use simul_box_oct_m
  use varinfo_oct_m
  use space_oct_m
  use states_mxll_oct_m

  implicit none

  private
  public ::                    &
    bc_mxll_init,              &
    bc_mxll_end,               &
    bc_mxll_t,                 &
    inner_and_outer_points_mapping,  &
    surface_grid_points_mapping


  type pml_t
    FLOAT                :: width
    integer              :: points_number
    integer, allocatable :: points_map(:)
    integer, allocatable :: points_map_inv(:)
    FLOAT                :: power
    FLOAT                :: refl_error
    FLOAT,   allocatable :: kappa(:,:)
    FLOAT,   allocatable :: sigma_e(:,:)
    FLOAT,   allocatable :: sigma_m(:,:)
    CMPLX,   allocatable :: a(:,:)
    CMPLX,   allocatable :: b(:,:)
    FLOAT,   allocatable :: c(:,:)
    FLOAT,   allocatable :: mask(:)
    CMPLX,   allocatable :: aux_ep(:,:,:)
    CMPLX,   allocatable :: aux_mu(:,:,:)
    CMPLX,   allocatable :: conv_plus(:,:,:)
    CMPLX,   allocatable :: conv_minus(:,:,:)
    CMPLX,   allocatable :: conv_plus_old(:,:,:)
    CMPLX,   allocatable :: conv_minus_old(:,:,:)
  end type pml_t

  type plane_wave_t
    integer                          :: points_number  !< number of points of plane wave boundary
    integer,             allocatable :: points_map(:) !< points map for plane waves boundary
    integer                          :: number !< number of plane waves given by user
    integer,             allocatable :: modus(:) !< input file modus, either parser or Maxwell function
    character(len=1024), allocatable :: e_field_string(:,:) !< string in case of parser
    FLOAT,               allocatable :: k_vector(:,:) !< k vector for each plane wave
    FLOAT,               allocatable :: v_vector(:,:) !< velocity vector for each plane wave
    CMPLX,               allocatable :: e_field(:,:) !< field amplitude for each plane wave
    type(mxf_t),         allocatable :: mx_function(:) !< Maxwell function for each plane wave
  end type plane_wave_t

  type bc_mxll_t
    integer              :: bc_type(MAX_DIM)
    integer              :: bc_ab_type(MAX_DIM)
    FLOAT                :: bc_bounds(2, MAX_DIM)
    logical              :: ab_user_def
    FLOAT,   allocatable :: ab_ufn(:)

    FLOAT                :: ab_width
    FLOAT                :: mask_width
    integer              :: mask_points_number(MAX_DIM)
    integer, allocatable :: mask_points_map(:,:)
    FLOAT,   allocatable :: mask(:,:)

    integer              :: der_bndry_mask_points_number
    integer, allocatable :: der_bndry_mask_points_map(:)
    FLOAT,   allocatable :: der_bndry_mask(:)

    type(pml_t)          :: pml       !< attributes of PML absorbing boundaries
    type(single_medium_box_t)   :: medium(3)    !< attributes of linear medium boundaries

    integer              :: constant_points_number
    integer, allocatable :: constant_points_map(:)
    CMPLX,   allocatable :: constant_rs_state(:,:)

    integer              :: mirror_points_number(3)
    integer, allocatable :: mirror_points_map(:,:)

    logical              :: do_plane_waves = .false.
    type(plane_wave_t)   :: plane_wave

    FLOAT                :: zero_width
    integer              :: zero_points_number(MAX_DIM)
    integer, allocatable :: zero_points_map(:,:)
    FLOAT,   allocatable :: zero(:,:)
  end type bc_mxll_t

  integer, public, parameter ::   &
    MXLL_BC_ZERO          = 0,    &
    MXLL_BC_CONSTANT      = 1,    &
    MXLL_BC_MIRROR_PEC    = 2,    &
    MXLL_BC_MIRROR_PMC    = 3,    &
    MXLL_BC_PLANE_WAVES   = 4,    &
    MXLL_BC_PERIODIC      = 5,    &
    MXLL_BC_MEDIUM        = 6

  integer, public, parameter ::   &
    MXLL_AB_NOT_ABSORBING = 0,    &
    MXLL_AB_MASK          = 1,    &
    MXLL_AB_CPML          = 2,    &
    MXLL_AB_MASK_ZERO     = 7

contains

  ! ---------------------------------------------------------
  subroutine bc_mxll_init(bc, namespace, space, gr, st, sb, dt)
    type(bc_mxll_t),          intent(inout) :: bc
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(inout) :: st
    type(simul_box_t),        intent(in)    :: sb
    FLOAT, optional,          intent(in)    :: dt

    integer             :: idim, nlines, icol, ncols, ab_shape_dim
    FLOAT               :: bounds(1:2,MAX_DIM), ab_bounds(1:2,MAX_DIM)
    type(block_t)       :: blk
    character(len=1024) :: string
    character(len=50)   :: str, ab_type_str
    logical             :: plane_waves_check = .false., ab_mask_check = .false., ab_pml_check = .false.
    logical             :: constant_check = .false., zero_check = .false.
    type(profile_t), save :: prof
    FLOAT :: ep_factor, mu_factor, sigma_e_factor, sigma_m_factor

    PUSH_SUB(bc_mxll_init)

    call profiling_in(prof, 'BC_MXLL_INIT')

    bc%ab_user_def = .false.
    bc%bc_ab_type(:) = MXLL_AB_NOT_ABSORBING ! default option

    !%Variable MaxwellAbsorbingBoundaries
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Type of absorbing boundaries used for Maxwell propagation in each direction.
    !%
    !% Example:
    !%
    !% <tt>%MaxwellAbsorbingBoundaries
    !% <br>&nbsp;&nbsp;   cpml | cpml | cpml
    !% <br>%</tt>
    !%
    !%Option not_absorbing 0
    !% No absorbing boundaries.
    !%Option mask 1
    !% A mask equal to the wavefunctions mask is applied to the Maxwell states at the boundaries
    !%Option cpml 2
    !% Perfectly matched layer absorbing boundary
    !%Option mask_zero 7
    !% Absorbing boundary region is set to zero
    !%End

    if(parse_block(namespace, 'MaxwellAbsorbingBoundaries', blk) == 0) then
      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      if (nlines /= 1) then
        message(1) = 'MaxwellAbsorbingBounaries has to consist of one line!'
        call messages_fatal(1, namespace=namespace)
      end if

      ncols = parse_block_cols(blk, 0)
      if (ncols /= 3) then
        message(1) = 'MaxwellAbsorbingBoundaries has to consist of three columns!'
        call messages_fatal(1, namespace=namespace)
      end if

      do icol = 1, ncols
        call parse_block_integer(blk, 0, icol-1, bc%bc_ab_type(icol))
      end do

    call parse_block_end(blk)
    end if

    do idim = 1, 3
      if (bc%bc_ab_type(idim) == MXLL_AB_MASK) ab_mask_check = .true.
      if (bc%bc_ab_type(idim) == MXLL_AB_CPML) ab_pml_check = .true.
      if (bc%bc_type(idim) == MXLL_BC_CONSTANT) constant_check = .true.
      if (bc%bc_type(idim) == MXLL_BC_ZERO) zero_check = .true.
    end do

    if (ab_mask_check .or. ab_pml_check) then

      write(str, '(a,i5)') 'Maxwell Absorbing Boundaries'
      call messages_print_stress(stdout, trim(str), namespace=namespace)

    end if

    do idim = 1, st%dim
      select case (bc%bc_type(idim))

      case (MXLL_BC_ZERO, MXLL_BC_MIRROR_PEC, MXLL_BC_MIRROR_PMC)

        bounds(1, idim) = (gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim))*gr%mesh%spacing(idim)
        bounds(2, idim) = (gr%mesh%idx%nr(2, idim)) * gr%mesh%spacing(idim)

      case (MXLL_BC_CONSTANT, MXLL_BC_PERIODIC)

        bounds(1, idim) = (gr%mesh%idx%nr(2, idim) - 2*gr%mesh%idx%enlarge(idim))*gr%mesh%spacing(idim)
        bounds(2, idim) = (gr%mesh%idx%nr(2, idim)) * gr%mesh%spacing(idim)

        if (bc%bc_type(idim) == MXLL_BC_PERIODIC .and. sb%latt%nonorthogonal) then
          message(1) = "Maxwell propagation does not work for non-orthogonal cells in periodic systems"
          call messages_fatal(1, namespace=namespace)
        end if

      case (MXLL_BC_PLANE_WAVES)

        bounds(1, idim) = (gr%mesh%idx%nr(2, idim) - 2*gr%mesh%idx%enlarge(idim))*gr%mesh%spacing(idim)
        bounds(2, idim) = (gr%mesh%idx%nr(2, idim)) * gr%mesh%spacing(idim)
        plane_waves_check = .true.
        bc%do_plane_waves = .true.

      case (MXLL_BC_MEDIUM)
        call bc_mxll_medium_init(gr, namespace, bounds, idim, ep_factor, mu_factor, sigma_e_factor, sigma_m_factor)
        call maxwell_medium_points_mapping(bc, gr%mesh, st, bounds)
        call bc_mxll_generate_medium(bc, gr, bounds, ep_factor, mu_factor, sigma_e_factor, sigma_m_factor)

      end select

      select type (box => gr%sb%box)
      type is (box_sphere_t)
        ab_shape_dim = 1
        if (space%is_periodic()) then
          message(1) = "Sphere box shape can only work for non-periodic systems"
          call messages_fatal(1, namespace=namespace)
        end if
      type is (box_parallelepiped_t)
        ab_shape_dim = space%dim
        ab_bounds(1, idim) = bounds(1, idim)
        ab_bounds(2, idim) = bounds(1, idim)
      class default
        message(1) = "Box shape for Maxwell propagation not supported yet"
        call messages_fatal(1, namespace=namespace)
      end select

      if (bc%bc_ab_type(idim) /= MXLL_AB_NOT_ABSORBING) then

        call messages_print_var_option(stdout, "MaxwellAbsorbingBoundaries", bc%bc_ab_type(idim))

        select case (bc%bc_ab_type(idim))
        case (MXLL_AB_MASK_ZERO)
          if (bc%bc_type(idim) == MXLL_BC_PERIODIC) then
            message(1) = "Zero absorbing boundary conditions do not work in periodic directions"
            call messages_fatal(1, namespace=namespace)
          end if

          call bc_mxll_ab_bounds_init(bc, gr, namespace, bounds, ab_bounds, idim)
          bc%zero_width = bc%ab_width

        case (MXLL_AB_MASK)
          call bc_mxll_ab_bounds_init(bc, gr, namespace, bounds, ab_bounds, idim)
          bc%mask_width = bc%ab_width

        case (MXLL_AB_CPML)
           call bc_mxll_pml_init(bc, gr, namespace, bounds, ab_bounds, idim)

        case default
          message(1) = "Absorbing boundary type not implemented for Maxwell propagation"
          call messages_fatal(1, namespace=namespace)
        end select

      end if

      select case (bc%bc_ab_type(idim))
      case (MXLL_AB_MASK, MXLL_AB_CPML, MXLL_AB_MASK_ZERO)
        bounds(1, idim) = ab_bounds(1, idim)
        bounds(2, idim) = bounds(2, idim)
        bc%bc_bounds(:, idim) = bounds(:, idim)
      case default
        bc%bc_bounds(:, idim) = bounds(:, idim)
      end select

      select type (box => gr%sb%box)
      type is (box_parallelepiped_t)
        select case (bc%bc_ab_type(idim))
        case(MXLL_AB_CPML)
          ab_type_str = "PML"
        case (MXLL_AB_MASK)
          ab_type_str = "Mask"
        case (MXLL_AB_MASK_ZERO)
          ab_type_str = "Zero"
        case default
          ab_type_str = ""
        end select

        if (bc%bc_ab_type(idim) == MXLL_AB_CPML .or. bc%bc_ab_type(idim) == MXLL_AB_MASK .or. &
             bc%bc_ab_type(idim) == MXLL_AB_MASK_ZERO) then
          string = trim(ab_type_str)//" Lower bound = "
          write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
               units_from_atomic(units_inp%length, ab_bounds(1, idim) ), ' [', &
               trim(units_abbrev(units_inp%length)), ']'
          write(message(1),'(a)') trim(string)

          string = trim(ab_type_str)//" Upper bound = "
          write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
               units_from_atomic(units_inp%length, ab_bounds(2, idim) ), ' [', &
               trim(units_abbrev(units_inp%length)), ']'

          write(message(2),'(a)') trim(string)
          call messages_info(2)
        end if

      class default

         write(message(1),'(a,es10.3,3a)') &
          "  Lower bound = ", units_from_atomic(units_inp%length, ab_bounds(1, idim) ),&
          ' [', trim(units_abbrev(units_inp%length)), ']'
        write(message(2),'(a,es10.3,3a)') & 
          "  Upper bound = ", units_from_atomic(units_inp%length, ab_bounds(2, idim) ),&
          ' [', trim(units_abbrev(units_inp%length)), ']'
        call messages_info(2)

      end select

    end do

    ! initialization of surfaces
    call maxwell_surfaces_init(gr%mesh, st, bounds)

    ! mapping of mask boundary points
    if (ab_mask_check) then
      call maxwell_mask_points_mapping(bc, gr%mesh, ab_bounds)
    end if

    ! mapping of pml boundary points
    if (ab_pml_check) then
      call maxwell_pml_points_mapping(bc, gr%mesh, ab_bounds)
    end if

    ! mapping of constant boundary points
    if (constant_check) then
      call maxwell_constant_points_mapping(bc, gr%mesh, bounds)
    end if

    ! mapping of plane waves boundary points
    if (plane_waves_check) then
      call maxwell_plane_waves_points_mapping(bc, gr%mesh, bounds)
      call maxwell_plane_waves_boundaries_init(bc, namespace)
    end if

    ! mapping of zero points
    if (zero_check) then
      call maxwell_zero_points_mapping(bc, gr%mesh, bounds)
    end if

    if (ab_mask_check) then
      call bc_mxll_generate_mask(bc, gr%mesh, ab_bounds)
    end if

    if (ab_pml_check) then
      call bc_mxll_generate_pml(bc%pml, gr, ab_bounds, dt)
    end if

    !call bc_generate_zero(bc, gr%mesh, ab_bounds)

    if(debug%info) call bc_mxll_write_info(bc, gr%mesh, namespace, space)

    if (ab_mask_check .or. ab_pml_check) then
      call messages_print_stress(stdout, namespace=namespace)
    end if

    call profiling_out(prof)

    POP_SUB(bc_mxll_init)
  end subroutine bc_mxll_init

  ! ---------------------------------------------------------
  subroutine bc_mxll_end(bc)
    type(bc_mxll_t),   intent(inout) :: bc

    integer :: idim

    PUSH_SUB(bc_mxll_end)

    SAFE_DEALLOCATE_A(bc%ab_ufn)

    SAFE_DEALLOCATE_A(bc%mask_points_map)
    SAFE_DEALLOCATE_A(bc%mask)

    SAFE_DEALLOCATE_A(bc%der_bndry_mask)
    SAFE_DEALLOCATE_A(bc%der_bndry_mask_points_map)

    call pml_end(bc%pml)
    do idim = 1, 3
      call single_medium_box_end(bc%medium(idim))
    end do

    SAFE_DEALLOCATE_A(bc%constant_points_map)
    SAFE_DEALLOCATE_A(bc%constant_rs_state)

    SAFE_DEALLOCATE_A(bc%mirror_points_map)

    call plane_wave_end(bc%plane_wave)

    SAFE_DEALLOCATE_A(bc%zero_points_map)
    SAFE_DEALLOCATE_A(bc%zero)

    POP_SUB(bc_mxll_end)
  end subroutine bc_mxll_end

  ! ---------------------------------------------------------
  subroutine pml_end(pml)
    type(pml_t),   intent(inout) :: pml

    PUSH_SUB(pml_end)

    SAFE_DEALLOCATE_A(pml%points_map)
    SAFE_DEALLOCATE_A(pml%points_map_inv)
    SAFE_DEALLOCATE_A(pml%kappa)
    SAFE_DEALLOCATE_A(pml%sigma_e)
    SAFE_DEALLOCATE_A(pml%sigma_m)
    SAFE_DEALLOCATE_A(pml%a)
    SAFE_DEALLOCATE_A(pml%b)
    SAFE_DEALLOCATE_A(pml%c)
    SAFE_DEALLOCATE_A(pml%mask)
    SAFE_DEALLOCATE_A(pml%aux_ep)
    SAFE_DEALLOCATE_A(pml%aux_mu)
    SAFE_DEALLOCATE_A(pml%conv_plus)
    SAFE_DEALLOCATE_A(pml%conv_minus)
    SAFE_DEALLOCATE_A(pml%conv_plus_old)
    SAFE_DEALLOCATE_A(pml%conv_minus_old)

    POP_SUB(pml_end)
  end subroutine pml_end

  ! ---------------------------------------------------------
  subroutine plane_wave_end(plane_wave)
    type(plane_wave_t),   intent(inout) :: plane_wave

    PUSH_SUB(plane_wave_end)

    SAFE_DEALLOCATE_A(plane_wave%points_map)
    SAFE_DEALLOCATE_A(plane_wave%modus)
    SAFE_DEALLOCATE_A(plane_wave%e_field_string)
    SAFE_DEALLOCATE_A(plane_wave%k_vector)
    SAFE_DEALLOCATE_A(plane_wave%v_vector)
    SAFE_DEALLOCATE_A(plane_wave%e_field)
    SAFE_DEALLOCATE_A(plane_wave%mx_function)

    POP_SUB(plane_wave_end)
  end subroutine plane_wave_end

  ! ---------------------------------------------------------
  subroutine bc_mxll_medium_init(gr, namespace, bounds, idim, ep_factor, mu_factor, sigma_e_factor, sigma_m_factor)
    type(grid_t),        intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    FLOAT,               intent(inout) :: bounds(:,:)
    integer,             intent(in)    :: idim
    FLOAT,               intent(out)   :: ep_factor
    FLOAT,               intent(out)   :: mu_factor
    FLOAT,               intent(out)   :: sigma_e_factor
    FLOAT,               intent(out)   :: sigma_m_factor

    FLOAT :: width
    type(profile_t), save :: prof

    PUSH_SUB(bc_mxll_medium_init)

    call profiling_in(prof, 'BC_MXLL_MEDIUM_INIT')

    !%Variable MediumWidth
    !%Type float
    !%Default 0.0 a.u.
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Width of the boundary region with medium
    !%End
    call parse_variable(namespace, 'MediumWidth', M_ZERO, width, units_inp%length)
    bounds(1,idim) = ( gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
    bounds(1,idim) = bounds(1,idim) - width
    bounds(2,idim) = ( gr%mesh%idx%nr(2, idim) ) * gr%mesh%spacing(idim)

    !%Variable MediumEpsilonFactor
    !%Type float
    !%Default 1.0.
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Linear medium electric susceptibility.
    !%End
    call parse_variable(namespace, 'MediumEpsilonFactor', M_ONE, ep_factor, unit_one)

    !%Variable MediumMuFactor
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Linear medium magnetic susceptibility.
    !%End
    call parse_variable(namespace, 'MediumMuFactor', M_ONE, mu_factor, unit_one)

    !%Variable MediumElectricSigma
    !%Type float
    !%Default 0.
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Electric conductivity of the linear medium.
    !%End

    call parse_variable(namespace, 'MediumElectricSigma', M_ZERO, sigma_e_factor, unit_one)
    !%Variable MediumMagneticSigma
    !%Type float
    !%Default 0.
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Magnetic conductivity of the linear medium.
    !%End
    call parse_variable(namespace, 'MediumMagneticSigma', M_ZERO, sigma_m_factor, unit_one)

    call profiling_out(prof)

    POP_SUB(bc_mxll_medium_init)
  end subroutine bc_mxll_medium_init

  ! ---------------------------------------------------------
  subroutine bc_mxll_ab_bounds_init(bc, gr, namespace, bounds, ab_bounds, idim)
    type(bc_mxll_t),     intent(inout) :: bc
    type(grid_t),        intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    FLOAT,               intent(inout) :: bounds(:,:), ab_bounds(:,:)
    integer,             intent(in)    :: idim

    FLOAT               :: default_width

    PUSH_SUB(bc_mxll_ab_bounds_init)

    !%Variable MaxwellABWidth
    !%Type float
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Width of the region used to apply the absorbing boundaries. The default value is twice
    !% the derivative order.
    !%End

    default_width = M_TWO * gr%der%order * maxval(gr%mesh%spacing(1:3))
    call parse_variable(namespace, 'MaxwellABWidth', default_width, bc%ab_width, units_inp%length)

    if (bc%ab_width < default_width) then
      bc%ab_width = default_width
      write(message(1),'(a)') 'Absorbing boundary width has to be larger or equal than twice the derivatives order times spacing!'
      write(message(2),'(a,es10.3)') 'Absorbing boundary width is set to: ', default_width
      call messages_info(2)
    end if

    ab_bounds(1, idim) = ab_bounds(2, idim) - bc%ab_width

    POP_SUB(bc_mxll_ab_bounds_init)
  end subroutine bc_mxll_ab_bounds_init

  ! ---------------------------------------------------------
  subroutine bc_mxll_pml_init(bc, gr, namespace, bounds, ab_bounds, idim)
    type(bc_mxll_t),     intent(inout) :: bc
    type(grid_t),        intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    FLOAT,               intent(inout) :: bounds(:,:), ab_bounds(:,:)
    integer,             intent(in)    :: idim

    PUSH_SUB(bc_mxll_pml_init)

    call bc_mxll_ab_bounds_init(bc, gr, namespace, bounds, ab_bounds, idim)
    bc%pml%width = bc%ab_width

    !%Variable MaxwellABPMLPower
    !%Type float
    !%Default 3.5
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Exponential of the polynomial profile for the non-physical conductivity of the PML.
    !% Should be between 2 and 4
    !%End
    call parse_variable(namespace, 'MaxwellABPMLPower', CNST(3.5), bc%pml%power, unit_one)

    !%Variable MaxwellABPMLReflectionError
    !%Type float
    !%Default 1.0e-16
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Tolerated reflection error for the PML
    !%End
    call parse_variable(namespace, 'MaxwellABPMLReflectionError', CNST(1.0e-16), bc%pml%refl_error, unit_one)

    POP_SUB(bc_mxll_pml_init)
  end subroutine bc_mxll_pml_init

  ! ---------------------------------------------------------
  subroutine bc_mxll_write_info(bc, mesh, namespace, space)
    type(bc_mxll_t),       intent(in) :: bc
    type(mesh_t),          intent(in) :: mesh
    type(namespace_t),     intent(in) :: namespace
    type(space_t),         intent(in) :: space

    integer :: err, idim, idim2
    FLOAT, allocatable :: tmp(:)
    logical :: mask_check, pml_check, medium_check
    character(1) :: dim_label(3)
    type(profile_t), save :: prof

    PUSH_SUB(bc_mxll_write_info)

    call profiling_in(prof, 'BC_MXLL_WRITE_INFO')

    mask_check = .false.
    pml_check = .false.
    medium_check = .false.

    dim_label = (/'x', 'y', 'z'/)

    do idim = 1, 3
      if (bc%bc_ab_type(idim) == MXLL_AB_MASK) then
        mask_check = .true.
      end if
      if (bc%bc_ab_type(idim) == MXLL_AB_CPML) then
        pml_check = .true.
      end if
      if (bc%bc_ab_type(idim) == MXLL_BC_MEDIUM) then
        medium_check = .true.
      end if
    end do

    if (mask_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      do idim = 1, 3
        tmp(:) = M_ZERO
        call get_mask_io_function(bc%mask, bc, tmp, idim)
        call write_files("maxwell_mask", tmp)
      end do
      SAFE_DEALLOCATE_A(tmp)
    else if (pml_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      do idim = 1, 3 
        ! sigma for electric field dim = idim
        tmp(:) = M_ONE
        call get_pml_io_function(bc%pml%sigma_e(:, idim), bc, tmp)
        call write_files("maxwell_sigma_e-"//dim_label(idim), tmp)

        ! sigma for magnetic dim = idim
        tmp(:) = M_ZERO
        call get_pml_io_function(bc%pml%sigma_m(:, 1), bc, tmp)
        call write_files("maxwell_sigma_m-"//dim_label(idim), tmp)

        ! pml_a for electric field dim = idim
        tmp(:) = M_ZERO
        call get_pml_io_function(TOFLOAT(bc%pml%a(:, idim)), bc, tmp)
        call write_files("maxwell_sigma_pml_a_e-"//dim_label(idim), tmp)

        ! pml_a for magnetic field dim = 1
        tmp(:) = M_ZERO
        call get_pml_io_function(aimag(bc%pml%a(:, idim)), bc, tmp)
        call write_files("maxwell_sigma_pml_a_m-"//dim_label(idim), tmp)
      end do
      SAFE_DEALLOCATE_A(tmp)
    end if
   
    if (medium_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      do idim = 1, 3
        ! medium epsilon
        tmp(:) = P_ep
        call get_medium_io_function(bc%medium(idim)%ep, bc, tmp, idim)
        call write_files("maxwell_ep"//dim_label(idim), tmp)
        ! medium mu
        tmp(:) = P_mu
        call get_medium_io_function(bc%medium(idim)%mu, bc, tmp, idim)
        call write_files("maxwell_mu"//dim_label(idim), tmp)
        ! medium epsilon
        tmp(:) = P_c
        call get_medium_io_function(bc%medium(idim)%c, bc, tmp, idim)
        call write_files("maxwell_c"//dim_label(idim), tmp)
        do idim2 = 1, 3
          ! medium epsilon aux field dim = idim
          tmp(:) = M_ZERO
          call get_medium_io_function(bc%medium(idim)%aux_ep(:, idim2), bc, tmp, idim)
          call write_files("maxwell_aux_ep-"//dim_label(idim)//"-"//dim_label(idim2), tmp)

          ! medium mu aux field dim = idim
          tmp(:) = M_ZERO
          call get_medium_io_function(bc%medium(idim)%aux_mu(:, idim2), bc, tmp, idim)
          call write_files("maxwell_aux_mu-"//dim_label(idim)//"-"//dim_label(idim2), tmp)
        end do
      end do

      SAFE_DEALLOCATE_A(tmp)
    end if

    call profiling_out(prof)

    POP_SUB(bc_mxll_write_info)
  contains

    subroutine get_pml_io_function(pml_func, bc, io_func)
      FLOAT,              intent(in)    :: pml_func(:)
      type(bc_mxll_t),    intent(in)    :: bc
      FLOAT,              intent(inout) :: io_func(:)

      integer :: ip, ip_in

      do ip_in = 1, bc%pml%points_number
        ip          = bc%pml%points_map(ip_in)
        io_func(ip) = pml_func(ip_in)
      end do

    end subroutine get_pml_io_function

    subroutine get_mask_io_function(mask_func, bc, io_func, idim)
      FLOAT,              intent(in)    :: mask_func(:,:)
      type(bc_mxll_t),    intent(in)    :: bc
      FLOAT,              intent(inout) :: io_func(:)
      integer,            intent(in)    :: idim

      integer :: ip, ip_in

      do ip_in = 1, bc%mask_points_number(idim)
        ip          = bc%mask_points_map(ip_in, idim)
        io_func(ip) = mask_func(ip_in, idim)
      end do

    end subroutine get_mask_io_function

    subroutine get_medium_io_function(medium_func, bc, io_func, idim)
      FLOAT,              intent(in)    :: medium_func(:)
      type(bc_mxll_t),    intent(in)    :: bc
      FLOAT,              intent(inout) :: io_func(:)
      integer,            intent(in)    :: idim

      integer :: ip, ip_in

      do ip_in = 1, bc%medium(idim)%points_number
        ip          = bc%medium(idim)%points_map(ip_in)
        io_func(ip) = medium_func(ip_in)
      end do

    end subroutine get_medium_io_function

    subroutine write_files(filename, tmp)
      character(len=*), intent(in) :: filename
      FLOAT,            intent(in) :: tmp(:)

      call dio_function_output(io_function_fill_how("VTK"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("AxisX"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("AxisY"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("AxisZ"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneX"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneY"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
      call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", trim(filename), namespace, space, mesh, tmp, &
        unit_one, err)
    end subroutine write_files

  end subroutine bc_mxll_write_info

  ! ---------------------------------------------------------
  subroutine maxwell_mask_points_mapping(bc, mesh, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, ip_in_max, point_info, idim
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_mask_points_mapping)

    call profiling_in(prof, 'MAXWELL_MASK_POINTS_MAPPING')
    ip_in_max = 1
    do idim = 1, 3
      if (bc%bc_ab_type(idim) == MXLL_AB_MASK) then
        ! allocate mask points map
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
          end if
        end do
        if (ip_in > ip_in_max) ip_in_max = ip_in
        bc%mask_points_number(idim) = ip_in
      end if
    end do
    SAFE_ALLOCATE(bc%mask(1:ip_in_max, idim))
    SAFE_ALLOCATE(bc%mask_points_map(1:ip_in_max, idim))

    do idim = 1, 3
      if (bc%bc_ab_type(idim) == MXLL_AB_MASK) then
        ! mask points mapping
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%mask_points_map(ip_in, idim) = ip
          end if
        end do
      end if
    end do

    call profiling_out(prof)
    POP_SUB(maxwell_mask_points_mapping)
  end subroutine maxwell_mask_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_pml_points_mapping(bc, mesh, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, point_info
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_pml_points_mapping)

    call profiling_in(prof, 'MAXWELL_PML_POINTS_MAPPING')

    ! allocate pml points map
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%pml%points_number = ip_in
    SAFE_ALLOCATE(bc%pml%points_map(1:ip_in))
    SAFE_ALLOCATE(bc%pml%points_map_inv(mesh%np))

    bc%pml%points_map_inv = 0
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%pml%points_map(ip_in) = ip
        bc%pml%points_map_inv(ip) = ip_in
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_pml_points_mapping)
  end subroutine maxwell_pml_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_constant_points_mapping(bc, mesh, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, point_info
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_constant_points_mapping)

    call profiling_in(prof, 'MAXWELL_CONSTANT_POINTS_MAP')

    ! allocate constant points map
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%constant_points_number = ip_in
    SAFE_ALLOCATE(bc%constant_points_map(1:ip_in))
    SAFE_ALLOCATE(bc%constant_rs_state(1:ip_in, 3))

    ! zero constant mapping
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%constant_points_map(ip_in) = ip
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_constant_points_mapping)
  end subroutine maxwell_constant_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_plane_waves_points_mapping(bc, mesh, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, point_info
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_plane_waves_points_mapping)

    call profiling_in(prof, 'MXLL_PLANE_WAVES_POINTS_MAP')

    ! allocate zero points map
    ip_in = 0
    do ip = 1, mesh%np
     call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%plane_wave%points_number = ip_in
    SAFE_ALLOCATE(bc%plane_wave%points_map(1:ip_in))

    ! zero points mapping
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%plane_wave%points_map(ip_in) = ip
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_plane_waves_points_mapping)
  end subroutine maxwell_plane_waves_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_zero_points_mapping(bc, mesh, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, ip_in_max, point_info, idim
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_zero_points_mapping)

    call profiling_in(prof, 'MXLL_ZERO_POINTS_MAPPING')

    ip_in_max = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_ZERO) then
        ! allocate zero points map
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
          end if
        end do
        if (ip_in > ip_in_max) ip_in_max = ip_in
      end if
    end do
    bc%zero_points_number = ip_in
    SAFE_ALLOCATE(bc%zero(1:ip_in_max,3))
    SAFE_ALLOCATE(bc%zero_points_map(1:ip_in_max, 3))

    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_ZERO) then
        ! zero points mapping
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%zero_points_map(ip_in, idim) = ip
          end if
        end do
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_zero_points_mapping)
  end subroutine maxwell_zero_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_medium_points_mapping(bc, mesh, st, bounds)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t), intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, ip_in_max, ip_bd, ip_bd_max, point_info, boundary_info, idim
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_medium_points_mapping)

    call profiling_in(prof, 'MXLL_MEDIUM_POINTS_MAPPING')

    ip_in_max = 0
    ip_bd_max = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_MEDIUM) then
        ! allocate pml points map
        ip_in = 0
        ip_bd = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          call maxwell_boundary_point_info(mesh, ip, bounds, boundary_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
          end if
          if ((boundary_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_bd = ip_bd + 1
          end if
        end do
        bc%medium(idim)%points_number = ip_in
        bc%medium(idim)%bdry_number = ip_bd
      end if
    end do
    do idim = 1, 3
      SAFE_ALLOCATE(bc%medium(idim)%aux_ep(1:ip_in, st%dim))
      SAFE_ALLOCATE(bc%medium(idim)%aux_mu(1:ip_in, st%dim))
      SAFE_ALLOCATE(bc%medium(idim)%points_map(1:ip_in_max))
      SAFE_ALLOCATE(bc%medium(idim)%bdry_map(1:ip_bd_max))
    end do

    ip_in = 0
    ip_bd = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == MXLL_BC_MEDIUM) then
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, point_info)
          call maxwell_boundary_point_info(mesh, ip, bounds, boundary_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%medium(idim)%points_map(ip_in) = ip
          end if
          if ((boundary_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_bd = ip_bd + 1
           bc%medium(idim)%bdry_map(ip_bd) = ip
          end if
        end do
      end if
    end do

    call profiling_out(prof)

    POP_SUB(maxwell_medium_points_mapping)
  end subroutine maxwell_medium_points_mapping

  ! ---------------------------------------------------------
  subroutine bc_mxll_generate_pml(pml, gr, bounds, dt)
    type(pml_t),        intent(inout) :: pml
    type(grid_t),       intent(in)    :: gr
    FLOAT,              intent(in)    :: bounds(:,:)
    FLOAT, optional,    intent(in)    :: dt

    integer :: ip, ip_in, idim
    FLOAT   :: width(3), ddv(3), ss_e, ss_m, ss_max, aa_e, aa_m, bb_e, bb_m, gg, hh, kk, ll_e, ll_m
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(bc_mxll_generate_pml)

    call profiling_in(prof, 'BC_MXLL_GENERATE_PML')

    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np, 1:gr%sb%dim))

    SAFE_ALLOCATE(pml%kappa(1:pml%points_number, 1:gr%sb%dim))
    SAFE_ALLOCATE(pml%sigma_e(1:pml%points_number, 1:gr%sb%dim))
    SAFE_ALLOCATE(pml%sigma_m(1:pml%points_number, 1:gr%sb%dim))
    SAFE_ALLOCATE(pml%a(1:pml%points_number, 1:gr%sb%dim))
    SAFE_ALLOCATE(pml%b(1:pml%points_number, 1:gr%sb%dim))
    SAFE_ALLOCATE(pml%c(1:pml%points_number, 1:3))
    SAFE_ALLOCATE(pml%mask(1:pml%points_number))
    SAFE_ALLOCATE(pml%conv_plus(1:pml%points_number, 1:3, 1:3))
    SAFE_ALLOCATE(pml%conv_minus(1:pml%points_number, 1:3, 1:3))
    SAFE_ALLOCATE(pml%conv_plus_old(1:pml%points_number, 1:3, 1:3))
    SAFE_ALLOCATE(pml%conv_minus_old(1:pml%points_number, 1:3, 1:3))
    SAFE_ALLOCATE(pml%aux_ep(1:pml%points_number, 1:3, 1:3))
    SAFE_ALLOCATE(pml%aux_mu(1:pml%points_number, 1:3, 1:3))

    pml%kappa                 = M_ONE
    pml%sigma_e               = M_ZERO
    pml%sigma_m               = M_ZERO
    pml%a                     = M_z0
    pml%b                     = M_z0
    pml%c                     = M_ZERO
    pml%mask                  = M_ONE
    pml%conv_plus             = M_z0
    pml%conv_minus            = M_z0
    pml%conv_plus_old         = M_z0
    pml%conv_minus_old        = M_z0

    width(1:3) = bounds(2, 1:3) - bounds(1, 1:3)

    ! PML variables for all boundary points
    do ip_in = 1, pml%points_number
      ip = pml%points_map(ip_in)
      ddv(1:3) = abs(gr%mesh%x(ip, 1:3)) - bounds(1, 1:3)
      do idim = 1, gr%sb%dim
        if (ddv(idim) >= M_ZERO) then
          gg     = (ddv(idim)/pml%width)**pml%power
          hh     = (M_ONE-ddv(idim)/pml%width)**pml%power
          kk     = M_ONE ! + (pml%kappa_max - M_ONE) * gg
          ss_max = -(pml%power + M_ONE)*P_c*P_ep*log(pml%refl_error)/(M_TWO * pml%width)
          ss_e   = gg * ss_max
          ss_m   = gg * ss_max ! * P_mu/P_ep
          ll_e   = ss_e*kk ! + kk**2*pml%alpha_max*hh
          ll_m   = ss_m*kk ! + kk**2*pml%alpha_max*hh
          bb_e   = exp(-(ss_e/(P_ep))*dt)
          bb_m   = exp(-(ss_m/(P_ep))*dt)
!          aa_e   = (ss_e/ll_e)*(bb_e - 1)
!          aa_m   = (ss_m/ll_m)*(bb_m - 1)
          aa_e   = (bb_e - 1)
          aa_m   = (bb_m - 1)
          if (ll_e == M_ZERO) aa_e = M_ZERO
          if (ll_m == M_ZERO) aa_m = M_ZERO
          pml%sigma_e(ip_in, idim) = ss_e
          pml%sigma_m(ip_in, idim) = ss_m
          pml%a(ip_in, idim)       = aa_e + M_zI * aa_m
          pml%b(ip_in, idim)       = bb_e + M_zI * bb_m
          pml%kappa(ip_in, idim)   = kk
          pml%mask(ip_in)          = pml%mask(ip_in) * (M_ONE - sin(ddv(idim)*M_PI/(M_TWO*(width(idim))))**2)
        else
          pml%kappa(ip_in, idim)   = M_ONE
          pml%sigma_e(ip_in, idim) = M_ZERO
          pml%sigma_m(ip_in, idiM) = M_ZERO
          pml%a(ip_in, idim)       = M_z0
          pml%b(ip_in, idim)       = M_z0
          pml%mask(ip_in)          = M_ONE
        end if
      end do
    end do

    ! PML auxiliary epsilon for all boundary points
    do idim = 1, gr%sb%dim
      tmp = P_ep
      do ip_in = 1, pml%points_number
        ip = pml%points_map(ip_in)
        tmp(ip) = P_ep / pml%kappa(ip_in, idim)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, pml%points_number
        ip = pml%points_map(ip_in)
        pml%aux_ep(ip_in, :, idim) = tmp_grad(ip, :)/(M_FOUR*P_ep*pml%kappa(ip_in, idim))
      end do
    end do

    ! PML auxiliary mu
    do idim = 1, gr%sb%dim
      tmp = P_mu
      do ip_in = 1, pml%points_number
        ip = pml%points_map(ip_in)
        tmp(ip) = P_mu / pml%kappa(ip_in, idim)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, pml%points_number
        ip = pml%points_map(ip_in)
        pml%aux_mu(ip_in, :, idim) = tmp_grad(ip, :)/(M_FOUR*P_mu*pml%kappa(ip_in, idim))
      end do
    end do

    ! PML auxiliary c for all boundary points
    do idim = 1, gr%sb%dim
      do ip_in = 1, pml%points_number
        pml%c(ip_in, idim) = P_c/pml%kappa(ip_in, idim)
      end do
    end do

    call profiling_out(prof)

    POP_SUB(bc_mxll_generate_pml)
  end subroutine bc_mxll_generate_pml

  ! ---------------------------------------------------------
  subroutine bc_mxll_generate_mask(bc, mesh, bounds)
    type(bc_mxll_t),    intent(inout) :: bc
    type(mesh_t),       intent(in)    :: mesh
    FLOAT,              intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, idim, ip_in_max
    FLOAT   :: ddv(3), tmp(3), width(3)
    FLOAT, allocatable :: mask(:)
    type(profile_t), save :: prof

    PUSH_SUB(bc_mxll_generate_mask)

    call profiling_in(prof, 'BC_MXLL_GENERATE_MASK')

    ip_in_max = maxval(bc%mask_points_number(:))

    SAFE_ALLOCATE(mask(mesh%np))

    mask(:) = M_ONE

    width(1:3) = bounds(2, 1:3) - bounds(1, 1:3)
    tmp(:)   = M_ZERO

    do ip = 1, mesh%np
      tmp = M_ONE
      mask(ip) = M_ONE
      ddv(1:3) = abs(mesh%x(ip, 1:3)) - bounds(1, 1:3)
      do idim = 1, mesh%sb%dim
        if(ddv(idim) >= M_ZERO ) then
          if (ddv(idim)  <=  width(idim)) then
            tmp(idim) = M_ONE - sin(ddv(idim) * M_PI / (M_TWO * (width(idim)) ))**2
          else
            tmp(idim) = M_ONE
          end if
        end if
        mask(ip) = mask(ip) * tmp(idim)
      end do
    end do

    do idim = 1, mesh%sb%dim
      do ip_in = 1, bc%mask_points_number(idim)
        ip = bc%mask_points_map(ip_in, idim)
        bc%mask(ip_in,idim) = mask(ip)
      end do
    end do

    SAFE_DEALLOCATE_A(mask)

    call profiling_out(prof)

    POP_SUB(bc_mxll_generate_mask)
  end subroutine bc_mxll_generate_mask

  ! ---------------------------------------------------------
  subroutine bc_mxll_generate_medium(bc, gr, bounds, ep_factor, mu_factor, sigma_e_factor, sigma_m_factor)
    type(bc_mxll_t),         intent(inout) :: bc
    type(grid_t),            intent(in)    :: gr
    FLOAT,                   intent(in)    :: bounds(:,:)
    FLOAT,                   intent(in)    :: ep_factor
    FLOAT,                   intent(in)    :: mu_factor
    FLOAT,                   intent(in)    :: sigma_e_factor
    FLOAT,                   intent(in)    :: sigma_m_factor
    
    integer :: ip, ipp, ip_in, ip_in_max, ip_bd, idim, point_info
    FLOAT   :: dd, dd_min, dd_max, xx(3), xxp(3)
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    type(profile_t), save :: prof

    PUSH_SUB(bc_mxll_generate_medium)

    call profiling_in(prof, 'BC_MXLL_GENERATE_MEDIUM')

    ip_in_max = max(bc%medium(1)%points_number, bc%medium(2)%points_number, bc%medium(3)%points_number)
    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    do idim = 1, 3
      SAFE_ALLOCATE(bc%medium(idim)%aux_ep(ip_in_max,gr%sb%dim))
      SAFE_ALLOCATE(bc%medium(idim)%aux_mu(ip_in_max,gr%sb%dim))
      SAFE_ALLOCATE(bc%medium(idim)%ep(ip_in_max))
      SAFE_ALLOCATE(bc%medium(idim)%mu(ip_in_max))
      SAFE_ALLOCATE(bc%medium(idim)%sigma_e(ip_in_max))
      SAFE_ALLOCATE(bc%medium(idim)%sigma_m(ip_in_max))
      SAFE_ALLOCATE(bc%medium(idim)%c(ip_in_max))
      SAFE_ALLOCATE(tmp(gr%mesh%np_part))
      SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%sb%dim))
      bc%medium(idim)%aux_ep = M_ZERO
      bc%medium(idim)%aux_mu = M_ZERO
      bc%medium(idim)%c = P_c

      tmp = P_ep
      do  ip = 1, gr%mesh%np_part
        call maxwell_box_point_info(bc, gr%mesh, ip, bounds, point_info)
        if ((point_info /= 0) .and. (abs(gr%mesh%x(ip, idim)) <= bounds(1, idim))) then
          xx(:) = gr%mesh%x(ip, :)
          dd_min = M_HUGE
          do ip_bd = 1, bc%medium(idim)%bdry_number
            ipp = bc%medium(idim)%bdry_map(ip_bd)
            xxp(:) = gr%mesh%x(ipp, :)
            dd = sqrt(sum((xx(1:3) - xxp(1:3))**2))
            if (dd < dd_min) dd_min = dd
          end do
          tmp(ip) = P_ep * (M_ONE + ep_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min-M_TWO*dd_max))))
        end if
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%medium(idim)%points_number
        ip = bc%medium(idim)%points_map(ip_in)
        bc%medium(idim)%aux_ep(ip_in, :) = &
          tmp_grad(ip, :)/(M_FOUR*P_ep*ep_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max-dd)))
      end do
    end do

    do idim = 1, 3
      tmp = P_mu
      do ip = 1, gr%mesh%np_part
        call maxwell_box_point_info(bc, gr%mesh, ip, bounds, point_info)
        if ((point_info == 1) .and. (abs(gr%mesh%x(ip, idim)) <= bounds(1, idim))) then
          xx(:) = gr%mesh%x(ip, :)
          dd_min = M_HUGE
          do ip_bd = 1, bc%medium(idim)%bdry_number
            ipp = bc%medium(idim)%bdry_map(ip_bd)
            xxp(:) = gr%mesh%x(ipp,:)
            dd = sqrt(sum((xx(1:3) - xxp(1:3))**2))
            if (dd < dd_min) dd_min = dd
          end do
          tmp(ip) = P_mu * (M_ONE + mu_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
        end if
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%medium(idim)%points_number
        ip = bc%medium(idim)%points_map(ip_in)
        bc%medium(idim)%aux_mu(ip_in, :) = &
          tmp_grad(ip, :)/(M_FOUR*P_mu*mu_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max-dd)))
      end do
    end do

    do idim = 1, 3
      do ip_in = 1, bc%medium(idim)%points_number
        ip = bc%medium(idim)%points_map(ip_in)
        xx(:) = gr%mesh%x(ip, :)
        dd_min = M_HUGE
        do ip_bd = 1, bc%medium(idim)%bdry_number
          ipp = bc%medium(idim)%bdry_map(ip_bd)
          xxp(:) = gr%mesh%x(ipp, :)
          dd = sqrt(sum((xx(1:3) - xxp(1:3))**2))
          if (dd < dd_min) dd_min = dd
        end do
        bc%medium(idim)%ep(ip_in) = P_ep * (M_ONE + ep_factor &
             * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium(idim)%mu(ip_in) = P_mu * (M_ONE + mu_factor &
             * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium(idim)%sigma_e(ip_in) = (M_ONE + sigma_e_factor &
             * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium(idim)%sigma_m(ip_in) = (M_ONE + sigma_m_factor &
             * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium(idim)%c(ip_in) = M_ONE/sqrt(bc%medium(idim)%ep(ip_in)*bc%medium(idim)%mu(ip_in))
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

    call profiling_out(prof)

    POP_SUB(bc_mxll_generate_medium)
  end subroutine bc_mxll_generate_medium

  ! ---------------------------------------------------------
  subroutine maxwell_plane_waves_boundaries_init(bc, namespace)
    type(bc_mxll_t),        intent(inout) :: bc
    type(namespace_t),      intent(in)    :: namespace

    type(block_t)        :: blk
    integer              :: il, nlines, ncols, ierr
    FLOAT                :: k_vector(3), vv(3), xx(3), rr, dummy(3), test, test_limit!, angle, sigma
    CMPLX                :: e_field(3)
    character(len=1024)  :: k_string(3)
    character(len=1024)  :: mxf_expression
    type(profile_t), save :: prof

    PUSH_SUB(maxwell_plane_waves_boundaries_init)

    call profiling_in(prof, 'MXLL_PLANE_WAVES_BOUND_INI')

    test_limit = CNST(10.0e-9)

    !%Variable MaxwellIncidentWaves
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% The initial electromagnetic fields can be set by the user 
    !% with the <tt>MaxwellIncidentWaves</tt> block variable.
    !% The electromagnetic fields have to fulfill the 
    !% Maxwells equations in vacuum.
    !%
    !% Example:
    !%
    !% <tt>%MaxwellIncidentWaves
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k1x" | "k1y" | "k1z" | "E1x" | "E1z" | "E1x"
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k2x" | "k2y" | "k2z" | "E2x" | "E2y" | "E2z"
    !% <br>&nbsp;&nbsp;   plane_wave_gauss       | "k3x" | "k3y" | "k3z" | "E3x" | "E3y" | "E3z" | "width" | "shift"
    !% <br>&nbsp;&nbsp;   plane_wave_mx_function | "E4x" | "E4y" | "E4z" | mx_envelope_name
    !% <br>%</tt>
    !%
    !% Description about MaxwellIncidentWaves follows
    !%
    !%Option plane_wave_parser 0
    !% Parser input modus
    !%Option plane_wave_mx_function 1
    !% The incident wave envelope is defined by an mx_function
    !%End

    if(parse_block(namespace, 'MaxwellIncidentWaves', blk) == 0) then

      call messages_print_stress(stdout, trim('Substitution of the electromagnetic incident waves'), namespace=namespace)

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      bc%plane_wave%number = nlines
      SAFE_ALLOCATE(bc%plane_wave%modus(nlines))
      SAFE_ALLOCATE(bc%plane_wave%e_field_string(3, nlines))
      SAFE_ALLOCATE(bc%plane_wave%e_field(3, nlines))
      SAFE_ALLOCATE(bc%plane_wave%k_vector(3, nlines))
      SAFE_ALLOCATE(bc%plane_wave%v_vector(3, nlines))
      SAFE_ALLOCATE(bc%plane_wave%mx_function(nlines))

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if ((ncols /= 5) .and. (ncols /= 7) .and. (ncols /= 9)) then
          message(1) = 'Each line in the MaxwellIncidentWaves block must have five, seven or nine columns.'
          call messages_fatal(1, namespace=namespace)
        end if

        ! check input modus e.g. parser of defined functions
        call parse_block_integer(blk, il - 1, 0, bc%plane_wave%modus(il))

        ! parse formula string
        if (bc%plane_wave%modus(il) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER) then

          call parse_block_string( blk, il - 1, 1, k_string(1))
          call parse_block_string( blk, il - 1, 2, k_string(2))
          call parse_block_string( blk, il - 1, 3, k_string(3))
          call parse_block_string( blk, il - 1, 4, bc%plane_wave%e_field_string(1, il))
          call parse_block_string( blk, il - 1, 5, bc%plane_wave%e_field_string(2, il))
          call parse_block_string( blk, il - 1, 6, bc%plane_wave%e_field_string(3, il))

          write(message(1), '(a,i2,a) ') 'Substituting electromagnetic incident wave ', il, ' with the expressions: '
          call messages_info(1)
          write(message(1), '(6a)')     '  Wave vector k(x)   = ', trim(k_string(1))
          write(message(2), '(2a)')     '  Wave vector k(y)   = ', trim(k_string(2))
          write(message(3), '(2a)')     '  Wave vector k(z)   = ', trim(k_string(3))
          write(message(4), '(2a)')     '  E-field(x) for t_0 = ', trim(bc%plane_wave%e_field_string(1, il))
          write(message(5), '(2a)')     '  E-field(y) for t_0 = ', trim(bc%plane_wave%e_field_string(2, il))
          write(message(6), '(2a)')     '  E-field(z) for t_0 = ', trim(bc%plane_wave%e_field_string(3, il))
          call messages_info(6)

          call conv_to_C_string(k_string(1))
          call conv_to_C_string(k_string(2))
          call conv_to_C_string(k_string(3))
          call conv_to_C_string(bc%plane_wave%e_field_string(1, il))
          call conv_to_C_string(bc%plane_wave%e_field_string(2, il))
          call conv_to_C_string(bc%plane_wave%e_field_string(3, il))

          xx(:) = M_ZERO
          rr    = M_ZERO
          call parse_expression(k_vector(1), dummy(1), 1, xx, rr, M_ZERO, k_string(1))
          call parse_expression(k_vector(2), dummy(2), 2, xx, rr, M_ZERO, k_string(2))
          call parse_expression(k_vector(3), dummy(3), 3, xx, rr, M_ZERO, k_string(3))
          k_vector(1) = units_to_atomic(unit_one/units_inp%length, k_vector(1))
          k_vector(2) = units_to_atomic(unit_one/units_inp%length, k_vector(2))
          k_vector(3) = units_to_atomic(unit_one/units_inp%length, k_vector(3))

          vv(:)    = k_vector(:) / sqrt(sum(k_vector(:)**2)) * P_c
          bc%plane_wave%k_vector(:,il) = k_vector(:)
          bc%plane_wave%v_vector(:,il) = vv(:)

        else if (bc%plane_wave%modus(il) == OPTION__MAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
          call parse_block_cmplx( blk, il - 1, 1, e_field(1))
          call parse_block_cmplx( blk, il - 1, 2, e_field(2))
          call parse_block_cmplx( blk, il - 1, 3, e_field(3))
          call parse_block_string( blk, il - 1, 4, mxf_expression)

          write(message(1), '(a,i2) ') 'Substituting electromagnetic incident wave ', il
          write(message(3), '(a)'    ) 'with the expression: '
          call messages_info(2)
          write(message(1), '(a,f9.4,sp,f9.4,"i")') '  E-field(x) complex amplitude  = ', real(e_field(1)), aimag(e_field(1))
          write(message(2), '(a,f9.4,sp,f9.4,"i")') '  E-field(y) complex amplitude  = ', real(e_field(2)), aimag(e_field(2))
          write(message(3), '(a,f9.4,sp,f9.4,"i")') '  E-field(z) complex amplitude  = ', real(e_field(3)), aimag(e_field(3))
          write(message(4), '(2a)'    )      '  Maxwell wave function name = ', trim(mxf_expression)
          call messages_info(4)
          call mxf_read(bc%plane_wave%mx_function(il), namespace, trim(mxf_expression), ierr)
          if (ierr /= 0) then            
            write(message(1),'(3A)') 'Error in the ""', trim(mxf_expression), &
              '"" field defined in the MaxwellIncidentWaves block'
            call messages_fatal(1, namespace=namespace)
          end if
          e_field  = units_to_atomic(units_inp%energy/units_inp%length, e_field)
          k_vector(1:3) = bc%plane_wave%mx_function(il)%k_vector(1:3)

          test = TOFLOAT(dot_product(k_vector(1:3), e_field(1:3)))
          if (abs(test) > test_limit) then
            message(1) = 'The wave vector k(:) or its electric field E-field(:) '
            message(2) = 'is not perpendicular enough.'
            call messages_fatal(2, namespace=namespace)
          end if
          if (Sqrt(sum(k_vector(:)**2)) < 1e-10) then
            message(1) = 'The k vector is not defined correctly.'
            call messages_fatal(1, namespace=namespace)
          end if

          bc%plane_wave%e_field(:,il)  = e_field(:)
          bc%plane_wave%k_vector(:,il) = k_vector(:)
          bc%plane_wave%v_vector(:,il) = k_vector(:) / Sqrt(sum(k_vector(:)**2)) * P_c

         end if
      end do

      call parse_block_end(blk)

      call messages_print_stress(stdout, namespace=namespace)

    end if

    call profiling_out(prof)

    POP_SUB(maxwell_plane_waves_boundaries_init)
  end subroutine maxwell_plane_waves_boundaries_init

  ! ---------------------------------------------------------
  subroutine maxwell_surfaces_init(mesh, st, bounds)
    type(mesh_t),             intent(in)    :: mesh
    type(states_mxll_t),      intent(inout) :: st
    FLOAT,                    intent(in)    :: bounds(:,:)

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_surfaces_init)

    call profiling_in(prof, 'MAXWELL_SURFACES_INIT')

    ! y-z surface at -x boundary
    st%surface(1, 1)%spacing   = M_HALF*(mesh%spacing(2) + mesh%spacing(3))
    st%surface(1, 1)%origin(:) = M_ZERO
    st%surface(1, 1)%origin(1) = -bounds(1, 1)
    st%surface(1, 1)%n(:) = M_ZERO
    st%surface(1, 1)%n(1) = -M_ONE
    st%surface(1, 1)%u(:) = M_ZERO
    st%surface(1, 1)%u(2) = -M_ONE
    st%surface(1, 1)%v(:) = M_ZERO
    st%surface(1, 1)%v(3) = M_ONE
    st%surface(1, 1)%nu   = -int(bounds(1, 2)/mesh%spacing(2))
    st%surface(1, 1)%mu   =  int(bounds(1, 2)/mesh%spacing(2))
    st%surface(1, 1)%nv   = -int(bounds(1, 3)/mesh%spacing(3))
    st%surface(1, 1)%mv   =  int(bounds(1, 3)/mesh%spacing(3))

    ! y-z surface at +x boundary
    st%surface(2, 1)%spacing   = M_HALF*(mesh%spacing(2) + mesh%spacing(3))
    st%surface(2, 1)%origin(:) = M_ZERO
    st%surface(2, 1)%origin(1) = bounds(1, 1)
    st%surface(2, 1)%n(:) = M_ZERO
    st%surface(2, 1)%n(1) = M_ONE
    st%surface(2, 1)%u(:) = M_ZERO
    st%surface(2, 1)%u(2) = M_ONE
    st%surface(2, 1)%v(:) = M_ZERO
    st%surface(2, 1)%v(3) = M_ONE
    st%surface(2, 1)%nu   = -int(bounds(1, 2)/mesh%spacing(2))
    st%surface(2, 1)%mu   =  int(bounds(1, 2)/mesh%spacing(2))
    st%surface(2, 1)%nv   = -int(bounds(1, 3)/mesh%spacing(3))
    st%surface(2, 1)%mv   =  int(bounds(1, 3)/mesh%spacing(3))

    ! x-z surface at -y boundary
    st%surface(1, 2)%spacing   = M_HALF*(mesh%spacing(1) + mesh%spacing(3))
    st%surface(1, 2)%origin(:) = M_ZERO
    st%surface(1, 2)%origin(2) = -bounds(1, 2)
    st%surface(1, 2)%n(:) = M_ZERO
    st%surface(1, 2)%n(2) = -M_ONE
    st%surface(1, 2)%u(:) = M_ZERO
    st%surface(1, 2)%u(1) = M_ONE
    st%surface(1, 2)%v(:) = M_ZERO
    st%surface(1, 2)%v(3) = M_ONE
    st%surface(1, 2)%nu   = -int(bounds(1, 1)/mesh%spacing(1))
    st%surface(1, 2)%mu   =  int(bounds(1, 1)/mesh%spacing(1))
    st%surface(1, 2)%nv   = -int(bounds(1, 3)/mesh%spacing(3))
    st%surface(1, 2)%mv   =  int(bounds(1, 3)/mesh%spacing(3))

    ! x-z surface at +y boundary
    st%surface(2, 2)%spacing   = M_HALF*(mesh%spacing(1) + mesh%spacing(3))
    st%surface(2, 2)%origin(:) = M_ZERO
    st%surface(2, 2)%origin(2) = bounds(1, 2)
    st%surface(2, 2)%n(:) = M_ZERO
    st%surface(2, 2)%n(2) = M_ONE
    st%surface(2, 2)%u(:) = M_ZERO
    st%surface(2, 2)%u(1) = M_ONE
    st%surface(2, 2)%v(:) = M_ZERO
    st%surface(2, 2)%v(3) = -M_ONE
    st%surface(2, 2)%nu   = -int(bounds(1, 1)/mesh%spacing(1))
    st%surface(2, 2)%mu   =  int(bounds(1, 1)/mesh%spacing(1))
    st%surface(2, 2)%nv   = -int(bounds(1, 3)/mesh%spacing(3))
    st%surface(2, 2)%mv   =  int(bounds(1, 3)/mesh%spacing(3))

    ! x-y surface at -z boundary
    st%surface(1, 3)%spacing   = M_HALF*(mesh%spacing(1) + mesh%spacing(2))
    st%surface(1, 3)%origin(:) = M_ZERO
    st%surface(1, 3)%origin(3) = -bounds(1, 3)
    st%surface(1, 3)%n(:) = M_ZERO
    st%surface(1, 3)%n(3) = -M_ONE
    st%surface(1, 3)%u(:) = M_ZERO
    st%surface(1, 3)%u(1) = M_ONE
    st%surface(1, 3)%v(:) = M_ZERO
    st%surface(1, 3)%v(2) = -M_ONE
    st%surface(1, 3)%nu   = -int(bounds(1, 1)/mesh%spacing(1))
    st%surface(1, 3)%mu   =  int(bounds(1, 1)/mesh%spacing(1))
    st%surface(1, 3)%nv   = -int(bounds(1, 2)/mesh%spacing(2))
    st%surface(1, 3)%mv   =  int(bounds(1, 2)/mesh%spacing(2))

    ! x-y surface at +z boundary
    st%surface(2, 3)%spacing   = M_HALF*(mesh%spacing(1) + mesh%spacing(2))
    st%surface(2, 3)%origin(:) = M_ZERO
    st%surface(2, 3)%origin(3) = bounds(1, 3)
    st%surface(2, 3)%n(:) = M_ZERO
    st%surface(2, 3)%n(3) = M_ONE
    st%surface(2, 3)%u(:) = M_ZERO
    st%surface(2, 3)%u(1) = M_ONE
    st%surface(2, 3)%v(:) = M_ZERO
    st%surface(2, 3)%v(2) = M_ONE
    st%surface(2, 3)%nu   = -int(bounds(1, 1)/mesh%spacing(1))
    st%surface(2, 3)%mu   =  int(bounds(1, 1)/mesh%spacing(1))
    st%surface(2, 3)%nv   = -int(bounds(1, 2)/mesh%spacing(2))
    st%surface(2, 3)%mv   =  int(bounds(1, 2)/mesh%spacing(2))

    call profiling_out(prof)

    POP_SUB(maxwell_surfaces_init)
  end subroutine maxwell_surfaces_init

  ! ---------------------------------------------------------
  subroutine maxwell_box_point_info(bc, mesh, ip, bounds, point_info) 
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(in)    :: ip
    FLOAT,               intent(in)    :: bounds(:,:)
    integer,             intent(out)   :: point_info

    FLOAT   :: rr, dd, xx(3), width(3)

    point_info = 0

    width(1:3) = bounds(2, 1:3) - bounds(1, 1:3)
    xx = M_ZERO
    xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)

    if (bc%ab_user_def) then

      dd = bc%ab_ufn(ip) - bounds(1, 1)
      if (dd > M_ZERO) then
        if (bc%ab_ufn(ip) < bounds(2, 1)) then
          point_info = 1
        end if
      end if

    else ! bc%ab_user_def == .false.

      select type (box => mesh%sb%box)
      type is (box_sphere_t)
        rr = norm2(xx -  box%center)
        dd = rr -  bounds(1, 1)
        if (dd > M_ZERO ) then
          if (dd  <  width(1)) then
            point_info = 1
          end if
        end if

      type is (box_parallelepiped_t)
        ! Limits of boundary region
        if (all(abs(xx(1:3) - box%center) <= bounds(2, 1:3))) then
          if (any(abs(xx(1:3) - box%center) > bounds(1, 1:3))) then
            point_info = 1
          else
            point_info = 0
          end if
        else
          point_info = -1
        end if

      class default
        ! Other box shapes are not supported.
        ASSERT(.false.)
      end select
    end if

  end subroutine maxwell_box_point_info

  ! ---------------------------------------------------------
  subroutine maxwell_boundary_point_info(mesh, ip, bounds, boundary_info)
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(in)    :: ip
    FLOAT,               intent(in)    :: bounds(:,:)
    integer,             intent(out)   :: boundary_info
 
    FLOAT   :: xx(3)

    boundary_info = 0

    xx = M_ZERO
    xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
    if (abs(xx(1)) == bounds(1, 1) .and. (all(abs(xx(2:3)) <= bounds(1, 2:3))) .or. &
        abs(xx(2)) == bounds(1, 2) .and. (all(abs(xx(1:3:2)) <= bounds(1, 1:3:2))) .or. &
        abs(xx(3)) == bounds(1, 3) .and. (all(abs(xx(1:2)) <= bounds(1, 1:2))) ) then
      boundary_info = 1
    else
      boundary_info = 0
    end if

  end subroutine maxwell_boundary_point_info

  ! ---------------------------------------------------------
  subroutine inner_and_outer_points_mapping(mesh, st, bounds)
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t), intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ip, ip_in, ip_bd, point_info
    FLOAT   :: xx(mesh%sb%dim)
    type(profile_t), save :: prof

    PUSH_SUB(inner_and_outer_points_mapping)

    call profiling_in(prof, 'INNER_AND_OUTER_POINTS_MAP')

    ! allocate inner and boundary points points map
    ip_in = 0
    ip_bd = 0
    do ip = 1, mesh%np
      xx = mesh%x(ip, :)
      if (all(abs(xx) <= bounds(2,1:mesh%sb%dim))) then
        if (any(abs(xx) > bounds(1,1:mesh%sb%dim))) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 0) then
        ip_in = ip_in + 1
      else
        ip_bd = ip_bd + 1
      end if
    end do
    st%inner_points_number = ip_in
    SAFE_ALLOCATE(st%inner_points_map(1:ip_in))
    SAFE_ALLOCATE(st%inner_points_mask(1:mesh%np))
    st%boundary_points_number = ip_bd
    SAFE_ALLOCATE(st%boundary_points_map(1:ip_bd))
    SAFE_ALLOCATE(st%boundary_points_mask(1:mesh%np))
    st%inner_points_mask = .false.
    st%boundary_points_mask = .false.

    ! inner and boundary points mapping
    ip_in = 0
    ip_bd = 0
    do ip = 1, mesh%np
      xx = mesh%x(ip, :)
      if (all(abs(xx) <= bounds(2,1:mesh%sb%dim))) then
        if (any(abs(xx) > bounds(1,1:mesh%sb%dim))) then
          point_info = 1
        else
          point_info = 0
        end if
      else
        point_info = -1
      end if
      if (point_info == 0) then
        ip_in = ip_in + 1
        st%inner_points_map(ip_in) = ip
        st%inner_points_mask(ip) = .true.
      else
        ip_bd = ip_bd + 1
        st%boundary_points_map(ip_bd) = ip
        st%boundary_points_mask(ip) = .true.
      end if
    end do

    call profiling_out(prof)

    POP_SUB(inner_and_outer_points_mapping)
  end subroutine inner_and_outer_points_mapping

  ! ---------------------------------------------------------
  subroutine surface_grid_points_mapping(mesh, st, bounds)
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t), intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)

    integer :: ix, ix_max, iix, iy, iy_max, iiy, iz, iz_max, iiz, idx1, idx2, ip_global, nn_max
    integer, allocatable :: nn(:,:,:,:)
    FLOAT   :: rr(3), delta(3), vec(2), min_1(3), max_1(3), min_2(3), max_2(3)
    type(profile_t), save :: prof

    PUSH_SUB(surface_grid_points_mapping)

    call profiling_in(prof, 'SURFACE_GRID_POINTS_MAPPING')

    st%surface_grid_rows_number(1) = 3
    ix_max  = st%surface_grid_rows_number(1)
    st%surface_grid_rows_number(2) = 3
    iy_max  = st%surface_grid_rows_number(2)
    st%surface_grid_rows_number(3) = 3
    iz_max  = st%surface_grid_rows_number(3)

    delta(1) = M_TWO * abs(bounds(1,1)) / float(ix_max)
    delta(2) = M_TWO * abs(bounds(1,2)) / float(iy_max)
    delta(3) = M_TWO * abs(bounds(1,3)) / float(iz_max)

    st%surface_grid_element(1) = delta(2) * delta(3)
    st%surface_grid_element(2) = delta(1) * delta(3)
    st%surface_grid_element(3) = delta(1) * delta(2)

    SAFE_ALLOCATE(nn(1:2, 1:3, 1:3, 1:3))

    st%surface_grid_center(1, 1, :, :) = -int(bounds(1,1))
    do iy = 1, iy_max
      do iz = 1, iz_max
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(1, 2, iy, iz) = int(rr(2))
        st%surface_grid_center(1, 3, iy, iz) = int(rr(3))
      end do
    end do
    st%surface_grid_center(2, 1, :, :) = int(bounds(1,1))
    do iy = 1, iy_max
      do iz = 1, iz_max
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(2, 2, iy, iz) = int(rr(2))
        st%surface_grid_center(2, 3, iy, iz) = int(rr(3))
      end do
    end do

    st%surface_grid_center(1, 2, :, :) = -int(bounds(1,2))
    do ix = 1, ix_max
      do iz = 1, iz_max
        rr(1) = -bounds(1,1) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(1, 1, ix, iz) = int(rr(1))
        st%surface_grid_center(1, 3, ix, iz) = int(rr(3))
      end do
    end do
    st%surface_grid_center(2, 2, :, :) = int(bounds(1,2))
    do ix = 1, ix_max
      do iz = 1, iz_max
        rr(1) = -bounds(1,2) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(3) = -bounds(1,3) + delta(3)/M_TWO + (iz-1) * delta(3)
        st%surface_grid_center(2, 1, ix, iz) = int(rr(1))
        st%surface_grid_center(2, 3, ix, iz) = int(rr(3))
      end do
    end do

    st%surface_grid_center(1, 3, :, :) = -int(bounds(1,3))
    do ix = 1, ix_max
      do iy = 1, iy_max
        rr(1) = -bounds(1,1) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        st%surface_grid_center(1, 1, ix, iy) = int(rr(1))
        st%surface_grid_center(1, 2, ix, iy) = int(rr(2))
      end do
    end do
    st%surface_grid_center(2, 3, :, :) = int(bounds(1,3))
    do ix = 1, ix_max
      do iy = 1, iy_max
        rr(1) = -bounds(1,2) + delta(1)/M_TWO + (ix-1) * delta(1)
        rr(2) = -bounds(1,2) + delta(2)/M_TWO + (iy-1) * delta(2)
        st%surface_grid_center(2, 1, ix, iy) = int(rr(1))
        st%surface_grid_center(2, 2, ix, iy) = int(rr(2))
      end do
    end do

    st%surface_grid_points_number(:,:,:) = 0

    nn_max = 0

    do iy = 1, iy_max
      do iz = 1, iz_max
        min_1(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_1(iy) = -bounds(1,2) + iy * delta(2)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iiy * mesh%spacing(2)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          st%surface_grid_points_number(1, idx1, idx2) = st%surface_grid_points_number(1, idx1, idx2) + 1
          if (nn_max < st%surface_grid_points_number(1, idx1, idx2)) then
            nn_max = st%surface_grid_points_number(1, idx1, idx2)
          end if
        end if
      end do
    end do

    do ix = 1, ix_max
      do iz = 1, iz_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iix = mesh%idx%nr(1, 1), mesh%idx%nr(2, 1)
      do iiz = mesh%idx%nr(1, 3), mesh%idx%nr(2, 3)
        vec(1)     = iix * mesh%spacing(1)
        vec(2)     = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          st%surface_grid_points_number(2,idx1,idx2) = st%surface_grid_points_number(2,idx1,idx2)+1
          if (nn_max < st%surface_grid_points_number(2, idx1, idx2)) then
            nn_max = st%surface_grid_points_number(2, idx1, idx2)
          end if
        end if
      end do
    end do

    do ix = 1, ix_max
      do iy = 1, iy_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_2(iy) = -bounds(1,2) + iy * delta(2)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        vec(1) = iix * mesh%spacing(1)
        vec(2) = iiy * mesh%spacing(2)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          st%surface_grid_points_number(3, idx1, idx2) = st%surface_grid_points_number(3, idx1, idx2) + 1
          if (nn_max < st%surface_grid_points_number(3, idx1, idx2)) then
            nn_max = st%surface_grid_points_number(3, idx1, idx2)
          end if
        end if
      end do
    end do

    ! originally there were three allocated of the same pointer here
    SAFE_ALLOCATE(st%surface_grid_points_map(1:2, 1:st%dim, 1:ix_max, 1:iy_max, 1:nn_max))

    nn(:,:,:,:) = 0

    do iy = 1, iy_max
      do iz = 1, iz_max
        min_1(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_1(iy) = -bounds(1,2) + iy * delta(2)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iiy * mesh%spacing(2)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          nn(1, 1, idx1, idx2) = nn(1, 1, idx1, idx2) + 1
          rr(1) = -bounds(1, 1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = iiz * mesh%spacing(3)
          iix = int(-bounds(1,1)/mesh%spacing(1))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(1, 1, idx1, idx2, nn(1, 1, idx1, idx2)) = ip_global
          nn(2, 1, idx1, idx2) = nn(2, 1, idx1, idx2) + 1
          rr(1) = bounds(1,1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = iiz * mesh%spacing(3)
          iix = int(bounds(1,1)/mesh%spacing(1))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(2, 1, idx1, idx2, nn(2, 1, idx1, idx2)) = ip_global
        end if
      end do
    end do

    do ix = 1, ix_max
      do iz = 1, iz_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iz) = -bounds(1,3) + (iz-1) * delta(3)
        max_2(iz) = -bounds(1,3) + iz * delta(3)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        vec(1) = iix * mesh%spacing(1)
        vec(2) = iiz * mesh%spacing(3)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          nn(1, 2, idx1, idx2) = nn(1, 2, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = -bounds(1, 2)
          rr(3) = iiz * mesh%spacing(3)
          iiy = int(-bounds(1,2)/mesh%spacing(2))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(1, 2, idx1, idx2, nn(1, 2, idx1, idx2)) = ip_global
          nn(2, 2, idx1, idx2) = nn(2, 2, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = bounds(1,2)
          rr(3) = iiz * mesh%spacing(3)
          iiy = int(bounds(1,2)/mesh%spacing(2))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(2, 2, idx1, idx2, nn(2, 2, idx1, idx2)) = ip_global
        end if
      end do
    end do

    do ix = 1, ix_max
      do iy = 1, iy_max
        min_1(ix) = -bounds(1,1) + (ix-1) * delta(1)
        max_1(ix) = -bounds(1,1) + ix * delta(1)
        min_2(iy) = -bounds(1,2) + (iy-1) * delta(2)
        max_2(iy) = -bounds(1,2) + iy * delta(2)
      end do
    end do
    do iix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
      do iiy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        vec(1) = iix * mesh%spacing(1)
        vec(2) = iiy * mesh%spacing(2)
        call get_surface_indices(vec, min_1, max_1, min_2, max_2, idx1, idx2)
        if (idx1 /= 0 .and. idx2 /= 0) then
          nn(1, 3, idx1, idx2) = nn(1, 3, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = -bounds(1,3)
          iiz = int(-bounds(1,3)/mesh%spacing(3))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(1, 3, idx1, idx2, nn(1, 3, idx1, idx2)) = ip_global
          nn(2, 3, idx1, idx2) = nn(2, 3, idx1, idx2) + 1
          rr(1) = iix * mesh%spacing(1)
          rr(2) = iiy * mesh%spacing(2)
          rr(3) = bounds(1,3)
          iiz = int(bounds(1,3)/mesh%spacing(3))
          ip_global = mesh_global_index_from_coords(mesh, [iix, iiy, iiz])
          st%surface_grid_points_map(2, 3, idx1, idx2, nn(2, 3, idx1, idx2)) = ip_global
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(nn)

    call profiling_out(prof)

    POP_SUB(surface_grid_points_mapping)
  contains

    subroutine get_surface_indices(vec, min_1, max_1, min_2, max_2, index_1, index_2)
      FLOAT,   intent(in)  :: vec(:)
      FLOAT,   intent(in)  :: min_1(:)
      FLOAT,   intent(in)  :: max_1(:)
      FLOAT,   intent(in)  :: min_2(:)
      FLOAT,   intent(in)  :: max_2(:)
      integer, intent(out) :: index_1
      integer, intent(out) :: index_2

      if (vec(1) >= min_1(1) .and. vec(1) <= max_1(1) .and. vec(2) >= min_2(1) .and. vec(2) <= max_2(1)) then
        index_1 = 1
        index_2 = 1
      else if (vec(1) >= min_1(2) .and. vec(1) <= max_1(2) .and. vec(2) >= min_2(1) .and. vec(2) <= max_2(1)) then
        index_1 = 2
        index_2 = 1
      else if (vec(1) >= min_1(3) .and. vec(1) <= max_1(3) .and. vec(2) >= min_2(1) .and. vec(2) <= max_2(1)) then
        index_1 = 3
        index_2 = 1
      else if (vec(1) >= min_1(1) .and. vec(1) <= max_1(1) .and. vec(2) >= min_2(2) .and. vec(2) <= max_2(2)) then
        index_1 = 1
        index_2 = 2
      else if (vec(1) >= min_1(2) .and. vec(1) <= max_1(2) .and. vec(2) >= min_2(2) .and. vec(2) <= max_2(2)) then
        index_1 = 2
        index_2 = 2
      else if (vec(1) >= min_1(3) .and. vec(1) <= max_1(3) .and. vec(2) >= min_2(2) .and. vec(2) <= max_2(2)) then
        index_1 = 3
        index_2 = 2
      else if (vec(1) >= min_1(1) .and. vec(1) <= max_1(1) .and. vec(2) >= min_2(3) .and. vec(2) <= max_2(3)) then
        index_1 = 1
        index_2 = 3
      else if (vec(1) >= min_1(2) .and. vec(1) <= max_1(2) .and. vec(2) >= min_2(3) .and. vec(2) <= max_2(3)) then
        index_1 = 2
        index_2 = 3
      else if (vec(1) >= min_1(3) .and. vec(1) <= max_1(3) .and. vec(2) >= min_2(3) .and. vec(2) <= max_2(3)) then
        index_1 = 3
        index_2 = 3
      else
        index_1 = 0
        index_2 = 0
      end if

    end subroutine get_surface_indices

  end subroutine surface_grid_points_mapping

end module maxwell_boundary_op_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
