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
  use derivatives_oct_m
  use io_oct_m
  use io_function_oct_m
  use io_oct_m
  use cube_function_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
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
  use unit_oct_m
  use unit_system_oct_m
  use simul_box_oct_m
  use varinfo_oct_m
  use states_mxll_oct_m
  
  implicit none

  private
  public ::                    &
    bc_mxll_init,              &
    bc_mxll_end,               &   
    bc_mxll_write_info,        &
    bc_mxll_t

  type bc_mxll_t
    integer           :: bc_type(MAX_DIM)
    integer           :: bc_ab_type(MAX_DIM)
    FLOAT             :: bc_bounds(2,MAX_DIM)
    logical           :: ab_user_def
    FLOAT,   pointer  :: ab_ufn(:)

    FLOAT             :: mask_width
    FLOAT             :: mask_alpha
    integer           :: mask_points_number(MAX_DIM)
    integer, allocatable  :: mask_points_map(:,:)
    FLOAT,   allocatable  :: mask(:,:)

    integer           :: der_bndry_mask_points_number
    integer, pointer  :: der_bndry_mask_points_map(:)
    FLOAT,   pointer  :: der_bndry_mask(:)

    FLOAT             :: pml_width
    integer           :: pml_points_number
    integer, pointer  :: pml_points_map(:)
    integer, pointer  :: pml_points_map_inv(:)
    FLOAT             :: pml_kappa_max
    FLOAT             :: pml_alpha_max
    FLOAT             :: pml_power
    FLOAT             :: pml_refl_error
    FLOAT,   pointer  :: pml_kappa(:,:)
    FLOAT,   pointer  :: pml_sigma_e(:,:)
    FLOAT,   pointer  :: pml_sigma_m(:,:)
    CMPLX,   pointer  :: pml_a(:,:)
    CMPLX,   pointer  :: pml_b(:,:)
    FLOAT,   pointer  :: pml_c(:,:)
    FLOAT,   pointer  :: pml_mask(:)
    CMPLX,   pointer  :: pml_aux_ep(:,:,:)
    CMPLX,   pointer  :: pml_aux_mu(:,:,:)
    CMPLX,   pointer  :: pml_conv_plus(:,:,:)
    CMPLX,   pointer  :: pml_conv_minus(:,:,:)
    CMPLX,   pointer  :: pml_conv_plus_old(:,:,:)
    CMPLX,   pointer  :: pml_conv_minus_old(:,:,:)

    integer           :: constant_points_number
    integer, pointer  :: constant_points_map(:)
    CMPLX,   pointer  :: constant_rs_state(:,:)

    integer           :: mirror_points_number(3)
    integer, pointer  :: mirror_points_map(:,:)

    logical           :: do_plane_waves
    integer           :: plane_waves_points_number
    integer,  pointer :: plane_waves_points_map(:)
    integer           :: plane_waves_number
    integer,  pointer :: plane_waves_modus(:)
    integer,  pointer :: plane_waves(:)
    character(len=1024), pointer :: plane_waves_e_field_string(:,:)
    FLOAT,    pointer :: plane_waves_k_vector(:,:)
    FLOAT,    pointer :: plane_waves_v_vector(:,:)
    FLOAT,    pointer :: plane_waves_e_field(:,:)
    type(mxf_t), pointer :: plane_waves_mx_function(:)
    type(mxf_t), pointer :: plane_waves_mx_phase(:)
    integer,  pointer :: plane_waves_oam(:)
    integer,  pointer :: plane_waves_sam(:)

    logical           :: bessel_beam = .false.

    FLOAT             :: medium_width
    FLOAT             :: medium_ep_factor
    FLOAT             :: medium_mu_factor
    FLOAT             :: medium_sigma_e_factor
    FLOAT             :: medium_sigma_m_factor
    FLOAT,   pointer  :: medium_ep(:,:)
    FLOAT,   pointer  :: medium_mu(:,:)
    FLOAT,   pointer  :: medium_sigma_e(:,:)
    FLOAT,   pointer  :: medium_sigma_m(:,:)
    FLOAT,   pointer  :: medium_c(:,:)
    integer           :: medium_points_number(MAX_DIM)
    integer, pointer  :: medium_points_map(:,:)
    FLOAT,   pointer  :: medium_aux_ep(:,:,:)
    FLOAT,   pointer  :: medium_aux_mu(:,:,:)
    integer           :: medium_bdry_number(MAX_DIM)
    integer, pointer  :: medium_bdry_map(:,:)

    FLOAT             :: zero_width
    integer           :: zero_points_number(MAX_DIM)
    integer, pointer  :: zero_points_map(:,:)
    FLOAT,   pointer  :: zero(:,:)

!    logical           :: maxwell_vacuum_fluctuation = .false.
  end type bc_mxll_t

  integer, public, parameter ::       &
    AB_NOT_ABSORBING        = 0,         &
    AB_MASK                 = 1,         &
    AB_MAXWELL_MASK         = 2,         &
    AB_UPML                 = 3,         &
    AB_CPML                 = 4,         &
    AB_MASK_ZERO            = 7

contains

  ! ---------------------------------------------------------
  subroutine bc_mxll_init(bc, namespace, gr, st, sb, geo, dt)
    type(bc_mxll_t),          intent(inout) :: bc
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(states_mxll_t),      intent(inout) :: st
    type(simul_box_t),        intent(in)    :: sb
    type(geometry_t),         intent(in)    :: geo
    FLOAT, optional,          intent(in)    :: dt

    integer             :: idim, ab_shape_dim, nlines, icol, ncols 
    FLOAT               :: bounds(1:2,MAX_DIM), ab_bounds(1:2,MAX_DIM)
    FLOAT               :: mask_width, pml_width, zero_width
    type(block_t)       :: blk
    character(len=1024) :: string
    character(len=50)   :: str
    logical             :: plane_waves_check = .false., ab_mask_check = .false., ab_pml_check = .false.
    logical             :: constant_check = .false., zero_check = .false.

    PUSH_SUB(bc_mxll_init)

    bc%ab_user_def = .false.

    !%Variable MaxwellBoundaryConditions
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellIncidentWaves
    !% <br>&nbsp;&nbsp;   maxwell_zero | maxwell_mirror_pec | maxwell_consant 
    !% <br>%</tt>
    !%
    !% Description follows
    !%
    !%Option maxwell_zero 0
    !% follows ...
    !%Option maxwell_constant 2
    !% follows ...
    !%Option maxwell_mirror_pec 3
    !% follows ...
    !%Option maxwell_mirror_pmc 4
    !% follows ...
    !%Option maxwell_plane_waves 5
    !% follows ...
    !%Option maxwell_periodic 6
    !% follows ...
    !%Option maxwell_medium 7
    !% follows ...
    !%Option maxwell_lossy_layer 8
    !% follows ...
    !%End
    if(parse_block(namespace, 'MaxwellBoundaryConditions', blk) == 0) then ! memleak in parse_block

      call messages_print_stress(stdout, trim('Maxwell boundary conditions:'), namespace=namespace)

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      if (nlines /= 1) then
        message(1) = 'MaxwellBoundaryConditions has to consist of one line!'
        call messages_fatal(1, namespace=namespace)
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 3) then
        message(1) = 'MaxwellBoundaryConditions has to consist of three columns!'
        call messages_fatal(1, namespace=namespace)
      end if
      do icol = 1, ncols
        call parse_block_integer(blk, 0, icol-1, bc%bc_type(icol))
        select case (bc%bc_type(icol))
          case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO)
          string = 'Zero'
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_CONSTANT)
          string = 'Constant'
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PEC)
          string = 'PEC Mirror'
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PMC)
          string = 'PMC Mirror'
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PERIODIC)
          string = 'Periodic'
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PLANE_WAVES)
          string = 'Plane waves'
          if (.not. (parse_is_defined(namespace, 'UserDefinedMaxwellIncidentWaves')) ) then
            write(message(1),'(a)') 'Input: Maxwell boundary condition option is set to "maxwell_plane_waves".'
            write(message(2),'(a)') 'Input: User defined Maxwell plane waves have to be defined!'
            call messages_fatal(2, namespace=namespace)
          end if
        case (OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM)
          string = 'Medium boundary'
        end select
        write(message(1),'(a,I1,a,a)') 'Maxwell boundary condition in direction ', icol, ': ', trim(string)
        call messages_info(1)
      end do

      call messages_print_stress(stdout, namespace=namespace)
    end if
    
    !%Variable MaxwellAbsorbingBoundaries
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Follows
    !%
    !% Example:
    !%
    !% <tt>%MaxwellAbsorbingBoundaries
    !% <br>&nbsp;&nbsp;   cpml | cpml | cpml
    !% <br>%</tt>
    !%
    !% Description follows
    !%
    !%Option not_absorbing 0
    !% No absorbing boundaries.
    !%Option mask 1
    !% A mask the same as for the matter wavefunctions is applied to the Maxwell states at the boundaries.
    !%Option maxwell_mask 2
    !% A different mask than for the matter to apply on Maxwell states
    !%Option upml 3
    !% UPML
    !%Option cpml 4
    !% CPML
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
    end if
    do icol=1, ncols
      call parse_block_integer(blk, 0, icol-1, bc%bc_ab_type(icol))
    end do


    do idim = 1, 3
      if (bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__MASK) then
        ab_mask_check = .true.
      end if
      if (bc%bc_ab_type(idim) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then
        ab_pml_check = .true.
      end if
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_CONSTANT) then
        constant_check = .true.
      end if
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO) then
        zero_check = .true.
      end if
    end do

    if (ab_mask_check .or. ab_pml_check) then
      write(str, '(a,i5)') 'Maxwell Absorbing Boundaries'
      call messages_print_stress(stdout, trim(str), namespace=namespace)
    end if

    do idim = 1, st%d%dim
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                     ) * gr%mesh%spacing(idim)
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_CONSTANT) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - 2 * gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                         ) * gr%mesh%spacing(idim)
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PEC) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                     ) * gr%mesh%spacing(idim)
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MIRROR_PMC) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                     ) * gr%mesh%spacing(idim)
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PERIODIC) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - 2 * gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                         ) * gr%mesh%spacing(idim)
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_PLANE_WAVES) then
        bounds(1, idim) = ( gr%mesh%idx%nr(2, idim) - 2 * gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(2, idim) = ( gr%mesh%idx%nr(2, idim)                                         ) * gr%mesh%spacing(idim)
        plane_waves_check = .true.
        bc%do_plane_waves = .true.
      else if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) then
        !%Variable MaxwellMediumWidth
        !%Type float
        !%Default 0.0 a.u.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% Width of the boundary region with medium
        !%End
        call parse_variable(namespace, 'MaxwellMediumWidth', M_ZERO, bc%medium_width, units_inp%length)
        bounds(1,idim) = ( gr%mesh%idx%nr(2, idim) - gr%mesh%idx%enlarge(idim) ) * gr%mesh%spacing(idim)
        bounds(1,idim) = bounds(1,idim) - bc%medium_width
        bounds(2,idim) = ( gr%mesh%idx%nr(2, idim) ) * gr%mesh%spacing(idim)
        !%Variable MaxwellEpsilonFactor
        !%Type float
        !%Default 1.0.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% Maxwell epsilon factor
        !%End
        call parse_variable(namespace, 'MaxwellEpsilonFactor', M_ONE, bc%medium_ep_factor, unit_one)
        !%Variable MaxwellMuFactor
        !%Type float
        !%Default 1.0.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% Maxwell mu factor
        !%End
        call parse_variable(namespace, 'MaxwellMuFactor', M_ZERO, bc%medium_mu_factor, unit_one)
        !%Variable MaxwellElectricSigma
        !%Type float
        !%Default 1.0.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% Maxwell electric sigma
        !%End
        call parse_variable(namespace, 'MaxwellElectricSigma', M_ZERO, bc%medium_sigma_e_factor, unit_one)
        !%Variable MaxwellMagneticSigma
        !%Type float
        !%Default 1.0.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% Maxwell magnetic sigma
        !%End
        call parse_variable(namespace, 'MaxwellMagneticSigma', M_ONE, bc%medium_sigma_m_factor, unit_one)
        call maxwell_medium_points_mapping(bc, gr%mesh, st, bounds, geo)
        call bc_mxll_generate_medium(bc, gr, bounds, geo)
      end if

      if(bc%bc_ab_type(idim) /= AB_NOT_ABSORBING) then

        call messages_print_var_option(stdout, "MaxwellAbsorbingBoundaries", bc%bc_ab_type(idim))

        if (bc%bc_ab_type(idim) == AB_MASK_ZERO) then
          !%Variable MaxwellABZeroWidth
          !%Type float
          !%Default 0.4 a.u.
          !%Section Time-Dependent::Absorbing Boundaries
          !%Description
          !% Width of the region used to apply the absorbing boundaries.
          !%End
          zero_width = ab_bounds(2,idim)-ab_bounds(1,idim)
          call parse_variable(namespace, 'MaxwellABZeroWidth', zero_width, zero_width, units_inp%length)
          bc%zero_width = zero_width

          if ( zero_width < (gr%der%order*gr%mesh%spacing(1)) .or. &
               zero_width < (gr%der%order*gr%mesh%spacing(2)) .or. &
               zero_width < (gr%der%order*gr%mesh%spacing(3)) ) then
            zero_width = max(gr%der%order*gr%mesh%spacing(1), &
                                gr%der%order*gr%mesh%spacing(2), &
                                gr%der%order*gr%mesh%spacing(2))
            write(message(1),'(a)') 'Zero absorbing width has to be larger or equal than derivatives order times spacing!'
            write(message(2),'(a,es10.3)') 'Zero absorbing width is set to: ', zero_width
            call messages_info(2)
          end if
        end if

        if (bc%bc_ab_type(idim) == AB_MASK) then

          if (gr%mesh%sb%box_shape == SPHERE) then
            ab_shape_dim = 1
          else if(gr%mesh%sb%box_shape == PARALLELEPIPED) then
            ab_shape_dim = sb%dim
            ab_bounds(1, idim) = bounds(1, idim)
            ab_bounds(2, idim) = bounds(1, idim)
          end if

          !%Variable MaxwellABMaskWidth
          !%Type float
          !%Default 0.4 a.u.
          !%Section Time-Dependent::Absorbing Boundaries
          !%Description
          !% Width of the region used to apply the absorbing boundaries.
          !%End
          mask_width = ab_bounds(2, idim) - ab_bounds(1, idim)
          call parse_variable(namespace, 'MaxwellABMaskWidth', mask_width, mask_width, units_inp%length)
          bc%mask_width = mask_width
          if ( mask_width < (gr%der%order*gr%mesh%spacing(1)) .or. &
               mask_width < (gr%der%order*gr%mesh%spacing(2)) .or. &
               mask_width < (gr%der%order*gr%mesh%spacing(3)) ) then
            mask_width = max(gr%der%order*gr%mesh%spacing(1), &
                                gr%der%order*gr%mesh%spacing(2), &
                                gr%der%order*gr%mesh%spacing(2))
            write(message(1),'(a)') 'Mask absorbing width has to be larger or equal than derivatives order times spacing!'
            write(message(2),'(a,es10.3)') 'Mask absorbing width is set to: ', mask_width
            call messages_info(2)
          end if
          !%Variable MaxwellABMaskAlpha
          !%Type float
          !%Default 0.4 a.u.
          !%Section Time-Dependent::Absorbing Boundaries
          !%Description
          !% Width of the region used to apply the absorbing boundaries.
          !%End
          call parse_variable(namespace, 'MaxwellABMaskWidth', M_TWO, bc%mask_alpha, units_inp%length)
          ab_bounds(1, idim) = ab_bounds(2, idim) - mask_width
        end if

        if (bc%bc_ab_type(idim) == AB_CPML) then

          if (gr%mesh%sb%box_shape == SPHERE) then
            ab_shape_dim = 1
          else if(gr%mesh%sb%box_shape == PARALLELEPIPED) then
            ab_shape_dim = sb%dim
            ab_bounds(1, idim) = bounds(1, idim)
            ab_bounds(2, idim) = bounds(1, idim)
          end if

          !%Variable MaxwellABPMLWidth
          !%Type float
          !%Default 0.4 a.u.
          !%Section Time-Dependent::Absorbing Boundaries
          !%Description
          !% Width of the region used to apply the absorbing boundaries.
          !%End
          pml_width = ab_bounds(2,idim) - ab_bounds(1,idim)
          call parse_variable(namespace, 'MaxwellABPMLWidth', pml_width, pml_width, units_inp%length)
          bc%pml_width = pml_width
          if (parse_is_defined(namespace, 'UserDefinedMaxwellIncidentWaves')) then
            if ( pml_width < (gr%der%order*gr%mesh%spacing(1)) .or. &
                 pml_width < (gr%der%order*gr%mesh%spacing(2)) .or. &
                 pml_width < (gr%der%order*gr%mesh%spacing(3)) ) then
              pml_width = max(gr%der%order*gr%mesh%spacing(1), &
                                  gr%der%order*gr%mesh%spacing(2), &
                                  gr%der%order*gr%mesh%spacing(2))
              write(message(1),'(a)') 'PML absorbing width has to be larger or equal than derivatives order times spacing!'
              write(message(2),'(a,es10.3)') 'PML absorbing width is set to: ', pml_width
              call messages_info(2)
            end if
          end if
          ab_bounds(1,idim) = ab_bounds(2,idim) - pml_width
        end if
      end if

      if (bc%bc_ab_type(idim) == AB_MASK) then
        bounds(1, idim) = ab_bounds(1, idim)
        bounds(2, idim) = bounds(2, idim)
        bc%bc_bounds(:,idim) = bounds(:,idim)
      else if (bc%bc_ab_type(idim) == AB_CPML) then
        bounds(1, idim) = ab_bounds(1, idim)
        bounds(2, idim) = bounds(2, idim)
        bc%bc_bounds(:, idim) = bounds(:, idim)
      else if (bc%bc_ab_type(idim) == AB_MASK_ZERO) then
        bounds(1, idim) = ab_bounds(1, idim)
        bounds(2, idim) = bounds(2, idim)
        bc%bc_bounds(:, idim) = bounds(:, idim)
      else
        bc%bc_bounds(:, idim) = bounds(:, idim)
      end if

      if (gr%mesh%sb%box_shape == PARALLELEPIPED) then
        if (bc%bc_ab_type(idim) == AB_CPML) then
          string = "  PML Lower bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(1, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(1),'(a)') trim(string)
          string = "  PML Upper bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(2, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(2),'(a)') trim(string) 
          call messages_info(2)
        else if (bc%bc_ab_type(idim) == AB_MASK) then
          string = "  Mask Lower bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(1, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(1),'(a)') trim(string)
          string = "  Mask Upper bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(2, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(2),'(a)') trim(string) 
          call messages_info(2)
        else if (bc%bc_ab_type(idim) == AB_MASK_ZERO) then
          string = "  Zero Lower bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(1, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(1),'(a)') trim(string)
          string = "  Zero Upper bound = "
            write(string,'(a,a,i1,a,es10.3,3a)') trim(string), "  dim ", idim, ":",&
                  units_from_atomic(units_inp%length, ab_bounds(2, idim) ), ' [', trim(units_abbrev(units_inp%length)),   ']'
          write(message(2),'(a)') trim(string) 
          call messages_info(2)
        end if
      else
        write(message(1),'(a,es10.3,3a)') & 
          "  Lower bound = ", units_from_atomic(units_inp%length, ab_bounds(1, idim) ),&
          ' [', trim(units_abbrev(units_inp%length)), ']'
        write(message(2),'(a,es10.3,3a)') & 
          "  Upper bound = ", units_from_atomic(units_inp%length, ab_bounds(2, idim) ),&
          ' [', trim(units_abbrev(units_inp%length)), ']'
        call messages_info(2)
      end if

    end do

    ! initialization of surfaces
    call maxwell_surfaces_init(gr%mesh, st, bounds)

    !%Variable MaxwellABPMLKappaMax
    !%Type float
    !%Default 0.4
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Follwos
    !%End
    call parse_variable(namespace, 'MaxwellABPMLKappaMax', CNST(2.0), bc%pml_kappa_max, unit_one)

    !%Variable MaxwellABPMLAlphaMax
    !%Type float
    !%Default 0.4
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Follwos
    !%End
    call parse_variable(namespace, 'MaxwellABPMLAlphaMax', CNST(1.0), bc%pml_alpha_max, unit_one)

    !%Variable MaxwellABPMLPower
    !%Type float
    !%Default 0.4
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Follwos
    !%End
    call parse_variable(namespace, 'MaxwellABPMLPower', CNST(3.5), bc%pml_power, unit_one)

    !%Variable MaxwellABPMLReflectionError
    !%Type float
    !%Default 0.4
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% Follwos
    !%End
    call parse_variable(namespace, 'MaxwellABPMLReflectionError', CNST(0.1), bc%pml_refl_error, unit_one)


    ! mapping of mask boundary points
    if (ab_mask_check) then
      call maxwell_mask_points_mapping(bc, gr%mesh, ab_bounds, geo)
    end if

    ! mapping of pml boundary points
    if (ab_pml_check) then
      call maxwell_pml_points_mapping(bc, gr%mesh, ab_bounds, geo)
    end if

    ! mapping of constant boundary points
    if (constant_check) then
      call maxwell_constant_points_mapping(bc, gr%mesh, bounds, geo)
    end if

    ! mapping of plane waves boundary points
    if (plane_waves_check) then
      call maxwell_plane_waves_points_mapping(bc, gr%mesh, bounds, geo)
      call maxwell_plane_waves_boundaries_init(bc, namespace)
    end if

    ! mapping of zero points
    if (zero_check) then
      call maxwell_zero_points_mapping(bc, gr%mesh, bounds, geo)
    end if

    if (ab_mask_check) then
      call bc_mxll_generate_mask(bc, gr%mesh, ab_bounds)
    end if

    if (ab_pml_check) then
      call bc_mxll_generate_pml(bc, gr, ab_bounds, dt)
    end if

    !call bc_generate_zero(bc, gr%mesh, ab_bounds)

    if(debug%info) call bc_mxll_write_info(bc, gr%mesh, namespace)

    if (ab_mask_check .or. ab_pml_check) then
      call messages_print_stress(stdout, namespace=namespace)
    end if

    POP_SUB(bc_mxll_init)
  end subroutine bc_mxll_init

  ! ---------------------------------------------------------
  subroutine bc_mxll_end(bc)
    type(bc_mxll_t),   intent(inout) :: bc
    PUSH_SUB(bc_mxll_end)

    if (allocated(bc%mask)) then
      SAFE_DEALLOCATE_A(bc%mask)
    end if

    if (bc%do_plane_waves) then
      SAFE_DEALLOCATE_P(bc%plane_waves_modus)
      SAFE_DEALLOCATE_P(bc%plane_waves_e_field_string)
      SAFE_DEALLOCATE_P(bc%plane_waves_e_field)
      SAFE_DEALLOCATE_P(bc%plane_waves_k_vector)
      SAFE_DEALLOCATE_P(bc%plane_waves_v_vector)
      SAFE_DEALLOCATE_P(bc%plane_waves_mx_function)
      SAFE_DEALLOCATE_P(bc%plane_waves_mx_phase)
      SAFE_DEALLOCATE_P(bc%plane_waves)
      SAFE_DEALLOCATE_P(bc%plane_waves_oam)
      SAFE_DEALLOCATE_P(bc%plane_waves_sam)
    end if

    POP_SUB(bc_mxll_end)
  end subroutine bc_mxll_end

  ! ---------------------------------------------------------
  subroutine bc_mxll_write_info(bc, mesh, namespace)
    type(bc_mxll_t),       intent(in) :: bc
    type(mesh_t),          intent(in) :: mesh
    type(namespace_t),     intent(in) :: namespace

    integer :: err, idim
    FLOAT, allocatable :: tmp(:)
    logical :: mask_check, pml_check, medium_check

    PUSH_SUB(bc_mxll_write_info)

    mask_check = .false.
    pml_check = .false.
    medium_check = .false.

    do idim=1, 3
      if (bc%bc_ab_type(idim) == AB_MASK) then
        mask_check = .true.
      end if
      if (bc%bc_ab_type(idim) == AB_CPML) then
        pml_check = .true.
      end if
      if (bc%bc_ab_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) then
        medium_check = .true.
      end if
    end do

    if (mask_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      call get_mask_io_function(bc%mask, bc, tmp)
      call write_files("maxwell_mask", tmp)
      SAFE_DEALLOCATE_A(tmp)
    else if (pml_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      ! sigma for electric field dim = 1
      tmp(:) = M_ONE
      call get_pml_io_function(bc%pml_sigma_e(:, 1), bc, tmp)
      call write_files("maxwell_sigma_e-x", tmp)
      ! sigma for electric field dim = 2
      tmp(:) = M_ONE
      call get_pml_io_function(bc%pml_sigma_e(:, 2), bc, tmp)
      call write_files("maxwell_sigma_e-y", tmp)
      ! sigma for electric field dim = 3
      tmp(:) = M_ONE
      call get_pml_io_function(bc%pml_sigma_e(:, 3), bc, tmp)
      call write_files("maxwell_sigma_e-z", tmp)
      ! sigma for magnetic field dim = 1
      tmp(:) = M_ZERO
      call get_pml_io_function(bc%pml_sigma_m(:, 1), bc, tmp)
      call write_files("maxwell_sigma_m-x", tmp)
      ! sigma for magnetic dim = 2
      tmp(:) = M_ZERO
      call get_pml_io_function(bc%pml_sigma_m(:, 2), bc, tmp)
      call write_files("maxwell_sigma_m-y", tmp)
      ! sigma for magnetic = 3
      tmp(:) = M_ZERO
      call get_pml_io_function(bc%pml_sigma_m(:, 3), bc, tmp)
      call write_files("maxwell_sigma_m-z", tmp)
      ! pml_a for electric field dim = 1
      tmp(:) = M_ZERO
      call get_pml_io_function(TOFLOAT(bc%pml_a(:, 1)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_e-x", tmp)
      ! pml_a for electric field dim = 2
      tmp(:) = M_ZERO
      call get_pml_io_function(TOFLOAT(bc%pml_a(:, 2)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_e-y", tmp)
      ! pml_a for electric field dim = 3
      tmp(:) = M_ZERO
      call get_pml_io_function(TOFLOAT(bc%pml_a(:, 3)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_e-z", tmp)
      ! pml_a for magnetic field dim = 1
      tmp(:) = M_ZERO
      call get_pml_io_function(aimag(bc%pml_a(:, 1)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_m-x", tmp)
      ! pml_a for magnetic field dim = 2
      tmp(:) = M_ZERO
      call get_pml_io_function(aimag(bc%pml_a(:, 2)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_m-y", tmp)
      ! pml_a for magnetic field dim = 3
      tmp(:) = M_ZERO
      call get_pml_io_function(aimag(bc%pml_a(:, 3)), bc, tmp)
      call write_files("maxwell_sigma_pml_a_m-z", tmp)
      SAFE_DEALLOCATE_A(tmp)
    end if
    if (medium_check) then
      SAFE_ALLOCATE(tmp(mesh%np))
      ! medium epsilon
      tmp(:) = P_ep
      call get_medium_io_function(bc%medium_ep, bc, tmp)
      call write_files("maxwell_ep", tmp)
      ! medium mu
      tmp(:) = P_mu
      call get_medium_io_function(bc%medium_mu, bc, tmp)
      call write_files("maxwell_mu", tmp)
      ! medium epsilon
      tmp(:) = P_c
      call get_medium_io_function(bc%medium_c, bc, tmp)
      call write_files("maxwell_c", tmp)
      ! medium epsilon aux field dim = 1
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_ep(:, 1, :), bc, tmp)
      call write_files("maxwell_aux_ep-x", tmp)
      ! medium epsilon aux field dim = 2
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_ep(:, 2, :), bc, tmp)
      call write_files("maxwell_aux_ep-y", tmp)
      ! medium epsilon aux field dim = 3
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_ep(:, 3, :), bc, tmp)
      call write_files("maxwell_aux_ep-z", tmp)
      ! medium mu aux field dim = 1
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_mu(:, 1, :), bc, tmp)
      call write_files("maxwell_aux_mu-x", tmp)
      ! medium mu aux field dim = 2
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_mu(:, 2, :), bc, tmp)
      call write_files("maxwell_aux_mu-y", tmp)
      ! medium mu aux field dim = 3
      tmp(:) = M_ZERO
      call get_medium_io_function(bc%medium_aux_mu(:, 3, :), bc, tmp)
      call write_files("maxwell_aux_mu-z", tmp)
      SAFE_DEALLOCATE_A(tmp)
    end if

    POP_SUB(bc_write_info)

    contains

      subroutine get_pml_io_function(pml_func, bc, io_func)
        FLOAT,              intent(in)    :: pml_func(:)
        type(bc_mxll_t),    intent(in)    :: bc
        FLOAT,              intent(inout) :: io_func(:)

        integer :: ip, ip_in

        do ip_in = 1, bc%pml_points_number
          ip          = bc%pml_points_map(ip_in)
          io_func(ip) = pml_func(ip_in)
        end do

      end subroutine get_pml_io_function

      subroutine get_mask_io_function(mask_func, bc, io_func)
        FLOAT,              intent(in)    :: mask_func(:,:)
        type(bc_mxll_t),    intent(in)    :: bc
        FLOAT,              intent(inout) :: io_func(:)

        integer :: ip, ip_in, idim

        do ip_in = 1, bc%mask_points_number(idim)
          ip          = bc%mask_points_map(ip_in, idim)
          io_func(ip) = mask_func(ip_in, idim)
        end do

      end subroutine get_mask_io_function

      subroutine get_medium_io_function(medium_func, bc, io_func)
        FLOAT,              intent(in)    :: medium_func(:,:)
        type(bc_mxll_t),    intent(in)    :: bc
        FLOAT,              intent(inout) :: io_func(:)

        integer :: ip, ip_in, idim

        do idim = 1, 3
          do ip_in = 1, bc%medium_points_number(idim)
            ip          = bc%medium_points_map(ip_in, idim)
            io_func(ip) = medium_func(ip_in, idim)
          end do
        end do

      end subroutine get_medium_io_function


      subroutine write_files(filename, tmp)
        character(len=*), intent(in) :: filename
        FLOAT,            intent(in) :: tmp(:)

        call dio_function_output(io_function_fill_how("VTK"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("AxisX"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("AxisY"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("AxisZ"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("PlaneX"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("PlaneY"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
        call dio_function_output(io_function_fill_how("PlaneZ"), "./td.general", trim(filename), &
                     namespace, mesh, tmp, unit_one, err)
      end subroutine write_files


  end subroutine bc_mxll_write_info

  ! ---------------------------------------------------------
  subroutine maxwell_mask_points_mapping(bc, mesh, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, ip_in_max, point_info, idim

    PUSH_SUB(maxwell_mask_points_mapping)

    ip_in_max = 1
    do idim = 1, 3
      if (bc%bc_ab_type(idim) == AB_MASK) then
        ! allocate mask points map
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
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

    do idim = 1,3
      if (bc%bc_ab_type(idim) == AB_MASK) then
        ! mask points mapping
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%mask_points_map(ip_in, idim) = ip
          end if
        end do
      end if
    end do

    POP_SUB(maxwell_mask_points_mapping)
  end subroutine maxwell_mask_points_mapping


  ! ---------------------------------------------------------
  subroutine maxwell_pml_points_mapping(bc, mesh, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, point_info

    PUSH_SUB(maxwell_pml_points_mapping)

    ! allocate pml points map
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%pml_points_number = ip_in
    SAFE_ALLOCATE(bc%pml_points_map(1:ip_in))
    SAFE_ALLOCATE(bc%pml_points_map_inv(mesh%np))

    bc%pml_points_map_inv = 0
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%pml_points_map(ip_in) = ip
        bc%pml_points_map_inv(ip) = ip_in
      end if
    end do

    POP_SUB(maxwell_pml_points_mapping)
  end subroutine maxwell_pml_points_mapping


  ! ---------------------------------------------------------
  subroutine maxwell_constant_points_mapping(bc, mesh, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, point_info

    PUSH_SUB(maxwell_constant_points_mapping)

    ! allocate constant points map
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
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
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%constant_points_map(ip_in) = ip
      end if
    end do

    POP_SUB(maxwell_constant_points_mapping)
  end subroutine maxwell_constant_points_mapping

  ! ---------------------------------------------------------
  subroutine maxwell_plane_waves_points_mapping(bc, mesh, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, point_info

    PUSH_SUB(maxwell_plane_waves_points_mapping)

    ! allocate zero points map
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
      end if
    end do
    bc%plane_waves_points_number = ip_in
    SAFE_ALLOCATE(bc%plane_waves_points_map(1:ip_in)) ! memleak

    ! zero points mapping
    ip_in = 0
    do ip = 1, mesh%np
      call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
      if (point_info == 1) then
        ip_in = ip_in + 1
        bc%plane_waves_points_map(ip_in) = ip
      end if
    end do

    POP_SUB(maxwell_plane_waves_points_mapping)
  end subroutine maxwell_plane_waves_points_mapping


  ! ---------------------------------------------------------
  subroutine maxwell_zero_points_mapping(bc, mesh, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, ip_in_max, point_info, idim

    PUSH_SUB(maxwell_zero_points_mapping)

    ip_in_max = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO) then
        ! allocate zero points map
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
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
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_ZERO) then
        ! zero points mapping
        ip_in = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%zero_points_map(ip_in, idim) = ip
          end if
        end do
      end if
    end do

    POP_SUB(maxwell_zero_points_mapping)
  end subroutine maxwell_zero_points_mapping


  ! ---------------------------------------------------------
  subroutine maxwell_medium_points_mapping(bc, mesh, st, bounds, geo)
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    type(states_mxll_t), intent(inout) :: st
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo

    integer :: ip, ip_in, ip_in_max, ip_bd, ip_bd_max, point_info, boundary_info, idim

    PUSH_SUB(maxwell_medium_points_mapping)

    ip_in_max = 0
    ip_bd_max = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) then
        ! allocate pml points map
        ip_in = 0
        ip_bd = 0
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
          call maxwell_boundary_point_info(mesh, ip, bounds, boundary_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
          end if
          if ((boundary_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_bd = ip_bd + 1
          end if
        end do
        bc%medium_points_number(idim) = ip_in
        bc%medium_bdry_number(idim) = ip_bd
      end if
    end do
    SAFE_ALLOCATE(bc%medium_aux_ep(1:ip_in, 1:st%d%dim, 3))
    SAFE_ALLOCATE(bc%medium_aux_mu(1:ip_in, 1:st%d%dim, 3))
    SAFE_ALLOCATE(bc%medium_points_map(1:ip_in_max, 3))
    SAFE_ALLOCATE(bc%medium_bdry_map(1:ip_bd_max, 3))

    ip_in = 0
    ip_bd = 0
    do idim = 1, 3
      if (bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) then
        do ip = 1, mesh%np
          call maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info)
          call maxwell_boundary_point_info(mesh, ip, bounds, boundary_info)
          if ((point_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_in = ip_in + 1
            bc%medium_points_map(ip_in,idim) = ip
          end if
          if ((boundary_info == 1) .and. (abs(mesh%x(ip, idim)) >= bounds(1, idim))) then
            ip_bd = ip_bd + 1
           bc%medium_bdry_map(ip_bd, idim) = ip
          end if
        end do
      end if
    end do

    POP_SUB(maxwell_medium_points_mapping)
  end subroutine maxwell_medium_points_mapping


  ! ---------------------------------------------------------
  subroutine bc_mxll_generate_pml(bc, gr, bounds, dt)
    type(bc_mxll_t),    intent(inout) :: bc
    type(grid_t),       intent(in)    :: gr
    FLOAT,              intent(in)    :: bounds(:,:)
    FLOAT, optional,    intent(in)    :: dt

    integer :: ip, ip_in, idim
    FLOAT   :: width(3), ddv(3), ss_e, ss_m, ss_max, aa_e, aa_m, bb_e, bb_m, gg, hh, kk, ll_e, ll_m
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)

    PUSH_SUB(bc_mxll_generate_pml)

    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np, 1:gr%mesh%sb%dim))

    SAFE_ALLOCATE(bc%pml_kappa(1:bc%pml_points_number, 1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(bc%pml_sigma_e(1:bc%pml_points_number, 1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(bc%pml_sigma_m(1:bc%pml_points_number, 1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(bc%pml_a(1:bc%pml_points_number, 1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(bc%pml_b(1:bc%pml_points_number, 1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(bc%pml_c(1:bc%pml_points_number, 1:3))
    SAFE_ALLOCATE(bc%pml_mask(1:bc%pml_points_number))
    SAFE_ALLOCATE(bc%pml_conv_plus(1:bc%pml_points_number, 1:3, 1:3))
    SAFE_ALLOCATE(bc%pml_conv_minus(1:bc%pml_points_number, 1:3, 1:3))
    SAFE_ALLOCATE(bc%pml_conv_plus_old(1:bc%pml_points_number, 1:3, 1:3))
    SAFE_ALLOCATE(bc%pml_conv_minus_old(1:bc%pml_points_number, 1:3, 1:3))
    SAFE_ALLOCATE(bc%pml_aux_ep(1:bc%pml_points_number, 1:3, 1:3))
    SAFE_ALLOCATE(bc%pml_aux_mu(1:bc%pml_points_number, 1:3, 1:3))

    bc%pml_kappa                 = M_ONE
    bc%pml_sigma_e               = M_ZERO
    bc%pml_sigma_m               = M_ZERO
    bc%pml_a                     = M_z0
    bc%pml_b                     = M_z0
    bc%pml_c                     = M_ZERO
    bc%pml_mask                  = M_ONE
    bc%pml_conv_plus             = M_z0
    bc%pml_conv_minus            = M_z0
    bc%pml_conv_plus_old         = M_z0
    bc%pml_conv_minus_old        = M_z0

    width(:) = bounds(2, :) - bounds(1, :)

    ! PML variables for all boundary points
    do ip_in = 1, bc%pml_points_number
      ip = bc%pml_points_map(ip_in)
      ddv(:) = abs(gr%mesh%x(ip, :)) - bounds(1, :)
      do idim = 1, gr%mesh%sb%dim
        if (ddv(idim) >= M_ZERO) then
          gg     = (ddv(idim)/bc%pml_width)**bc%pml_power
          hh     = (M_ONE-ddv(idim)/bc%pml_width)**bc%pml_power
          kk     = M_ONE ! + (bc%pml_kappa_max - M_ONE) * gg
          ss_max = -(bc%pml_power + M_ONE)*P_c*P_ep*log(bc%pml_refl_error)/(M_TWO * bc%pml_width)
          ss_e   = gg * ss_max
          ss_m   = gg * ss_max ! * P_mu/P_ep
          ll_e   = ss_e*kk ! + kk**2*bc%pml_alpha_max*hh
          ll_m   = ss_m*kk ! + kk**2*bc%pml_alpha_max*hh
          bb_e   = exp(-(ss_e/(P_ep))*dt)
          bb_m   = exp(-(ss_m/(P_ep))*dt)
!          aa_e   = (ss_e/ll_e)*(bb_e - 1)
!          aa_m   = (ss_m/ll_m)*(bb_m - 1)
          aa_e   = (bb_e - 1)
          aa_m   = (bb_m - 1)
          if (ll_e == M_ZERO) aa_e = M_ZERO
          if (ll_m == M_ZERO) aa_m = M_ZERO
          bc%pml_sigma_e(ip_in, idim) = ss_e
          bc%pml_sigma_m(ip_in, idim) = ss_m
          bc%pml_a(ip_in, idim)       = aa_e + M_zI * aa_m
          bc%pml_b(ip_in, idim)       = bb_e + M_zI * bb_m
          bc%pml_kappa(ip_in, idim)   = kk
          bc%pml_mask(ip_in)          = bc%pml_mask(ip_in) * (M_ONE - sin(ddv(idim)*M_PI/(M_TWO*(width(idim))))**2)
        else
          bc%pml_kappa(ip_in, idim)   = M_ONE
          bc%pml_sigma_e(ip_in, idim) = M_ZERO
          bc%pml_sigma_m(ip_in, idiM) = M_ZERO
          bc%pml_a(ip_in, idim)       = M_z0
          bc%pml_b(ip_in, idim)       = M_z0
          bc%pml_mask(ip_in)          = M_ONE
        end if
      end do
    end do

    ! PML auxiliary epsilon for all boundary points
    do idim = 1, gr%mesh%sb%dim
      tmp = P_ep
      do ip_in = 1, bc%pml_points_number
        ip = bc%pml_points_map(ip_in)
        tmp(ip) = P_ep / bc%pml_kappa(ip_in, idim)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%pml_points_number
        ip = bc%pml_points_map(ip_in)
        bc%pml_aux_ep(ip_in, :, idim) = tmp_grad(ip, :)/(M_FOUR*P_ep*bc%pml_kappa(ip_in, idim))
      end do
    end do

    ! PML auxiliary mu
    do idim = 1, gr%mesh%sb%dim
      tmp = P_mu
      do ip_in = 1, bc%pml_points_number
        ip = bc%pml_points_map(ip_in)
        tmp(ip) = P_mu / bc%pml_kappa(ip_in, idim)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%pml_points_number
        ip = bc%pml_points_map(ip_in)
        bc%pml_aux_mu(ip_in, :, idim) = tmp_grad(ip, :)/(M_FOUR*P_mu*bc%pml_kappa(ip_in, idim))
      end do
    end do

    ! PML auxiliary c for all boundary points
    do idim = 1, gr%mesh%sb%dim
      do ip_in = 1, bc%pml_points_number
        bc%pml_c(ip_in, idim) = P_c/bc%pml_kappa(ip_in, idim)
      end do
    end do

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

    PUSH_SUB(bc_mxll_generate_mask)

    ip_in_max = maxval(bc%mask_points_number(:))

    SAFE_ALLOCATE(bc%mask(ip_in_max, 3))
    SAFE_ALLOCATE(mask(mesh%np))

    mask(:) = M_ONE

    width(:) = bounds(2, :) - bounds(1, :)
    tmp(:)   = M_ZERO

    do ip = 1, mesh%np
      tmp = M_ONE
      mask(ip) = M_ONE
      ddv(:) = abs(mesh%x(ip, :)) - bounds(1, :)
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

    POP_SUB(bc_mxll_generate_mask)
  end subroutine bc_mxll_generate_mask


  ! ---------------------------------------------------------
  subroutine bc_mxll_generate_medium(bc, gr, bounds, geo)
    type(bc_mxll_t),         intent(inout) :: bc
    type(grid_t),            intent(in)    :: gr
    FLOAT,                   intent(in)    :: bounds(:,:)
    type(geometry_t),        intent(in)    :: geo

    integer :: ip, ipp, ip_in, ip_in_max, ip_bd, idim, point_info
    FLOAT   :: dd, dd_min, dd_max, xx(3), xxp(3)
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)

    PUSH_SUB(bc_mxll_generate_medium)

    ip_in_max = maxval(bc%medium_points_number(:))

    SAFE_ALLOCATE(bc%medium_aux_ep(ip_in_max,gr%mesh%sb%dim, 3))
    SAFE_ALLOCATE(bc%medium_aux_mu(ip_in_max,gr%mesh%sb%dim, 3))
    SAFE_ALLOCATE(bc%medium_ep(ip_in_max, 3))
    SAFE_ALLOCATE(bc%medium_mu(ip_in_max, 3))
    SAFE_ALLOCATE(bc%medium_sigma_e(ip_in_max, 3))
    SAFE_ALLOCATE(bc%medium_sigma_m(ip_in_max, 3))
    SAFE_ALLOCATE(bc%medium_c(ip_in_max, 3))
    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%mesh%sb%dim))
    bc%medium_aux_ep = M_ZERO
    bc%medium_aux_mu = M_ZERO
    bc%medium_c = P_c

    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    do idim = 1, 3
      tmp = P_ep
      do  ip = 1, gr%mesh%np_part
        call maxwell_box_point_info(bc, gr%mesh, ip, bounds, geo, point_info)
        if ((point_info /= 0) .and. (abs(gr%mesh%x(ip, idim)) <= bounds(1, idim))) then
          xx(:) = gr%mesh%x(ip, :)
          dd_min = M_HUGE
          do ip_bd=1, bc%medium_bdry_number(idim)
            ipp = bc%medium_bdry_map(ip_bd, idim)
            xxp(:) = gr%mesh%x(ipp, :)
            dd = sqrt((xx(1) - xxp(1))**2 + (xx(2) - xxp(2))**2 + (xx(3) - xxp(3))**2)
            if (dd < dd_min) dd_min = dd
          end do
          tmp(ip) = P_ep * (M_ONE + bc%medium_ep_factor &
                  * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min-M_TWO*dd_max)) ) )
        end if
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%medium_points_number(idim)
        ip = bc%medium_points_map(ip_in, idim)
        bc%medium_aux_ep(ip_in, :, idim) = &
          tmp_grad(ip, :)/(M_FOUR*P_ep*bc%medium_ep_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max-dd)))
      end do
    end do

    do idim = 1, 3
      tmp = P_mu
      do  ip = 1, gr%mesh%np_part
        call maxwell_box_point_info(bc, gr%mesh, ip, bounds, geo, point_info)
        if ((point_info == 1) .and. (abs(gr%mesh%x(ip, idim)) <= bounds(1, idim))) then
          xx(:) = gr%mesh%x(ip, :)
          dd_min = M_HUGE
          do ip_bd = 1, bc%medium_bdry_number(idim)
            ipp = bc%medium_bdry_map(ip_bd, idim)
            xxp(:) = gr%mesh%x(ipp,:)
            dd = sqrt((xx(1) - xxp(1))**2 + (xx(2) - xxp(2))**2 + (xx(3) - xxp(3))**2)
            if (dd < dd_min) dd_min = dd
          end do
          tmp(ip) = P_mu * (M_ONE + bc%medium_mu_factor &
                  * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        end if
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, bc%medium_points_number(idim)
        ip = bc%medium_points_map(ip_in, idim)
        bc%medium_aux_mu(ip_in, :, idim) = & 
          tmp_grad(ip, :)/(M_FOUR*P_mu*bc%medium_mu_factor * M_ONE/(M_ONE + exp(-M_FIVE/dd_max-dd)))
      end do
    end do

    do idim = 1, 3
      do ip_in=1, bc%medium_points_number(idim)
        ip = bc%medium_points_map(ip_in,idim)
        xx(:) = gr%mesh%x(ip,:)
        dd_min = M_HUGE
        do ip_bd = 1, bc%medium_bdry_number(idim)
          ipp = bc%medium_bdry_map(ip_bd, idim)
          xxp(:) = gr%mesh%x(ipp, :)
          dd = sqrt((xx(1) - xxp(1))**2 + (xx(2) - xxp(2))**2 + (xx(3) - xxp(3))**2)
          if (dd < dd_min) dd_min = dd
        end do
        bc%medium_ep(ip_in, idim) = P_ep * (M_ONE+bc%medium_ep_factor & 
                                          * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium_mu(ip_in, idim) = P_mu * (M_ONE+bc%medium_mu_factor &
                                          * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium_sigma_e(ip_in, idim) = (M_ONE+bc%medium_sigma_e_factor & 
                                               * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium_sigma_m(ip_in, idim) = (M_ONE+bc%medium_sigma_m_factor &
                                               * M_ONE/(M_ONE + exp( -M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) ) )
        bc%medium_c(ip_in, idim) = M_ONE/sqrt(bc%medium_ep(ip_in, idim)*bc%medium_mu(ip_in, idim))
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)
 
    POP_SUB(bc_mxll_generate_medium)
  end subroutine bc_mxll_generate_medium


  ! ---------------------------------------------------------
  subroutine maxwell_plane_waves_boundaries_init(bc, namespace)
    type(bc_mxll_t),        intent(inout) :: bc
    type(namespace_t),      intent(in)    :: namespace

    type(block_t)        :: blk
    integer              :: il, nlines, ncols, ierr
    integer              :: oam, sam
    FLOAT                :: k_vector(MAX_DIM), e_field(MAX_DIM), vv(MAX_DIM), xx(MAX_DIM), rr, dummy(MAX_DIM), test, test_limit!, angle, sigma
    character(len=1024)  :: k_string(MAX_DIM)
    character(len=1024)  :: mxf_expression

    PUSH_SUB(maxwell_plane_waves_boundaries_init)

    test_limit = CNST(10.0e-9)

    !%Variable UserDefinedMaxwellIncidentWaves
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% The initial electromagnetic fields can be set by the user 
    !% with the <tt>UserDefinedMaxwellIncidentWaves</tt> block variable.
    !% The electromagnetic fields have to fulfill the 
    !% Maxwells equations in vacuum.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellIncidentWaves
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k1x" | "k1y" | "k1z" | "E1x" | "E1z" | "E1x" | plane_wave
    !% <br>&nbsp;&nbsp;   plane_wave_parser      | "k2x" | "k2y" | "k2z" | "E2x" | "E2y" | "E2z" | no_plane_wave
    !% <br>&nbsp;&nbsp;   plane_wave_gauss       | "k3x" | "k3y" | "k3z" | "E3x" | "E3y" | "E3z" | "width"           | "shift" | plane_wave
    !% <br>&nbsp;&nbsp;   plane_wave_mx_function | "k4x" | "k4y" | "k4z" | "E4x" | "E4y" | "E4z" | mx_envelope_name  | phase   | plane_wave
    !% <br>&nbsp;&nbsp;   bessel_mx_function      | "k5x" | "k5y" | "k5z" | "E5x" | "E5y" | "E5z" | "angle" | "orbital_momentum" | "sigma" | mx_envelope_name | phase
    !% <br>%</tt>
    !%
    !% Description about UserDefinedMaxwellIncidentWaves follows
    !%
    !%Option plane_wave_parser 0
    !% Parser input modus
    !%Option plane_wave_mx_function 1
    !% The incident wave envelope is defined by an mx_function
    !%Option bessel_mx_function 4
    !% Follows!
    !%End

    if(parse_block(namespace, 'UserDefinedMaxwellIncidentWaves', blk) == 0) then ! memleak in parse_block

      call messages_print_stress(stdout, trim('Substitution of the electromagnetic incident waves'), namespace=namespace)

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      bc%plane_waves_number = nlines
      SAFE_ALLOCATE(bc%plane_waves_modus(nlines))
      SAFE_ALLOCATE(bc%plane_waves_e_field_string(MAX_DIM, nlines))
      SAFE_ALLOCATE(bc%plane_waves_e_field(MAX_DIM, nlines))
      SAFE_ALLOCATE(bc%plane_waves_k_vector(MAX_DIM, nlines))
      SAFE_ALLOCATE(bc%plane_waves_v_vector(MAX_DIM, nlines))
      SAFE_ALLOCATE(bc%plane_waves_mx_function(nlines))
      SAFE_ALLOCATE(bc%plane_waves_mx_phase(nlines))
      SAFE_ALLOCATE(bc%plane_waves(nlines))
      SAFE_ALLOCATE(bc%plane_waves_oam(nlines))
      SAFE_ALLOCATE(bc%plane_waves_sam(nlines))

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if ((ncols /= 6) .and. (ncols /= 8) .and. (ncols /= 10) .and. (ncols /= 9)) then
          message(1) = 'Each line in the UserDefinedMaxwellIncidentWaves block must have'
          message(2) = 'six, eight, nine or ten columns.'
          call messages_fatal(2, namespace=namespace)
        end if

        ! check input modus e.g. parser of defined functions
        call parse_block_integer(blk, il - 1, 0, bc%plane_waves_modus(il))

        ! parse formula string
        if (bc%plane_waves_modus(il) == OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__PLANE_WAVE_PARSER) then

          call parse_block_string( blk, il - 1, 1, k_string(1))
          call parse_block_string( blk, il - 1, 2, k_string(2))
          call parse_block_string( blk, il - 1, 3, k_string(3))
          call parse_block_string( blk, il - 1, 4, bc%plane_waves_e_field_string(1, il))
          call parse_block_string( blk, il - 1, 5, bc%plane_waves_e_field_string(2, il))
          call parse_block_string( blk, il - 1, 6, bc%plane_waves_e_field_string(3, il))
          call parse_block_integer(blk, il - 1, 7, bc%plane_waves(il))

          write(message(1), '(a,i2,a) ') 'Substituting electromagnetic incident wave ', il, ' with the expressions: '
          call messages_info(1)
          write(message(1), '(6a)')     '  Wave vector k(x)   = ', trim(k_string(1))
          write(message(2), '(2a)')     '  Wave vector k(y)   = ', trim(k_string(2))
          write(message(3), '(2a)')     '  Wave vector k(z)   = ', trim(k_string(3))
          write(message(4), '(2a)')     '  E-field(x) for t_0 = ', trim(bc%plane_waves_e_field_string(1, il))
          write(message(5), '(2a)')     '  E-field(y) for t_0 = ', trim(bc%plane_waves_e_field_string(2, il))
          write(message(6), '(2a)')     '  E-field(z) for t_0 = ', trim(bc%plane_waves_e_field_string(3, il))
          call messages_info(6)

          call conv_to_C_string(k_string(1))
          call conv_to_C_string(k_string(2))
          call conv_to_C_string(k_string(3))
          call conv_to_C_string(bc%plane_waves_e_field_string(1, il))
          call conv_to_C_string(bc%plane_waves_e_field_string(2, il))
          call conv_to_C_string(bc%plane_waves_e_field_string(3, il))

          xx(:) = M_ZERO
          rr    = M_ZERO
          call parse_expression(k_vector(1), dummy(1), 1, xx, rr, M_ZERO, k_string(1))
          call parse_expression(k_vector(2), dummy(2), 2, xx, rr, M_ZERO, k_string(2))
          call parse_expression(k_vector(3), dummy(3), 3, xx, rr, M_ZERO, k_string(3))
          k_vector(1) = units_to_atomic(unit_one/units_inp%length, k_vector(1))
          k_vector(2) = units_to_atomic(unit_one/units_inp%length, k_vector(2))
          k_vector(3) = units_to_atomic(unit_one/units_inp%length, k_vector(3))

          vv(:)    = k_vector(:) / sqrt(sum(k_vector(:)**2)) * P_c
          bc%plane_waves_k_vector(:,il) = k_vector(:)
          bc%plane_waves_v_vector(:,il) = vv(:)

        else if (bc%plane_waves_modus(il) == OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__PLANE_WAVE_MX_FUNCTION) then
          call parse_block_float( blk, il - 1, 1, e_field(1))
          call parse_block_float( blk, il - 1, 2, e_field(2))
          call parse_block_float( blk, il - 1, 3, e_field(3))
          call parse_block_string( blk, il - 1, 4, mxf_expression)
          call parse_block_integer(blk, il - 1, 5, bc%plane_waves(il))

          write(message(1), '(a,i2) ') 'Substituting electromagnetic incident wave ', il
          write(message(3), '(a)'    ) 'with the expression: '
          call messages_info(2)
          write(message(1), '(a,es9.2)')     '  E-field(x) amplitude       = ', e_field(1)
          write(message(2), '(a,es9.2)')     '  E-field(y) amplitude       = ', e_field(2)
          write(message(3), '(a,es9.2)')     '  E-field(z) amplitude       = ', e_field(3)
          write(message(4), '(2a)'    )      '  Maxwell wave function name = ', trim(mxf_expression)
          call messages_info(4)
          call mxf_read(bc%plane_waves_mx_function(il), namespace, trim(mxf_expression), ierr)
          if (ierr /= 0) then            
            write(message(1),'(3A)') 'Error in the ""', trim(mxf_expression), &
              '"" field defined in the UserDefinedMaxwellIncidentWaves block'
            call messages_fatal(1, namespace=namespace)
          end if
          e_field  = units_to_atomic(units_inp%energy/units_inp%length, e_field)
          k_vector(:) = bc%plane_waves_mx_function(il)%k_vector(:)

          test = ddot_product(k_vector, e_field)
          if (abs(test) > test_limit) then
            message(1) = 'The wave vector k(:) or its electric field E-field(:) '
            message(2) = 'is not perpendicular enough.'
            call messages_fatal(2, namespace=namespace)
          end if
          if (Sqrt(sum(k_vector(:)**2)) < 1e-10) then
            message(1) = 'The k vector is not defined correctly.'
            call messages_fatal(1, namespace=namespace)
          end if

          bc%plane_waves_e_field(:,il)  = e_field(:)
          bc%plane_waves_k_vector(:,il) = k_vector(:)
          bc%plane_waves_v_vector(:,il) = k_vector(:) / Sqrt(sum(k_vector(:)**2)) * P_c

        else if (bc%plane_waves_modus(il) == OPTION__USERDEFINEDMAXWELLINCIDENTWAVES__BESSEL_MX_FUNCTION) then

          bc%bessel_beam = .true.

          call parse_block_float(blk, il - 1, 1, e_field(1))
          call parse_block_float(blk, il - 1, 2, e_field(2))
          call parse_block_float(blk, il - 1, 3, e_field(3))
          call parse_block_float(blk, il - 1, 4, k_vector(1))
          call parse_block_float(blk, il - 1, 5, k_vector(2))
          call parse_block_float(blk, il - 1, 6, k_vector(3))
          call parse_block_integer(blk, il - 1,  7, oam)
          call parse_block_integer(blk, il - 1,  8, sam)

          write(message(1), '(a,i2) ') 'Substituting electromagnetic incident wave ', il
          write(message(3), '(a)'    ) 'with the expression: '
          call messages_info(2)
          write(message(1), '(a,es9.2)')     '  E-field(x) amplitude       = ', e_field(1)
          write(message(2), '(a,es9.2)')     '  E-field(y) amplitude       = ', e_field(2)
          write(message(3), '(a,es9.2)')     '  E-field(z) amplitude       = ', e_field(3)
          !write(message(4), '(2a)'    )      '  Maxwell wave function name = ', trim(mxf_expression)
          call messages_info(3)
          !call mxf_read(bc%plane_waves_mx_function(il), trim(mxf_expression), ierr)
          !if (ierr /= 0) then            
          !  write(message(1),'(3A)') 'Error in the "', trim(mxf_expression), &
          !    '" field defined in the UserDefinedMaxwellIncidentWaves block'
          !  call messages_fatal(1, namespace=namespace)
          !end if
          e_field  = units_to_atomic(units_inp%energy/units_inp%length, e_field)

          bc%plane_waves_e_field(:, il)  = e_field(:)
          bc%plane_waves_k_vector(:, il) = k_vector(:)

          bc%plane_waves_v_vector(:, il) = M_ZERO 
          bc%plane_waves_v_vector(3, il) = sqrt(k_vector(1)**2 + k_vector(2)**2)*P_c/k_vector(2)

          bc%plane_waves_oam(il) = oam
          bc%plane_waves_sam(il) = sam
         end if
      end do

      call messages_print_stress(stdout, namespace=namespace)

    end if

    POP_SUB(maxwell_plane_waves_boundaries_init)
  end subroutine maxwell_plane_waves_boundaries_init


  ! ---------------------------------------------------------
  subroutine maxwell_surfaces_init(mesh, st, bounds)
    type(mesh_t),             intent(in)    :: mesh
    type(states_mxll_t),      intent(inout) :: st
    FLOAT,                    intent(in)    :: bounds(:,:)

    PUSH_SUB(maxwell_surfaces_init)

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

    POP_SUB(maxwell_surfaces_init)
  end subroutine maxwell_surfaces_init


  ! ---------------------------------------------------------
  subroutine maxwell_box_point_info(bc, mesh, ip, bounds, geo, point_info) 
    type(bc_mxll_t),     intent(inout) :: bc
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(in)    :: ip
    FLOAT,               intent(in)    :: bounds(:,:)
    type(geometry_t),    intent(in)    :: geo
    integer,             intent(out)   :: point_info
 
    FLOAT   :: rr, dd, xx(MAX_DIM), width(MAX_DIM)
    
    point_info = 0
    
    width(:) = bounds(2, :) - bounds(1, :)
    xx = M_ZERO
    xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
    rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))

    if(bc%ab_user_def) then

      dd = bc%ab_ufn(ip) - bounds(1, 1)
      if(dd > M_ZERO) then
        if(bc%ab_ufn(ip) < bounds(2, 1) ) then
           point_info = 1
        end if
      end if

    else ! bc%ab_user_def == .false.

      if(mesh%sb%box_shape == SPHERE) then

        dd = rr -  bounds(1, 1)
        if(dd > M_ZERO ) then
          if (dd  <  width(1)) then
            point_info = 1
          end if
        end if

      else if (mesh%sb%box_shape == PARALLELEPIPED) then

        ! Limits of boundary region
        if ( (abs(xx(1)) <= bounds(2, 1)) .and. (abs(xx(2)) <= bounds(2, 2)) .and. (abs(xx(3)) <= bounds(2, 3)) ) then
          if ( (abs(xx(1)) > bounds(1, 1)) .or. (abs(xx(2)) > bounds(1, 2)) .or. (abs(xx(3)) > bounds(1, 3)) ) then
            point_info = 1
          else
            point_info = 0
          end if
        else
          point_info = -1
        end if

      else

        if(mesh_inborder(mesh, geo, ip, dd, width(1))) then
          point_info = 1
        end if

      end if
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
    if ( abs(xx(1)) == bounds(1, 1) .and. (abs(xx(2)) <= bounds(1, 2) .and. abs(xx(3)) <= bounds(1, 3)) .or. &
         abs(xx(2)) == bounds(1, 2) .and. (abs(xx(1)) <= bounds(1, 1) .and. abs(xx(3)) <= bounds(1, 3)) .or. &
         abs(xx(3)) == bounds(1, 3) .and. (abs(xx(1)) <= bounds(1, 1) .and. abs(xx(2)) <= bounds(1, 2)) ) then
      boundary_info = 1
    else
      boundary_info = 0
    end if

  end subroutine maxwell_boundary_point_info


end module maxwell_boundary_op_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
