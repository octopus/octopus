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

module grid_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use double_grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use mesh_oct_m
  use mesh_init_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use simul_box_oct_m
  use stencil_oct_m
  use stencil_cube_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                &
    grid_t,                &
    grid_init_stage_0,     &
    grid_init_stage_1,     &
    grid_init_stage_2,     &
    grid_end,              &
    grid_write_info,       &
    grid_create_largergrid

  type grid_t
    type(simul_box_t)           :: sb
    type(mesh_t)                :: mesh
    type(derivatives_t)         :: der
    type(double_grid_t)         :: dgrid
    type(stencil_t)             :: stencil
  end type grid_t


contains

  !-------------------------------------------------------------------
  !>
  !! "Zero-th" stage of grid initialization. It initializes the simulation box.
  subroutine grid_init_stage_0(gr, geo, space)
    type(grid_t),          intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(space_t),         intent(in)    :: space

    PUSH_SUB(grid_init_stage_0)

    call simul_box_init(gr%sb, geo, space)
      
    POP_SUB(grid_init_stage_0)
  end subroutine grid_init_stage_0


  !-------------------------------------------------------------------
  subroutine grid_init_stage_1(gr, geo)
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    type(stencil_t) :: cube
    integer :: enlarge(1:MAX_DIM)
    type(block_t) :: blk
    integer :: idir
    FLOAT :: def_h, def_rsize
    FLOAT :: grid_spacing(1:MAX_DIM)

    PUSH_SUB(grid_init_stage_1)

    call geometry_grid_defaults(geo, def_h, def_rsize)
    call geometry_grid_defaults_info(geo)
    
    ! initialize to -1
    grid_spacing = -M_ONE

    !%Variable Spacing
    !%Type float
    !%Section Mesh
    !%Description
    !% The spacing between the points in the mesh. This controls the
    !% quality of the discretization: smaller spacing gives more
    !% precise results but increased computational cost.
    !%
    !% In periodic directions, your spacing may be slightly different
    !% than what you request here, since the box size must be an
    !% integer multiple of the spacing.
    !%
    !% The default value is defined by the species if only default pseudopotentials are used.
    !% Otherwise, there is no default.
    !%
    !% It is possible to have a different spacing in each one of the Cartesian directions
    !% if we define <tt>Spacing</tt> as block of the form
    !%
    !% <tt>%Spacing
    !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
    !% <br>%</tt>
    !%End

    if(parse_block('Spacing', blk) == 0) then
      if(parse_block_cols(blk,0) < gr%sb%dim) call messages_input_error('Spacing')
      do idir = 1, gr%sb%dim
        call parse_block_float(blk, 0, idir - 1, grid_spacing(idir), units_inp%length)
        if(def_h > M_ZERO) call messages_check_def(grid_spacing(idir), .true., def_h, 'Spacing', units_out%length)
      end do
      call parse_block_end(blk)
    else
      call parse_variable('Spacing', -M_ONE, grid_spacing(1), units_inp%length)
      grid_spacing(1:gr%sb%dim) = grid_spacing(1)
      if(def_h > M_ZERO) call messages_check_def(grid_spacing(1), .true., def_h, 'Spacing', units_out%length)
    end if

    do idir = 1, gr%sb%dim
      if(grid_spacing(idir) < M_EPSILON) then
        if(def_h > M_ZERO .and. def_h < huge(def_h)) then
          grid_spacing(idir) = def_h
          write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", idir, &
            ") [", trim(units_abbrev(units_out%length)), "] = ",                        &
            units_from_atomic(units_out%length, grid_spacing(idir))
          call messages_info(1)
        ! Note: the default automatically matches the 'recommended' value compared by messages_check_def above.
        else
          message(1) = 'Either:'
          message(2) = "   *) variable 'Spacing' is not defined and"
          message(3) = "      I can't find a suitable default"
          message(4) = "   *) your input for 'Spacing' is negative or zero"
          call messages_fatal(4)
        end if
      end if
    end do

    ! initialize derivatives
    call derivatives_init(gr%der, gr%sb)

    call double_grid_init(gr%dgrid, gr%sb)

    enlarge = 0
    enlarge(1:gr%sb%dim) = 2
    enlarge = max(enlarge, double_grid_enlarge(gr%dgrid))
    enlarge = max(enlarge, gr%der%n_ghost)

    ! now we generate the mesh and the derivatives
    call mesh_init_stage_1(gr%mesh, gr%sb, grid_spacing, enlarge)

    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian.
    call stencil_cube_get_lapl(cube, gr%sb%dim, order = 2)
    call stencil_union(gr%sb%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    call mesh_init_stage_2(gr%mesh, gr%sb, geo, gr%stencil)

    POP_SUB(grid_init_stage_1)

  end subroutine grid_init_stage_1

  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, mc, geo)
    type(grid_t), target, intent(inout) :: gr
    type(multicomm_t),    intent(in)    :: mc
    type(geometry_t),     intent(in)    :: geo

    PUSH_SUB(grid_init_stage_2)

    call mesh_init_stage_3(gr%mesh, gr%stencil, mc)

    call nl_operator_global_init()
    call derivatives_build(gr%der, gr%mesh)
    call mesh_check_symmetries(gr%mesh, gr%mesh%sb)

    ! print info concerning the grid
    call grid_write_info(gr, geo, stdout)

    POP_SUB(grid_init_stage_2)
  end subroutine grid_init_stage_2

  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    PUSH_SUB(grid_end)

    call nl_operator_global_end()
    call double_grid_end(gr%dgrid)
    call derivatives_end(gr%der)
    call mesh_end(gr%mesh)
    call stencil_end(gr%stencil)

    POP_SUB(grid_end)
  end subroutine grid_end

  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, geo, iunit)
    type(grid_t),     intent(in) :: gr
    type(geometry_t), intent(in) :: geo
    integer,          intent(in) :: iunit

    PUSH_SUB(grid_write_info)

    if(.not.mpi_grp_is_root(mpi_world)) then
      if(debug%info) call messages_debug_newlines(6)
      POP_SUB(grid_write_info)
      return
    end if

    call messages_print_stress(iunit, "Grid")
    call simul_box_write_info(gr%sb, geo, iunit)

    message(1) = "Main mesh:"
    call messages_info(1, iunit)
    call mesh_write_info(gr%mesh, iunit)
    call messages_print_stress(iunit)

    POP_SUB(grid_write_info)
  end subroutine grid_write_info

  !-------------------------------------------------------------------
  subroutine grid_create_largergrid(grin, geo, mc, grout)
    type(grid_t),      intent(in)  :: grin
    type(geometry_t),  intent(in)  :: geo
    type(multicomm_t), intent(in)  :: mc
    type(grid_t),      intent(out) :: grout

    FLOAT :: avg

    PUSH_SUB(grid_create_largergrid)

    call simul_box_copy(grout%sb, grin%sb)

    ! Modification of the simulation box.
    select case(grout%sb%box_shape)
    case(CYLINDER)
      grout%sb%box_shape = CYLINDER
      avg = M_HALF * (grin%sb%rsize + grin%sb%xsize)
      if(grin%sb%rsize <= grin%sb%xsize) then
        grout%sb%rsize = avg
        grout%sb%xsize = grin%sb%xsize
      else
        grout%sb%xsize = avg
        grout%sb%rsize = grin%sb%rsize
      end if
      grout%sb%lsize(1)        = grout%sb%xsize
      grout%sb%lsize(2:grout%sb%dim) = grout%sb%rsize
    case(PARALLELEPIPED, MINIMUM)
      grout%sb%box_shape = SPHERE
      grout%sb%rsize = sqrt( sum(grout%sb%lsize(:)**2) )
      grout%sb%lsize(1:grout%sb%dim) = grout%sb%rsize
    case default
      write(message(1),'(a)') 'grid_create_largergrid: Internal octopus error -- unsupported box shape for this routine.'
      call messages_fatal(1)
    end select

    call stencil_copy(grin%stencil, grout%stencil)
    call grid_init_stage_1(grout, geo)
    call mesh_init_stage_3(grout%mesh, grout%stencil, mc)
    call derivatives_build(grout%der, grout%mesh)

    POP_SUB(grid_create_largergrid)
  end subroutine grid_create_largergrid

end module grid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
