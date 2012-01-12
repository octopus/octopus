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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module grid_m
  use cube_m
  use curvilinear_m
  use datasets_m
  use derivatives_m
  use double_grid_m
  use geometry_m
  use global_m
  use ob_interface_m
  use ob_grid_m
  use mesh_m
  use mesh_init_m
  use messages_m
  use mpi_m
  use multicomm_m
  use multigrid_m
  use nl_operator_m
  use parser_m
  use profiling_m
  use space_m
  use simul_box_m
  use stencil_m
  use stencil_cube_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                &
    grid_t,                &
    grid_init_stage_0,     &
    grid_init_stage_1,     &
    grid_init_stage_2,     &
    grid_end,              &
    grid_write_info,       &
    grid_create_multigrid, &
    grid_create_largergrid

  type grid_t
    type(simul_box_t)           :: sb
    type(mesh_t)                :: mesh
    type(interface_t)           :: intf(2*MAX_DIM)
    type(multigrid_level_t)     :: fine
    type(derivatives_t)         :: der
    type(curvilinear_t)         :: cv
    type(multigrid_t), pointer  :: mgrid
    type(ob_grid_t)             :: ob_grid
    type(double_grid_t)         :: dgrid
    logical                     :: have_fine_mesh
    type(stencil_t)             :: stencil
  end type grid_t


contains

  !-------------------------------------------------------------------
  !>
  !! "Zero-th" stage of grid initialization. It initializes the open boundaries stuff
  !! (if necessary), and the simulation box.
  subroutine grid_init_stage_0(gr, geo, space)
    type(grid_t),          intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(space_t),         intent(in)    :: space

    PUSH_SUB(grid_init_stage_0)

    call ob_grid_init(gr%ob_grid, geo)
    if(gr%ob_grid%open_boundaries) then
      call simul_box_init(gr%sb, geo, space, gr%ob_grid%transport_mode, &
        gr%ob_grid%lead(:)%sb, gr%ob_grid%lead(:)%info)
    else
      call simul_box_init(gr%sb, geo, space)
    end if

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

    !%Variable UseFineMesh
    !%Type logical
    !%Default no
    !%Section Mesh
    !%Description
    !% If enabled, <tt>Octopus</tt> will use a finer mesh for the calculation
    !% of the forces or other sensitive quantities. The default is no.
    !%End
    if (gr%sb%dim == 3) then 
      call parse_logical(datasets_check('UseFineMesh'), .false., gr%have_fine_mesh)
    else
      gr%have_fine_mesh = .false.
    end if

    if(gr%have_fine_mesh) call messages_experimental("UseFineMesh")

    call geometry_grid_defaults(geo, def_h, def_rsize)

    ! initialize to -1
    grid_spacing = -M_ONE

#if defined(HAVE_GDLIB)
    if(gr%sb%box_shape == BOX_IMAGE) then 
      ! grid_spacing is determined from lsize and the size of the image
      grid_spacing(1:2) = M_TWO*gr%sb%lsize(1:2)/real(gr%sb%image_size(1:2), REAL_PRECISION)
    else
#endif

    !%Variable Spacing
    !%Type float
    !%Section Mesh
    !%Description
    !% The spacing between the points in the mesh. If using curvilinear
    !% coordinates, this is a canonical spacing that will be changed locally by the
    !% transformation. In periodic directions, your spacing may be slightly larger than
    !% what you request here, since the box size must be an integer multiple of the spacing.
    !%
    !% It is possible to have a different spacing in each one of the Cartesian directions
    !% if we define <tt>Spacing</tt> as block of the form
    !%
    !% <tt>%Spacing
    !% <br>&nbsp;&nbsp;spacing_x | spacing_y | spacing_z
    !% <br>%</tt>
    !%End

    if(parse_block(datasets_check('Spacing'), blk) == 0) then
      if(parse_block_cols(blk,0) < gr%sb%dim) call input_error('Spacing')
      do idir = 1, gr%sb%dim
        call parse_block_float(blk, 0, idir - 1, grid_spacing(idir), units_inp%length)
      end do
      call parse_block_end(blk)
    else
      call parse_float(datasets_check('Spacing'), -M_ONE, grid_spacing(1), units_inp%length)
      grid_spacing(1:gr%sb%dim) = grid_spacing(1)
    end if

    do idir = 1, gr%sb%dim
      if(grid_spacing(idir) < M_ZERO) then
        if(def_h > M_ZERO .and. def_h < huge(def_h)) then
          grid_spacing(idir) = def_h
          write(message(1), '(a,i1,3a,f6.3)') "Info: Using default spacing(", idir, &
            ") [", trim(units_abbrev(units_out%length)), "] = ",                        &
            units_from_atomic(units_out%length, grid_spacing(idir))
          call messages_info(1)
        else
          message(1) = 'Either:'
          message(2) = "   *) variable 'Spacing' is not defined and"
          message(4) = "      I can't find a suitable default"
          message(3) = "   *) your input for 'Spacing' is negative"
          call messages_fatal(4)
        end if
      end if
      if(def_rsize > M_ZERO) call messages_check_def(grid_spacing(idir), def_rsize, 'Spacing')
    end do

#if defined(HAVE_GDLIB)
  end if
#endif

    ! initialize curvilinear coordinates
    call curvilinear_init(gr%cv, gr%sb, geo, grid_spacing)

    ! initialize derivatives
    call derivatives_init(gr%der, gr%sb, gr%cv%method /= CURV_METHOD_UNIFORM)

    call double_grid_init(gr%dgrid, gr%sb)

    enlarge = 0
    enlarge(1:gr%sb%dim) = 2
    enlarge = max(enlarge, double_grid_enlarge(gr%dgrid))
    enlarge = max(enlarge, gr%der%n_ghost)

    ! now we generate the mesh and the derivatives
    call mesh_init_stage_1(gr%mesh, gr%sb, gr%cv, grid_spacing, enlarge, gr%ob_grid)

    if(gr%ob_grid%open_boundaries) call mesh_read_lead(gr%ob_grid, gr%mesh)

    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the Laplacian.
    call stencil_cube_get_lapl(cube, gr%sb%dim, order = 2)
    call stencil_union(gr%sb%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)

    call mesh_init_stage_2(gr%mesh, gr%sb, geo, gr%cv, gr%stencil)

    POP_SUB(grid_init_stage_1)

  end subroutine grid_init_stage_1


  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, mc, geo)
    type(grid_t), target, intent(inout) :: gr
    type(multicomm_t),    intent(in)    :: mc
    type(geometry_t),     intent(in)    :: geo

    integer :: il
    type(mpi_grp_t) :: grp

    PUSH_SUB(grid_init_stage_2)

    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)) then
      call mpi_grp_init(grp, mc%group_comm(P_STRATEGY_DOMAINS))
      call mesh_init_stage_3(gr%mesh, gr%stencil, grp)
    else
      call mesh_init_stage_3(gr%mesh)
    end if

    call nl_operator_global_init()
    call derivatives_build(gr%der, gr%mesh)

    ! we need the derivative for the interface, therefore do the initialization here
    if(gr%ob_grid%open_boundaries) then
      do il = 1, NLEADS
        call interface_init(gr%der, gr%intf(il), il, gr%ob_grid%lead(il)%sb%lsize)
      end do
    end if

    ! initialize a finer mesh to hold the density, for this we use the
    ! multigrid routines
    
    if(gr%have_fine_mesh) then

      if(gr%mesh%parallel_in_domains) then
        message(1) = 'UseFineMesh does not work with domain parallelization.'
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(gr%fine%mesh)
      SAFE_ALLOCATE(gr%fine%der)
      
      call multigrid_mesh_double(geo, gr%cv, gr%mesh, gr%fine%mesh, gr%stencil)
      
      call derivatives_init(gr%fine%der, gr%mesh%sb, gr%cv%method .ne. CURV_METHOD_UNIFORM)
      
      if(gr%mesh%parallel_in_domains) then
        call mesh_init_stage_3(gr%fine%mesh, gr%stencil, gr%mesh%mpi_grp)
      else
        call mesh_init_stage_3(gr%fine%mesh)
      end if
      
      call multigrid_get_transfer_tables(gr%fine%tt, gr%fine%mesh, gr%mesh)
      
      call derivatives_build(gr%fine%der, gr%fine%mesh)

      gr%fine%der%coarser => gr%der
      gr%der%finer =>  gr%fine%der
      gr%fine%der%to_coarser => gr%fine%tt
      gr%der%to_finer => gr%fine%tt

    else
      gr%fine%mesh => gr%mesh
      gr%fine%der => gr%der
    end if

    ! multigrids are not initialized by default
    nullify(gr%mgrid)

    ! print info concerning the grid
    call grid_write_info(gr, geo, stdout)
    if(gr%ob_grid%open_boundaries) call ob_grid_write_info(gr%ob_grid, stdout)

    POP_SUB(grid_init_stage_2)
  end subroutine grid_init_stage_2


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    integer :: il

    PUSH_SUB(grid_end)

    call nl_operator_global_end()

    if(gr%have_fine_mesh) then
      call derivatives_end(gr%fine%der)
      call mesh_end(gr%fine%mesh)
      SAFE_DEALLOCATE_P(gr%fine%mesh)
      SAFE_DEALLOCATE_P(gr%fine%der)
      SAFE_DEALLOCATE_P(gr%fine%tt%to_coarse)
      SAFE_DEALLOCATE_P(gr%fine%tt%to_fine1)
      SAFE_DEALLOCATE_P(gr%fine%tt%to_fine2)
      SAFE_DEALLOCATE_P(gr%fine%tt%to_fine4)
      SAFE_DEALLOCATE_P(gr%fine%tt%to_fine8)
      SAFE_DEALLOCATE_P(gr%fine%tt%fine_i)
    end if

    call double_grid_end(gr%dgrid)

    call derivatives_end(gr%der)
    call curvilinear_end(gr%cv)
    call mesh_end(gr%mesh)

    if(associated(gr%mgrid)) then
      call multigrid_end(gr%mgrid)
      SAFE_DEALLOCATE_P(gr%mgrid)
    end if

    if(gr%ob_grid%open_boundaries) then
      do il=1, NLEADS
        call interface_end(gr%intf(il))
        ! As usual with open boundaries code, this object is
        ! half-created so we cannot call mesh_end and we have to
        ! deallocate by hand. Please fix this.
        SAFE_DEALLOCATE_P(gr%ob_grid%lead(il)%mesh%idx%lxyz)
        SAFE_DEALLOCATE_P(gr%ob_grid%lead(il)%mesh%idx%lxyz_inv)        
      end do
    end if

    call ob_grid_end(gr%ob_grid)

    call stencil_end(gr%stencil)

    POP_SUB(grid_end)
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, geo, iunit)
    type(grid_t),     intent(in) :: gr
    type(geometry_t), intent(in) :: geo
    integer,          intent(in) :: iunit

    integer :: il

    PUSH_SUB(grid_write_info)

    if(.not.mpi_grp_is_root(mpi_world)) then
      if(in_debug_mode) call messages_debug_newlines(6)
      POP_SUB(grid_write_info)
      return
    end if

    call messages_print_stress(iunit, "Grid")
    call simul_box_write_info(gr%sb, geo, iunit)

    if(gr%have_fine_mesh) then
      message(1) = "Wave-functions mesh:"
      call messages_info(1, iunit)
      call mesh_write_info(gr%mesh, iunit)
      message(1) = "Density mesh:"
    else
      message(1) = "Main mesh:"
    end if
    call messages_info(1, iunit)
    call mesh_write_info(gr%fine%mesh, iunit)

    if(gr%ob_grid%open_boundaries) then
      do il = 1, NLEADS
        call interface_write_info(gr%intf(il), iunit)
      end do
    end if
    if (gr%mesh%use_curvilinear) then
      call curvilinear_write_info(gr%cv, iunit)
    end if
    call messages_print_stress(iunit)

    POP_SUB(grid_write_info)
  end subroutine grid_write_info


  !-------------------------------------------------------------------
  subroutine grid_create_multigrid(gr, geo)
    type(grid_t), intent(inout) :: gr
    type(geometry_t), intent(in) :: geo

    PUSH_SUB(grid_create_multigrid)

    SAFE_ALLOCATE(gr%mgrid)
    call multigrid_init(gr%mgrid, geo, gr%cv, gr%mesh, gr%der, gr%stencil)

    POP_SUB(grid_create_multigrid)
  end subroutine grid_create_multigrid


  !-------------------------------------------------------------------
  subroutine grid_create_largergrid(grin, geo, grout)
    type(grid_t), intent(in)     :: grin
    type(geometry_t), intent(in) :: geo
    type(grid_t), intent(out)    :: grout

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
    case(PARALLELEPIPED, MINIMUM, BOX_IMAGE, BOX_USDEF)
      grout%sb%box_shape = SPHERE
      grout%sb%rsize = sqrt( sum(grout%sb%lsize(:)**2) )
      grout%sb%lsize(1:grout%sb%dim) = grout%sb%rsize
    case default
      write(message(1),'(a)') 'Internal octopus error.'
      call messages_fatal(1)
    end select

    call stencil_copy(grin%stencil, grout%stencil)

    call grid_init_stage_1(grout, geo)

    if(grin%mesh%parallel_in_domains) then
      call mesh_init_stage_3(grout%mesh, grout%stencil, grin%mesh%mpi_grp)
    else
      call mesh_init_stage_3(grout%mesh)
    end if

    call derivatives_build(grout%der, grout%mesh)

    ! multigrids are not initialized by default
    nullify(grout%mgrid)

    POP_SUB(grid_create_largergrid)
  end subroutine grid_create_largergrid

end module grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
