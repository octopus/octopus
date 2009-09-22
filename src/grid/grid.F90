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
  use curvilinear_m
  use datasets_m
  use derivatives_m
  use double_grid_m
  use geometry_m
  use global_m
  use ob_interface_m
  use mesh_m
  use mesh_init_m
  use messages_m
  use mpi_m
  use multicomm_m
  use multigrid_m
  use nl_operator_m
  use loct_parser_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use stencil_cube_m

  implicit none

  private
  public ::                &
    grid_t,                &
    grid_init_stage_1,     &
    grid_init_stage_2,     &
    grid_end,              &
    grid_write_info,       &
    grid_create_multigrid, &
    grid_create_largergrid

  type grid_t
    type(simul_box_t)           :: sb
    type(mesh_t)                :: mesh
    type(interface_t), pointer  :: intf(:)
    type(multigrid_level_t)     :: fine
    type(derivatives_t)         :: der
    type(curvilinear_t)         :: cv
    type(multigrid_t), pointer  :: mgrid
    type(double_grid_t)         :: dgrid
    logical                     :: have_fine_mesh
    type(stencil_t)             :: stencil
  end type grid_t


contains

  !-------------------------------------------------------------------
  subroutine grid_init_stage_1(gr, geo)
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    type(stencil_t) :: cube

    call push_sub('grid.grid_init_stage_1')

    !%Variable UseFineMesh
    !%Type logical
    !%Default no
    !%Section Mesh
    !%Description
    !% If enabled, octopus will use a finer mesh for the calculation
    !% of the forces or other sensitive quantities. The default is
    !% disable.
    !%End
    if (gr%sb%dim == 3) then 
      call loct_parse_logical(datasets_check('UseFineMesh'), .false., gr%have_fine_mesh)
    else
      gr%have_fine_mesh = .false.
    end if

    ! initialize curvilinear coordinates
    call curvilinear_init(gr%sb, geo, gr%cv)

    ! initilize derivatives
    call derivatives_init(gr%der, gr%sb, gr%cv%method /= CURV_METHOD_UNIFORM)

    call double_grid_init(gr%dgrid, gr%sb)

    ! now we generate the mesh and the derivatives
    call mesh_init_stage_1(gr%mesh, gr%sb, gr%cv, &
         enlarge = max(gr%der%n_ghost, double_grid_enlarge(gr%dgrid)))

    ! the stencil used to generate the grid is a union of a cube (for
    ! multigrid) and the laplacian.
    call stencil_cube_get_lapl(cube, gr%sb%dim, order = 2)
    call stencil_union(gr%sb%dim, cube, gr%der%lapl%stencil, gr%stencil)
    call stencil_end(cube)
    
    call mesh_init_stage_2(gr%mesh, gr%sb, geo, gr%cv, gr%stencil)

    call pop_sub()

  end subroutine grid_init_stage_1


  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, mc, geo)
    type(grid_t), target, intent(inout) :: gr
    type(multicomm_t),    intent(in)    :: mc
    type(geometry_t),     intent(in)    :: geo

    integer :: il
    type(mpi_grp_t) :: grp

    call push_sub('grid.grid_init_stage_2')

    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)) then
      call mpi_grp_init(grp, mc%group_comm(P_STRATEGY_DOMAINS))
      call mesh_init_stage_3(gr%mesh, gr%stencil, grp)
    else
      call mesh_init_stage_3(gr%mesh)
    end if

    if(gr%sb%open_boundaries) then
      SAFE_ALLOCATE(gr%intf(1:NLEADS))
      do il = 1, NLEADS
        call interface_init(gr%mesh, gr%sb, gr%der, gr%intf(il), il)
      end do
    else
      nullify(gr%intf)
    end if

    call nl_operator_global_init()
    call derivatives_build(gr%der, gr%mesh)

    ! initialize a finer mesh to hold the density, for this we use the
    ! multigrid routines
    
    if(gr%have_fine_mesh) then
      
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
      
      call mesh_write_info(gr%fine%mesh, stdout)

    else
      gr%fine%mesh => gr%mesh
      gr%fine%der => gr%der
    end if

    ! multigrids are not initialized by default
    nullify(gr%mgrid)

    ! print info concerning the grid
    call grid_write_info(gr, geo, stdout)

    call pop_sub()
  end subroutine grid_init_stage_2


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    integer :: il

    call push_sub('grid.grid_end')

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
    call mesh_end(gr%mesh)

    if(associated(gr%mgrid)) then
      call multigrid_end(gr%mgrid)
      SAFE_DEALLOCATE_P(gr%mgrid)
    end if

    if(gr%sb%open_boundaries) then
      do il=1, NLEADS
        call interface_end(gr%intf(il))
      end do
      SAFE_DEALLOCATE_P(gr%intf)
    end if

    call stencil_end(gr%stencil)

    call pop_sub()
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, geo, iunit)
    type(grid_t),     intent(in) :: gr
    type(geometry_t), intent(in) :: geo
    integer,          intent(in) :: iunit

    integer :: il

    call push_sub('grid.grid_write_info')

    if(.not.mpi_grp_is_root(mpi_world)) then
      if(in_debug_mode) call write_debug_newlines(6)
      call pop_sub(); return
    end if

    call messages_print_stress(iunit, "Grid")
    call simul_box_write_info(gr%sb, geo, iunit)
    call mesh_write_info(gr%mesh, iunit)
    if(gr%sb%open_boundaries) then
      do il = 1, NLEADS
        call interface_write_info(gr%intf(il), iunit)
      end do
    end if
    if (gr%mesh%use_curvilinear) then
      call curvilinear_write_info(gr%cv, iunit)
    end if
    call messages_print_stress(iunit)

    call pop_sub()
  end subroutine grid_write_info


  !-------------------------------------------------------------------
  subroutine grid_create_multigrid(gr, geo)
    type(grid_t), intent(inout) :: gr
    type(geometry_t), intent(in) :: geo

    call push_sub('grid.grid_create_multigrid')

    SAFE_ALLOCATE(gr%mgrid)
    call multigrid_init(gr%mgrid, geo, gr%cv, gr%mesh, gr%der, gr%stencil)

    call pop_sub()
  end subroutine grid_create_multigrid


  !-------------------------------------------------------------------
  subroutine grid_create_largergrid(grin, geo, grout)
    type(grid_t), intent(in)     :: grin
    type(geometry_t), intent(in) :: geo
    type(grid_t), intent(out)    :: grout

    FLOAT :: avg

    call push_sub('grid.grid_create_largergrid')

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
      call write_fatal(1)
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

    call pop_sub()
  end subroutine grid_create_largergrid

end module grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
