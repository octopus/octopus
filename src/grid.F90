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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module grid_m
  use global_m
  use lib_oct_parser_m
  use messages_m
  use datasets_m
  use mesh_m
  use simul_box_m
  use geometry_m
  use functions_m
  use curvlinear_m
  use multigrid_m
  use par_vec_m
  use mpi_m
  use multicomm_m

  implicit none

  private
  public ::                &
    grid_t,                &
    grid_init_stage_1,     &
    grid_init_stage_2,     &
    grid_end,              &
    grid_write_info,       &
    grid_create_multigrid

  type grid_t
    type(simul_box_t)           :: sb
    type(mesh_t)                :: m
    type(f_der_t)               :: f_der
    type(curvlinear_t)          :: cv
    type(multigrid_t), pointer  :: mgrid
  end type grid_t


contains

  !-------------------------------------------------------------------
  subroutine grid_init_stage_1(gr, geo)
    type(grid_t),      intent(inout) :: gr
    type(geometry_t),  intent(in)    :: geo

    call push_sub('grid.grid_init_stage_1')

    ! initialize curvlinear coordinates
    call curvlinear_init(gr%sb, gr%cv)

    ! initilize derivatives
    call f_der_init(gr%f_der, gr%sb, gr%cv%method.ne.CURV_METHOD_UNIFORM)

    ! now we generate create the mesh and the derivatives
    call mesh_init_stage_1(gr%sb, gr%m, geo, gr%cv, gr%f_der%n_ghost)
    call mesh_init_stage_2(gr%sb, gr%m, geo, gr%cv)

    call pop_sub()
  end subroutine grid_init_stage_1


  !-------------------------------------------------------------------
  subroutine grid_init_stage_2(gr, mc, geo)
    type(grid_t),      intent(inout) :: gr
    type(multicomm_t), intent(in)    :: mc
    type(geometry_t),  intent(in)    :: geo

    call push_sub('grid.grid_init_stage_2')

    call mesh_init_stage_3(gr%m, geo, gr%cv,  &
      gr%f_der%der_discr%lapl%stencil,        &
      gr%f_der%der_discr%lapl%n, mc)

    call f_der_build(gr%sb, gr%m, gr%f_der)

    ! multigrids are not initialized by default
    nullify(gr%mgrid)

    ! print info concerning the grid
    call grid_write_info(gr, stdout)

    call pop_sub()
  end subroutine grid_init_stage_2


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_t), intent(inout) :: gr

    call push_sub('grid.grid_end')

    call f_der_end(gr%f_der)
    call mesh_end(gr%m)

    if(associated(gr%mgrid)) then
      call multigrid_end(gr%mgrid)
      deallocate(gr%mgrid); nullify(gr%mgrid)
    end if

    call pop_sub()
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, iunit)
    type(grid_t), intent(in) :: gr
    integer,      intent(in) :: iunit

    if(.not.mpi_grp_is_root(mpi_world)) then
      if(in_debug_mode) call write_debug_newlines(6)
      return
    end if

    call messages_print_stress(iunit, "Grid")
    call simul_box_write_info(gr%sb, iunit)
    call mesh_write_info(gr%m, iunit)
    if (gr%m%use_curvlinear) then
      call curvlinear_write_info(gr%cv, iunit)
    end if
    call messages_print_stress(iunit)

  end subroutine grid_write_info


  !-------------------------------------------------------------------
  subroutine grid_create_multigrid(gr, geo)
    type(grid_t), intent(inout) :: gr
    type(geometry_t), intent(in) :: geo

    ALLOCATE(gr%mgrid, 1)
    call multigrid_init(geo, gr%cv, gr%m, gr%f_der, gr%mgrid)

  end subroutine grid_create_multigrid

end module grid_m
