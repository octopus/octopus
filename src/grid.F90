!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module grid
  use mesh
  use simul_box
  use geometry
  use functions
  use curvlinear
  use multigrid
  use par_vec
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  use mpi_mod
#endif

  implicit none

  private
  public :: grid_type, &
            grid_init, &
            grid_end,  &
            grid_write_info, &
            grid_create_multigrid

  type grid_type
    type(simul_box_type)       :: sb
    type(geometry_type)        :: geo
    type(mesh_type)            :: m
    type(f_der_type)           :: f_der

    type(curvlinear_type) :: cv
    type(multigrid_type), pointer  :: mgrid

  end type grid_type

contains

  !-------------------------------------------------------------------
  subroutine grid_init(gr, domain_comm_of_node)
    type(grid_type),      intent(out) :: gr
    integer, intent(in) :: domain_comm_of_node(0:mpiv%numprocs-1)

    logical :: filter
    integer :: comm

    call push_sub('grid.grid_init')

    ! initialize curvlinear coordinates
    call curvlinear_init(gr%sb, gr%cv)

    ! initilize derivatives
    call f_der_init(gr%f_der, gr%sb, gr%cv%method.ne.CURV_METHOD_UNIFORM)

    ! now we generate create the mesh and the derivatives
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    comm = MPI_COMM_WORLD
#endif
    call mesh_init(gr%sb, gr%m, gr%geo, gr%cv, gr%f_der%n_ghost, &
         gr%f_der%der_discr%lapl%stencil,                        &
         gr%f_der%der_discr%lapl%n, comm)

    call f_der_build(gr%sb, gr%m, gr%f_der)

    ! do we want to filter out the external potentials, or not.
    call loct_parse_logical(check_inp('FilterPotentials'), .false., filter)
    if(filter) call geometry_filter(gr%geo, mesh_gcutoff(gr%m))

    ! Now that we are really done with initializing the geometry, print debugging information.
    if(in_debug_mode) call geometry_debug(gr%geo, 'debug')

    ! multigrids are not initialized by default
    nullify(gr%mgrid)

    ! print info concerning the grid
    call grid_write_info(gr, stdout)

    call pop_sub()
  end subroutine grid_init


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_type), intent(inout) :: gr

    call push_sub('grid.grid_end')

    call f_der_end(gr%f_der)
    call mesh_end(gr%m)

    if(associated(gr%mgrid)) then
      call multigrid_end(gr%mgrid)
      deallocate(gr%mgrid); nullify(gr%mgrid)
    end if

    call geometry_end(gr%geo)

    call pop_sub()
  end subroutine grid_end


  !-------------------------------------------------------------------
  subroutine grid_write_info(gr, iunit)
    type(grid_type), intent(in) :: gr
    integer,         intent(in) :: iunit

    if(mpiv%node .ne. 0) then
       if(in_debug_mode) call write_debug_newlines(4)
       return
    endif

    write(iunit,'(/,a)') stars
    call simul_box_write_info(gr%sb, iunit)
    call mesh_write_info(gr%m, iunit)
    write(iunit,'(a,/)') stars

  end subroutine grid_write_info


  subroutine grid_create_multigrid(gr)
    type(grid_type), intent(inout) :: gr

    allocate(gr%mgrid)
    call multigrid_init(gr%geo, gr%cv, gr%m, gr%f_der, gr%mgrid)

  end subroutine grid_create_multigrid

end module grid
