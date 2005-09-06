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

  implicit none

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
  subroutine grid_init(gr)
    type(grid_type),     intent(out) :: gr

    logical :: filter

    call push_sub('grid.grid_init')

    ! initilize geometry
    call geometry_init_xyz(gr%geo)
    call geometry_init_species(gr%geo)
    call geometry_init_vel(gr%geo)

    ! initialize simulation box
    call simul_box_init(gr%sb, gr%geo)

    ! initialize curvlinear coordinates
    call curvlinear_init(gr%sb, gr%cv)

    ! initilize derivatives
    call f_der_init(gr%f_der, gr%sb, gr%cv%method.ne.CURV_METHOD_UNIFORM)

    ! now we generate create the mesh and the derivatives
    call mesh_init(gr%sb, gr%m, gr%geo, gr%cv, gr%f_der%n_ghost)
#ifdef HAVE_MPI
    ! Initialize parallel vectors. The points of the stencil are
    ! picked out of the laplacian.
    call vec_init_default(gr%m, gr%f_der%der_discr%lapl%stencil, &
                          gr%f_der%der_discr%lapl%n, gr%m%vp)
    ! Fill in local point numbers.
    call mesh_par_adj(gr%m)
#endif
    call f_der_build(gr%sb, gr%m, gr%f_der)

    ! do we want to filter out the external potentials, or not.
    call loct_parse_logical(check_inp('FilterPotentials'), .false., filter)
    if(filter) call geometry_filter(gr%geo, mesh_gcutoff(gr%m))

    ! Now that we are really done with initializing the geometry, print debugging information.
    if(conf%verbose>=VERBOSE_DEBUG) then
      call geometry_debug(gr%geo, 'debug')
    end if

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
#ifdef HAVE_MPI
    call vec_end(gr%m%vp)
#endif
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


#ifdef HAVE_MPI
    if(mpiv%node .ne. 0) return
#endif
    if(iunit==stdout.and.conf%verbose<VERBOSE_NORMAL) return

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
