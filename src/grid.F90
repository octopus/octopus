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

  implicit none
  
  type grid_type
    type(simul_box_type)       :: sb
    type(mesh_type),  pointer  :: m
    type(f_der_type), pointer  :: f_der

    type(curvlinear_type) :: cv

    integer                    :: n_multi
    type(mesh_type),  pointer  :: multi_m(:)
    type(f_der_type), pointer  :: multi_f_der(:)
  end type grid_type

contains

  !-------------------------------------------------------------------
  subroutine grid_init(gr, geo)
    type(grid_type),     intent(out) :: gr
    type(geometry_type), intent(in)  :: geo

    call push_sub('grid_init')

    ! this will count the number of grids contained in our multigrid
    gr%n_multi = 0

    allocate(gr%m, gr%f_der)
    nullify(gr%multi_m, gr%multi_f_der)

    ! initialize simulation box
    call simul_box_init(gr%sb, geo)

    ! initialize curvlinear coordinates
    call curvlinear_init(gr%sb%lsize(:), gr%cv)

    ! initilize derivatives
    call f_der_init(gr%f_der, gr%cv%method.ne.CURV_METHOD_UNIFORM)

    ! now we generate create the mesh and the derivatives
    call mesh_init(gr%m, gr%sb, gr%cv, geo, gr%f_der%n_ghost(1))
    call f_der_build(gr%f_der, gr%m)

    call grid_write_info(gr, stdout)

    call pop_sub()
  end subroutine grid_init


  !-------------------------------------------------------------------
  subroutine grid_end(gr)
    type(grid_type), intent(inout) :: gr

    integer :: i

    call push_sub('grid_end')

    ASSERT(associated(gr%m))
    ASSERT(gr%n_multi >= 0)

    call f_der_end(gr%f_der)
    call mesh_end(gr%m)

    deallocate(gr%m, gr%f_der)
    nullify   (gr%m, gr%f_der)

    if(gr%n_multi > 0) then
      do i = 1, gr%n_multi
        call f_der_end( gr%multi_f_der(i))
        call mesh_end(gr%multi_m(i))
      end do

      deallocate(gr%multi_m, gr%multi_f_der)
      nullify(gr%multi_m, gr%multi_f_der)
    end if

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
end module grid
