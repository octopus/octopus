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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, gr, td)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    call push_sub('opt_control_parameters.parameters_init')

    ALLOCATE(cp%laser(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter))
    cp%laser(:, :) = M_ZERO

  end subroutine parameters_init


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, gr, td, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    type(epot_t), intent(in) :: ep

    integer :: jj, kk
    FLOAT   :: t
    FLOAT, allocatable :: field(:)

    call push_sub('opt_control_parameters.parameters_set')

    ALLOCATE(field(1:NDIM), NDIM)
    cp%laser = M_ZERO;

      do jj = 0, 2*td%max_iter
        t = td%dt*(jj-1)/M_TWO
        !i = int(abs(M_TWO*t/l(1)%dt) + M_HALF)
        do kk=1, ep%no_lasers
          call laser_field(gr%sb, ep%no_lasers, ep%lasers, t, field)
          cp%laser(:,jj) = cp%laser(:,jj) + field(:)
        end do
      enddo

    deallocate(field)
    call pop_sub()
  end subroutine parameters_set


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, gr, td, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    type(epot_t), intent(inout) :: ep

    call push_sub('opt_control_paramters.parameters_to_h')

    ! TODO Probably is this better if this was not a pointer?
    ep%lasers(1)%numerical => cp%laser

    call pop_sub()
  end subroutine parameters_to_h


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    call push_sub('opt_control_parameters.parameters_end')

    deallocate(cp%laser)
    nullify(cp%laser)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, steps, n_dim, las, dt)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: steps
    integer,          intent(in) :: n_dim
    FLOAT,            intent(in) :: las(1:n_dim,0:2*steps)
    FLOAT,            intent(in) :: dt
    integer :: i, iunit
    FLOAT   :: tgrid(0:2*steps)

    call push_sub('opt_control_parameters.parameters_write')
    
    call t_lookup(2*steps+1,dt/CNST(2.0),tgrid)
    iunit = io_open(filename, action='write')
    do i = 0, 2*steps, 2
       write(iunit, '(4es30.16e4)') tgrid(i), las(:, i)
    end do
    call io_close(iunit)

    call pop_sub()
  end subroutine parameters_write




  
