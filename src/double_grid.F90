!! Copyright (C) 2002-2006 X. Andrade
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
!! $Id: specie.F90 2711 2007-02-13 17:36:18Z xavier $

#include "global.h"

module double_grid_m
  
  use global_m
  use mesh_m
  use messages_m
  use lib_oct_gsl_spline_m

  implicit none
  
  private
  
  public ::              &
       double_grid_t,    &
       double_grid_init, &
       double_grid_end,  &
       double_grid_apply

  type double_grid_t
     private

     FLOAT   :: h_fine(MAX_DIM)
     FLOAT   :: h_coarse(MAX_DIM)
     integer :: nn
     integer :: dim

  end type double_grid_t
  
contains
  
  subroutine double_grid_init(this, m, dim)
    type(double_grid_t), intent(out) :: this
    type(mesh_t),        intent(in)  :: m
    integer,             intent(in)  :: dim

    this%nn = 3
    this%h_fine(1:MAX_DIM) = m%h(1:MAX_DIM)/dble(this%nn)
    this%h_coarse(1:MAX_DIM) = m%h(1:MAX_DIM)
    this%dim = dim
    
  end subroutine double_grid_init

  subroutine double_grid_end(this)
    type(double_grid_t), intent(in) :: this
 
  end subroutine double_grid_end

  FLOAT function double_grid_apply(this, ps_spline, x_atom, x_grid) result(ww)
    type(double_grid_t), intent(in) :: this
    type(loct_spline_t), intent(in)  :: ps_spline
    FLOAT,               intent(in)  :: x_atom(1:MAX_DIM)
    FLOAT,               intent(in)  :: x_grid(1:MAX_DIM)
    
    FLOAT :: x_fine_grid(1:MAX_DIM), r
   
    integer :: kk, idim
    
    r = sqrt(sum( (x_grid(1:this%dim)-x_atom(1:this%dim))**2 ))
    ww = loct_splint(ps_spline, r)

  end function double_grid_apply

end module double_grid_m
