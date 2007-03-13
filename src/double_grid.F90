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

  use datasets_m
  use global_m
  use mesh_m
  use messages_m
  use lib_oct_gsl_spline_m
  use lib_oct_parser_m
  use specie_m

  implicit none
  
  private
  
  public ::              &
       double_grid_t,    &
       double_grid_init, &
       double_grid_end,  &
       double_grid_apply, &
       dg_add_localization_density

  type double_grid_t

     FLOAT   :: h_fine(MAX_DIM)
     FLOAT   :: h_coarse(MAX_DIM)
     integer :: nn
     integer :: dim
     integer :: loc_function
     type(loct_spline_t) :: rho_corr
     type(loct_spline_t) :: pot_corr

  end type double_grid_t
  
  integer, parameter :: &
       F_NONE = 0,      &
       F_GAUSSIAN = 1 
  
  
contains
  
  subroutine double_grid_init(this, m, dim)
    type(double_grid_t), intent(out) :: this
    type(mesh_t),        intent(in)  :: m
    integer,             intent(in)  :: dim

    this%nn = 3
    this%h_fine(1:MAX_DIM) = m%h(1:MAX_DIM)/dble(this%nn)
    this%h_coarse(1:MAX_DIM) = m%h(1:MAX_DIM)
    this%dim = dim
    
    !%Variable LocalizationDensity
    !%Type integer
    !%Default none
    !%Section Mesh
    !%Description
    !% When a localization density is used, the local part of the
    !% pseudopotential is separated in a localized part and a long
    !% range term. The long range term is calculated by solving the
    !% poisson equation of a localized charge density. This separation
    !% is required to apply to use the double grid technique or to
    !% apply a filter to the local part of the pseudopotential. This
    !% variable selects which kind of function is used for the
    !% localized charge density. By default this technique is not
    !% used.
    !%Option none 0 
    !% The potential is not separated. This is the default.
    !%Option gaussian 3
    !% A gaussian charge distribution is used, the correspoding
    !% potential is erf(r)/r.
    !%End

    call loct_parse_int(check_inp('LocalizationDensity'), F_GAUSSIAN, this%loc_function)

    call functions_init()

    contains

      subroutine functions_init()
        integer, parameter :: nn = 3000
        FLOAT,   parameter :: rmax = CNST(30.0), sigma = CNST(0.625)
        
        integer :: ii
        FLOAT :: dr, r(1:nn), rho(1:nn), pot(1:nn)

        dr = rmax/(nn-M_ONE)
        
        r(1) = M_ZERO
        pot(1) = M_TWO/(sqrt(M_TWO*M_PI)*sigma)
        do ii = 1, nn

          rho(ii) = M_ONE/(sigma*sqrt(M_TWO*M_PI))**3*exp(-M_HALF*r(ii)**2/sigma**2)

          if ( ii /= 1 ) pot(ii) = erf(r(ii)/(sigma*sqrt(M_TWO)))/r(ii)

          if ( ii /= nn ) r(ii+1) = r(ii) + dr
        end do

        call loct_spline_init(this%rho_corr)
        call loct_spline_init(this%pot_corr)
        
        call loct_spline_fit(nn, r, rho, this%rho_corr)
        call loct_spline_fit(nn, r, pot, this%pot_corr)

      end subroutine functions_init

  end subroutine double_grid_init

  subroutine double_grid_end(this)
    type(double_grid_t), intent(inout) :: this
 
        call loct_spline_end(this%rho_corr)
        call loct_spline_end(this%pot_corr)

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

  logical function dg_add_localization_density(this) 
    type(double_grid_t), intent(in) :: this
    
    dg_add_localization_density = ( this%loc_function /= F_NONE )
    
  end function dg_add_localization_density

end module double_grid_m
