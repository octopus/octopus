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

module spline_filter_m
  use global_m
  use splines_m
  use loct_math_m
  use messages_m
  use c_pointer_m

  implicit none

  private

  public ::               &
    spline_filter_ft,     &
    spline_filter_bessel, &  
    filter_mask       

  integer, parameter :: mask_n = 201
  FLOAT :: mask_x(mask_n), mask_y(mask_n)

contains

  !     The function spline_filter_ft permits to filter out
  !     high-values of a given spline function, either in real or in
  !     Fourier space.  If the optional argument fs is supplied, a
  !     filter in Fourier space will be done by moving to the Fourier
  !     representation and then calling spline_cut with cutoff = fs(1)
  !     and beta = fs(2). If the optional argument rs is supplied, a
  !     filter in real space will be done by calling spline_cut with
  !     cutoff = rs(1) and beta = rs(2). If both arguments are
  !     supplied, the Fourier filter will be applied *before*.
  !

  subroutine spline_filter_ft(spl, fs, rs)
    type(spline_t),      intent(inout) :: spl
    FLOAT, optional,     intent(in)    :: fs(2)
    FLOAT, optional,     intent(in)    :: rs(2)

    type(spline_t) :: splw

    if(present(fs)) then
      call spline_init(splw)
      call spline_3dft(spl, splw, CNST(2.0)*fs(1))
      call spline_cut(splw, fs(1), fs(2))
      call spline_3dft(splw, spl)
      call spline_times(TOFLOAT(CNST(1.0)/(CNST(2.0)*M_PI)**3), spl)
      call spline_end(splw)
    end if
    if(present(rs)) then
      call spline_cut(spl, rs(1), rs(2))
    end if

  end subroutine spline_filter_ft

  subroutine spline_filter_bessel(spl, l, fs, rs)
    type(spline_t), intent(inout) :: spl
    integer, intent(in) :: l
    FLOAT, intent(in), optional :: fs(2)
    FLOAT, intent(in), optional :: rs(2)

    type(spline_t) :: splw

    if(present(fs)) then
      call spline_init(splw)
      call spline_besselft(spl, splw, l, CNST(2.0)*fs(1))
      call spline_cut(splw, fs(1), fs(2))
      call spline_besselft(splw, spl, l)
      call spline_end(splw)
    end if
    if(present(rs)) then
      call spline_cut(spl, rs(1), rs(2))
    end if

  end subroutine spline_filter_bessel

  subroutine filter_mask(spl, l,rnlmax, qmax, alpha, gamma)
    type(spline_t), intent(inout) :: spl
    integer,        intent(in)    :: l
    FLOAT,          intent(in)    :: rnlmax
    FLOAT,          intent(in)    :: qmax
    FLOAT,          intent(in)    :: alpha
    FLOAT,          intent(in)    :: gamma

    type(spline_t) :: mask, splw
    FLOAT :: rcut, beta

    call push_sub('filters.filter_mask')

    call mask_function_init()

    rcut = gamma*rnlmax

    ! we define the mask function as f(r/rcut)
    mask_x = mask_x/rcut
    call spline_fit(mask_n, mask_x, mask_y, mask)

    ! apply the mask function
    call spline_div(spl, mask)

    ! transform to bessel space
    call spline_init(splw)
    call spline_besselft(spl, splw, l, M_TWO*qmax)

    ! apply a cutoff
    beta = log(M_ONE/CNST(1e-5))/(M_ONE/alpha - 1)**2
    call spline_cut(spl, qmax/alpha, beta)

    ! transform back to real space
    call spline_besselft(splw, spl, l)
    call spline_end(splw)

    ! remove the mask function
    call spline_mult(spl, mask)

    call pop_sub()

  end subroutine filter_mask

#include "filters_data.F90"

end module spline_filter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
