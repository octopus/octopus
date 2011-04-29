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
  use io_m
  use splines_m
  use loct_math_m
  use messages_m
  use c_pointer_m

  implicit none

  private

  public ::               &
    spline_filter_ft,     &
    spline_filter_bessel, &  
    spline_filter_mask,   &
    spline_filter_mask_init

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
  !----------------------------------------------------------------------------
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

  
  !----------------------------------------------------------------------------
  subroutine spline_filter_bessel(spl, l, qmax, alpha, beta_fs, rcut, beta_rs)
    type(spline_t), intent(inout) :: spl
    integer, intent(in) :: l
    FLOAT, intent(in)   :: qmax, alpha, beta_fs, rcut, beta_rs

    type(spline_t) :: splw

    call spline_init(splw)
    call spline_besselft(spl, splw, l, M_FOUR*qmax)
    call spline_cut(splw, alpha*qmax, beta_fs)
    call spline_besselft(splw, spl, l)
    call spline_end(splw)
    call spline_cut(spl, rcut, beta_rs)

  end subroutine spline_filter_bessel


  !----------------------------------------------------------------------------
  subroutine spline_filter_mask_init()
    integer :: iunit, i

    PUSH_SUB(spline_filter_mask_init)

    iunit = io_open(trim(conf%share)//"/filter_mask.data", action='read', status='old', die=.true.)

    do i = 1, mask_n
      read(iunit, *) mask_x(i), mask_y(i)
    end do
    
    call io_close(iunit)

    POP_SUB(spline_filter_mask_init)
  end subroutine spline_filter_mask_init


  !----------------------------------------------------------------------------
  subroutine spline_filter_mask(spl, l, rmax, qmax, alpha, gamma)
    type(spline_t), intent(inout) :: spl
    integer,        intent(in)    :: l
    FLOAT,          intent(in)    :: rmax
    FLOAT,          intent(in)    :: qmax
    FLOAT,          intent(in)    :: alpha
    FLOAT,          intent(in)    :: gamma

    type(spline_t) :: mask, splw
    FLOAT :: local_mask_x(mask_n), rcut, beta

    PUSH_SUB(spline_filter_mask)

    rcut = gamma*rmax

    ! we define the mask function as f(r/rcut)
    local_mask_x = mask_x*rcut
    call spline_fit(mask_n, local_mask_x, mask_y, mask)

    ! apply the mask function
    call spline_div(spl, mask)

    ! transform to bessel space
    call spline_init(splw)
    call spline_besselft(spl, splw, l, M_FOUR*qmax)

    ! apply a cutoff
    beta = log(CNST(1e5))/(alpha - 1)**2
    call spline_cut(splw, qmax/alpha, beta)

    ! transform back to real space
    call spline_besselft(splw, spl, l)
    call spline_end(splw)

    ! remove the mask function
    call spline_mult(spl, mask)

    call spline_end(mask)

    POP_SUB(spline_filter_mask)

  end subroutine spline_filter_mask

end module spline_filter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
