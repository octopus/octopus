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
!! $Id: nl_operator.F90 2793 2007-03-26 00:20:03Z xavier $

#include "global.h"

module operate_m
  implicit none

  private
  public ::   &
    doperate, &
    zoperate

contains

  ! ---------------------------------------------------------
  subroutine doperate(np, nn, w, opi, fi, fo)
    integer, intent(in) :: np
    integer, intent(in) :: nn
    FLOAT,   intent(in) :: w(1:nn)
    integer, intent(in) :: opi(1:nn, 1:np)
    FLOAT,   intent(in) :: fi(1:np)
    FLOAT,   intent(out):: fo(1:np) 

    integer :: ii

    do ii = 1, np
      fo(ii) = sum(w(1:nn)  * fi(opi(1:nn, ii)))
    end do
  end subroutine doperate


  ! ---------------------------------------------------------
  subroutine zoperate(np, nn, w, opi, fi, fo)
    integer, intent(in) :: np
    integer, intent(in) :: nn
    FLOAT,   intent(in) :: w(1:nn)
    integer, intent(in) :: opi(1:nn, 1:np)
    CMPLX,   intent(in) :: fi(1:np)
    CMPLX,   intent(out):: fo(1:np) 

    integer :: ii

    do ii = 1, np
      fo(ii) = sum(w(1:nn)  * fi(opi(1:nn, ii)))
    end do
  end subroutine zoperate

end module operate_m
