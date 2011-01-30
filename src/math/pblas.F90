!! Copyright (C) 2011 X. Andrade
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
!! $Id: lapack.F90 5550 2009-06-03 20:53:01Z xavier $

#include "global.h"

module pblas_m
  implicit none

  private 

  public ::            &
    pblas_gemm

  interface pblas_gemm
    subroutine pdgemm ( transa, transb, m, n, k, alpha, & 
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc )
      
      character*1      :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      real(8)          :: alpha, beta
      real(8)          :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pdgemm

    subroutine pzgemm ( transa, transb, m, n, k, alpha, & 
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc )
      
      character*1      :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      complex(8)       :: alpha, beta
      complex(8)       :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pzgemm
  end interface

end module pblas_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
