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
    pblas_gemm,        &
    pblas_herk,        &
    pblas_trsm
  
  interface pblas_gemm
    subroutine pdgemm ( transa, transb, m, n, k, alpha, & 
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc )
      implicit none

      character*1      :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      real(8)          :: alpha, beta
      real(8)          :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pdgemm

    subroutine pzgemm ( transa, transb, m, n, k, alpha, & 
      a, ia, ja, desca, b, ib, jb, descb, &
      beta, c, ic, jc, descc )
      implicit none

      character*1      :: transa, transb
      integer          :: m, n, k, ia, ja, ib, jb, ic, jc
      complex(8)       :: alpha, beta
      complex(8)       :: a, b, c
      integer          :: desca, descb, descc
    end subroutine pzgemm
  end interface

  interface pblas_herk
    subroutine pdsyrk(uplo, trans, n, k, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
      implicit none
      
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      real(8),      intent(in)    :: alpha
      real(8),      intent(in)    :: a    
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(8),      intent(in)    :: beta
      real(8),      intent(inout) :: c
      integer,      intent(in)    :: ic
      integer,      intent(in)    :: jc
      integer,      intent(in)    :: descc
    end subroutine pdsyrk

    subroutine pzherk(uplo, trans, n, k, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
      implicit none

      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: trans
      integer,      intent(in)    :: n
      integer,      intent(in)    :: k
      complex(8),   intent(in)    :: alpha
      complex(8),   intent(in)    :: a    
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      complex(8),   intent(in)    :: beta
      complex(8),   intent(inout) :: c
      integer,      intent(in)    :: ic
      integer,      intent(in)    :: jc
      integer,      intent(in)    :: descc
    end subroutine pzherk
 
  end interface pblas_herk

  interface pblas_trsm
    subroutine pdtrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      implicit none

      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      real(8),      intent(in)    :: alpha
      real(8),      intent(in)    :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      real(8),      intent(inout) :: b
      integer,      intent(in)    :: ib
      integer,      intent(in)    :: jb
      integer,      intent(in)    :: descb
    end subroutine pdtrsm

    subroutine pztrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      implicit none

      character(1), intent(in)    :: side
      character(1), intent(in)    :: uplo
      character(1), intent(in)    :: transa
      character(1), intent(in)    :: diag
      integer,      intent(in)    :: m
      integer,      intent(in)    :: n
      complex(8),   intent(in)    :: alpha
      complex(8),   intent(in)    :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      complex(8),   intent(inout) :: b
      integer,      intent(in)    :: ib
      integer,      intent(in)    :: jb
      integer,      intent(in)    :: descb
    end subroutine pztrsm

  end interface pblas_trsm

end module pblas_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
