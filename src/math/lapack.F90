!! Copyright (C) 2009 X. Andrade
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

! -----------------------------------------------------------------------
! This module contains interfaces for LAPACK routines
! -----------------------------------------------------------------------

module lapack_m
  implicit none

  interface lapack_potrf
    subroutine spotrf(uplo, n, a, lda, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(4),      intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine spotrf

    subroutine dpotrf(uplo, n, a, lda, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(8),      intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine dpotrf

    subroutine cpotrf(uplo, n, a, lda, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(4),   intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine cpotrf

    subroutine zpotrf(uplo, n, a, lda, info)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(8),   intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine zpotrf
  end interface

end module lapack_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
