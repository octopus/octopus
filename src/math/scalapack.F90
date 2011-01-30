!! Copyright (C) 2011 D. Strubbe
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
! This module contains interfaces for ScaLAPACK and BLACS routines
! Interfaces are from http://www.netlib.org/scalapack/tools, double, complex16
! and from http://www.netlib.org/blacs/BLACS/QRef.html (entered manually...)
! -----------------------------------------------------------------------

module scalapack_m
  implicit none

  interface
    INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
      INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
    end FUNCTION NUMROC
  end interface

  interface
    INTEGER FUNCTION ICEIL( INUM, IDENOM )
      INTEGER            IDENOM, INUM
    end FUNCTION ICEIL
  end interface

  interface
    SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
    end SUBROUTINE DESCINIT
  end interface

  interface scalapack_geqrf
    SUBROUTINE PDGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO )
      INTEGER            IA, INFO, JA, LWORK, M, N
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A, TAU, WORK
    end SUBROUTINE PDGEQRF
    
    SUBROUTINE PZGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO )
      INTEGER            IA, INFO, JA, LWORK, M, N
      INTEGER            DESCA( * )
      COMPLEX*16         A, TAU, WORK
    end SUBROUTINE PZGEQRF
  end interface scalapack_geqrf

  interface scalapack_orgqr
    SUBROUTINE PDORGQR( M, N, K, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO ) 
      INTEGER            IA, INFO, JA, K, LWORK, M, N
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A, TAU, WORK
    end SUBROUTINE PDORGQR
    
    SUBROUTINE PZUNGQR( M, N, K, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO )
      INTEGER            IA, INFO, JA, K, LWORK, M, N
      INTEGER            DESCA( * )
      COMPLEX*16         A, TAU, WORK
    end SUBROUTINE PZUNGQR
  end interface scalapack_orgqr
  
  interface
    SUBROUTINE PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV
      DOUBLE PRECISION   A, B
    end SUBROUTINE PDGESV
  end interface

  interface
    SUBROUTINE PZGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV
      COMPLEX*16         A, B
    end SUBROUTINE PZGESV
  end interface

  interface
    SUBROUTINE PDSYEVX( JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, &
      VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL, &
      ICLUSTR, GAP, INFO )
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      DOUBLE PRECISION   ABSTOL, ORFAC, VL, VU
      INTEGER            DESCA( * ), DESCZ( * ), ICLUSTR, IFAIL, IWORK
      DOUBLE PRECISION   A, GAP, W, WORK, Z
    end SUBROUTINE PDSYEVX
  end interface

  interface
    subroutine pdsyev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info)
      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      real(8),          intent(inout) :: a
      integer,          intent(in)    :: ia
      integer,          intent(in)    :: ja
      integer,          intent(in)    :: desca
      real(8),          intent(out)   :: w
      real(8),          intent(out)   :: z
      integer,          intent(in)    :: iz
      integer,          intent(in)    :: jz
      integer,          intent(in)    :: descz
      real(8),          intent(out)   :: work
      integer,          intent(in)    :: lwork
      integer,          intent(out)   :: info
    end subroutine pdsyev

    subroutine pzheev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, rwork, lrwork, info)
      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      complex(8),       intent(inout) :: a
      integer,          intent(in)    :: ia
      integer,          intent(in)    :: ja
      integer,          intent(in)    :: desca
      real(8),          intent(out)   :: w
      complex(8),       intent(out)   :: z
      integer,          intent(in)    :: iz
      integer,          intent(in)    :: jz
      integer,          intent(in)    :: descz
      complex(8),       intent(out)   :: work
      integer,          intent(in)    :: lwork
      complex(8),       intent(out)   :: rwork
      integer,          intent(in)    :: lrwork
      integer,          intent(out)   :: info
    end subroutine pzheev
  end interface

end module scalapack_m

!! local Variables:
!! mode: f90
!! coding: utf-8
!! End:
