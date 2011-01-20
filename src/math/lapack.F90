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

  interface lapack_sygv
    subroutine ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      real(4),      intent(inout) :: a, b    ! a(lda,n), b(ldb,n)
      real(4),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine ssygv

    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      real(8),      intent(inout) :: a, b    ! a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine dsygv
  end interface

  interface lapack_hegv
    subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      complex(4),   intent(inout) :: a, b     ! a(lda,n), b(ldb,n)
      real(4),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(4),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine chegv

    subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      complex(8),   intent(inout) :: a, b     ! a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine zhegv
  end interface

  interface lapack_geev
    subroutine sgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      real(4),      intent(inout) :: a ! a(lda,n)
      real(4),      intent(out)   :: w, vl, vr ! w(n), vl(ldvl,n), vl(ldvr,n)
      real(4),      intent(out)   :: rwork ! rwork(max(1,2n))
      real(4),      intent(out)   :: work  ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine sgeev

    subroutine dgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      real(8),      intent(inout) :: a ! a(lda,n)
      real(8),      intent(out)   :: w, vl, vr ! w(n), vl(ldvl,n), vl(ldvr,n)
      real(8),      intent(out)   :: rwork ! rwork(max(1,2n))
      real(8),      intent(out)   :: work  ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine dgeev

    subroutine cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      complex(4),   intent(inout) :: a ! a(lda,n)
      complex(4),   intent(out)   :: w, vl, vr ! w(n), vl(ldvl,n), vl(ldvr,n)
      real(4),      intent(out)   :: rwork ! rwork(max(1,2n))
      complex(4),   intent(out)   :: work  ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine cgeev

    subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      complex(8),        intent(inout) :: a ! a(lda,n)
      complex(8),        intent(out)   :: w, vl, vr ! w(n), vl(ldvl,n), vl(ldvr,n)
      real(8),        intent(out)   :: rwork ! rwork(max(1,2n))
      complex(8),        intent(out)   :: work  ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine zgeev
  end interface


  interface lapack_syev
    subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      real(4),      intent(inout) :: a       ! a(lda,n)
      real(4),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine ssyev

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      real(8),      intent(inout) :: a       ! a(lda,n)
      real(8),      intent(out)   :: w, work ! w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine dsyev
  end interface

  interface lapack_heev
    subroutine cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      complex(4),   intent(inout) :: a        ! a(lda,n)
      real(4),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(4),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine cheev

    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      complex(8),   intent(inout) :: a        ! a(lda,n)
      real(8),      intent(out)   :: w, rwork ! w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     ! work(lwork)
      integer,      intent(out)   :: info
    end subroutine zheev
  end interface

  interface
    SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE DGEQRF
  end interface

  interface
    SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, M, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE ZGEQRF
  end interface

  interface
    SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE DORGQR
  end interface

  interface
    SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE ZUNGQR
  end interface

end module lapack_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
