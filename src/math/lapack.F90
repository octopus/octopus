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
!> This module contains interfaces for LAPACK routines
! -----------------------------------------------------------------------

module lapack_m
  implicit none

  public ! only interfaces in this module

  interface lapack_potrf
    subroutine spotrf(uplo, n, a, lda, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(4),      intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine spotrf

    subroutine dpotrf(uplo, n, a, lda, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(8),      intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine dpotrf

    subroutine cpotrf(uplo, n, a, lda, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(4),   intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine cpotrf

    subroutine zpotrf(uplo, n, a, lda, info)
      implicit none
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(8),   intent(inout) :: a
      integer,      intent(in)    :: lda
      integer,      intent(out)   :: info
    end subroutine zpotrf
  end interface lapack_potrf

  interface lapack_sygv
    subroutine ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      real(4),      intent(inout) :: a, b    !< a(lda,n), b(ldb,n)
      real(4),      intent(out)   :: w, work !< w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine ssygv

    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: itype, n, lda, ldb, lwork
      real(8),      intent(inout) :: a, b    !< a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, work !< w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine dsygv
  end interface lapack_sygv

  interface lapack_hegv
    subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      complex(4),   intent(inout) :: a, b     !< a(lda,n), b(ldb,n)
      real(4),      intent(out)   :: w, rwork !< w(n), rwork(max(1,3*n-2))
      complex(4),   intent(out)   :: work     !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine chegv

    subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, itype, lda, ldb, lwork
      complex(8),   intent(inout) :: a, b     !< a(lda,n), b(ldb,n)
      real(8),      intent(out)   :: w, rwork !< w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine zhegv
  end interface lapack_hegv

  interface lapack_geev
    subroutine sgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      real(4),      intent(inout) :: a !< a(lda,n)
      real(4),      intent(out)   :: w, vl, vr !< w(n), vl(ldvl,n), vl(ldvr,n)
      real(4),      intent(out)   :: rwork !< rwork(max(1,2n))
      real(4),      intent(out)   :: work  !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine sgeev

    subroutine dgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      real(8),      intent(inout) :: a !< a(lda,n)
      real(8),      intent(out)   :: w, vl, vr !< w(n), vl(ldvl,n), vl(ldvr,n)
      real(8),      intent(out)   :: rwork !< rwork(max(1,2n))
      real(8),      intent(out)   :: work  !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine dgeev

    subroutine cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      complex(4),   intent(inout) :: a !< a(lda,n)
      complex(4),   intent(out)   :: w, vl, vr !< w(n), vl(ldvl,n), vl(ldvr,n)
      real(4),      intent(out)   :: rwork !< rwork(max(1,2n))
      complex(4),   intent(out)   :: work  !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine cgeev

    subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      character(1), intent(in)    :: jobvl, jobvr
      integer,      intent(in)    :: n, lda, ldvl, ldvr, lwork
      complex(8),        intent(inout) :: a !< a(lda,n)
      complex(8),        intent(out)   :: w, vl, vr !< w(n), vl(ldvl,n), vl(ldvr,n)
      real(8),        intent(out)   :: rwork !< rwork(max(1,2n))
      complex(8),        intent(out)   :: work  !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine zgeev
  end interface lapack_geev


  interface lapack_syev
    subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      real(4),      intent(inout) :: a       !< a(lda,n)
      real(4),      intent(out)   :: w, work !< w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine ssyev

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      real(8),      intent(inout) :: a       !< a(lda,n)
      real(8),      intent(out)   :: w, work !< w(n), work(lwork)
      integer,      intent(out)   :: info
    end subroutine dsyev
  end interface lapack_syev

  interface lapack_heev
    subroutine cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      complex(4),   intent(inout) :: a        !< a(lda,n)
      real(4),      intent(out)   :: w, rwork !< w(n), rwork(max(1,3*n-2))
      complex(4),   intent(out)   :: work     !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine cheev

    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n, lda, lwork
      complex(8),   intent(inout) :: a        !< a(lda,n)
      real(8),      intent(out)   :: w, rwork !< w(n), rwork(max(1,3*n-2))
      complex(8),   intent(out)   :: work     !< work(lwork)
      integer,      intent(out)   :: info
    end subroutine zheev
  end interface lapack_heev

  interface lapack_geqrf
    subroutine sgeqrf( m, n, a, lda, tau, work, lwork, info )
      integer            info, lda, lwork, m, n
      real(4)            a, tau, work
    end subroutine sgeqrf

    subroutine dgeqrf( m, n, a, lda, tau, work, lwork, info )
      integer            info, lda, lwork, m, n
      double precision   a, tau, work
    end subroutine dgeqrf

    subroutine cgeqrf( m, n, a, lda, tau, work, lwork, info )
      integer            info, lda, lwork, m, n
      complex(4)         a, tau, work
    end subroutine cgeqrf
    
    subroutine zgeqrf( m, n, a, lda, tau, work, lwork, info )
      integer            info, lda, lwork, m, n
      complex(8)         a, tau, work
    end subroutine zgeqrf
  end interface lapack_geqrf
  
  interface lapack_orgqr 
    subroutine dorgqr( m, n, k, a, lda, tau, work, lwork, info )
      integer            info, k, lda, lwork, m, n
      double precision   a, tau, work
    end subroutine dorgqr
    
    subroutine zungqr( m, n, k, a, lda, tau, work, lwork, info )
      integer            info, k, lda, lwork, m, n
      complex(8)         a, tau, work
    end subroutine zungqr

    subroutine sorgqr( m, n, k, a, lda, tau, work, lwork, info )
      integer            info, k, lda, lwork, m, n
      real               a, tau, work
    end subroutine sorgqr

    subroutine cungqr( m, n, k, a, lda, tau, work, lwork, info )
      integer            info, k, lda, lwork, m, n
      complex(4)         a, tau, work
    end subroutine cungqr
  end interface lapack_orgqr

  interface lapack_sygvx
    subroutine dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &
      m, w, z, ldz, work, lwork, iwork, ifail, info)
      implicit none

      integer,             intent(in)    :: itype !< Specifies the problem: 1: A*x = l*B*x 2:  A*B*x = l*x 3: B*A*x = l*x
      character(len=1),    intent(in)    :: jobz  !< N: Compute eigenvalues only; V: Compute eigenvalues and eigenvectors.
      character(len=1),    intent(in)    :: range !< A: all eigenval V: all eigenval in (VL,VU] I: IL-th through IU-th eigenval
      character(len=1),    intent(in)    :: uplo  !< U: Upper triangle of A and B stored L: Lower triangle of A and B stored
      integer,             intent(in)    :: n     !< The order of the matrix pencil (A,B)
      real(8),             intent(inout) :: a     !< a(:) On entry, the symmetric matrix A. On exit, destroyed
      integer,             intent(in)    :: lda   !< The leading dimension of the array A
      real(8),             intent(inout) :: b     !< b(:)
      integer,             intent(in)    :: ldb
      real(8),             intent(in)    :: vl
      real(8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      real(8),             intent(out)   :: w
      real(8),             intent(out)   :: z
      integer,             intent(in)    :: ldz
      real(8),             intent(out)   :: work   !< work(:)
      integer,             intent(in)    :: lwork
      integer,             intent(out)   :: iwork  !< iwork(1:5*n)
      integer,             intent(out)   :: ifail  !< ifail(1:n)
      integer,             intent(out)   :: info
    end subroutine dsygvx

    subroutine ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &
      m, w, z, ldz, work, lwork, iwork, ifail, info)
      implicit none

      integer,             intent(in)    :: itype !< Specifies the problem: 1: A*x = l*B*x 2:  A*B*x = l*x 3: B*A*x = l*x
      character(len=1),    intent(in)    :: jobz  !< N: Compute eigenvalues only; V: Compute eigenvalues and eigenvectors.
      character(len=1),    intent(in)    :: range !< A: all eigenval V: all eigenval in (VL,VU] I: IL-th through IU-th eigenval
      character(len=1),    intent(in)    :: uplo  !< U: Upper triangle of A and B stored L: Lower triangle of A and B stored
      integer,             intent(in)    :: n     !< The order of the matrix pencil (A,B)
      real(4),             intent(inout) :: a     !< a(:) On entry, the symmetric matrix A. On exit, destroyed
      integer,             intent(in)    :: lda   !< The leading dimension of the array A
      real(4),             intent(inout) :: b     !< b(:)
      integer,             intent(in)    :: ldb
      real(4),             intent(in)    :: vl
      real(4),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(4),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      real(4),             intent(out)   :: w
      real(4),             intent(out)   :: z
      integer,             intent(in)    :: ldz
      real(4),             intent(out)   :: work   !< work(:)
      integer,             intent(in)    :: lwork
      integer,             intent(out)   :: iwork  !< iwork(1:5*n)
      integer,             intent(out)   :: ifail  !< ifail(1:n)
      integer,             intent(out)   :: info
    end subroutine ssygvx
  end interface lapack_sygvx

  interface lapack_hegvx
    subroutine zhegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &
      m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      implicit none

      integer,             intent(in)    :: itype
      character(len=1),    intent(in)    :: jobz
      character(len=1),    intent(in)    :: range
      character(len=1),    intent(in)    :: uplo
      integer,             intent(in)    :: n
      complex(8),          intent(inout) :: a
      integer,             intent(in)    :: lda
      complex(8),          intent(inout) :: b
      integer,             intent(in)    :: ldb
      real(8),             intent(in)    :: vl
      real(8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      real(8),             intent(out)   :: w
      complex(8),          intent(out)   :: z
      integer,             intent(in)    :: ldz
      complex(8),          intent(out)   :: work
      integer,             intent(in)    :: lwork
      real(8),             intent(out)   :: rwork !< rwork(1:7*n)
      integer,             intent(out)   :: iwork !< iwork(1:5*n)
      integer,             intent(out)   :: ifail !< ifail(1:n)
      integer,             intent(out)   :: info
    end subroutine zhegvx

    subroutine chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, &
      m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      implicit none

      integer,             intent(in)    :: itype
      character(len=1),    intent(in)    :: jobz
      character(len=1),    intent(in)    :: range
      character(len=1),    intent(in)    :: uplo
      integer,             intent(in)    :: n
      complex(4),          intent(inout) :: a
      integer,             intent(in)    :: lda
      complex(4),          intent(inout) :: b
      integer,             intent(in)    :: ldb
      real(4),             intent(in)    :: vl
      real(4),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(4),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      real(4),             intent(out)   :: w
      complex(4),          intent(out)   :: z
      integer,             intent(in)    :: ldz
      complex(4),          intent(out)   :: work
      integer,             intent(in)    :: lwork
      real(4),             intent(out)   :: rwork !< rwork(1:7*n)
      integer,             intent(out)   :: iwork !< iwork(1:5*n)
      integer,             intent(out)   :: ifail !< ifail(1:n)
      integer,             intent(out)   :: info
    end subroutine chegvx
  end interface lapack_hegvx

end module lapack_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
