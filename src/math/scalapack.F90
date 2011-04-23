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
!> This module contains interfaces for ScaLAPACK routines
!! Interfaces are from http://www.netlib.org/scalapack/tools, double, complex16
! -----------------------------------------------------------------------

module scalapack_m
  implicit none
#ifdef HAVE_SCALAPACK

  interface
    integer function iceil(inum, idenom)
      implicit none

      integer            idenom, inum
    end function iceil
  end interface

  interface
    subroutine descinit(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info)
      implicit none

      integer            icsrc, ictxt, info, irsrc, lld, m, mb, n, nb
      integer            desc
    end subroutine descinit
  end interface

  interface scalapack_geqrf
    subroutine pdgeqrf(m, n, a, ia, ja, desca, tau, work, lwork, info)
      implicit none

      integer            ia, info, ja, lwork, m, n
      integer            desca
      double precision   a, tau, work
    end subroutine pdgeqrf
    
    subroutine pzgeqrf(m, n, a, ia, ja, desca, tau, work, lwork, info)
      implicit none
      
      integer            ia, info, ja, lwork, m, n
      integer            desca
      complex*16         a, tau, work
    end subroutine pzgeqrf
  end interface scalapack_geqrf

  interface scalapack_orgqr
    subroutine pdorgqr(m, n, k, a, ia, ja, desca, tau, work, lwork, info) 
      implicit none

      integer            ia, info, ja, k, lwork, m, n
      integer            desca
      double precision   a, tau, work
    end subroutine pdorgqr
    
    subroutine pzungqr(m, n, k, a, ia, ja, desca, tau, work, lwork, info)
      implicit none
      
      integer            ia, info, ja, k, lwork, m, n
      integer            desca
      complex*16         a, tau, work
    end subroutine pzungqr
  end interface scalapack_orgqr
  
  interface
    subroutine pdgesv(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info)
      implicit none

      integer            ia, ib, info, ja, jb, n, nrhs
      integer            desca, descb, ipiv
      double precision   a, b
    end subroutine pdgesv
  end interface

  interface
    subroutine pzgesv(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info)
      implicit none

      integer            ia, ib, info, ja, jb, n, nrhs
      integer            desca, descb, ipiv
      complex*16         a, b
    end subroutine pzgesv
  end interface

  interface
    subroutine pdsyev(jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info)
      implicit none
      
      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n
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
      implicit none

      character,        intent(in)    :: jobz
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n
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

  interface scalapack_syev
    subroutine pdsyevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, &
         m, nz, w, orfac, z, iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info )
      
      character(1), intent(in)    :: jobz     
      character(1), intent(in)    :: range    
      character(1), intent(in)    :: uplo     
      integer,      intent(in)    :: n         
      real(8),      intent(inout) :: a         
      integer,      intent(in)    :: ia        
      integer,      intent(in)    :: ja        
      integer,      intent(in)    :: desca     
      real(8),      intent(in)    :: vl        
      real(8),      intent(in)    :: vu        
      integer,      intent(in)    :: il        
      integer,      intent(in)    :: iu        
      real(8) ,     intent(in)    :: abstol    
      integer,      intent(out)   :: m         
      integer,      intent(out)   :: nz       
      real(8),      intent(out)   :: w        
      real(8),      intent(in)    :: orfac   
      real(8),      intent(out)   :: z       
      integer,      intent(in)    :: iz      
      integer,      intent(in)    :: jz      
      integer,      intent(in)    :: descz   
      real(8),      intent(out)   :: work    
      integer,      intent(in)    :: lwork   
      integer,      intent(inout) :: iwork   
      integer,      intent(in)    :: liwork  
      integer,      intent(out)   :: ifail  
      integer,      intent(out)   :: iclustr
      real(8),      intent(out)   :: gap    
      integer,      intent(out)   :: info   

    end subroutine pdsyevx

    subroutine pzheevx( jobz, range, uplo, n, a, ia, ja, desca, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, &
         jz, descz, work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info )
      character(1), intent(in)    :: jobz 
      character(1), intent(in)    :: range
      character(1), intent(in)    :: uplo 
      integer,      intent(in)    :: n    
      complex(8),   intent(inout) :: a    
      integer,      intent(in)    :: ia   
      integer,      intent(in)    :: ja   
      integer,      intent(in)    :: desca
      real(8),      intent(in)    :: vl         
      real(8),      intent(in)    :: vu         
      integer,      intent(in)    :: il
      integer,      intent(in)    :: iu
      real(8),      intent(in)    :: abstol     
      integer,      intent(out)   :: m 
      integer,      intent(out)   :: nz
      real(8),      intent(out)   :: w          
      real(8),      intent(in)    :: orfac      
      complex(8),   intent(out)   :: z 
      integer,      intent(in)    :: iz         
      integer,      intent(in)    :: jz         
      integer,      intent(in)    :: descz      
      complex(8),   intent(out)   :: work       
      integer,      intent(in)    :: lwork      
      real(8),      intent(out)   :: rwork      
      integer,      intent(in)    :: lrwork    
      integer,      intent(inout) :: iwork     
      integer,      intent(in)    :: liwork    
      integer,      intent(out)   :: ifail     
      integer,      intent(out)   :: iclustr   
      real(8),      intent(out)   :: gap      
      integer,      intent(out)   :: info    
    end subroutine pzheevx
  end interface scalapack_syev
  
  interface scalapack_sygvx
    subroutine pdsygvx(ibtype, jobz, range, uplo, n, a, ia, ja,       &
      desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz,  &
      work, lwork, iwork, liwork, ifail, iclustr, gap, info)
      implicit none

      integer,             intent(in)    :: ibtype
      character,           intent(in)    :: jobz
      character,           intent(in)    :: range
      character,           intent(in)    :: uplo
      integer,             intent(in)    :: n
      real(8),             intent(inout) :: a
      integer,             intent(in)    :: ia
      integer,             intent(in)    :: ja
      integer,             intent(in)    :: desca
      real(8),             intent(inout) :: b
      integer,             intent(in)    :: ib
      integer,             intent(in)    :: jb
      integer,             intent(in)    :: descb
      real(8),             intent(in)    :: vl
      real(8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      integer,             intent(out)   :: nz
      real(8),             intent(in)    :: w
      real(8),             intent(in)    :: orfac
      real(8),             intent(inout) :: z
      integer,             intent(in)    :: iz
      integer,             intent(in)    :: jz
      integer,             intent(in)    :: descz
      real(8),             intent(out)   :: work
      integer,             intent(in)    :: lwork
      integer,             intent(out)   :: iwork
      integer,             intent(in)    :: liwork
      integer,             intent(out)   :: ifail
      integer,             intent(out)   :: iclustr
      real(8),             intent(out)   :: gap
      integer,             intent(out)   :: info
    end subroutine pdsygvx
  end interface scalapack_sygvx
  
  ! -------------------------------------------------------------
 
  interface scalapack_hegvx
    subroutine pzhegvx(ibtype, jobz, range, uplo, n, a, ia, ja,       &
      desca, b, ib, jb, descb, vl, vu, il, iu, abstol, m, nz, w, orfac, z, iz, jz, descz,  &
      work, lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
      implicit none

      integer,             intent(in)    :: ibtype
      character,           intent(in)    :: jobz
      character,           intent(in)    :: range
      character,           intent(in)    :: uplo
      integer,             intent(in)    :: n
      complex(8),          intent(inout) :: a
      integer,             intent(in)    :: ia
      integer,             intent(in)    :: ja
      integer,             intent(in)    :: desca
      complex(8),          intent(inout) :: b
      integer,             intent(in)    :: ib
      integer,             intent(in)    :: jb
      integer,             intent(in)    :: descb
      real(8),             intent(in)    :: vl
      real(8),             intent(in)    :: vu
      integer,             intent(in)    :: il
      integer,             intent(in)    :: iu
      real(8),             intent(in)    :: abstol
      integer,             intent(out)   :: m
      integer,             intent(out)   :: nz
      real(8),             intent(in)    :: w
      real(8),             intent(in)    :: orfac
      complex(8),          intent(inout) :: z
      integer,             intent(in)    :: iz
      integer,             intent(in)    :: jz
      integer,             intent(in)    :: descz
      complex(8),          intent(out)   :: work
      integer,             intent(in)    :: lwork
      real(8),             intent(out)   :: rwork
      integer,             intent(in)    :: lrwork
      integer,             intent(out)   :: iwork
      integer,             intent(in)    :: liwork
      integer,             intent(out)   :: ifail
      integer,             intent(out)   :: iclustr
      real(8),             intent(out)   :: gap
      integer,             intent(out)   :: info
    end subroutine pzhegvx
  end interface scalapack_hegvx

  ! -------------------------------------------------------------

  interface scalapack_potrf
    subroutine pdpotrf(uplo, n, a, ia, ja, desca, info)
      implicit none

      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      real(8),      intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      integer,      intent(out)   :: info
    end subroutine pdpotrf
    
    subroutine pzpotrf(uplo, n, a, ia, ja, desca, info)
      implicit none

      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n
      complex(8),   intent(inout) :: a
      integer,      intent(in)    :: ia
      integer,      intent(in)    :: ja
      integer,      intent(in)    :: desca
      integer,      intent(out)   :: info
    end subroutine pzpotrf
    
  end interface scalapack_potrf

#endif
end module scalapack_m

!! local Variables:
!! mode: f90
!! coding: utf-8
!! End:
