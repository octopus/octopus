!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

! This module is intended to contain "only mathematical" functions and        !
!	procedures.                                                           !

module math
  use global

  implicit none

  interface
     !--------------------------------------------------------------------------------!
     ! To improve performance only the first element of an array should be passed to  !
     ! the BLAS and LAPACK routines. The interfaces in this module enforce that. Real !
     ! rank and dimensions of the arrays are indicated in front of the corresponding  !
     ! scalar arguments.                                                              !
     !--------------------------------------------------------------------------------!

     !--- BLAS interfaces ---!

     subroutine dscal(n, da, dx, incx)
       ! scales a vector by a constant
       use global
       integer, intent(in) :: n, incx
       real(r8), intent(in) :: da
       real(r8), intent(inout) :: dx ! dx(n)
     end subroutine dscal

     subroutine zscal(n, za, zx, incx)
       ! scales a vector by a constant
       use global
       integer, intent(in) :: n, incx
       complex(r8), intent(in) :: za
       complex(r8), intent(inout) :: zx ! zx(n)
     end subroutine zscal

     subroutine daxpy(n, da, dx, incx, dy, incy)
       ! constant times a vector plus a vector
       use global
       integer, intent(in) :: n, incx, incy
       real(r8), intent(in) :: da, dx ! dx(n)
       real(r8), intent(inout) :: dy ! dy(n)
     end subroutine daxpy

     subroutine zaxpy(n, za, zx, incx, zy, incy)
       ! constant times a vector plus a vector
       use global
       integer, intent(in) :: n, incx, incy
       complex(r8), intent(in) :: za, zx ! zx(n)
       complex(r8), intent(inout) :: zy ! zy(n)
     end subroutine zaxpy

     function ddot(n, dx, incx, dy, incy)
       ! forms the dot product of two vectors
       use global
       real(r8) :: ddot
       integer, intent(in) :: n, incx, incy
       real(r8), intent(in) :: dx, dy ! dx(n), dy(n)
     end function ddot

     function zdotc(n, zx, incx, zy, incy)
       ! forms the dot product of two vectors
       use global
       complex(r8) :: zdotc
       integer, intent(in) :: n, incx, incy
       complex(r8), intent(in) :: zx, zy ! zx(n), zy(n)
     end function zdotc

     function dnrm2(n, dx, incx)
       ! returns the euclidean norm of a vector
       use global
       real(r8) :: dnmr2
       integer, intent(in) :: n, incx
       real(r8), intent(in) :: dx ! dx(n)
     end function dnrm2

     function dznrm2(n, zx, incx)
       ! returns the euclidean norm of a vector
       use global
       real(r8) :: dznrm2
       integer, intent(in) :: n, incx
       complex(r8), intent(in) :: zx ! zx(n)
     end function dznrm2

     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       ! matrix-matrix multiplication
       use global
       character(1), intent(in) :: transa, transb
       integer, intent(in) :: m, n, k, lda, ldb, ldc
       real(r8), intent(in) :: alpha, beta
       real(r8), intent(in) :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
       real(r8), intent(in) :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
       real(r8), intent(inout) :: c ! c(ldc,n) 
     end subroutine dgemm

     subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       ! matrix-matrix multiplication
       use global
       character(1), intent(in) :: transa, transb
       integer, intent(in) :: m, n, k, lda, ldb, ldc
       complex(r8), intent(in) :: alpha, beta
       complex(r8), intent(in) :: a ! a(lda,ka)    ka=k if transa='N' or 'n'; m otherwise
       complex(r8), intent(in) :: b ! b(ldb,kb)    kb=k if transa='N' or 'n'; m otherwise
       complex(r8), intent(inout) :: c ! c(ldc,n)
     end subroutine zgemm

     subroutine dcopy(n, dx, incx, dy, incy)
       ! copies a vector, x, to a vector, y
       use global
       integer, intent(in) :: n, incx, incy
       real(r8), intent(in) :: dx ! dx(n)
       real(r8), intent(out) :: dy ! dy(n)
     end subroutine dcopy

     subroutine zcopy(n, zx, incx, zy, incy)
       ! copies a vector, x, to a vector, y
       use global
       integer, intent(in) :: n, incx, incy
       complex(r8), intent(in) :: zx ! dz(n)
       complex(r8), intent(out) :: zy ! zy(n)
     end subroutine zcopy


     !--- LAPACK interfaces ---!

     subroutine dgetrf(m, n, a, lda, ipiv, info)
       ! computes an LU factorization of a general M-by-N matrix A
       ! using partial pivoting with row interchanges.
       use global
       integer, intent(in) :: m, n, lda
       real(r8), intent(inout) :: a ! a(lda,n)
       integer, intent(out) :: ipiv ! ipiv(min(m,n))
       integer, intent(out) :: info
     end subroutine dgetrf

     subroutine zgetrf(m, n, a, lda, ipiv, info)
       ! computes an LU factorization of a general M-by-N matrix A
       ! using partial pivoting with row interchanges.
       use global
       integer, intent(in) :: m, n, lda
       complex(r8), intent(inout) :: a ! a(lda,n)
       integer, intent(out) :: ipiv ! ipiv(min(m,n))
       integer, intent(out) :: info
     end subroutine zgetrf

     subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
     ! computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
       use global
       character(1), intent(in) :: jobz, uplo
       integer, intent(in) :: n, lda, lwork
       real(r8), intent(inout) :: a ! a(lda,n)
       real(r8), intent(out) :: w, work ! w(n), work(lwork)
       integer, intent(out) :: info 
     end subroutine dsyev

     subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
     ! computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
       use global
       character(1), intent(in) :: jobz, uplo
       integer, intent(in) :: n, lda, lwork
       complex(r8), intent(inout) :: a ! a(lda,n)
       real(r8), intent(out) :: w, rwork ! w(n), rwork(max(1,3*n-2))
       complex(r8), intent(out) :: work ! work(lwork)
       integer, intent(out) :: info
     end subroutine zheev

     subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
       ! computes all the eigenvalues, and optionally, the eigenvectors of a real
       ! generalized symmetric-definite eigenproblem, of the form  A*x=(lambda)*B*x,
       ! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
       ! Here A and B are assumed to be symmetric and B is also positive definite.
       use global
       character(1), intent(in) :: jobz, uplo
       integer, intent(in) :: itype, n, lda, ldb, lwork
       real(r8), intent(inout) :: a, b ! a(lda,n), b(ldb,n)
       real(r8), intent(out) :: w, work ! w(n), work(lwork)
       integer, intent(out) :: info
     end subroutine dsygv

     subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
       ! computes all the eigenvalues, and optionally, the eigenvectors of a complex 
       ! generalized Hermitian-definite eigenproblem, of the form  A*x=(lambda)*B*x,  
       ! A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
       ! Here A and B are assumed to be Hermitian and B is also positive definite.
       use global
       character(1), intent(in) :: jobz, uplo
       integer, intent(in) :: n, itype, lda, ldb, lwork
       complex(r8), intent(inout) :: a, b ! a(lda,n), b(ldb,n)
       real(r8), intent(out) :: w, rwork ! w(n), rwork(max(1,3*n-2))
       complex(r8), intent(out) :: work ! work(lwork)
       integer, intent(out) :: info 
     end subroutine zhegv

  end interface

contains

! a simple congruent random number generator
subroutine quickrnd(iseed, rnd)
  integer, intent(inout) :: iseed
  real(r8), intent(inout) :: rnd

  integer, parameter :: im=6075, ia=106, ic=1283

  iseed = mod(iseed*ia + ic, im)
  rnd = real(iseed, r8)/real(im, r8)
  
end subroutine quickrnd

! Step function, needed for definition of fermi function.
function stepf(x)
  real(r8), intent(in) ::  x
  real(r8) :: stepf

  if (x.gt.100.0_r8) then
     stepf = 0.0_r8
  elseif (x.lt.-100.0_r8) then
     stepf = 2.0_r8
  else
     stepf = 2.0_r8 / ( 1.0_r8 + exp(x) )
  endif

end function stepf

! This is a Numerical Recipes based subroutine
! computes real spherical harmonics ylm in the direction of vector r:
!    ylm = c * plm( cos(theta) ) * sin(m*phi)   for   m <  0
!    ylm = c * plm( cos(theta) ) * cos(m*phi)   for   m >= 0
! with (theta,phi) the polar angles of r, c a positive normalization
subroutine grylmr(x, y, z, li, mi, ylm, grylm)
  integer, intent(in) :: li, mi
  real(r8), intent(in) :: x, y, z
  real(r8), intent(out) :: ylm, grylm(3)

  integer, parameter :: lmaxd = 20
  real(r8), parameter :: tiny=1.e-30_r8, zero=0.0_r8, half=0.5_r8, one=1.0_r8, &
       two=2.0_r8, three=3.0_r8, six=6.0_r8
  integer :: i, ilm0, l, m, mabs
  integer, save :: lmax = -1

  real(r8) :: cmi, cosm, cosmm1, cosphi, dphase, dplg, fac, &
       fourpi, plgndr, phase, pll, pmm, pmmp1, sinm, &
       sinmm1, sinphi, r2, rsize, Rx, Ry, Rz, xysize
  real(r8), save :: c(0:(lmaxd+1)*(lmaxd+1))
  

! evaluate normalization constants once and for all
  if (li.gt.lmax) then
    fourpi = two**4*datan(one)
    do l = 0, li
      ilm0 = l*l + l
      do m = 0, l
        fac = (2*l+1)/fourpi
        do i = l - m + 1, l + m
          fac = fac/i
        end do
        c(ilm0 + m) = sqrt(fac)
        ! next line because ylm are real combinations of m and -m
        if(m.ne.0) c(ilm0 + m) = c(ilm0 + m)*sqrt(two)
        c(ilm0 - m) = c(ilm0 + m)
      end do
    end do
    lmax = li
  end if

  ! if l=0, no calculations are required
  if (li.eq.0) then
    ylm = c(0)
    grylm(:) = zero
    return
  end if

  ! if r=0, direction is undefined => make ylm=0 except for l=0
  r2 = x**2 + y**2 + z**2
  if(r2.lt.tiny) then
    ylm = zero
    grylm(:) = zero
    return
  endif
  rsize = sqrt(r2)

  Rx = x/rsize
  Ry = y/rsize
  Rz = z/rsize

  ! explicit formulas for l=1 and l=2
  if(li.eq.1) then
    select case(mi)
    case(-1)
      ylm = (-c(1))*Ry
      grylm(1) = c(1)*Rx*Ry/rsize
      grylm(2) = (-c(1))*(one - Ry*Ry)/rsize
      grylm(3) = c(1)*Rz*Ry/rsize 
    case(0)
      ylm = c(2)*Rz
      grylm(1) = (-c(2))*Rx*Rz/rsize
      grylm(2) = (-c(2))*Ry*Rz/rsize
      grylm(3) = c(2)*(one - Rz*Rz)/rsize
    case(1)
      ylm = (-c(3))*Rx
      grylm(1) = (-c(3))*(one - Rx*Rx)/rsize
      grylm(2) = c(3)*Ry*Rx/rsize
      grylm(3) = c(3)*Rz*Rx/rsize
    end select
    return
  end if
         
  if(li.eq.2) then
    select case(mi)
    case(-2)
      ylm = c(4)*six*Rx*Ry
      grylm(1) = (-c(4))*six*(two*Rx*Rx*Ry - Ry)/rsize
      grylm(2) = (-c(4))*six*(two*Ry*Rx*Ry - Rx)/rsize
      grylm(3) = (-c(4))*six*(two*Rz*Rx*Ry)/rsize
    case(-1)
      ylm = (-c(5))*three*Ry*Rz
      grylm(1) = c(5)*three*(two*Rx*Ry*Rz)/rsize
      grylm(2) = c(5)*three*(two*Ry*Ry*Rz - Rz)/rsize
      grylm(3) = c(5)*three*(two*Rz*Ry*Rz - Ry)/rsize
    case(0)
      ylm = c(6)*half*(three*Rz*Rz - one)
      grylm(1) = (-c(6))*three*(Rx*Rz*Rz)/rsize
      grylm(2) = (-c(6))*three*(Ry*Rz*Rz)/rsize
      grylm(3) = (-c(6))*three*(Rz*Rz - one)*Rz/rsize
    case(1)
      ylm = (-c(7))*three*Rx*Rz
      grylm(1) = c(7)*three*(two*Rx*Rx*Rz - Rz)/rsize
      grylm(2) = c(7)*three*(two*Ry*Rx*Rz)/rsize
      grylm(3) = c(7)*three*(two*Rz*Rx*Rz - Rx)/rsize
    case(2)
      ylm = c(8)*three*(Rx*Rx - Ry*Ry)
      grylm(1) = (-c(8))*six*(Rx*Rx - Ry*Ry - one)*Rx/rsize
      grylm(2) = (-c(8))*six*(Rx*Rx - Ry*Ry + one)*Ry/rsize
      grylm(3) = (-c(8))*six*(Rx*Rx - Ry*Ry)*Rz/rsize
    end select
    return
  end if

! general algorithm based on routine plgndr of numerical recipes
  mabs = abs(mi)
  xysize = sqrt(max(Rx*Rx + Ry*Ry, tiny))
  cosphi = Rx/xysize
  sinphi = Ry/xysize
  cosm = one
  sinm = zero
  do m = 1, mabs
    cosmm1 = cosm
    sinmm1 = sinm
    cosm = cosmm1*cosphi - sinmm1*sinphi
    sinm = cosmm1*sinphi + sinmm1*cosphi
  end do
     
  if(mi.lt.0) then
    phase = sinm
    dphase = mabs*cosm
  else
    phase = cosm
    dphase = (-mabs)*sinm
  end if

  pmm = one
  fac = one

  if(mabs.gt.zero) then
    do i = 1, mabs
      pmm = (-pmm)*fac*xysize
      fac = fac + two
    end do
  end if

  if(li.eq.mabs) then
    plgndr = pmm
    dplg = (-li)*Rz*pmm/(xysize**2)
  else
    pmmp1 = Rz*(2*mabs + 1)*pmm
    if(li.eq.mabs + 1) then
      plgndr = pmmp1
      dplg = -((li*Rz*pmmp1 - (mabs + li)*pmm)/(xysize**2))
    else
      do l = mabs + 2, li
        pll = (Rz*(2*l - 1)*pmmp1 - (l + mabs - 1)*pmm)/(l - mabs)
        pmm = pmmp1
        pmmp1 = pll
      end do
      plgndr = pll
      dplg = -((li*Rz*pll - (l + mabs - 1)*pmm)/(xysize**2))
    end if
  end if   

  ilm0 = li*li + li
  cmi = c(ilm0 + mi)
  ylm = cmi*plgndr*phase
  grylm(1) = (-cmi)*dplg*Rx*Rz*phase/rsize     &
       -cmi*plgndr*dphase*Ry/(rsize*xysize**2)
  grylm(2) = (-cmi)*dplg*Ry*Rz*phase/rsize     &
       +cmi*plgndr*dphase*Rx/(rsize*xysize**2)
  grylm(3)= cmi*dplg*(one - Rz*Rz)*phase/rsize
   
  return
end subroutine grylmr

! Compute the Weights for finite-difference calculations:
!
!  N -> highest order fo the derivative to be approximated
!  M -> number of grid points to be used in the approsimation.
!
!  c(j,k,i) -> ith order derivative at kth-order approximation
!              j=0,k: the coefficients acting of each point
subroutine weights(N, M, cc)
  integer, intent(in) :: N, M
  real(r8), intent(out) :: cc(0:M, 0:M, 0:N)

  integer :: i, j, k, mn
  real(r8) :: c1, c2, c3, c4, c5, xi
  real(r8) :: x(0:M)

  ! grid-points for one-side finite-difference formulas on an equi.spaced grid
  ! x(:) = (/(i,i=0,M)/) 

  ! grid-points for centered finite-difference formulas on an equi.spaced grid
  mn = M/2
  x(:) = (/0,(-i,i,i=1,mn)/)

  xi = 0.0_r8  ! point at which the approx. are to be accurate

  cc = 0.0_r8
  cc(0,0,0) = 1.0_r8

  c1 = 1.0_r8
  c4 = x(0) - xi
       
  do j = 1, M
    mn = min(j,N)
    c2 = 1.0_r8
    c5 = c4
    c4 = x(j) - xi
    
    do k = 0, j - 1
      c3 = x(j) - x(k)
      c2 = c2*c3
      
      if (j <= N) cc(k, j - 1, j) = 0.0_r8
      cc(k, j, 0) = c4*cc(k, j - 1, 0)/c3
      
      do i = 1, mn
        cc(k, j, i) = (c4*cc(k, j - 1, i) - i*cc(k, j - 1, i - 1))/c3
      end do
      
      cc(j, j, 0) = -c1*c5*cc(j - 1, j - 1, 0) / c2         
    end do
    
    do i = 1, mn
      cc(j, j, i) = c1*(i*cc(j - 1, j - 1, i - 1) - c5*cc(j - 1, j - 1, i))/c2
    end do

    c1 = c2
  end do

end subroutine weights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ddet and zdet return the determinant of a square matrix a of dimensions (n,n)
! (ddet for real and zdet for complex numbers)
real(r8) function ddet(a, n)
  integer, intent(in) :: n
  real(r8) :: a(n, n)

  integer :: i, info, ipiv(n)

  call dgetrf(n, n, a(1,1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1),'(a,i4)') 'In det, LAPACK dgetrf returned error code ',info
    call write_fatal(1)
  endif
  ddet = M_ONE
  do i = 1, n
     ddet = ddet*a(i, i)
     if(ipiv(i).ne.i) then
       ipiv(ipiv(i)) = ipiv(i)
       ddet = -ddet
     endif
  enddo

end function ddet

complex(r8) function zdet(a, n)
  integer, intent(in) :: n
  complex(r8) :: a(n, n)

  integer :: info, i, ipiv(n)

  call zgetrf(n, n, a(1,1), n, ipiv(1), info)
  if(info < 0) then
    write(message(1),'(a,i4)') 'In zdet, LAPACK zgetrf returned error code ',info
    call write_fatal(1)
  endif
  zdet = M_z1
  do i = 1, n
     zdet = zdet*a(i, i)
     if(ipiv(i).ne.i) then
       ipiv(ipiv(i)) = ipiv(i)
       zdet = -zdet
     endif
  enddo

end function zdet

!!! These routines diagonalize a matrix using LAPACK
subroutine diagonalise(dim, a, b, e)
  integer  :: dim
  real(r8) :: a(dim, dim), b(dim, dim)
  real(r8) :: e(dim)
  
  integer :: info, lwork
  real(r8), allocatable :: work(:)
  
  lwork = 6*dim
  allocate(work(lwork))
  b = a
  call dsyev('V', 'U', dim, b(1,1), dim, e(1), work(1), lwork, info)
  if(info.ne.0) then
    write(message(1),'(a,i5)') 'dsyev returned error message ', info
    call write_fatal(1)
  endif
  deallocate(work)
end subroutine diagonalise

subroutine ziagonalise(dim, a, b, e)
  integer     :: dim
  complex(r8) :: a(dim, dim), b(dim, dim)
  real(r8)    :: e(dim)
  
  integer :: info, lwork
  complex(r8), allocatable :: work(:)
  real(r8), allocatable :: rwork(:)
  
  lwork = 6*dim
  allocate(work(lwork), rwork(max(1,3*dim-2)))
  b = a
  call zheev('V','U', dim, b(1,1), dim, e(1), work(1), lwork, rwork(1), info)
  if(info.ne.0) then
    write(message(1),'(a,i5)') 'dsyev returned error message ', info
    call write_fatal(1)
  endif
  deallocate(work, rwork)
end subroutine ziagonalise

#include "undef.F90"
#include "complex.F90"
#include "math_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "math_inc.F90"

end module math
