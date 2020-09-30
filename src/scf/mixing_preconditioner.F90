!! Copyright (C) 2020 N. Tancogne-Dejean
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module mixing_preconditioner_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                               &
    mixing_preconditioner_t,              &
    mixing_preconditioner_init,           &
    mixing_preconditioner_end

  integer, public, parameter ::     &
    MIX_PRE_NONE  = 0,              &
    MIX_KERKER    = 1,              &
    MIX_TKERKER   = 2,              &
    MIX_RESTA     = 3

  type mixing_preconditioner_t
    private
    integer :: type

    !Parameter of different preconditioners
    FLOAT :: alpha0
    FLOAT :: kf
    FLOAT :: eps0

    integer, public :: tt !> Order of the fiting expansion 

    FLOAT, allocatable, public :: coeff(:) !> Coefficient of the nonlinear fit of the preconditioner

  end type mixing_preconditioner_t


contains

  ! ---------------------------------------------------------
  subroutine mixing_preconditioner_init(this)
    type(mixing_preconditioner_t), intent(inout) :: this

    integer :: ip, np, ncoeff, iter
    FLOAT, allocatable :: qq(:), Pq(:), coeff(:), alpha(:,:), covar(:,:)
    FLOAT :: alamda, chisq
    FLOAT, parameter :: tol=CNST(1e-2) 

    PUSH_SUB(mixing_preconditioner_init)

    this%alpha0 = CNST(0.17)
    this%kf = CNST(0.529) 
    this%eps0 = M_ZERO

    this%tt = 2

    !We evaluate the preconditioner on a grid ranging from 0 to 20 Bohr^-1
    !with a spacing of 0.1 Bohr^-1
    np = 200
    SAFE_ALLOCATE(qq(1:np))
    SAFE_ALLOCATE(Pq(1:np))

    do ip = 1, np
     qq(ip) = CNST(0.1)*(ip-1) + CNST(1e-10)
     Pq(ip) = max(this%alpha0, qq(ip)**2/(this%kf**2 + qq(ip)**2))
    end do

    ncoeff = 2*this%tt+1
    SAFE_ALLOCATE(coeff(1:ncoeff))
    SAFE_ALLOCATE(covar(1:ncoeff, 1:ncoeff))
    SAFE_ALLOCATE(alpha(1:ncoeff, 1:ncoeff))

    SAFE_ALLOCATE(this%coeff(1:ncoeff))
    
    alamda = -M_ONE !Initialization
    !Initial guess
    do ip = 1, ncoeff
      coeff(ip) = ip
    end do
    chisq = M_ONE
    iter = 1
    do while (chisq > tol .and. iter < 100) 
      call mrqmin(this, qq, Pq, coeff, covar, alpha, chisq, precondtioner_model, alamda) 
      iter = iter + 1
    end do
 
    !We store the result of the nonlinear fit
    this%coeff = coeff

    !We free the memory
    alamda = M_ZERO
    call mrqmin(this, qq, Pq, coeff, covar, alpha, chisq, precondtioner_model, alamda)

    SAFE_DEALLOCATE_A(alpha)
    SAFE_DEALLOCATE_A(covar)
    SAFE_DEALLOCATE_A(qq)
    SAFE_DEALLOCATE_A(Pq)
    SAFE_DEALLOCATE_A(coeff)

    POP_SUB(mixing_preconditioner_init)
  end subroutine mixing_preconditioner_init

  ! ---------------------------------------------------------
  subroutine mixing_preconditioner_end(this)
    type(mixing_preconditioner_t), intent(inout) :: this

    PUSH_SUB(mixing_preconditioner_end)

    SAFE_DEALLOCATE_A(this%coeff)

    POP_SUB(mixing_preconditioner_end)
  end subroutine mixing_preconditioner_end

  ! ---------------------------------------------------------
  !This is the model used to fit the Preconditioner in Fourier space
  ! See Eq. (4) of arXiv:1910:07071
  subroutine precondtioner_model(this, x, a, yfit, dyda)
    type(mixing_preconditioner_t), intent(in) :: this
    FLOAT,                         intent(in) :: x(:), a(:)
    FLOAT,                        intent(out) :: yfit(:), dyda(:,:)

    integer :: ip, np, na, ia

    np = size(x)
    na = size(a)

    !Evaluate the fit function and its derivatives
    do ip = 1, np
      yfit(ip) = a(1)
      dyda(ip, 1) = M_ONE
      do ia = 2, na, 2
        yfit(ip) = yfit(ip) + (a(ia) * x(ip)**2)/(x(ip)**2 + a(ia+1))
        dyda(ip, ia)   = (x(ip)**2)/(x(ip)**2 + a(ia+1))
        dyda(ip, ia+1) = -(a(ia) * x(ip)**2)/((x(ip)**2 + a(ia+1))**2)
      end do
    end do

  end subroutine precondtioner_model


  ! ---------------------------------------------------------
  ! In order to fit the Fourier space preconditioner by a sum of fractional functions
  ! we use the Levenberg-Marquardt method, as provided by Numerical Recipes
  ! This is adapted for the need of Octopus
  subroutine mrqmin(this, x, y, a, covar, alpha, chisq, funcs, alamda)
    type(mixing_preconditioner_t), intent(in)    :: this
    FLOAT,                         intent(in)    :: x(:)
    FLOAT,                         intent(in)    :: y(:)
    FLOAT,                         intent(inout) :: a(:)
    FLOAT,                         intent(inout) :: covar(:,:), alpha(:,:)
    FLOAT,                         intent(out)   :: chisq
    FLOAT,                         intent(inout) :: alamda

    interface
      subroutine funcs(this, x, a, yfit, dyda)
        import mixing_preconditioner_t
        type(mixing_preconditioner_t), intent(in) :: this
        FLOAT,                         intent(in) :: x(:), a(:)
        FLOAT,                        intent(out) :: yfit(:), dyda(:,:)
      end subroutine funcs
    end interface

    integer :: ma,ndata

    call mrqmin_private
 
   contains
 
    subroutine mrqmin_private()
      FLOAT, save :: ochisq
      FLOAT, allocatable, save :: atry(:), beta(:)
      FLOAT, allocatable, save :: da(:,:)

      integer :: ia

      ASSERT(size(x) == size(y))
      ndata = size(x)

      ASSERT(size(a) == size(covar,1))
      ASSERT(size(a) == size(covar,2))
      ASSERT(size(a) == size(alpha,1))
      ASSERT(size(a) == size(alpha,2))
      ma = size(a)

      if (alamda < M_ZERO) then
        allocate(atry(ma), beta(ma), da(ma,1))
        alamda = CNST(0.001)
        call mrqcof(this, a, alpha, beta)
        ochisq = chisq
        atry = a
      end if

      covar(1:ma, 1:ma) = alpha(1:ma, 1:ma)
      do ia = 1, ma
        covar(ia, ia) = covar(ia, ia) * (M_ONE+alamda)
      end do
      da(1:ma, 1) = beta(1:ma)

      call gaussj(covar(1:ma, 1:ma), da(1:ma, 1:1))

      if (alamda == M_ZERO) then 
        deallocate(atry, beta, da)
        return
      end if
 
      atry = a + da(:,1) 

      call mrqcof(this, atry, covar, da(1:ma, 1))

      if (chisq < ochisq) then
        alamda = CNST(0.1) * alamda
        ochisq = chisq
        alpha(1:ma, 1:ma) = covar(1:ma, 1:ma)
        beta(1:ma) = da(1:ma, 1)
        a = atry
      else
        alamda = CNST(10.0) * alamda
        chisq = ochisq
      end if
   end subroutine mrqmin_private
 
   subroutine mrqcof(this, a, alpha, beta)
     type(mixing_preconditioner_t), intent(in)  :: this
     FLOAT,                         intent(in)  :: a(:)
     FLOAT,                         intent(out) :: alpha(:,:), beta(:)

     integer :: k,l,m
     FLOAT :: dyda(size(x),size(a))
     FLOAT, dimension(size(x)) :: dy, ymod

     call funcs(this, x, a, ymod, dyda)

     ! Loop over all the data.
     dy = y - ymod
     do l = 1, ma
       do m = 1, l
         alpha(l,m) = dot_product(dyda(:, l), dyda(:,m))
         alpha(m,l) = alpha(l,m) !Symmetric side
       end do
       beta(l) = dot_product(dy, dyda(:, l))
     end do
     !Find chi2
     chisq = dot_product(dy, dy)
 
   end subroutine mrqcof
 end subroutine mrqmin

 ! ---------------------------------------------------------
 !Linear equation solution by Gauss-Jordan elimination,from Numerical recipes
 ! On output,a is replaced by its matrix inverse, and b is replaced by the corresponding 
 ! set of solution vectors.
 SUBROUTINE gaussj(a, b)
   FLOAT, intent(inout) :: a(:,:), b(:,:)

   integer, dimension(size(a,1)) :: ipiv,indxr,indxc
   FLOAT :: pivinv, dum, big
   FLOAT, dimension(size(a,1)) :: dumc
   integer :: i, l, j, k, n, mm
   integer :: irow, icol

   ASSERT(size(a,1) == size(a,2))
   ASSERT(size(a,1) == size(b,1))
   n = size(a,1)
   mm = size(b, 2)

   ipiv=0

   do i = 1, n

     big = M_ZERO
     do j =1, n 
       if (ipiv(j) /= 1) then
         do k = 1, n
           if (ipiv(k) == 0) then
             if (abs(a(j, k)) >= big) then
               big = abs(a(j, k))
               irow = j
               icol = k
             end if
           else if(ipiv(k) > 1) then
             message(1) = 'gaussj: singular matrix (1)'
             call messages_fatal(1)
           end if
         end do
       end if
     end do
     ipiv(icol) = ipiv(icol) + 1

     if (irow /= icol) then
       do l = 1, n
         dum = a(irow, l)
         a(irow, l) = a(icol, l)
         a(icol, l) = dum
       end do
       do l = 1, mm
         dum = b(irow, l)
         b(irow, l) = b(icol, l)
         b(icol, l) = dum
       end do
     end if
 
     indxr(i) = irow
     indxc(i) = icol

     if (a(icol, icol) == M_ZERO) then
       message(1) = 'gaussj: singular matrix (2)'
       call messages_fatal(1)
     end if
 
     pivinv        = M_ONE / a(icol, icol)
     a(icol,icol)  = M_ONE
     a(icol,:)     = a(icol,:) * pivinv
     b(icol,:)     = b(icol,:) * pivinv
     dumc          = a(:,icol)
     a(:,icol)     = M_ZERO
     a(icol,icol)  = pivinv
     a(1:icol-1,:) = a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
     b(1:icol-1,:) = b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
     a(icol+1:,:)  = a(icol+1:,:) -outerprod(dumc(icol+1:),a(icol,:))
     b(icol+1:,:)  = b(icol+1:,:) -outerprod(dumc(icol+1:),b(icol,:))

   end do
 
   do l = n, 1, -1
     if(indxr(l) /= indxc(l)) then
       do k = 1, n
         dum = a(k, indxr(l))
         a(k, indxr(l)) = a(k, indxc(l))
         a(k, indxc(l)) = dum
       end do
     end if
   end do

   contains

     function outerprod(a, b)
       FLOAT, intent(in) :: a(:), b(:)
       FLOAT :: outerprod(size(a), size(b))

       outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
     end function

 end subroutine gaussj

end module mixing_preconditioner_oct_m
  
