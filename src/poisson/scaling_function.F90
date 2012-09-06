!! Copyright (C) 2011 X. Andrade
!! Copyright (C) Luigi Genovese, Thierry Deutsch, CEA Grenoble, 2006
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
!! $Id: scaling_function.F90 2544 2006-11-03 17:41:04Z xavier $

#include <global.h>

module scaling_function_m
  use blas_m
  use global_m
  use messages_m

  implicit none

  private

  public ::            &
    scaling_function,  &
    scf_recursion
    
contains

  !!****h* BigDFT/scaling_function
  !! NAME
  !!   scaling_function
  !!
  !! FUNCTION
  !!   Calculate the values of a scaling function in real uniform grid
  !!
  !! SOURCE
  !!
  subroutine scaling_function(itype,nd,nrange,a,x)
    integer, intent(in) :: itype
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(kind=8), dimension(0:nd), intent(out) :: a,x

    real(kind=8), dimension(:), allocatable :: y
    integer :: i,nt,ni

    !Only itype=8,14,16,20,24,30,40,50,60,100
    select case(itype)
    case(8)
      !O.K.
    case default
      message(1) = "Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100."
      call messages_fatal(1)
    end select
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

    !Give the range of the scaling function
    !from -itype to itype
    ni=2*itype
    nrange = ni
    allocate(y(0:nd))

    ! plot scaling function
    x = 0.0_8
    y = 0.0_8

    nt=ni
    x(nt/2-1)=1.d0
    loop1: do
      nt=2*nt
      !	write(6,*) 'nd,nt',nd,nt
      select case(itype)
      case(8)
        call back_trans_8(nd,nt,x,y)
      end select
      call blas_copy(nt, y(0), 1, x(0) ,1)
      if (nt.eq.nd) then
        exit loop1
      end if
    end do loop1

    !open (unit=1,file='scfunction',status='unknown')
    do i=0,nd
      a(i) = 1.d0*i*ni/nd-(.5d0*ni-1.d0)
      !write(1,*) 1.d0*i*ni/nd-(.5d0*ni-1.d0),x(i)
    end do
    !close(1)

    deallocate(y)
  end subroutine scaling_function
  !!***

  !!****h* BigDFT/scf_recursion
  !! NAME
  !!   scf_recursion
  !!
  !! FUNCTION
  !!   Do iterations to go from p0gauss to pgauss
  !!   order interpolating scaling function
  !!
  !! SOURCE
  !!
  subroutine scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
    integer, intent(in) :: itype,n_iter,n_range
    real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
    real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)

    !Only itype=8,14,16,20,24,30,40,50,60,100
    select case(itype)
    case(8)
      !O.K.
    case default
      message(1) = "Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100."
      call messages_fatal(1)
    end select

    select case(itype)
    case(8)
      call scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
    end select

  end subroutine scf_recursion
  !!***

  !!****h* BigDFT/for_trans_8
  !! NAME
  !!   for_trans_8
  !!
  !! FUNCTION
  !!   forward wavelet transform
  !!   nd: length of data set
  !!   nt length of data in data set to be transformed
  !!   m filter length (m has to be even!)
  !!   x input data, y output data
  !!
  !! SOURCE
  !!
  subroutine for_trans_8(nd,nt,x,y)
    integer, intent(in) :: nd,nt
    real(kind=8), intent(in) :: x(0:nd-1)
    real(kind=8), intent(out) :: y(0:nd-1)

    !Local variables
    integer :: i,j,ind

#include "lazy_8_inc.F90"

    do i=0,nt/2-1
      y(     i)=0.d0
      y(nt/2+i)=0.d0

      do j=-m+1,m

        ! periodically wrap index if necessary
        ind=j+2*i
        loop99: do
          if (ind.lt.0) then 
            ind=ind+nt
            cycle loop99
          end if
          if (ind.ge.nt) then 
            ind=ind-nt
            cycle loop99
          end if
          exit loop99
        end do loop99

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
      end do

    end do

  end subroutine for_trans_8
  !!***


  !!****h* BigDFT/back_trans_8
  !! NAME
  !!   back_trans_8
  !!
  !! FUNCTION
  !!
  !! SOURCE
  !!
  ! backward wavelet transform
  ! nd: length of data set
  ! nt length of data in data set to be transformed
  ! m filter length (m has to be even!)
  ! x input data, y output data
  subroutine back_trans_8(nd,nt,x,y)
    integer, intent(in) :: nd,nt
    real(kind=8), intent(in) :: x(0:nd-1)
    real(kind=8), intent(out) :: y(0:nd-1)

    integer :: i,j,ind

#include "lazy_8_inc.F90"

    do i=0,nt/2-1
      y(2*i+0)=0.d0
      y(2*i+1)=0.d0

      do j=-m/2,m/2-1

        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
          if (ind.lt.0) then 
            ind=ind+nt/2
            cycle loop99
          end if
          if (ind.ge.nt/2) then 
            ind=ind-nt/2
            cycle loop99
          end if
          exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
      end do

    end do

  end subroutine back_trans_8
  !!***
  
  !!****h* BigDFT/scf_recursion_8
  !! NAME
  !!   scf_recursion_8
  !!
  !! FUNCTION
  !!   Do iterations to go from p0gauss to pgauss
  !!   8th-order interpolating scaling function
  !!
  !! SOURCE
  !!
  subroutine scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
    integer, intent(in) :: n_iter,n_range
    real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
    real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)

    real(kind=8) :: kern,kern_tot
    integer :: i_iter,i,j,ind

#include "lazy_8_inc.F90"

    !Start the iteration to go from p0gauss to pgauss
    loop_iter_scf: do i_iter=1,n_iter
      kern_1_scf(:) = kernel_scf(:)
      kernel_scf(:) = 0.d0
      loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
          ind = 2*i-j
          if (abs(ind) > n_range) then
            kern = 0.d0
          else
            kern = kern_1_scf(ind)
          end if
          kern_tot = kern_tot + ch(j)*kern
        end do
        if (kern_tot == 0.d0) then
          !zero after (be sure because strictly == 0.d0)
          exit loop_iter_i
        else
          kernel_scf( i) = 0.5d0*kern_tot
          kernel_scf(-i) = kernel_scf(i)
        end if
      end do loop_iter_i
    end do loop_iter_scf
  end subroutine scf_recursion_8
  !!***
  
end module scaling_function_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
