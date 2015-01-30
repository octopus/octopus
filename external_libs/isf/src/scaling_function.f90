!> @file
!!  Routines which define and use scaling functions
!! @author
!! Copyright (C) 2002-2007 BigDFT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS 


!> Calculate the values of a scaling function in real uniform grid
subroutine scaling_function(itype,nd,nrange,a,x)

  use Poisson_Solver, only: dp
  use memory_profiling
  use dynamic_memory
  implicit none
  !Arguments
  !Type of interpolating functions
  integer, intent(in) :: itype
  !Number of points: must be 2**nex
  integer, intent(in) :: nd
  integer, intent(out) :: nrange
  real(kind=8), dimension(0:nd), intent(out) :: a,x
  !Local variables
  character(len=*), parameter :: subname='scaling_function'
  real(kind=8), dimension(:), allocatable :: y
  integer :: i,nt,ni
  
  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select
!!  write(unit=*,fmt="(1x,a,i0,a)") &
!!       "Use interpolating scaling functions of ",itype," order"

  !Give the range of the scaling function
  !from -itype to itype
  ni=2*itype
  nrange = ni
  
  y = f_malloc(0.to.nd,id='y')
  
  ! plot scaling function
  call zero(nd+1,x)
  call zero(nd+1,y)
  nt=ni
  x(nt/2-1)=1.d0
  loop1: do
     nt=2*nt
     ! write(6,*) 'nd,nt',nd,nt
     select case(itype)
     case(8)
        call back_trans_8(nd,nt,x,y)
     case(14)
        call back_trans_14(nd,nt,x,y)
     case(16)
        call back_trans_16(nd,nt,x,y)
     case(20)
        call back_trans_20(nd,nt,x,y)
     case(24)
        call back_trans_24(nd,nt,x,y)
     case(30)
        call back_trans_30(nd,nt,x,y)
     case(40)
        call back_trans_40(nd,nt,x,y)
     case(50)
        call back_trans_50(nd,nt,x,y)
     case(60)
        call back_trans_60(nd,nt,x,y)
     case(100)
        call back_trans_100(nd,nt,x,y)
     end select
     !call dcopy(nt,y,1,x,1)
     do i=0,nt-1
        x(i)=y(i)
     end do
     if (nt.eq.nd) then
        exit loop1
     end if
  end do loop1

!  open (unit=1,file='scfunction',status='unknown')
  do i=0,nd
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
!  write(1,*) a(i),x(i)
  end do
!  close(1)

  call f_free(y)
END SUBROUTINE scaling_function


!> Calculate the values of the wavelet function in a real uniform mesh.
subroutine wavelet_function(itype,nd,a,x)

  use Poisson_Solver, only: dp
  use memory_profiling
  use dynamic_memory
  implicit none
  !Arguments
  !Type of the interpolating scaling function
  integer, intent(in) :: itype
  !must be 2**nex
  integer, intent(in) :: nd
  real(kind=8), dimension(0:nd), intent(out) :: a,x
  !Local variables
  character(len=*), parameter :: subname='wavelet_function'
  real(kind=8), dimension(:), allocatable :: y
  integer :: i,nt,ni

  !Only itype=8,14,16,20,24,30,40,50,60,100
  Select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select

  !Give the range of the scaling function
  !from -itype to itype
  ni=2*itype

  y = f_malloc(0.to.nd,id='y')
  
  ! plot wavelet 
  call zero(nd+1,x)
  call zero(nd+1,y)
  nt=ni
  x(nt+nt/2-1)=1.d0
  loop3: do
     nt=2*nt
     !write(6,*) 'nd,nt',nd,nt
     select case(itype)
     case(8)
        call back_trans_8(nd,nt,x,y)
     case(14)
        call back_trans_14(nd,nt,x,y)
     case(16)
        call back_trans_16(nd,nt,x,y)
     case(20)
        call back_trans_20(nd,nt,x,y)
     case(24)
        call back_trans_24(nd,nt,x,y)
     case(30)
        call back_trans_30(nd,nt,x,y)
     case(40)
        call back_trans_40(nd,nt,x,y)
     case(50)
        call back_trans_50(nd,nt,x,y)
     case(60)
        call back_trans_60(nd,nt,x,y)
     case(100)
        call back_trans_100(nd,nt,x,y)
     end select
     !call dcopy(nt,y,1,x,1)
     do i=0,nt-1
        x(i)=y(i)
     end do
     if (nt.eq.nd) then
        exit loop3
     end if
  end do loop3

  !open (unit=1,file='wavelet',status='unknown')
  do i=0,nd-1
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-.5d0)
     !write(1,*) a(i),x(i)
  end do
  !close(1)

  call f_free(y)
 
END SUBROUTINE wavelet_function


!> Do iterations to go from p0gauss to pgauss
!! order interpolating scaling function
subroutine scf_recursion(itype,n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: itype,n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables

  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select

  select case(itype)
  case(8)
     call scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
  case(14)
     call scf_recursion_14(n_iter,n_range,kernel_scf,kern_1_scf)
  case(16)
     call scf_recursion_16(n_iter,n_range,kernel_scf,kern_1_scf)
  case(20)
     call scf_recursion_20(n_iter,n_range,kernel_scf,kern_1_scf)
  case(24)
     call scf_recursion_24(n_iter,n_range,kernel_scf,kern_1_scf)
  case(30)
     call scf_recursion_30(n_iter,n_range,kernel_scf,kern_1_scf)
  case(40)
     call scf_recursion_40(n_iter,n_range,kernel_scf,kern_1_scf)
  case(50)
     call scf_recursion_50(n_iter,n_range,kernel_scf,kern_1_scf)
  case(60)
     call scf_recursion_60(n_iter,n_range,kernel_scf,kern_1_scf)
  case(100)
     call scf_recursion_100(n_iter,n_range,kernel_scf,kern_1_scf)
  end select

END SUBROUTINE scf_recursion


subroutine zero(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n)
  !Local variables
  integer :: i
  do i=1,n
     x(i)=0.d0
  end do
END SUBROUTINE zero


!> Forward wavelet transform
!!   m filter length (m has to be even!)
subroutine for_trans_8(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_8.inc'

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

END SUBROUTINE for_trans_8


!> Backward wavelet transform
subroutine back_trans_8(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_8.inc'
  
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
        
END SUBROUTINE back_trans_8

!> Tests the 4 orthogonality relations of the filters
subroutine ftest_8
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  integer :: i,j,l
  real(kind=8) :: t1,t2,t3,t4,eps

  include 'lazy_8.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_8


!> Do iterations to go from p0gauss to pgauss
!! 8th-order interpolating scaling function
subroutine scf_recursion_8(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_8.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_8


!> Forward wavelet transform
subroutine for_trans_14(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_14.inc'

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
  
END SUBROUTINE for_trans_14


!> Backward wavelet transform
subroutine back_trans_14(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_14.inc'
  
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
        
END SUBROUTINE back_trans_14


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_14
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_14.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_14


!> Do iterations to go from p0gauss to pgauss
!! 14th-order interpolating scaling function
subroutine scf_recursion_14(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_14.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_14


!> Forward wavelet transform
subroutine for_trans_16(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_16.inc'

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
  
END SUBROUTINE for_trans_16


!> Backward wavelet transform
subroutine back_trans_16(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_16.inc'
  
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
        
END SUBROUTINE back_trans_16


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_16
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_16.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_16


!> Do iterations to go from p0gauss to pgauss
!! 16th-order interpolating scaling function
subroutine scf_recursion_16(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_16.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_16


!> Forward wavelet transform
subroutine for_trans_20(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_20.inc'

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
  
END SUBROUTINE for_trans_20


!> Backward wavelet transform
subroutine back_trans_20(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_20.inc'
  
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
        
END SUBROUTINE back_trans_20



!> Tests the 4 orthogonality relations of the filters
subroutine ftest_20
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_20.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_20


!> Do iterations to go from p0gauss to pgauss
!! 20th-order interpolating scaling function
subroutine scf_recursion_20(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_20.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_20


!> Forward wavelet transform
subroutine for_trans_24(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_24.inc'

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
  
END SUBROUTINE for_trans_24


!> Backward wavelet transform
subroutine back_trans_24(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_24.inc'
  
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
        
END SUBROUTINE back_trans_24


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_24
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_24.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_24


!> Do iterations to go from p0gauss to pgauss
!! 24th-order interpolating scaling function
subroutine scf_recursion_24(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_24.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_24


!> Forward wavelet transform
subroutine for_trans_30(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_30.inc'

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
  
END SUBROUTINE for_trans_30


!> Backward wavelet transform
subroutine back_trans_30(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_30.inc'
  
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
        
END SUBROUTINE back_trans_30


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_30()
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_30.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_30


!> Do iterations to go from p0gauss to pgauss
!! 30th-order interpolating scaling function
subroutine scf_recursion_30(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_30.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_30


!> Forward wavelet transform
subroutine for_trans_40(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_40.inc'

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
  
END SUBROUTINE for_trans_40


!> Backward wavelet transform
subroutine back_trans_40(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_40.inc'
  
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
        
END SUBROUTINE back_trans_40


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_40
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_40.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_40



!> Do iterations to go from p0gauss to pgauss
!! 40th-order interpolating scaling function
subroutine scf_recursion_40(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_40.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_40


!> Forward wavelet transform
subroutine for_trans_50(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_50.inc'

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
  
END SUBROUTINE for_trans_50


!> Backward wavelet transform
subroutine back_trans_50(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_50.inc'
  
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
        
END SUBROUTINE back_trans_50


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_50
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_50.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_50


!> Do iterations to go from p0gauss to pgauss
!! 50th-order interpolating scaling function
subroutine scf_recursion_50(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_50.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_50


!> Forward wavelet transform
subroutine for_trans_60(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_60.inc'

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
  
END SUBROUTINE for_trans_60


!> Backward wavelet transform
subroutine back_trans_60(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_60.inc'
  
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
        
END SUBROUTINE back_trans_60


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_60
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_60.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_60


!> Do iterations to go from p0gauss to pgauss
!! 60th-order interpolating scaling function
subroutine scf_recursion_60(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_60.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_60


!> Forward wavelet transform
subroutine for_trans_100(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_100.inc'

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
  
END SUBROUTINE for_trans_100


!> Backward wavelet transform
subroutine back_trans_100(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd !< Length of data set
  integer, intent(in) :: nt !< Length of data in data set to be transformed
  real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
  real(kind=8), intent(out) :: y(0:nd-1) !< Output data
  !Local variables
  integer :: i,j,ind

  include 'lazy_100.inc'
  
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
        
END SUBROUTINE back_trans_100


!> Tests the 4 orthogonality relations of the filters
subroutine ftest_100
  implicit none
  !Arguments
  !Local variables
  character(len=*), parameter :: fmt22 = "(a,i3,i4,4(e17.10))"
  real(kind=8) :: t1,t2,t3,t4,eps
  integer :: i,j,l

  include 'lazy_100.inc'
  
  ! do i=-m,m
  ! write(6,*) i,ch(i),cg(i)
  ! end do
  
  do i=-m,m
     do j=-m,m
        t1=0.d0
        t2=0.d0
        t3=0.d0
        t4=0.d0
        do l=-3*m,3*m
           if ( l-2*i.ge.-m .and. l-2*i.le.m  .and. &
                l-2*j.ge.-m .and. l-2*j.le.m ) then
              t1=t1+ch(l-2*i)*cht(l-2*j)
              t2=t2+cg(l-2*i)*cgt(l-2*j)
              t3=t3+ch(l-2*i)*cgt(l-2*j)
              t4=t4+cht(l-2*i)*cg(l-2*j)
           end if
        end do
        eps=1.d-10
        if (i.eq.j) then
           if (abs(t1-1.d0).gt.eps .or. abs(t2-1.d0).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then 
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        else
           if (abs(t1).gt.eps .or. abs(t2).gt.eps .or. &
             & abs(t3).gt.eps  .or. abs(t4).gt.eps ) then
              write(6,fmt22) 'Orthogonality ERROR', i,j,t1,t2,t3,t4
           end if
        end if
     end do
  end do
  
  write(6,*) 'FILTER TEST PASSED'
  
END SUBROUTINE ftest_100


!> Do iterations to go from p0gauss to pgauss
!! 100th-order interpolating scaling function
subroutine scf_recursion_100(n_iter,n_range,kernel_scf,kern_1_scf)
  implicit none
  !Arguments
  integer, intent(in) :: n_iter,n_range
  real(kind=8), intent(inout) :: kernel_scf(-n_range:n_range)
  real(kind=8), intent(out) :: kern_1_scf(-n_range:n_range)
  !Local variables
  real(kind=8) :: kern,kern_tot
  integer :: i_iter,i,j,ind

  include "lazy_100.inc"

  !Start the iteration to go from p0gauss to pgauss
  loop_iter_scf: do i_iter=1,n_iter
     !kern_1_scf(:) = kernel_scf(:)
     !kernel_scf(:) = 0.d0
     call dcopy(2*n_range+1,kernel_scf(-n_range),1,kern_1_scf(-n_range),1)
     call dscal(2*n_range+1,0.0d0,kernel_scf(-n_range),1)
     loop_iter_i: do i=0,n_range
        kern_tot = 0.d0
        do j=-m,m
           ind = 2*i-j
           !if (abs(ind) > n_range) then
           if (ind > n_range .or. ind < -n_range ) then
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
           kernel_scf(-i) = 0.5d0*kern_tot
        end if
     end do loop_iter_i
  end do loop_iter_scf
END SUBROUTINE scf_recursion_100
