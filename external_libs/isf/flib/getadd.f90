!> @file
!! Routines to deal with the address of objects or external functions
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Module used by the module to manage the memory allocations
!! needed to pass to the C routines the correct address
!! in order to take the address of the metadata
module metadata_interfaces
  implicit none

  private

  integer, parameter :: longsize=20              !<could be lower
  character(len=*), parameter :: fmtlong='(i20)' !< conversion of long integer

  interface
     subroutine geti1(array,iadd)
       implicit none
       integer, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti1

     subroutine geti2(array,iadd)
       implicit none
       integer, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti2

     subroutine geti3(array,iadd)
       implicit none
       integer, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti3

     subroutine geti4(array,iadd)
       implicit none
       integer, dimension(:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti4

     !character templates, not the length is added
     subroutine getc1(length,array,iadd)
       implicit none
       integer, intent(in) :: length
       character(len=length), dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getc1

     subroutine getl1(array,iadd)
       implicit none
       logical, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl1

     subroutine getl2(array,iadd)
       implicit none
       logical, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl2

     subroutine getl3(array,iadd)
       implicit none
       logical, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl3

     subroutine getr1(array,iadd)
       implicit none
       real, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getr1

     subroutine getr2(array,iadd)
       implicit none
       real, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getr2

     subroutine getr3(array,iadd)
       implicit none
       real, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getr3

     subroutine getdp1(array,iadd)
       implicit none
       double precision, dimension(:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp1

     subroutine getdp2(array,iadd)
       implicit none
       double precision, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp2

     subroutine getdp3(array,iadd)
       implicit none
       double precision, dimension(:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp3

     subroutine getdp4(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp4

     subroutine getdp5(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp5

     subroutine getdp6(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp6

     subroutine getdp7(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:,:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp7

     subroutine getz2(array,iadd)
       implicit none
       double complex, dimension(:,:), allocatable, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getz2

     subroutine getdp1ptr(array,iadd)
       implicit none
       double precision, dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp1ptr

     subroutine getdp2ptr(array,iadd)
       implicit none
       double precision, dimension(:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp2ptr

     subroutine getdp3ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp3ptr

     subroutine getdp4ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp4ptr

     subroutine getdp5ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp5ptr

     subroutine getdp6ptr(array,iadd)
       implicit none
       double precision, dimension(:,:,:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getdp6ptr

     subroutine geti1ptr(array,iadd)
       implicit none
       integer, dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti1ptr

     subroutine geti2ptr(array,iadd)
       implicit none
       integer, dimension(:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti2ptr

     subroutine geti3ptr(array,iadd)
       implicit none
       integer, dimension(:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti3ptr

     subroutine geti4ptr(array,iadd)
       implicit none
       integer, dimension(:,:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine geti4ptr

     subroutine getl2ptr(array,iadd)
       implicit none
       logical, dimension(:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl2ptr

     subroutine getl3ptr(array,iadd)
       implicit none
       logical, dimension(:,:,:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getl3ptr

     subroutine getc1ptr(length,array,iadd)
       implicit none
       integer, intent(in) :: length
       character(len=length), dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getc1ptr

     subroutine getz1ptr(array,iadd)
       implicit none
       double complex, dimension(:), pointer, intent(in) :: array
       integer(kind=8), intent(out) :: iadd
     end subroutine getz1ptr


  end interface

interface pad_array
  module procedure pad_i1,pad_i2,pad_i3,pad_i4
  module procedure pad_c1
  module procedure pad_l1,pad_l2,pad_l3
  module procedure pad_r1,pad_r2,pad_r3
  module procedure pad_dp1,pad_dp2,pad_dp3,pad_dp4,pad_dp5,pad_dp6,pad_dp7
  module procedure pad_z1, pad_z2
end interface

!>procedure to retrieve the value of the address of the first element of an array
!!this is made necessary from the fact that the treatment of character arrays is 
!!not standard
interface loc_arr
   module procedure la_i1,la_i2,la_i3,la_i4
   module procedure la_r1,la_r2,la_r3
   module procedure la_d1,la_d2,la_d3,la_d4,la_d5,la_d6,la_d7
   module procedure la_l1,la_l2,la_l3
   module procedure la_z1,la_z2
   module procedure la_c1
end interface

public :: pad_array,geti1,geti2,geti3,geti4
public :: getc1
public :: getl1,getl2,getl3
public :: getr1,getr2,getr3
public :: getdp1,getdp2,getdp3,getdp4,getdp5,getdp6,getdp7!,getlongaddress
public :: getz2
public :: getdp1ptr,getdp2ptr,getdp3ptr,getdp4ptr,getdp5ptr,getdp6ptr
public :: geti1ptr,geti2ptr,geti3ptr,geti4ptr
public :: getl2ptr, getl3ptr
public :: getz1ptr
public :: getc1ptr
public :: address_toi,long_toa,loc_arr

contains

  subroutine pad_i1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    integer, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_i1

  subroutine pad_i2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    integer, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_i2

  subroutine pad_i3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    integer, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_i3

  subroutine pad_i4(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(4), intent(in) :: shp
    integer, dimension(shp(1),shp(2),shp(3),shp(4)+ndebug), intent(out) :: array
    
    call pad_integer(array,init_to_zero,product(shp),product(shp(1:3))*(shp(4)+ndebug))

  end subroutine pad_i4

  subroutine pad_c1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    character(len=*), dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_character(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_c1

  subroutine pad_l1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    logical, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_logical(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_l1

  subroutine pad_l2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    logical, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_logical(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_l2

  subroutine pad_l3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    logical, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array

    call pad_logical(array,init_to_zero,product(shp),&
         product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_l3

  subroutine pad_r1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    real, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_simple(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_r1

  subroutine pad_r2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    real, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_simple(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_r2

  subroutine pad_r3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    real, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array

    call pad_simple(array,init_to_zero,product(shp),&
         product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_r3

  subroutine pad_dp1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    double precision, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_dp1

  subroutine pad_dp2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    double precision, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_dp2

  subroutine pad_dp3(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(3), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:2))*(shp(3)+ndebug))

  end subroutine pad_dp3

  subroutine pad_dp4(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(4), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:3))*(shp(4)+ndebug))

  end subroutine pad_dp4

  subroutine pad_dp5(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(5), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4),shp(5)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:4))*(shp(5)+ndebug))

  end subroutine pad_dp5

  subroutine pad_dp6(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(6), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:5))*(shp(6)+ndebug))

  end subroutine pad_dp6

  subroutine pad_dp7(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(7), intent(in) :: shp
    double precision, dimension(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)+ndebug), intent(out) :: array
    
    call pad_double(array,init_to_zero,product(shp),product(shp(1:6))*(shp(7)+ndebug))

  end subroutine pad_dp7

  subroutine pad_z1(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(1), intent(in) :: shp
    double complex, dimension(shp(1)+ndebug), intent(out) :: array
    
    call pad_double_complex(array,init_to_zero,shp(1),shp(1)+ndebug)

  end subroutine pad_z1

  subroutine pad_z2(array,init_to_zero,shp,ndebug)
    implicit none
    logical, intent(in) :: init_to_zero
    integer, intent(in) :: ndebug
    integer, dimension(2), intent(in) :: shp
    double complex, dimension(shp(1),shp(2)+ndebug), intent(out) :: array
    
    call pad_double_complex(array,init_to_zero,product(shp),product(shp(1:1))*(shp(2)+ndebug))

  end subroutine pad_z2

  subroutine pad_double(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    double precision, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) call razero(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)=d_nan()
    end do
  end subroutine pad_double

  subroutine pad_double_complex(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    double complex, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) call razero(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)=(1.d0,1.d0)*d_nan()
    end do
  end subroutine pad_double_complex

  subroutine pad_simple(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    real, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i,i_nan
    real :: r_nan1
    equivalence (r_nan1,i_nan)

    if (init) call razero_simple(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)=r_nan()
    end do
  end subroutine pad_simple

  subroutine pad_logical(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    logical, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) then
       do i=1,ndim_tot
          array(i)=.false.
       end do
    end if
    do i=ndim_tot+1,ndim_extra
       array(i)=.true.
    end do
  end subroutine pad_logical

  subroutine pad_integer(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    integer, dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) call razero_integer(ndim_tot,array)
    do i=ndim_tot+1,ndim_extra
       array(i)= 2147483647 !i_nan
    end do
  end subroutine pad_integer

  subroutine pad_character(array,init,ndim_tot,ndim_extra)
    implicit none
    logical, intent(in) :: init
    integer, intent(in) :: ndim_tot, ndim_extra
    character(len=*), dimension(ndim_extra), intent(out) :: array
    !local variables
    integer :: i

    if (init) then
       do i=1,ndim_tot
          array(i)=repeat(' ',len(array(1)))
       end do
    end if
    do i=ndim_tot+1,ndim_extra
       array(i)=repeat('X',len(array(1)))
    end do
  end subroutine pad_character

  !> Function which specify NaN according to IEEE specifications
  function d_nan()
   implicit none
   double precision :: d_nan
   !local variables
   double precision :: dnan
   integer, dimension(2) :: inan
   equivalence (dnan, inan)
   ! This first assignment is for big-endian machines
   inan(1) = 2147483647
   ! The second assignment is for little-endian machines
   inan(2) = 2147483647
   d_nan = dnan
  end function d_nan

  !> Function which specify NaN according to IEEE specifications
  function r_nan()
   implicit none
   real :: r_nan
   !local variables
   real :: rnan
   integer :: inan
   equivalence (rnan, inan)
   inan = 2147483647
   r_nan = rnan
  end function r_nan

  function exa_toi(a)
    character(len=1), intent(in) :: a
    integer :: exa_toi
    select case(a)
    case('0')
       exa_toi=0
    case('1')
       exa_toi=1
    case('2')
       exa_toi=2
    case('3')
       exa_toi=3
    case('4')
       exa_toi=4
    case('5')
       exa_toi=5
    case('6')
       exa_toi=6
    case('7')
       exa_toi=7
    case('8')
       exa_toi=8
    case('9')
       exa_toi=9
    case('a')
       exa_toi=10
    case('b')
       exa_toi=11
    case('c')
       exa_toi=12
    case('d')
       exa_toi=13
    case('e')
       exa_toi=14
    case('f')
       exa_toi=15
    case default
       !raise an error
       print *,'a: ',a
       stop 'undefined value'
    end select
    
  end function exa_toi

  function address_toi(address)
    character(len=*), intent(in) ::  address
    integer(kind=8) :: address_toi
    !local variables
    integer :: i,l
    integer(kind=8) :: j
    character(len=1) :: a

    l=len_trim(address)
    address_toi=0
    do i=l-2,1,-1
       a=address(i+2:i+2)
       j=int(16**(l-2-i),kind=8)*int(exa_toi(a),kind=8)
       !print *,'i,a',i,a,exa_toi(a)
       address_toi=address_toi+j
    end do
    
  end function address_toi

  pure function long_toa(iadd)
    use yaml_strings
    implicit none 
    integer(kind=8), intent(in) :: iadd
    character(len=longsize) :: long_toa
    
    long_toa=adjustl(yaml_toa(iadd,fmt=fmtlong))

  end function long_toa

  !loc array functions
  function la_i1(array) result(la)
    implicit none
    integer, dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_i1
  function la_i2(array) result(la)
    implicit none
    integer, dimension(:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1))
  end function la_i2
  function la_i3(array) result(la)
    implicit none
    integer, dimension(:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1))
  end function la_i3
  function la_i4(array) result(la)
    implicit none
    integer, dimension(:,:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1,1))
  end function la_i4
  function la_r1(array) result(la)
    implicit none
    real, dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_r1
  function la_r2(array) result(la)
    implicit none
    real, dimension(:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1))
  end function la_r2
  function la_r3(array) result(la)
    implicit none
    real, dimension(:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1))
  end function la_r3
  function la_d1(array) result(la)
    implicit none
    double precision, dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_d1
  function la_d2(array) result(la)
    implicit none
    double precision, dimension(:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1))
  end function la_d2
  function la_d3(array) result(la)
    implicit none
    double precision, dimension(:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1))
  end function la_d3
  function la_d4(array) result(la)
    implicit none
    double precision, dimension(:,:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1,1))
  end function la_d4
  function la_d5(array) result(la)
    implicit none
    double precision, dimension(:,:,:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1,1,1))
  end function la_d5
  function la_d6(array) result(la)
    implicit none
    double precision, dimension(:,:,:,:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1,1,1,1))
  end function la_d6
  function la_d7(array) result(la)
    implicit none
    double precision, dimension(:,:,:,:,:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1,1,1,1,1))
  end function la_d7
  function la_l1(array) result(la)
    implicit none
    logical, dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_l1
  function la_l2(array) result(la)
    implicit none
    logical, dimension(:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1))
  end function la_l2
  function la_l3(array) result(la)
    implicit none
    logical, dimension(:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1))
  end function la_l3
  function la_z1(array) result(la)
    implicit none
    double complex, dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_z1
  function la_z2(array) result(la)
    implicit none
    double complex, dimension(:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1))
  end function la_z2
  function la_z3(array) result(la)
    implicit none
    double complex, dimension(:,:,:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1,1,1))
  end function la_z3
  function la_c1(array) result(la)
    implicit none
    character(len=*), dimension(:), intent(in) :: array
    include 'getadd-c-inc.f90' 
    la=f_loc(array(1))
  end function la_c1


end module metadata_interfaces


!> Routine to call an external routine with an integer argument
subroutine call_external(routine,args)
  implicit none
  external :: routine                  !< Routine to be called
  integer(kind=8), intent(in) :: args  !< Argument of the called routine


  print *,'calling external, args',args
  
  if (args==0) then
     call routine()
  else
     call routine(args)
  end if
end subroutine call_external


!> Call the external routine with no argument
!! to be generalized to the case where extra arguments are needed
recursive subroutine call_external_f(routine)!,args)
  implicit none
  external :: routine                  !< Routine to be called
!  integer(kind=8), intent(in) :: args


!  print *,'calling external, args',args
  
!  if (args==0) then

     call routine()
!  else
!     call routine(args)
!  end if

end subroutine call_external_f


!> Function which identifies the address of the scalar object
!! associated to a unknown quantity
function f_loc(routine)
  implicit none
  external :: routine       !< Object
  integer(kind=8) :: f_loc  !< Address of the object routine

  call getlongaddress(routine,f_loc)

end function f_loc

