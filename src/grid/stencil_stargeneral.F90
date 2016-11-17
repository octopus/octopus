!! Copyright (C) 2016 U. De Giovannini, H Huebener
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
!! $Id: stencil_stargeneral.F90 10978 2013-07-11 15:28:46Z micael $

#include "global.h"

module stencil_stargeneral_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use stencil_oct_m

  implicit none

  private
  public ::                     &
    stencil_stargeneral_size_lapl, &
    stencil_stargeneral_extent,    &
    stencil_stargeneral_get_lapl,  &
    stencil_stargeneral_pol_lapl,  &
    stencil_stargeneral_get_arms



contains
  
  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_get_arms(this, sb, spacing)
    type(stencil_t),     intent(inout) :: this 
    type(simul_box_t),      intent(in) :: sb
    FLOAT, intent(in)   :: spacing(sb%dim)
    integer :: idim, dim
    FLOAT   :: vec1(1:3), vec2(1:3), theta, arm(1:3), ratio
    integer :: n1, n2, n3
    
    PUSH_SUB(stencil_stargeneral_get_arms)  

    dim = sb%dim    
       
    vec1(:) = M_ZERO
    vec2(:) = M_ZERO
    
    this%stargeneral%narms = 0

    if (dim == 1 ) then 
      !we are done 
      POP_SUB(stencil_stargeneral_get_arms)      
      return 
    end if   
    
    vec1(1:dim)=sb%rlattice_primitive(1:dim, 1)
    vec2(1:dim)=sb%rlattice_primitive(1:dim, 2)
    !get the angle between the primitive vectors
    theta = acos(dot_product(vec1(1:dim),vec2(1:dim)))

! find a reasonable direction for arm 1
    ratio = (M_TWO*sb%lsize(1)/spacing(1)) / (M_TWO*sb%lsize(2)/spacing(2))
    if(ratio > M_ONE)then
      ratio = ratio + M_HALF
      n1 = aint(ratio); n2 = 1
    else
      ratio = M_ONE/ratio + M_HALF
      n1 = 1; n2 = aint(ratio)
    end if

    if (theta < M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      arm(1:3) = (/n1,-n2,0/)
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = arm(1:dim)
    else if(theta > M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      arm(1:3) = (/n1,+n2,0/)
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = arm(1:dim)
    end if
    !if theta == pi/2 we do not need additional arms
    
! NB: you have supposed the axis of the 2D system is along z.
    if (dim == 2 ) then 
      !we are done 
      POP_SUB(stencil_stargeneral_get_arms)      
      return 
    end if   

    ! dim>2
    
    vec1(1:dim)=sb%rlattice_primitive(1:dim, 2)
    vec2(1:dim)=sb%rlattice_primitive(1:dim, 3)
    !get the angle between the primitive vectors
    theta = acos(dot_product(vec1(1:dim),vec2(1:dim)))

! find a reasonable direction for arm 2
    ratio = (M_TWO*sb%lsize(2)/spacing(2)) / (M_TWO*sb%lsize(3)/spacing(3))
    if(ratio > M_ONE)then
      ratio = ratio + M_HALF
      n2 = aint(ratio); n3 = 1
    else
      ratio = M_ONE/ratio + M_HALF
      n2 = 1; n3 = aint(ratio)
    end if

    if (theta < M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = (/0,n2,-n3/)
    else if(theta > M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = (/0,n2,+n3/)
    end if
    !if theta == pi/2 we do not need additional arms
      
    vec1(1:dim)=sb%rlattice_primitive(1:dim, 3)
    vec2(1:dim)=sb%rlattice_primitive(1:dim, 1)
    !get the angle between the primitive vectors
    theta = acos(dot_product(vec1(1:dim),vec2(1:dim)))

! find a reasonable direction for arm 2
    ratio = (M_TWO*sb%lsize(3)/spacing(3)) / (M_TWO*sb%lsize(1)/spacing(1))
    if(ratio > M_ONE)then
      ratio = ratio + M_HALF
      n3 = aint(ratio); n1 = 1
    else
      ratio = M_ONE/ratio + M_HALF
      n3 = 1; n1 = aint(ratio)
    end if

    if (theta < M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = (/-n1,0,n3/)
    else if(theta > M_PI*M_HALF) then
      this%stargeneral%narms = this%stargeneral%narms + 1
      this%stargeneral%arms(this%stargeneral%narms, 1:dim) = (/+n1,0,n3/)
    end if
    !if theta == pi/2 we do not need additional arms

      
      
    POP_SUB(stencil_stargeneral_get_arms)      
  end subroutine stencil_stargeneral_get_arms


  ! ---------------------------------------------------------
  integer function stencil_stargeneral_size_lapl(this, dim, order) result(n)
    type(stencil_t),     intent(inout) :: this 
    integer, intent(in) :: dim
    integer, intent(in) :: order

    PUSH_SUB(stencil_stargeneral_size_lapl)

    !normal star
    n = 2*dim*order + 1

    ! star general 
    n = n + 2 * order * this%stargeneral%narms 
    

    POP_SUB(stencil_stargeneral_size_lapl)
  end function stencil_stargeneral_size_lapl


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given discretization order.
  integer function stencil_stargeneral_extent(dir, order)
    integer, intent(in) :: dir
    integer, intent(in) :: order

    integer :: extent

    PUSH_SUB(stencil_stargeneral_extent)

    extent = 0
    if(dir >= 1.or.dir <= 3) then
      if(order <= 2) then
        extent = 2
      else
        extent = order
      end if
    end if
    stencil_stargeneral_extent = extent

    POP_SUB(stencil_stargeneral_extent)
  end function stencil_stargeneral_extent



  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_get_lapl(this, dim, order)
    type(stencil_t), intent(out) :: this
    integer,         intent(in)  :: dim
    integer,         intent(in)  :: order

    integer :: i, j, n
    logical :: got_center

    PUSH_SUB(stencil_stargeneral_get_lapl)

    call stencil_allocate(this, stencil_stargeneral_size_lapl(this, dim, order))

    n = 1
    select case(dim)
    case(1)
      n = 1
      do i = 1, dim
        do j = -order, order
          if(j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
    case(2)
      n = 1
      do i = 1, dim
        do j = -order, order
          if(j == 0) cycle
          n = n + 1
          this%points(i, n) = j
        end do
      end do
      
      do j = -order, order
        if(j == 0) cycle
        do i = 1, this%stargeneral%narms 
          n = n + 1
          this%points(1:2, n) = this%stargeneral%arms(i, 1:2)*j
        end do 
      end do
      
    case(3)
      got_center = .false.
      
      n = 0
      do i = 1, dim
        do j = -order, order
          
          ! count center only once
          if(j == 0) then
            if(got_center) then
              cycle
            else
              got_center = .true.
            end if

          end if
          n = n + 1
          this%points(i, n) = j
        end do
      end do

      do j = -order, order
        if(j == 0) cycle
        do i = 1, this%stargeneral%narms 
          n = n + 1
          this%points(1:3, n) = this%stargeneral%arms(i, 1:3)*j
        end do 
      end do

    end select

    call stencil_init_center(this)

    POP_SUB(stencil_stargeneral_get_lapl)
  end subroutine stencil_stargeneral_get_lapl




  ! ---------------------------------------------------------
  subroutine stencil_stargeneral_pol_lapl(this, dim, order, pol)
    type(stencil_t), intent(in) :: this
    integer, intent(in)          :: dim
    integer, intent(in)          :: order
    integer, intent(out)         :: pol(:,:) !< pol(dim, order)

    integer :: i, j, n, j1, j2

    PUSH_SUB(stencil_stargeneral_pol_lapl)

    n = 1
    select case(dim)
    case(1)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do
    case(2)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do

      do j = 1, 2*order
        do i = 1, this%stargeneral%narms 
          n = n + 1
          if (sum(this%stargeneral%arms(i,1:dim))==0 )then          
            pol(1:2, n) = (/j,1/)          
          else  
            pol(1:2, n) = (/1,j/)          
          end if  
        end do
      end do
      
    case(3)
      n = 1
      pol(:,:) = 0
      do i = 1, dim
        do j = 1, 2*order
          n = n + 1
          pol(i, n) = j
        end do
      end do

      do i = 1, this%stargeneral%narms 

        if (this%stargeneral%arms(i,1)==0) then
          ! sum(this%stargeneral%arms(i,1:dim))==0 just checks whether we have a -1 in the arm vector or not
          if(this%stargeneral%arms(i,2)*this%stargeneral%arms(i,3) < 0 )then
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/0,j1,j2/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/0,j2,j1/)
              end do
            end do
          end if
        end if

        if (this%stargeneral%arms(i,2)==0) then
          ! sum(this%stargeneral%arms(i,1:dim))==0 just checks whether we have a -1 in the arm vector or not
          if(this%stargeneral%arms(i,3)*this%stargeneral%arms(i,1) < 0 )then
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/j2,0,j1/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/j1,0,j2/)
              end do
            end do
          end if
        end if

        if (this%stargeneral%arms(i,3)==0) then
          ! sum(this%stargeneral%arms(i,1:dim))==0 just checks whether we have a -1 in the arm vector or not
          if(this%stargeneral%arms(i,1)*this%stargeneral%arms(i,2) < 0 )then
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/j1,j2,0/)
              end do
            end do
          else
            do j1 = 1, 2*order
              do j2 = 1, 2*order-j1
                n = n + 1
                pol(1:3, n) = (/j2,j1,0/)
              end do
            end do
          end if
        end if
      end do ! end loop on number of arms

!       !FCC
!       do j = 1, 2*order
!         n = n + 1
!         pol(1:3, n) = (/j,1,0/)
!         n = n + 1
!         pol(1:3, n) = (/1,0,j/)
!         n = n + 1
!         pol(1:3, n) = (/0,j,1/)
!       end do
      
!       !HEX
!       do j = 1, 2*order
!         n = n + 1
!         pol(1:3, n) = (/j, 1, 0/)
!       end do
      

    end select

    POP_SUB(stencil_stargeneral_pol_lapl)
  end subroutine stencil_stargeneral_pol_lapl


end module stencil_stargeneral_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
