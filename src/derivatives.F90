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

module derivatives
  use mesh

  implicit none

  type der_discr_type
    type(mesh_type), pointer :: m             ! pointer to the underlying mesh
    integer                  :: order         ! order of the discretization (value depends on stencil)
    integer                  :: stencil_type  ! type of discretization
    integer                  :: stencil_np    ! number of points in stencil
    integer,         pointer :: stencil(:,:)  ! the stencil
  end type der_discr_type

  integer, parameter :: &
     DER_STAR = 0

contains
  subroutine derivatives_init(der)
    type(der_discr_type), intent(out) :: der

    call push_sub('derivatives_init')

    call loct_parse_int('DerivativesStencil', DER_STAR, der%stencil_type)
    call loct_parse_int('DerivativesOrder', 4, der%order)

    call derivatives_get_stencil(der)

    call pop_sub()
  end subroutine derivatives_init


  subroutine derivatives_end(der)
    type(der_discr_type), intent(inout) :: der
    call push_sub('derivatives_init')

    ASSERT(associated(der%stencil))

    deallocate(der%stencil); nullify(der%stencil)

    call pop_sub()
  end subroutine derivatives_end


  subroutine derivatives_get_stencil(der)
    type(der_discr_type), intent(inout) :: der
    
    call push_sub('der_get_stencil')
    
    ASSERT(.not.associated(der%stencil))

    select case(der%stencil_type)
    case(DER_STAR)
      call star()
    end select

    call pop_sub()

  contains
    subroutine star
      integer :: i, j, n
      
      der%stencil_np = 2*conf%dim*der%order + 1
      allocate(der%stencil(3, der%stencil_np))
      der%stencil(:,:) = 0
      n = 1
      do i = 1, conf%dim
        do j = -der%order, der%order
          if(j == 0) cycle
          n = n + 1
          der%stencil(i, n) = j
        end do
      end do

    end subroutine star
  end subroutine derivatives_get_stencil

  subroutine derivatives_build(der)
    type(der_discr_type), intent(inout) :: der

  end subroutine derivatives_build

end module derivatives
