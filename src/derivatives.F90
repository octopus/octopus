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
  use math, only: weights
  use mesh
  use nl_operator

  implicit none

  type der_discr_type
    type(mesh_type), pointer :: m             ! pointer to the underlying mesh
    integer                  :: order         ! order of the discretization (value depends on stencil)
    integer                  :: stencil_type  ! type of discretization

    integer                  :: boundaries(3) ! bounday conditions
    logical                  :: zero_bc

    type(nl_operator_type), pointer :: lapl
    type(nl_operator_type), pointer :: grad(:)
  end type der_discr_type

  integer, parameter ::    &
     DER_BC_ZERO_F  = 0,   &  ! function si zero at the boundaries
     DER_BC_ZERO_DF = 1,   &  ! first derivative of the function is zero
     DER_BC_PERIOD  = 2       ! boundary is periodic

  integer, parameter ::    &
     DER_STAR        = 0,  &
     DER_VARIATIONAL = 1

contains

  ! ---------------------------------------------------------
  subroutine derivatives_init(m, der, n_ghost)
    type(mesh_type),      pointer     :: m
    type(der_discr_type), intent(out) :: der
    integer,              intent(out) :: n_ghost(:)

    integer :: i

    call push_sub('derivatives_init')

    der%m => m    ! make a pointer to the underlying mesh

    call loct_parse_int('DerivativesStencil', DER_STAR, der%stencil_type)
    call loct_parse_int('DerivativesOrder', 4, der%order)

    ! construct lapl and grad structures
    call derivatives_get_stencil_lapl(der)
    call derivatives_get_stencil_grad(der)

    ! find out the bounday conditions
    call loct_parse_int('DerivativesBoundaries', DER_BC_ZERO_F, i)
    if((i.ne.DER_BC_ZERO_F).and.(i.ne.DER_BC_ZERO_DF)) then
      write(message(1), '(a,i2,a)') 'DerivativesBoundaries = "', i, '" is unknown to octopus'
      call write_fatal(1)
    end if
    
    der%zero_bc = (i == DER_BC_ZERO_F)
    der%boundaries(:) = i
    der%boundaries(1:conf%periodic_dim) = DER_BC_PERIOD

    ! find out how many ghost points we need in which dimension
    n_ghost(:) = 0
    do i = 1, conf%dim
      if(der%boundaries(i) == DER_BC_ZERO_F) then
!        n_ghost(i) = max(maxval(abs(der%grad(i)%stencil(i,:))), maxval(abs(der%lapl%stencil(i,:))))
        n_ghost(i) = maxval(abs(der%lapl%stencil(i,:)))
      end if
    end do

    call pop_sub()
  end subroutine derivatives_init


  ! ---------------------------------------------------------
  subroutine derivatives_end(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i

    call push_sub('derivatives_init')

    ASSERT(associated(der%lapl))
    ASSERT(associated(der%grad))

    call nl_operator_end(der%lapl)
    do i = 1, conf%dim
      call nl_operator_end(der%grad(i))
    end do

    deallocate(der%lapl, der%grad)
    nullify   (der%lapl, der%grad, der%m)

    call pop_sub()
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der)
    type(der_discr_type), intent(inout) :: der

    integer :: n
    
    ASSERT(.not.associated(der%lapl))

    call push_sub('der_get_stencil_lapl')
    
    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      n = 2*conf%dim*der%order + 1
    end select
    
    ! initialize nl operator
    allocate(der%lapl)
    call nl_operator_init(der%lapl, n)
    
    ! create stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      call star(der%lapl%stencil)
    end select

    call pop_sub()

  contains
    subroutine star(stencil)
      integer, intent(out) :: stencil(:,:)
      integer :: i, j, n
      
      stencil(:,:) = 0
      n = 1
      do i = 1, conf%dim
        do j = -der%order, der%order
          if(j == 0) cycle
          n = n + 1
          stencil(i, n) = j
        end do
      end do

    end subroutine star
  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_grad(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i, n
    
    ASSERT(.not.associated(der%grad))

    call push_sub('der_get_stencil_grad')
    
    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      n = 2*der%order + 1
    end select
    
    ! initialize nl operator
    allocate(der%grad(conf%dim))
    do i = 1, conf%dim
      call nl_operator_init(der%grad(i), n)

      ! create stencil
      select case(der%stencil_type)
      case(DER_STAR, DER_VARIATIONAL)
        call star(i, der%grad(i)%stencil)
      end select
    end do

    call pop_sub()

  contains
    subroutine star(i, stencil)
      integer, intent(in)  :: i
      integer, intent(out) :: stencil(:,:)
      integer :: j, n
      
      stencil(:,:) = 0
      n = 1
      do j = -der%order, der%order
        stencil(i, n) = j
        n = n + 1
      end do

    end subroutine star
  end subroutine derivatives_get_stencil_grad


  ! ---------------------------------------------------------
  subroutine derivatives_build(der)
    type(der_discr_type), intent(inout) :: der

    call derivatives_build_lapl(der)
    call derivatives_build_grad(der)
    
  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine derivatives_build_lapl(der)
    type(der_discr_type), intent(inout) :: der

    ASSERT(associated(der%lapl))

    select case(der%stencil_type)
    case(DER_STAR)  
      call nl_operator_build(der%m, der%lapl, der%m%np, .true.)
      call star(der%lapl)
    case(DER_VARIATIONAL)
      call nl_operator_build(der%m, der%lapl, der%m%np, .true.)
      call variational(der%lapl)
    end select

  contains
    subroutine star(lapl)
      type(nl_operator_type), intent(inout) :: lapl

      integer :: k, i, j, morder
      FLOAT, allocatable :: cc(:,:,:)

      ASSERT(lapl%n >= 3)

      morder = (lapl%n - 1)/conf%dim
      allocate(cc(0:morder, 0:morder, 0:2))
      call weights(2, morder, cc)
      lapl%w(1,1) = cc(0, morder, 2)*sum(1/der%m%h(1:conf%dim)**2)
      
      k = 1
      do i = 1, conf%dim
        do j = -der%order, -1
          k = k + 1
          lapl%w(k,1) = cc(-2*j-1, morder, 2) / der%m%h(i)**2
        end do

        do j = 1, der%order
          k = k + 1
          lapl%w(k,1) = cc(2*j,   morder, 2) / der%m%h(i)**2
        end do
      end do

      deallocate(cc)
    end subroutine star

    ! ----------------------------------------
    ! Implements the variational discretization of the laplacian
    ! as proposed by P. Maragakis, J. Soler, and E. Kaxiras, PRB 64, 193101 (2001)
    ! ----------------------------------------
    subroutine variational(lapl)
      type(nl_operator_type), intent(inout) :: lapl

      integer :: i, j, k
      FLOAT, parameter   :: pi2 = M_PI**2
      FLOAT, allocatable :: fp(:)

      allocate(fp(der%order+1))
      select case(der%order)
      case(1) 
        fp(1) = -pi2/M_TWO
        fp(2) =  pi2/M_FOUR
      case(2) 
        fp(1) = -M_HALF-M_THREE*pi2/M_EIGHT
        fp(2) =  pi2/M_FOUR
        fp(3) =  M_ONE/M_FOUR - pi2/CNST(16.0)
      case(3)
        fp(1) = -M_FIVE/M_SIX - M_FIVE*pi2/CNST(16.0)
        fp(2) =  M_ONE/CNST(12.0) + CNST(15.0)*pi2/CNST(64.0)
        fp(3) =  M_FIVE/CNST(12.0) - M_THREE*pi2/CNST(32.0)
        fp(4) = -M_ONE/CNST(12.0) + pi2/CNST(64.0)
      case(4)
        fp(1) = -CNST(77.0)/CNST(72.0) - CNST(35.0)*pi2/CNST(128.0)
        fp(2) =  M_EIGHT/CNST(45.0) + M_SEVEN*pi2/CNST(32.0)
        fp(3) =  CNST(23.0)/CNST(45.0) - M_SEVEN*pi2/CNST(64.0)
        fp(4) = -M_EIGHT/CNST(45.0) + pi2/CNST(32.0)
        fp(5) =  CNST(17.0)/CNST(720.0) - pi2/CNST(256.0)
      case(5)
        fp(1) = -CNST(449.0)/CNST(360.0) - CNST(63.0)*pi2/CNST(256.0)
        fp(2) =  M_FOUR/CNST(15.0) + CNST(105.0)*pi2/CNST(512.0)
        fp(3) =  CNST(59.0)/CNST(105.0) - CNST(15.0)*pi2/CNST(128.0)
        fp(4) = -CNST(82.0)/CNST(315.0) + CNST(45.0)*pi2/CNST(1024.0)
        fp(5) =  CNST(311.0)/CNST(5040.0) - M_FIVE*pi2/CNST(512.0)
        fp(6) = -M_TWO/CNST(315.0) + pi2/CNST(1024.0)
      case(6)
        fp(1) = -CNST(2497.0)/CNST(1800.0) - CNST(231.0)*pi2/CNST(1024.0)
        fp(2) =  CNST(26.0)/CNST(75.0) + CNST(99.0)*pi2/CNST(512.0)
        fp(3) =  CNST(493.0)/CNST(840.0) - CNST(495.0)*pi2/CNST(4096.0)
        fp(4) = -CNST(103.0)/CNST(315.0) + CNST(55.0)*pi2/CNST(1024.0)
        fp(5) =  CNST(2647.0)/CNST(25200.0) - CNST(33.0)*pi2/CNST(2048.0)
        fp(6) = -CNST(31.0)/CNST(1575.0) + M_THREE*pi2/CNST(1024.0)
        fp(7) =  M_ONE/CNST(600.0) - pi2/CNST(4096.0)
      end select

      lapl%w(1,1) = fp(1)*sum(1/der%m%h(1:conf%dim)**2)

      k = 1
      do i = 1, conf%dim
        do j = -der%order, -1
          k = k + 1
          lapl%w(k,1) = fp(-j+1) / der%m%h(i)**2
        end do

        do j = 1, der%order
          k = k + 1
          lapl%w(k,1) = fp( j+1) / der%m%h(i)**2
        end do
      end do
    
      deallocate(fp)
    end subroutine variational
  end subroutine derivatives_build_lapl

  ! ---------------------------------------------------------
  subroutine derivatives_build_grad(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i

    ASSERT(associated(der%grad))

    do i = 1, conf%dim
      select case(der%stencil_type)
      case(DER_STAR, DER_VARIATIONAL)  
        call nl_operator_build(der%m, der%grad(i), der%m%np, .true.)
        call star(der%grad(i))
      end select
    end do

  contains
    subroutine star(grad)
      type(nl_operator_type), intent(inout) :: grad

      integer :: i, j, k, morder
      FLOAT, allocatable :: cc(:,:,:)

      ASSERT(grad%n >= 3)

      morder = grad%n - 1
      allocate(cc(0:morder, 0:morder, 0:1))
      call weights(1, morder, cc)
      
      k = 1
      do j = -der%order, -1
        grad%w(k,1) = cc(-2*j-1, morder, 1) / der%m%h(i)
        k = k + 1
      end do

      grad%w(k,1) = cc(0, morder, 1) / der%m%h(i)

      do j = 1, der%order
        k = k + 1
        grad%w(k,1) = cc(2*j, morder, 1) / der%m%h(i)
      end do

      deallocate(cc)
    end subroutine star

  end subroutine derivatives_build_grad

#include "undef.F90"
#include "real.F90"
#include "derivatives_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "derivatives_inc.F90"

end module derivatives
