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
  use global
  use mesh
  use nl_operator
  use stencil_star
  use stencil_variational
  use stencil_cube

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
     DER_BC_ZERO_F  = 0,   &  ! function is zero at the boundaries
     DER_BC_ZERO_DF = 1,   &  ! first derivative of the function is zero
     DER_BC_PERIOD  = 2       ! boundary is periodic

  integer, parameter ::    &
     DER_STAR        = 1,  &
     DER_VARIATIONAL = 2,  &
     DER_CUBE        = 3

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
    if(der%stencil_type<DER_STAR.or.der%stencil_type>DER_CUBE) then
      write(message(1), '(a,i2,a)') 'DerivativesStencil = "', der%stencil_type, '" is unknown to octopus'
      call write_fatal(1)
    end if

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

    call push_sub('derivatives_end')

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
    
    call push_sub('derivatives_get_stencil_lapl')
    
    ASSERT(.not.associated(der%lapl))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR,DER_VARIATIONAL)
      n = stencil_star_size_lapl(der%order)
    case(DER_CUBE)
      n = stencil_cube_size_lapl(der%order)
    end select

    ! initialize nl operator
    allocate(der%lapl)
    call nl_operator_init(der%lapl, n)
    
    ! create stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      call stencil_star_get_lapl(der%order, der%lapl%stencil)
    case(DER_CUBE)
      call stencil_cube_get_lapl(der%order, der%lapl%stencil)
    end select

    call pop_sub()

  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_grad(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i, n
    
    call push_sub('derivatives_get_stencil_grad')
    
    ASSERT(.not.associated(der%grad))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      n = stencil_star_size_grad(der%order)
    case(DER_CUBE)
      n = stencil_cube_size_grad(der%order)
    end select
    
    ! initialize nl operator
    allocate(der%grad(conf%dim))
    do i = 1, conf%dim
      call nl_operator_init(der%grad(i), n)

      ! create stencil
      select case(der%stencil_type)
      case(DER_STAR, DER_VARIATIONAL)
        call stencil_star_get_grad(i, der%order, der%grad(i)%stencil)
      case(DER_CUBE)
        call stencil_cube_get_grad(i, der%order, der%grad(i)%stencil)
      end select
    end do

    call pop_sub()

  end subroutine derivatives_get_stencil_grad


  ! ---------------------------------------------------------
  subroutine derivatives_build(der)
    type(der_discr_type), intent(inout) :: der

    call push_sub('derivatives_build')
    
    call derivatives_build_lapl(der)
    call derivatives_build_grad(der)
    
    call pop_sub()
  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine derivatives_build_lapl(der)
    type(der_discr_type), intent(inout) :: der

    integer, allocatable :: polynomials(:,:), rhs(:)
    logical :: const_w

    call push_sub('derivatives_build_lapl')
    
    ASSERT(associated(der%lapl))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_CUBE)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL.and.der%m%use_curvlinear))

    if(der%m%use_curvlinear.or.der%stencil_type==DER_CUBE) then
      call nl_operator_build(der%m, der%lapl, der%m%np, .not.der%m%use_curvlinear)
      allocate(polynomials(conf%dim, der%lapl%n), rhs(der%lapl%n))

      select case(der%stencil_type)
      case(DER_STAR)
        call stencil_star_polynomials_lapl(polynomials, der%order)
      case(DER_CUBE)
        call stencil_cube_polynomials_lapl(polynomials, der%order)
      end select

      call get_rhs()
      call make_discretization(der%m, polynomials, rhs, der%lapl)
      deallocate(polynomials, rhs)
    else
      call nl_operator_build(der%m, der%lapl, der%m%np, .true.)
      select case(der%stencil_type)
      case(DER_STAR)  
        call stencil_star_coeff_lapl(der%m%h, der%order, der%lapl)
      case(DER_VARIATIONAL)
        call stencil_variational_coeff_lapl(der%m%h, der%order, der%lapl)
      end select
    end if

    call pop_sub()

  contains

    subroutine get_rhs()
      integer :: i, j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do i = 1, conf%dim
        do j = 1, der%lapl%n
          this_one = .true.
          do k = 1, conf%dim
            if(k == i .and. polynomials(k, j).ne.2) this_one = .false.
            if(k.ne.i .and. polynomials(k, j).ne.0) this_one = .false.
          end do
          if(this_one) rhs(j) = M_TWO
        end do
      end do

    end subroutine get_rhs

  end subroutine derivatives_build_lapl


  ! ---------------------------------------------------------
  subroutine derivatives_build_grad(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i
    integer, allocatable :: polynomials(:,:), rhs(:)

    call push_sub('derivatives_build_grad')

    ASSERT(associated(der%grad))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_CUBE)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL.and.der%m%use_curvlinear))

    do i = 1, conf%dim
      if(der%m%use_curvlinear.or.der%stencil_type==DER_CUBE) then
        call nl_operator_build(der%m, der%grad(i), der%m%np, .not.der%m%use_curvlinear)
        allocate(polynomials(conf%dim, der%grad(i)%n), rhs(der%grad(i)%n))
      
        select case(der%stencil_type)
        case(DER_STAR)
          call stencil_star_polynomials_grad(i, polynomials, der%order)
        case(DER_CUBE)
          call stencil_cube_polynomials_grad(i, polynomials, der%order)
        end select

        call get_rhs()
        call make_discretization(der%m, polynomials, rhs, der%grad(i))
        deallocate(polynomials, rhs)
      else
        call nl_operator_build(der%m, der%grad(i), der%m%np, .true.)
        select case(der%stencil_type)
        case(DER_STAR, DER_VARIATIONAL)
          call stencil_star_coeff_grad(der%m%h(i), der%order, der%grad(i))
        end select
      end if
    end do

    call pop_sub()

  contains
    subroutine get_rhs()
      integer :: j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do j = 1, der%grad(i)%n
        this_one = .true.
        do k = 1, conf%dim
          if(k == i .and. polynomials(k, j).ne.1) this_one = .false.
          if(k.ne.i .and. polynomials(k, j).ne.0) this_one = .false.
        end do
        if(this_one) rhs(j) = M_ONE
      end do

    end subroutine get_rhs

  end subroutine derivatives_build_grad


  ! ---------------------------------------------------------
  subroutine make_discretization(m, pol, rhs_, op)
    type(mesh_type),        intent(in)    :: m
    integer,                intent(in)    :: pol(:,:)  ! pol(conf%dim, op%n)
    integer,                intent(in)    :: rhs_(:)   ! rhs_(op%n, 1)
    type(nl_operator_type), intent(inout) :: op
    
    integer :: p, p_max, i, j, k
    FLOAT   :: x(conf%dim)
    FLOAT, allocatable :: mat(:,:), rhs(:,:), sol(:,:)
    
    call push_sub('make_discretization')

    allocate(mat(op%n, op%n), rhs(op%n, 1), sol(op%n, 1))
    rhs(:, 1) = rhs_(:)
    
    p_max = op%np
    if(op%const_w) p_max = 1

    do p = 1, p_max
      mat(1,:) = M_ONE
      do i = 1, op%n
        x(1:conf%dim) = m%x(p, 1:conf%dim) - m%x(op%i(i, p), 1:conf%dim)
        do j = 2, op%n
          mat(j, i) = product(x(1:conf%dim)**pol(1:conf%dim, j))
        end do
      end do
      
      call lalg_linsyssolve(op%n, 1, mat, rhs, sol)
      op%w(:, p) = sol(:, 1)
      
    end do

    deallocate(mat, rhs, sol)

    call pop_sub()
  end subroutine make_discretization

#include "undef.F90"
#include "real.F90"
#include "derivatives_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "derivatives_inc.F90"

end module derivatives
