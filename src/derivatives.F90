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
  use lib_adv_alg
  use stencil_star
  use stencil_variational
  use stencil_cube

  implicit none

  private
  public :: der_discr_type, &
            derivatives_init, &
            derivatives_end, &
            derivatives_build, &
            dderivatives_lapl, &
            zderivatives_lapl, &
            dderivatives_laplt, &
            zderivatives_laplt, &
            dderivatives_grad, &
            zderivatives_grad, &
            dderivatives_div, &
            zderivatives_div, &
            dderivatives_curl, &
            zderivatives_curl

  type der_discr_type
    type(mesh_type), pointer :: m             ! pointer to the underlying mesh
    integer                  :: order         ! order of the discretization (value depends on stencil)
    integer                  :: stencil_type  ! type of discretization

    integer                  :: boundaries(3) ! bounday conditions
    logical                  :: zero_bc

    FLOAT :: lapl_cutoff ! If the so-called variational discretization is used, this controls a
                         ! possible filter on the Laplacian.

    type(nl_operator_type), pointer :: op(:)  ! op(1:conf%dim) => gradient
                                              ! op(conf%dim+1) => laplacian
    type(nl_operator_type), pointer :: lapl   ! these are just shortcuts for op
    type(nl_operator_type), pointer :: grad(:)

    type(nl_operator_type) :: laplt ! The transponse of the Laplacian.
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

    call loct_parse_int(check_inp('DerivativesStencil'), DER_STAR, der%stencil_type)
    if(der%stencil_type<DER_STAR.or.der%stencil_type>DER_CUBE) then
      write(message(1), '(a,i2,a)') 'DerivativesStencil = "', der%stencil_type, '" is unknown to octopus'
      call write_fatal(1)
    end if
    if(der%stencil_type == DER_VARIATIONAL) then
       call loct_parse_float(check_inp('DerivativesLaplacianFilter'), M_ONE, der%lapl_cutoff)
    endif

    call loct_parse_int(check_inp('DerivativesOrder'), 4, der%order)

    ! construct lapl and grad structures
    allocate(der%op(conf%dim + 1))
    der%grad => der%op
    der%lapl => der%op(conf%dim + 1)

    call derivatives_get_stencil_lapl(der)
    call derivatives_get_stencil_grad(der)

    ! find out the bounday conditions
    call loct_parse_int(check_inp('DerivativesBoundaries'), DER_BC_ZERO_F, i)
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

    ASSERT(associated(der%op))

    do i = 1, conf%dim+1
      call nl_operator_end(der%op(i))
    end do

    deallocate(der%op)
    nullify   (der%op, der%lapl, der%grad, der%m)

    call pop_sub()
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der)
    type(der_discr_type), intent(inout) :: der

    integer :: n
    
    call push_sub('derivatives_get_stencil_lapl')
    
    ASSERT(associated(der%lapl))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR,DER_VARIATIONAL)
      n = stencil_star_size_lapl(der%order)
    case(DER_CUBE)
      n = stencil_cube_size_lapl(der%order)
    end select

    ! initialize nl operator
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
    
    ASSERT(associated(der%grad))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      n = stencil_star_size_grad(der%order)
    case(DER_CUBE)
      n = stencil_cube_size_grad(der%order)
    end select
    
    ! initialize nl operator
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

    integer, allocatable :: polynomials(:,:)
    FLOAT,   allocatable :: rhs(:,:)
    logical :: const_w    
    integer :: i

    call push_sub('derivatives_build')
    
    ASSERT(associated(der%op))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_CUBE)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL.and.der%m%use_curvlinear))

    ! build operators
    do i = 1, conf%dim+1
      call nl_operator_build(der%m, der%op(i), der%m%np, .not.der%m%use_curvlinear)
    end do

    if(der%m%use_curvlinear.or.der%stencil_type==DER_CUBE) then

      select case(der%stencil_type)
      case(DER_STAR) ! laplacian and gradient have different stencils
        do i = 1, conf%dim + 1
          allocate(polynomials(conf%dim, der%op(i)%n), rhs(der%op(i)%n,1))

          if(i <= conf%dim) then  ! gradient
            call stencil_star_polynomials_grad(i, polynomials, der%order)
            call get_rhs_grad(i, rhs(:,1))
          else                      ! laplacian
            call stencil_star_polynomials_lapl(polynomials, der%order)
            call get_rhs_lapl(rhs(:,1))
          end if

          call make_discretization(der%m, polynomials, rhs, 1, der%op(i:i))
          deallocate(polynomials, rhs)
        end do

      case(DER_CUBE) ! laplacian and gradient have similar stencils
        allocate(polynomials(conf%dim, der%op(1)%n), rhs(der%op(1)%n,conf%dim+1))
        call stencil_cube_polynomials_lapl(polynomials, der%order)

        do i = 1, conf%dim
          call get_rhs_grad(i, rhs(:,i))
        end do
        call get_rhs_lapl(rhs(:,conf%dim+1))

        call make_discretization(der%m, polynomials, rhs, conf%dim+1, der%op(:))

        deallocate(polynomials, rhs)
      end select


    else ! we have the explicit coefficients

      ! get laplacian
      select case(der%stencil_type)
      case(DER_STAR)  
         write(*,*) 'der%m%h',der%m%h
        call stencil_star_coeff_lapl(der%m%h, der%order, der%lapl)
      case(DER_VARIATIONAL)
        call stencil_variational_coeff_lapl(der%m%h, der%order, der%lapl, alpha = der%lapl_cutoff)
      end select

      ! get gradient (we use the same both for star and variational)
      do i = 1, conf%dim
        call stencil_star_coeff_grad(der%m%h(i), der%order, der%grad(i))
      end do
    end if

    if(der%m%use_curvlinear) then
       call nl_operator_transpose(der%lapl, der%laplt)
    else
       der%laplt = der%lapl
    endif

    call pop_sub()
  contains

    subroutine get_rhs_lapl(rhs)
      FLOAT, intent(out) :: rhs(:)

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

    end subroutine get_rhs_lapl

    subroutine get_rhs_grad(dir, rhs)
      integer, intent(in)  :: dir
      FLOAT,   intent(out) :: rhs(:)

      integer :: j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do j = 1, der%grad(dir)%n
        this_one = .true.
        do k = 1, conf%dim
          if(k == dir .and. polynomials(k, j).ne.1) this_one = .false.
          if(k.ne.dir .and. polynomials(k, j).ne.0) this_one = .false.
        end do
        if(this_one) rhs(j) = M_ONE
      end do

    end subroutine get_rhs_grad
  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine make_discretization(m, pol, rhs, n, op)
    type(mesh_type),        intent(in)    :: m
    integer,                intent(in)    :: pol(:,:)  ! pol(conf%dim, op%n)
    integer,                intent(in)    :: n
    FLOAT,                  intent(inout) :: rhs(:,:)   ! rhs_(op%n, n)
    type(nl_operator_type), intent(inout) :: op(:)
    
    integer :: p, p_max, i, j, k, pow_max
    FLOAT   :: x(conf%dim)
    FLOAT, allocatable :: mat(:,:), sol(:,:), powers(:,:)
    
    call push_sub('make_discretization')

    allocate(mat(op(1)%n, op(1)%n), sol(op(1)%n, n))
    
    ! use to generate power lookup table
    pow_max = maxval(pol)
    allocate(powers(conf%dim, 0:pow_max))
    powers(:,:) = M_ZERO
    powers(:,0) = M_ONE

    p_max = op(1)%np
    if(op(1)%const_w) p_max = 1

    do p = 1, p_max
      mat(1,:) = M_ONE
      do i = 1, op(1)%n
        x(1:conf%dim) = m%x(p, 1:conf%dim) - m%x(op(1)%i(i, p), 1:conf%dim)

        ! calculate powers
        do j = 1, conf%dim
          powers(j,1) = x(j)
          do k = 2, pow_max
            powers(j,k) = x(j)*powers(j,k-1)
          end do
        end do

        ! generate the matrix
        do j = 2, op(1)%n
          mat(j, i) = powers(1, pol(1, j))
          do k = 2, conf%dim
            mat(j, i) = mat(j, i)*powers(k, pol(k, j))
          end do
        end do
      end do
      
      call lalg_linsyssolve(op(1)%n, n, mat, rhs, sol)
      do i = 1, n
        op(i)%w(:, p) = sol(:, n)
      end do

    end do

    deallocate(mat, sol, powers)

    call pop_sub()
  end subroutine make_discretization

#include "undef.F90"
#include "real.F90"
#include "derivatives_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "derivatives_inc.F90"

end module derivatives
