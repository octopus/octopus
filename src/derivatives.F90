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
!!
!! $Id$

#include "global.h"

module derivatives
  use global
  use messages
  use syslabels
  use varinfo
  use lib_oct_parser
  use mesh
  use nl_operator
  use lib_adv_alg
  use stencil_star
  use stencil_variational
  use stencil_cube
  use stencil_starplus
  use simul_box

  implicit none

  private
  public ::              &
    der_discr_type,      &
    derivatives_init,    &
    derivatives_end,     &
    derivatives_build,   &
    dderivatives_lapl,   &
    dderivatives_laplt,  &
    dderivatives_grad,   &
    dderivatives_div,    &
    dderivatives_curl,   &
    zderivatives_lapl,   &
    zderivatives_laplt,  &
    zderivatives_grad,   &
    zderivatives_div,    &
    zderivatives_curl

  integer, parameter ::  &
    DER_BC_ZERO_F  = 0,  &  ! function is zero at the boundaries
    DER_BC_ZERO_DF = 1,  &  ! first derivative of the function is zero
    DER_BC_PERIOD  = 2      ! boundary is periodic

  integer, parameter ::  &
    DER_STAR        = 1, &
    DER_VARIATIONAL = 2, &
    DER_CUBE        = 3, &
    DER_STARPLUS    = 4

  type der_discr_type
    type(mesh_type), pointer :: m             ! pointer to the underlying mesh
    integer                  :: dim           ! dimensionality of the space (sb%dim)
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


contains

  ! ---------------------------------------------------------
  subroutine derivatives_init(sb, der, n_ghost, use_curvilinear)
    type(simul_box_type), intent(in)  :: sb
    type(der_discr_type), intent(out) :: der
    integer,              intent(out) :: n_ghost(:)
    logical,              intent(in)  :: use_curvilinear

    integer :: i

    call push_sub('derivatives.derivatives_init')

    ! copy this value to my structure
    der%dim = sb%dim

    !%Variable DerivativesStencil
    !%Type integer
    !%Default stencil_star
    !%Section Mesh::Derivatives
    !%Description
    !% Decides what kind of stencil is used, i.e. what points, around
    !% each point in the mesh, are the neighboring points used in the
    !% expression of the differential operator.
    !%
    !% If curvilinear coordinates are to be used, then only the "stencil_starplus"
    !% or the "stencil_cube" may be used. We only recommend the "stencil_starplus",
    !% the cube typically needing way too much memory resources.
    !%Option stencil_star 1
    !% A star around each point (i.e., only points in the axis).
    !%Option stencil_variational 2
    !% Same as the star, but with coefficients built in a different way.
    !%Option stencil_cube 3
    !% A cube of points around each point.
    !%Option stencil_starplus 4
    !% The star, plus a number of off-axis points.
    !%End
    if(use_curvilinear) then
      call loct_parse_int(check_inp('DerivativesStencil'), DER_STARPLUS, der%stencil_type)
    else
      call loct_parse_int(check_inp('DerivativesStencil'), DER_STAR, der%stencil_type)
    endif
    if(.not.varinfo_valid_option('DerivativesStencil', der%stencil_type)) call input_error('DerivativesStencil')
    if(use_curvilinear  .and.  der%stencil_type < DER_CUBE) call input_error('DerivativesStencil')
    if(der%stencil_type == DER_VARIATIONAL) then
      call loct_parse_float(check_inp('DerivativesLaplacianFilter'), M_ONE, der%lapl_cutoff)
    end if

    call loct_parse_int(check_inp('DerivativesOrder'), 4, der%order)

    ! construct lapl and grad structures
    allocate(der%op(der%dim + 1))
    der%grad => der%op
    der%lapl => der%op(der%dim + 1)

    call derivatives_get_stencil_lapl(der)
    call derivatives_get_stencil_grad(der)

    ! find out the bounday conditions
    call loct_parse_int(check_inp('DerivativesBoundaries'), DER_BC_ZERO_F, i)
    if((i < DER_BC_ZERO_F).or.(i > DER_BC_PERIOD)) then
      write(message(1), '(a,i2,a)') 'DerivativesBoundaries = "', i, '" is unknown to octopus'
      call write_fatal(1)
    end if

    der%zero_bc = (i == DER_BC_ZERO_F)
    der%boundaries(:) = i
    der%boundaries(1:sb%periodic_dim) = DER_BC_PERIOD

    ! find out how many ghost points we need in which dimension
    n_ghost(:) = 0
    do i = 1, der%dim
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

    call push_sub('derivatives.derivatives_end')

    ASSERT(associated(der%op))

    do i = 1, der%dim+1
      call nl_operator_end(der%op(i))
    end do

    deallocate(der%op)
    nullify   (der%op, der%lapl, der%grad)

    call pop_sub()
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der)
    type(der_discr_type), intent(inout) :: der

    integer :: n

    call push_sub('derivatives.derivatives_get_stencil_lapl')

    ASSERT(associated(der%lapl))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR,DER_VARIATIONAL)
      n = stencil_star_size_lapl(der%dim, der%order)
    case(DER_CUBE)
      n = stencil_cube_size_lapl(der%dim, der%order)
    case(DER_STARPLUS)
      n = stencil_starplus_size_lapl(der%dim, der%order)
    end select

    ! initialize nl operator
    call nl_operator_init(der%lapl, n)

    ! create stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      call stencil_star_get_lapl(der%dim, der%order, der%lapl%stencil)
    case(DER_CUBE)
      call stencil_cube_get_lapl(der%dim, der%order, der%lapl%stencil)
    case(DER_STARPLUS)
      call stencil_starplus_get_lapl(der%dim, der%order, der%lapl%stencil)
    end select

    call pop_sub()

  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_grad(der)
    type(der_discr_type), intent(inout) :: der

    integer :: i, n

    call push_sub('derivatives.derivatives_get_stencil_grad')

    ASSERT(associated(der%grad))

    ! get size of stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      n = stencil_star_size_grad(der%order)
    case(DER_CUBE)
      n = stencil_cube_size_grad(der%dim, der%order)
    case(DER_STARPLUS)
      n = stencil_starplus_size_grad(der%dim, der%order)
    end select

    ! initialize nl operator
    do i = 1, der%dim
      call nl_operator_init(der%grad(i), n)

      ! create stencil
      select case(der%stencil_type)
      case(DER_STAR, DER_VARIATIONAL)
        call stencil_star_get_grad(i, der%order, der%grad(i)%stencil)
      case(DER_CUBE)
        call stencil_cube_get_grad(der%dim, der%order, der%grad(i)%stencil)
      case(DER_STARPLUS)
        call stencil_starplus_get_grad(der%dim, i, der%order, der%grad(i)%stencil)
      end select
    end do

    call pop_sub()

  end subroutine derivatives_get_stencil_grad


  ! ---------------------------------------------------------
  subroutine derivatives_build(m, der)
    type(mesh_type), target, intent(in)    :: m
    type(der_discr_type),    intent(inout) :: der

    integer, allocatable :: polynomials(:,:)
    FLOAT,   allocatable :: rhs(:,:)
    integer :: i, j, k

    type(nl_operator_type) :: auxop

    call push_sub('derivatives.derivatives_build')

    ASSERT(associated(der%op))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_STARPLUS)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL.and.m%use_curvlinear))

    der%m => m    ! make a pointer to the underlying mesh

    ! build operators
    do i = 1, der%dim+1
      call nl_operator_build(m, der%op(i), der%m%np, .not.m%use_curvlinear)
    end do

    if(m%use_curvlinear.or.der%stencil_type==DER_CUBE) then

      select case(der%stencil_type)
      case(DER_STAR) ! laplacian and gradient have different stencils
        do i = 1, der%dim + 1
          allocate(polynomials(der%dim, der%op(i)%n), rhs(der%op(i)%n,1))

          if(i <= der%dim) then  ! gradient
            call stencil_star_polynomials_grad(i, der%order, polynomials)
            call get_rhs_grad(i, rhs(:,1))
          else                      ! laplacian
            call stencil_star_polynomials_lapl(der%dim, der%order, polynomials)
            call get_rhs_lapl(rhs(:,1))
          end if

          call make_discretization(der%dim, der%m, polynomials, rhs, 1, der%op(i:i))
          deallocate(polynomials, rhs)
        end do

      case(DER_CUBE) ! laplacian and gradient have similar stencils
        allocate(polynomials(der%dim, der%op(1)%n), rhs(der%op(1)%n, der%dim+1))
        call stencil_cube_polynomials_lapl(der%dim, der%order, polynomials)

        do i = 1, der%dim
          call get_rhs_grad(i, rhs(:,i))
        end do
        call get_rhs_lapl(rhs(:, der%dim+1))

        call make_discretization(der%dim, der%m, polynomials, rhs, der%dim+1, der%op(:))

        deallocate(polynomials, rhs)
      case(DER_STARPLUS)
        do i = 1, der%dim
          allocate(polynomials(der%dim, der%op(i)%n), rhs(der%op(i)%n, 1))
          call stencil_starplus_pol_grad(der%dim, i, der%order, polynomials)
          call get_rhs_grad(i, rhs(:, 1))
          call make_discretization(der%dim, der%m, polynomials, rhs, 1, der%op(i:i))
          deallocate(polynomials, rhs)
        end do
        allocate(polynomials(der%dim, der%op(der%dim+1)%n), rhs(der%op(i)%n, 1))
        call stencil_starplus_pol_lapl(der%dim, der%order, polynomials)
        call get_rhs_lapl(rhs(:, 1))
        call make_discretization(der%dim, der%m, polynomials, rhs, 1, der%op(der%dim+1:der%dim+1))
        deallocate(polynomials, rhs)
      end select

      ! Here the Laplacian is forced to be self-adjoint, and the gradient to be skew-selfadjoint
      if(m%use_curvlinear) then
        do i = 1, der%dim
          call nl_operator_init(auxop, der%grad(i)%n)
          auxop%stencil = der%grad(i)%stencil
          call nl_operator_build(m, auxop, der%m%np, const_w = .false.)
          call nl_operator_skewadjoint(der%grad(i), auxop, der%m)
          call nl_operator_equal(der%grad(i), auxop)
          call nl_operator_end(auxop)
        end do
        call nl_operator_init(auxop, der%lapl%n)
        auxop%stencil = der%lapl%stencil
        call nl_operator_build(m, auxop, der%m%np, const_w = .false.)
        call nl_operator_selfadjoint(der%lapl, auxop, der%m)
        call nl_operator_equal(der%lapl, auxop)
        call nl_operator_end(auxop)
      end if

      ! Here I will nullify all the coefficients that are outside the box (only for
      ! the case of non-constant weights == curvilinear coordinates).
      ! WARNING: Same thing should be done for the gradients. The subroutines in
      ! derivatives_inc.F90 should then be changed accordingly.
      if(m%use_curvlinear) then
        do i = 1, m%np
          do j = 1, der%lapl%n
            k = der%lapl%i(j, i)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
            if(k>m%vp%np_local(m%vp%partno) + m%vp%np_ghost(m%vp%partno)) then
#else
            if(k>m%np_global) then
#endif
              der%lapl%w_re(j, i) = M_ZERO
              der%lapl%i(j, i) = i
            end if
          end do
        end do
      end if

    else ! we have the explicit coefficients

      ! get laplacian
      select case(der%stencil_type)
      case(DER_STAR)
        call stencil_star_coeff_lapl(der%dim, der%order, m%h(1:der%dim), der%lapl)
      case(DER_VARIATIONAL)
        call stencil_variational_coeff_lapl(der%dim, der%order, m%h, der%lapl, alpha = der%lapl_cutoff)
      case(DER_STARPLUS)
        ! equivalent to normal stencil.
        der%stencil_type = DER_STAR
        call stencil_star_coeff_lapl(der%dim, der%order, m%h(1:der%dim), der%lapl)
      end select

      ! get gradient (we use the same both for star and variational)
      do i = 1, der%dim
        call stencil_star_coeff_grad(der%order, m%h(i), der%grad(i))
      end do
    end if

    call pop_sub()

  contains

    subroutine get_rhs_lapl(rhs)
      FLOAT, intent(out) :: rhs(:)

      integer :: i, j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do i = 1, der%dim
        do j = 1, der%lapl%n
          this_one = .true.
          do k = 1, der%dim
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
        do k = 1, der%dim
          if(k == dir .and. polynomials(k, j).ne.1) this_one = .false.
          if(k.ne.dir .and. polynomials(k, j).ne.0) this_one = .false.
        end do
        if(this_one) rhs(j) = M_ONE
      end do

    end subroutine get_rhs_grad

  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine make_discretization(dim, m, pol, rhs, n, op)
    integer,                intent(in)    :: dim
    type(mesh_type),        intent(in)    :: m
    integer,                intent(in)    :: pol(:,:)  ! pol(dim, op%n)
    integer,                intent(in)    :: n
    FLOAT,                  intent(inout) :: rhs(:,:)   ! rhs_(op%n, n)
    type(nl_operator_type), intent(inout) :: op(:)

    integer :: p, p_max, i, j, k, pow_max
    FLOAT   :: x(dim)
    FLOAT, allocatable :: mat(:,:), sol(:,:), powers(:,:)

    call push_sub('derivatives.make_discretization')

    allocate(mat(op(1)%n, op(1)%n), sol(op(1)%n, n))

    message(1) = 'Info: Generating weights for finite-difference discretization.'
    call write_info(1)

    ! use to generate power lookup table
    pow_max = maxval(pol)
    allocate(powers(dim, 0:pow_max))
    powers(:,:) = M_ZERO
    powers(:,0) = M_ONE

    p_max = op(1)%np
    if(op(1)%const_w) p_max = 1

    do p = 1, p_max
      mat(1,:) = M_ONE
      do i = 1, op(1)%n
        x(1:dim) = m%x(op(1)%i(i, p), 1:dim) - m%x(p, 1:dim)

        ! calculate powers
        do j = 1, dim
          powers(j,1) = x(j)
          do k = 2, pow_max
            powers(j,k) = x(j)*powers(j,k-1)
          end do
        end do

        ! generate the matrix
        do j = 2, op(1)%n
          mat(j, i) = powers(1, pol(1, j))
          do k = 2, dim
            mat(j, i) = mat(j, i)*powers(k, pol(k, j))
          end do
        end do
      end do

      call lalg_linsyssolve(op(1)%n, n, mat, rhs, sol)
      do i = 1, n
        op(i)%w_re(:, p) = sol(:, n)
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
