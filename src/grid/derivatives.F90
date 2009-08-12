!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module derivatives_m
  use batch_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use math_m
  use mesh_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use stencil_cube_m
  use stencil_star_m
  use stencil_starplus_m
  use stencil_variational_m
  use varinfo_m

  implicit none

  private
  public ::                             &
    derivatives_t,                      &
    derivatives_init,                   &
    derivatives_end,                    &
    derivatives_build,                  &
    dderivatives_lapl,                  &
    zderivatives_lapl,                  &
    derivatives_lapl_diag,              &
    dderivatives_laplt,                 &
    zderivatives_laplt,                 &
    dderivatives_grad,                  &
    zderivatives_grad,                  &
    dderivatives_oper,                  &
    zderivatives_oper,                  &
    dderivatives_oper_batch,            &
    zderivatives_oper_batch,            &
    dderivatives_div,                   &
    zderivatives_div,                   &
    dderivatives_curl,                  &
    zderivatives_curl,                  &
    dset_bc,                            &
    zset_bc,                            &
    dset_bc_batch,                      &
    zset_bc_batch,                      &
    stencil_extent,                     &
    derivatives_handle_batch_t,         &
    dderivatives_batch_start,           &
    dderivatives_batch_finish,          &
    zderivatives_batch_start,           &
    zderivatives_batch_finish,          &
    df_angular_momentum,                &
    zf_angular_momentum,                &
    df_l2, zf_l2

  integer, parameter ::     &
    DER_BC_ZERO_F    = 0,   &  ! function is zero at the boundaries
    DER_BC_ZERO_DF   = 1,   &  ! first derivative of the function is zero
    DER_BC_PERIOD    = 2       ! boundary is periodic

  integer, parameter ::     &
    DER_STAR         = 1,   &
    DER_VARIATIONAL  = 2,   &
    DER_CUBE         = 3,   &
    DER_STARPLUS     = 4


  type derivatives_t
    type(mesh_t), pointer :: mesh          ! pointer to the underlying mesh
    integer               :: dim           ! dimensionality of the space (sb%dim)
    integer               :: order         ! order of the discretization (value depends on stencil)
    integer               :: stencil_type  ! type of discretization

    FLOAT                 :: masses(MAX_DIM)     ! we can have different weights (masses) per space direction

    integer               :: boundaries(MAX_DIM) ! boundary conditions
    logical               :: zero_bc
    logical               :: periodic_bc

    ! If the so-called variational discretization is used, this controls a
    ! possible filter on the Laplacian.
    FLOAT :: lapl_cutoff   

    type(nl_operator_t), pointer :: op(:)  ! op(1:conf%dim) => gradient
    ! op(conf%dim+1) => laplacian
    type(nl_operator_t), pointer :: lapl   ! these are just shortcuts for op
    type(nl_operator_t), pointer :: grad(:)

    type(nl_operator_t) :: laplt ! The transpose of the Laplacian.
    integer             :: n_ghost(MAX_DIM)   ! ghost points to add in each dimension
#if defined(HAVE_MPI)
    integer             :: comm_method 
#endif
  end type derivatives_t


  type derivatives_handle_t
    private
#ifdef HAVE_MPI
    type(pv_handle_t) :: pv_h
    logical           :: parallel_in_domains
#endif
    FLOAT, pointer    :: df(:)
    CMPLX, pointer    :: zf(:)
    FLOAT, pointer    :: dlapl(:)
    CMPLX, pointer    :: zlapl(:)
    logical           :: ghost_update
  end type derivatives_handle_t


  type derivatives_handle_batch_t
    private
#ifdef HAVE_MPI
    type(pv_handle_batch_t)      :: pv_h
#endif
    type(derivatives_t), pointer :: der
    type(nl_operator_t), pointer :: op
    type(batch_t),       pointer :: ff
    type(batch_t),       pointer :: opff
    logical                      :: ghost_update
  end type derivatives_handle_batch_t

  type(profile_t), save :: set_bc_prof
#ifdef HAVE_MPI
  type(profile_t), save :: set_bc_comm_prof
#endif

contains

  ! ---------------------------------------------------------
  subroutine derivatives_init(der, sb, use_curvilinear)
    type(derivatives_t), intent(out) :: der
    type(simul_box_t),   intent(in)  :: sb
    logical,             intent(in)  :: use_curvilinear

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
      call loct_parse_int(datasets_check('DerivativesStencil'), DER_STARPLUS, der%stencil_type)
    else
      call loct_parse_int(datasets_check('DerivativesStencil'), DER_STAR, der%stencil_type)
    endif
    if(.not.varinfo_valid_option('DerivativesStencil', der%stencil_type)) call input_error('DerivativesStencil')
    call messages_print_var_option(stdout, "DerivativesStencil", der%stencil_type)

    if(use_curvilinear  .and.  der%stencil_type < DER_CUBE) call input_error('DerivativesStencil')
    if(der%stencil_type == DER_VARIATIONAL) then
      call loct_parse_float(datasets_check('DerivativesLaplacianFilter'), M_ONE, der%lapl_cutoff)
    end if

    !%Variable DerivativesOrder
    !%Type integer
    !%Default 4
    !%Section Mesh::Derivatives
    !%Description
    !% This variable gives the discretization order for the approximation of
    !% the differential operators. This means, basically, that
    !% <tt>DerivativesOrder</tt> points are used in each positive/negative
    !% spatial direction, e. g. <tt>DerivativesOrder = 1</tt> would give
    !% the well-known three-point formula in 1D.
    !% The number of points actually used for the Laplacian
    !% depends on the stencil used:
    !%
    !% <tt>stencil_star</tt>: 2*<tt>DerivativesOrder</tt>*dim+1
    !%
    !% <tt>stencil_cube</tt>: (2*<tt>DerivativesOrder</tt>+1)^dim
    !%
    !% <tt>stencil_starplus</tt>: 2*<tt>DerivativesOrder</tt>+1+n with n being 12
    !% in 2D and 44 in 3D.
    !%End
    call loct_parse_int(datasets_check('DerivativesOrder'), 4, der%order)

#ifdef HAVE_MPI
    !%Variable ParallelizationOfDerivatives
    !%Type integer
    !%Default non_blocking
    !%Section Execution::Parallelization
    !%Description
    !% This option selects how the communication required for the
    !% synchronization is performed. The default is non_blocking.
    !%Option blocking 1
    !% Blocking communication.
    !%Option non_blocking 2
    !% Communication is based on non blocking point to point communication.
    !%Option non_blocking_collective 3
    !% Non-blocking collective communication (requires libnbc).
    !%End
    
    call loct_parse_int(datasets_check('ParallelizationOfDerivatives'), NON_BLOCKING, der%comm_method)
    
    if(.not. varinfo_valid_option('ParallelizationOfDerivatives', der%comm_method)) then
      call input_error('ParallelizationOfDerivatives')
    end if

#ifndef HAVE_LIBNBC
    if(der%comm_method == NON_BLOCKING_COLLECTIVE) then
      message(1) = "Error: libnbc is not available. Check the ParallelizationOfDerivatives variable."
      call write_fatal(1)
    end if
#endif

    call obsolete_variable('OverlapDerivatives', 'ParallelizationOfDerivatives')
#endif

    ! if needed, der%masses should be initialized in modelmb_particles_init
    der%masses = M_ONE

    ! construct lapl and grad structures
    SAFE_ALLOCATE(der%op(1:der%dim + 1))
    der%grad => der%op
    der%lapl => der%op(der%dim + 1)

    call derivatives_get_stencil_lapl(der)
    call derivatives_get_stencil_grad(der)

    ! find out the bounday conditions
    call loct_parse_int(datasets_check('DerivativesBoundaries'), DER_BC_ZERO_F, i)
    if((i < DER_BC_ZERO_F).or.(i > DER_BC_PERIOD)) then
      write(message(1), '(a,i2,a)') 'DerivativesBoundaries = "', i, '" is unknown to octopus'
      call write_fatal(1)
    end if

    der%zero_bc = (i == DER_BC_ZERO_F .and. sb%periodic_dim < 3)
    der%periodic_bc = (sb%periodic_dim > 0)
    der%boundaries(:) = i
    der%boundaries(1:sb%periodic_dim) = DER_BC_PERIOD

    ! find out how many ghost points we need in which dimension
    der%n_ghost(:) = 0
    do i = 1, der%dim
      if(der%boundaries(i) == DER_BC_ZERO_F .or. der%boundaries(i) == DER_BC_PERIOD) then
        der%n_ghost(i) = maxval(abs(der%lapl%stencil%points(i,:)))
      end if
    end do

    call pop_sub()
  end subroutine derivatives_init


  ! ---------------------------------------------------------
  subroutine derivatives_end(der)
    type(derivatives_t), intent(inout) :: der

    integer :: i

    call push_sub('derivatives.derivatives_end')

    ASSERT(associated(der%op))

    do i = 1, der%dim+1
      call nl_operator_end(der%op(i))
    end do

    SAFE_DEALLOCATE_P(der%op)
    nullify   (der%op, der%lapl, der%grad)

    call pop_sub()
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  ! Returns maximum extension of the stencil in spatial direction
  ! dir = 1, 2, 3 for a given derivative der.
  integer function stencil_extent(der, dir)
    type(derivatives_t), intent(in) :: der
    integer,           intent(in) :: dir

    call push_sub('derivatives.stencil_extent')

    select case(der%stencil_type)
      case(DER_STAR)
        stencil_extent = stencil_star_extent(dir, der%order)
      case(DER_VARIATIONAL)
        stencil_extent = stencil_variational_extent(dir, der%order)
      case(DER_CUBE)
        stencil_extent = stencil_cube_extent(dir, der%order)
      case(DER_STARPLUS)
        stencil_extent = stencil_cube_extent(dir, der%order)
      end select
      
    call pop_sub()
  end function stencil_extent


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der)
    type(derivatives_t), intent(inout) :: der

    call push_sub('derivatives.derivatives_get_stencil_lapl')

    ASSERT(associated(der%lapl))

    ! initialize nl operator
    call nl_operator_init(der%lapl, "Laplacian")

    ! create stencil
    select case(der%stencil_type)
    case(DER_STAR, DER_VARIATIONAL)
      call stencil_star_get_lapl(der%lapl%stencil, der%dim, der%order)
    case(DER_CUBE)
      call stencil_cube_get_lapl(der%lapl%stencil, der%dim, der%order)
    case(DER_STARPLUS)
      call stencil_starplus_get_lapl(der%lapl%stencil, der%dim, der%order)
    end select

    call pop_sub()

  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  ! Returns the diagonal elements of the laplacian needed for preconditioning
  subroutine derivatives_lapl_diag(der, lapl)
    type(derivatives_t), intent(in)  :: der
    FLOAT,               intent(out) :: lapl(:)  ! lapl(mesh%np)

    call push_sub('derivatives.derivatives_lapl_diag')

    ASSERT(ubound(lapl, DIM=1) >= der%mesh%np)

    ! the laplacian is a real operator
    call dnl_operator_operate_diag(der%lapl, lapl)

    call pop_sub()

  end subroutine derivatives_lapl_diag


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_grad(der)
    type(derivatives_t), intent(inout) :: der

    integer  :: ii
    character :: dir_label, dir_labels(4) = (/'X', 'Y', 'Z', 'W' /)

    call push_sub('derivatives.derivatives_get_stencil_grad')

    ASSERT(associated(der%grad))

    ! initialize nl operator
    do ii = 1, der%dim
      dir_label = ' '
      if(ii < 5) dir_label = dir_labels(ii)

      call nl_operator_init(der%grad(ii), "Gradient "//dir_label)

      ! create stencil
      select case(der%stencil_type)
      case(DER_STAR, DER_VARIATIONAL)
        call stencil_star_get_grad(der%grad(ii)%stencil, ii, der%order)
      case(DER_CUBE)
        call stencil_cube_get_grad(der%grad(ii)%stencil, der%dim, der%order)
      case(DER_STARPLUS)
        call stencil_starplus_get_grad(der%grad(ii)%stencil, der%dim, ii, der%order)
      end select
    end do

    call pop_sub()

  end subroutine derivatives_get_stencil_grad


  ! ---------------------------------------------------------
  subroutine derivatives_build(der, mesh)
    type(derivatives_t),    intent(inout) :: der
    type(mesh_t),   target, intent(in)    :: mesh

    integer, allocatable :: polynomials(:,:)
    FLOAT,   allocatable :: rhs(:,:)
    integer :: i
    logical :: const_w_, cmplx_op_

    type(nl_operator_t) :: auxop

    call push_sub('derivatives.derivatives_build')

    ASSERT(associated(der%op))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_STARPLUS)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL .and. mesh%use_curvilinear))

    der%mesh => mesh    ! make a pointer to the underlying mesh

    const_w_  = .true.
    cmplx_op_ = .false.

    ! need non-constant weights for curvilinear and scattering meshes
    if(mesh%use_curvilinear) const_w_ = .false.

    ! build operators
    do i = 1, der%dim+1
      call nl_operator_build(mesh, der%op(i), der%mesh%np, const_w = const_w_, cmplx_op = cmplx_op_)
    end do

    select case(der%stencil_type)

    case(DER_STAR) ! laplacian and gradient have different stencils
      do i = 1, der%dim + 1
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%size))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))

        if(i <= der%dim) then  ! gradient
          call stencil_star_polynomials_grad(i, der%order, polynomials)
          call get_rhs_grad(i, rhs(:,1))
        else                      ! laplacian
          call stencil_star_polynomials_lapl(der%dim, der%order, polynomials)
          call get_rhs_lapl(rhs(:,1))
        end if

        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(i:i))
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do

    case(DER_CUBE) ! laplacian and gradient have similar stencils
      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(1)%stencil%size))
      SAFE_ALLOCATE(rhs(1:der%op(1)%stencil%size, 1:der%dim + 1))
      call stencil_cube_polynomials_lapl(der%dim, der%order, polynomials)

      do i = 1, der%dim
        call get_rhs_grad(i, rhs(:,i))
      end do
      call get_rhs_lapl(rhs(:, der%dim+1))

      call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, der%dim+1, der%op(:))

      SAFE_DEALLOCATE_A(polynomials)
      SAFE_DEALLOCATE_A(rhs)

    case(DER_STARPLUS)
      do i = 1, der%dim
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%size))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
        call stencil_starplus_pol_grad(der%dim, i, der%order, polynomials)
        call get_rhs_grad(i, rhs(:, 1))
        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(i:i))
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do
      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(der%dim+1)%stencil%size))
      SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
      call stencil_starplus_pol_lapl(der%dim, der%order, polynomials)
      call get_rhs_lapl(rhs(:, 1))
      call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(der%dim+1:der%dim+1))
      SAFE_DEALLOCATE_A(polynomials)
      SAFE_DEALLOCATE_A(rhs)

    case(DER_VARIATIONAL)
      ! we have the explicit coefficients
      call stencil_variational_coeff_lapl(der%dim, der%order, mesh%h, der%lapl, alpha = der%lapl_cutoff)

    end select

    ! Here the Laplacian is forced to be self-adjoint, and the gradient to be skew-selfadjoint
    if(mesh%use_curvilinear .and. (.not. der%mesh%sb%mr_flag)) then
      do i = 1, der%dim
        call nl_operator_init(auxop, "auxop")
        call nl_operator_skewadjoint(der%grad(i), auxop, der%mesh)

        call nl_operator_end(der%grad(i))
        call nl_operator_copy(der%grad(i), auxop)
        call nl_operator_end(auxop)
      end do
      call nl_operator_init(auxop, "auxop")
      call nl_operator_selfadjoint(der%lapl, auxop, der%mesh)

      call nl_operator_end(der%lapl)
      call nl_operator_copy(der%lapl, auxop)
      call nl_operator_end(auxop)
    end if

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine get_rhs_lapl(rhs)
      FLOAT, intent(out) :: rhs(:)

      integer :: i, j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do i = 1, der%dim
        do j = 1, der%lapl%stencil%size
          this_one = .true.
          do k = 1, der%dim
            if(k == i .and. polynomials(k, j).ne.2) this_one = .false.
            if(k.ne.i .and. polynomials(k, j).ne.0) this_one = .false.
          end do
          if(this_one) rhs(j) = M_TWO
        end do
      end do

    end subroutine get_rhs_lapl

    ! ---------------------------------------------------------
    subroutine get_rhs_grad(dir, rhs)
      integer, intent(in)  :: dir
      FLOAT,   intent(out) :: rhs(:)

      integer :: j, k
      logical :: this_one

      ! find right hand side for operator
      rhs(:) = M_ZERO
      do j = 1, der%grad(dir)%stencil%size
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
  subroutine derivatives_make_discretization(dim, mesh, masses, pol, rhs, n, op)
    integer,                intent(in)    :: dim
    type(mesh_t),           intent(in)    :: mesh
    FLOAT,                  intent(in)    :: masses(:)
    integer,                intent(in)    :: pol(:,:)
    integer,                intent(in)    :: n
    FLOAT,                  intent(inout) :: rhs(:,:)
    type(nl_operator_t),    intent(inout) :: op(:)

    integer :: p, p_max, i, j, k, pow_max
    FLOAT   :: x(MAX_DIM)
    FLOAT, allocatable :: mat(:,:), sol(:,:), powers(:,:)

    call push_sub('derivatives.derivatives_make_discretization')

    SAFE_ALLOCATE(mat(1:op(1)%stencil%size, 1:op(1)%stencil%size))
    SAFE_ALLOCATE(sol(1:op(1)%stencil%size, 1:n))

    message(1) = 'Info: Generating weights for finite-difference discretization.'
    call write_info(1)

    ! use to generate power lookup table
    pow_max = maxval(pol)
    SAFE_ALLOCATE(powers(1:dim, 0:pow_max))
    powers(:,:) = M_ZERO
    powers(:,0) = M_ONE

    p_max = op(1)%np
    if(op(1)%const_w) p_max = 1

    do p = 1, p_max
      mat(1,:) = M_ONE
      do i = 1, op(1)%stencil%size
        if(mesh%use_curvilinear) then
          forall(j = 1:dim) x(j) = mesh%x(p + op(1)%ri(i, op(1)%rimap(p)), j) - mesh%x(p, j)
        else
          forall(j = 1:dim) x(j) = real(op(1)%stencil%points(j, i), REAL_PRECISION)*mesh%h(j)
        end if

        forall(j = 1:dim) x(j) = x(j)*sqrt(masses(j))

        ! calculate powers
        do j = 1, dim
          powers(j, 1) = x(j)
          do k = 2, pow_max
            powers(j, k) = x(j)*powers(j, k-1)
          end do
        end do

        ! generate the matrix
        do j = 2, op(1)%stencil%size
          mat(j, i) = powers(1, pol(1, j))
          do k = 2, dim
            mat(j, i) = mat(j, i)*powers(k, pol(k, j))
          end do
        end do
      end do

      call lalg_linsyssolve(op(1)%stencil%size, n, mat, rhs, sol)
      do i = 1, n
        op(i)%w_re(:, p) = sol(:, n)
      end do

    end do

    SAFE_DEALLOCATE_A(mat)
    SAFE_DEALLOCATE_A(sol)
    SAFE_DEALLOCATE_A(powers)

    call pop_sub()
  end subroutine derivatives_make_discretization


  ! ---------------------------------------------------------
  subroutine derivatives_handle_init(this, der)
    type(derivatives_handle_t),  intent(out) :: this
    type(derivatives_t),         intent(in)  :: der

    call push_sub('derivatives.derivatives_handle_init')

    nullify(this%df, this%zf)
    nullify(this%dlapl, this%zlapl)
#ifdef HAVE_MPI
    this%parallel_in_domains = der%mesh%parallel_in_domains
    if(this%parallel_in_domains) call pv_handle_init(this%pv_h, der%mesh%vp, der%comm_method)
#endif

    call pop_sub()
  end subroutine derivatives_handle_init


  ! ---------------------------------------------------------
  subroutine derivatives_handle_end(this)
    type(derivatives_handle_t), intent(inout) :: this

    call push_sub('derivatives.derivatives_handle_end')

    nullify(this%df, this%zf)
    nullify(this%dlapl, this%zlapl)
#ifdef HAVE_MPI
    if(this%parallel_in_domains) call pv_handle_end(this%pv_h)
#endif

    call pop_sub()
  end subroutine derivatives_handle_end


#ifdef HAVE_MPI    
  ! ---------------------------------------------------------
  logical function derivatives_overlap(this) result(overlap)
    type(derivatives_t), intent(in) :: this

    call push_sub('derivatives.derivatives_overlap')

    overlap = this%comm_method /= BLOCKING  

    call pop_sub()
  end function derivatives_overlap
#endif
  
#include "undef.F90"
#include "real.F90"
#include "derivatives_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "derivatives_inc.F90"

end module derivatives_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
