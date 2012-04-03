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
  use boundaries_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use nl_operator_m
  use opencl_m
  use par_vec_m
  use parser_m
  use profiling_m
  use simul_box_m
  use stencil_cube_m
  use stencil_star_m
  use stencil_starplus_m
  use stencil_variational_m
  use transfer_table_m
  use utils_m
  use varinfo_m

  implicit none

  private
  public ::                             &
    derivatives_t,                      &
    derivatives_init,                   &
    derivatives_end,                    &
    derivatives_build,                  &
    derivatives_stencil_extent,         &
    derivatives_handle_batch_t,         &
    dderivatives_test,                  &
    zderivatives_test,                  &
    dderivatives_set_bc,                &
    zderivatives_set_bc,                &
    dderivatives_batch_set_bc,          &
    zderivatives_batch_set_bc,          &
    dderivatives_batch_start,           &
    zderivatives_batch_start,           &
    dderivatives_batch_finish,          &
    zderivatives_batch_finish,          &
    dderivatives_batch_perform,         &
    zderivatives_batch_perform,         &
    dderivatives_perform,               &
    zderivatives_perform,               &
    dderivatives_lapl,                  &
    zderivatives_lapl,                  &
    derivatives_lapl_diag,              &
    dderivatives_grad,                  &
    zderivatives_grad,                  &
    dderivatives_div,                   &
    zderivatives_div,                   &
    dderivatives_curl,                  &
    zderivatives_curl


  integer, parameter ::     &
    DER_BC_ZERO_F    = 0,   &  !< function is zero at the boundaries
    DER_BC_ZERO_DF   = 1,   &  !< first derivative of the function is zero
    DER_BC_PERIOD    = 2       !< boundary is periodic

  integer, parameter ::     &
    DER_STAR         = 1,   &
    DER_VARIATIONAL  = 2,   &
    DER_CUBE         = 3,   &
    DER_STARPLUS     = 4

  integer, parameter ::     &
    BLOCKING = 1,           &
    NON_BLOCKING = 2 

  type derivatives_t
    type(boundaries_t)    :: boundaries
    type(mesh_t), pointer :: mesh          !< pointer to the underlying mesh
    integer               :: dim           !< dimensionality of the space (sb%dim)
    integer               :: order         !< order of the discretization (value depends on stencil)
    integer               :: stencil_type  !< type of discretization

    FLOAT                 :: masses(MAX_DIM)     !< we can have different weights (masses) per space direction

    logical               :: zero_bc
    logical               :: periodic_bc

    !> If the so-called variational discretization is used, this controls a
    !! possible filter on the Laplacian.
    FLOAT :: lapl_cutoff   

    type(nl_operator_t), pointer :: op(:)  !< op(1:conf%dim) => gradient
    !! op(conf%dim+1) => Laplacian
    type(nl_operator_t), pointer :: lapl   !< these are just shortcuts for op
    type(nl_operator_t), pointer :: grad(:)

    integer                      :: n_ghost(MAX_DIM)   !< ghost points to add in each dimension
#if defined(HAVE_MPI)
    integer                      :: comm_method 
#endif
    type(derivatives_t),    pointer :: finer
    type(derivatives_t),    pointer :: coarser
    type(transfer_table_t), pointer :: to_finer
    type(transfer_table_t), pointer :: to_coarser
  end type derivatives_t

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
    logical                      :: factor_present
    FLOAT                        :: dfactor
    CMPLX                        :: zfactor
  end type derivatives_handle_batch_t

  type(profile_t), save :: gradient_prof, divergence_prof, curl_prof
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

    integer :: idir

    PUSH_SUB(derivatives_init)

    ! copy this value to my structure
    der%dim = sb%dim

    !%Variable DerivativesStencil
    !%Type integer
    !%Default stencil_star
    !%Section Mesh::Derivatives
    !%Description
    !% Decides what kind of stencil is used, <i>i.e.</i> which points, around
    !% each point in the mesh, are the neighboring points used in the
    !% expression of the differential operator.
    !%
    !% If curvilinear coordinates are to be used, then only the <tt>stencil_starplus</tt>
    !% or the <tt>stencil_cube</tt> may be used. We only recommend the <tt>stencil_starplus</tt>,
    !% since the cube typically needs way too much memory resources.
    !%Option stencil_star 1
    !% A star around each point (<i>i.e.</i>, only points on the axis).
    !%Option stencil_variational 2
    !% Same as the star, but with coefficients built in a different way.
    !%Option stencil_cube 3
    !% A cube of points around each point.
    !%Option stencil_starplus 4
    !% The star, plus a number of off-axis points.
    !%End
    if(use_curvilinear) then
      call parse_integer(datasets_check('DerivativesStencil'), DER_STARPLUS, der%stencil_type)
    else
      call parse_integer(datasets_check('DerivativesStencil'), DER_STAR, der%stencil_type)
    endif
    if(.not.varinfo_valid_option('DerivativesStencil', der%stencil_type)) call input_error('DerivativesStencil')
    call messages_print_var_option(stdout, "DerivativesStencil", der%stencil_type)

    if(use_curvilinear  .and.  der%stencil_type < DER_CUBE) call input_error('DerivativesStencil')
    if(der%stencil_type == DER_VARIATIONAL) then
      call parse_float(datasets_check('DerivativesLaplacianFilter'), M_ONE, der%lapl_cutoff)
    end if

    !%Variable DerivativesOrder
    !%Type integer
    !%Default 4
    !%Section Mesh::Derivatives
    !%Description
    !% This variable gives the discretization order for the approximation of
    !% the differential operators. This means, basically, that
    !% <tt>DerivativesOrder</tt> points are used in each positive/negative
    !% spatial direction, <i>e.g.</i> <tt>DerivativesOrder = 1</tt> would give
    !% the well-known three-point formula in 1D.
    !% The number of points actually used for the Laplacian
    !% depends on the stencil used:
    !%
    !% <tt>stencil_star</tt>: 2*<tt>DerivativesOrder</tt>*<i>dim</i>+1
    !%
    !% <tt>stencil_cube</tt>: (2*<tt>DerivativesOrder</tt>+1)^<i>dim</i>
    !%
    !% <tt>stencil_starplus</tt>: 2*<tt>DerivativesOrder</tt>+1+<i>n</i> with <i>n</i> being 12
    !% in 2D and 44 in 3D.
    !%End
    call parse_integer(datasets_check('DerivativesOrder'), 4, der%order)

#ifdef HAVE_MPI
    !%Variable ParallelizationOfDerivatives
    !%Type integer
    !%Default non_blocking
    !%Section Execution::Parallelization
    !%Description
    !% This option selects how the communication of mesh boundaries is performed.
    !%Option blocking 1
    !% Blocking communication.
    !%Option non_blocking 2
    !% Communication is based on non-blocking point-to-point communication. This is the default.
    !%End
    
    call parse_integer(datasets_check('ParallelizationOfDerivatives'), NON_BLOCKING, der%comm_method)
    
    if(.not. varinfo_valid_option('ParallelizationOfDerivatives', der%comm_method)) then
      call input_error('ParallelizationOfDerivatives')
    end if

    call messages_obsolete_variable('OverlapDerivatives', 'ParallelizationOfDerivatives')
#endif

    ! if needed, der%masses should be initialized in modelmb_particles_init
    der%masses = M_ONE

    ! construct lapl and grad structures
    SAFE_ALLOCATE(der%op(1:der%dim + 1))
    der%grad => der%op
    der%lapl => der%op(der%dim + 1)

    call derivatives_get_stencil_lapl(der)
    call derivatives_get_stencil_grad(der)

    der%zero_bc = (sb%periodic_dim < 3)
    der%periodic_bc = (sb%periodic_dim > 0)

    ! find out how many ghost points we need in each dimension
    der%n_ghost(:) = 0
    do idir = 1, der%dim
      der%n_ghost(idir) = maxval(abs(der%lapl%stencil%points(idir, :)))
    end do

    nullify(der%coarser)
    nullify(der%finer)
    nullify(der%to_coarser)
    nullify(der%to_finer)

    POP_SUB(derivatives_init)
  end subroutine derivatives_init


  ! ---------------------------------------------------------
  subroutine derivatives_end(der)
    type(derivatives_t), intent(inout) :: der

    integer :: idim

    PUSH_SUB(derivatives_end)

    ASSERT(associated(der%op))

    do idim = 1, der%dim+1
      call nl_operator_end(der%op(idim))
    end do

    SAFE_DEALLOCATE_P(der%op)
    nullify(der%lapl, der%grad)

    nullify(der%coarser)
    nullify(der%finer)
    nullify(der%to_coarser)
    nullify(der%to_finer)

    call boundaries_end(der%boundaries)

    POP_SUB(derivatives_end)
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  !> Returns maximum extension of the stencil in spatial direction
  !! dir = 1, 2, 3 for a given derivative der.
  integer function derivatives_stencil_extent(der, dir) result(extent)
    type(derivatives_t), intent(in) :: der
    integer,             intent(in) :: dir

    PUSH_SUB(stencil_extent)

    select case(der%stencil_type)
      case(DER_STAR)
        extent = stencil_star_extent(dir, der%order)
      case(DER_VARIATIONAL)
        extent = stencil_variational_extent(dir, der%order)
      case(DER_CUBE)
        extent = stencil_cube_extent(dir, der%order)
      case(DER_STARPLUS)
        extent = stencil_cube_extent(dir, der%order)
      end select
      
    POP_SUB(stencil_extent)
  end function derivatives_stencil_extent


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der)
    type(derivatives_t), intent(inout) :: der

    PUSH_SUB(derivatives_get_stencil_lapl)

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

    POP_SUB(derivatives_get_stencil_lapl)

  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  !> Returns the diagonal elements of the Laplacian, needed for preconditioning
  subroutine derivatives_lapl_diag(der, lapl)
    type(derivatives_t), intent(in)  :: der
    FLOAT,               intent(out) :: lapl(:)  ! lapl(mesh%np)

    PUSH_SUB(derivatives_lapl_diag)

    ASSERT(ubound(lapl, DIM=1) >= der%mesh%np)

    ! the Laplacian is a real operator
    call dnl_operator_operate_diag(der%lapl, lapl)

    POP_SUB(derivatives_lapl_diag)

  end subroutine derivatives_lapl_diag


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_grad(der)
    type(derivatives_t), intent(inout) :: der

    integer  :: ii
    character :: dir_label

    PUSH_SUB(derivatives_get_stencil_grad)

    ASSERT(associated(der%grad))

    ! initialize nl operator
    do ii = 1, der%dim
      dir_label = ' '
      if(ii < 5) dir_label = index2axis(ii)

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

    POP_SUB(derivatives_get_stencil_grad)

  end subroutine derivatives_get_stencil_grad


  ! ---------------------------------------------------------
  subroutine derivatives_build(der, mesh)
    type(derivatives_t),    intent(inout) :: der
    type(mesh_t),   target, intent(in)    :: mesh

    integer, allocatable :: polynomials(:,:)
    FLOAT,   allocatable :: rhs(:,:)
    integer :: i
    logical :: const_w_, cmplx_op_
    character(len=32) :: name
    type(nl_operator_t) :: auxop

    PUSH_SUB(derivatives_build)

    call boundaries_init(der%boundaries, mesh)

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

    case(DER_STAR) ! Laplacian and gradient have different stencils
      do i = 1, der%dim + 1
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%size))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))

        if(i <= der%dim) then  ! gradient
          call stencil_star_polynomials_grad(i, der%order, polynomials)
          call get_rhs_grad(i, rhs(:,1))
          name = index2axis(i) // "-gradient"
        else                      ! Laplacian
          call stencil_star_polynomials_lapl(der%dim, der%order, polynomials)
          call get_rhs_lapl(rhs(:,1))
          name = "Laplacian"
        end if

        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(i:i), name)
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do

    case(DER_CUBE) ! Laplacian and gradient have similar stencils
      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(1)%stencil%size))
      SAFE_ALLOCATE(rhs(1:der%op(1)%stencil%size, 1:der%dim + 1))
      call stencil_cube_polynomials_lapl(der%dim, der%order, polynomials)

      do i = 1, der%dim
        call get_rhs_grad(i, rhs(:,i))
      end do
      call get_rhs_lapl(rhs(:, der%dim+1))

      name = "derivatives"
      call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, der%dim+1, der%op(:), name)

      SAFE_DEALLOCATE_A(polynomials)
      SAFE_DEALLOCATE_A(rhs)

    case(DER_STARPLUS)
      do i = 1, der%dim
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%size))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
        call stencil_starplus_pol_grad(der%dim, i, der%order, polynomials)
        call get_rhs_grad(i, rhs(:, 1))
        name = index2axis(i) // "-gradient"
        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(i:i), name)
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do
      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(der%dim+1)%stencil%size))
      SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
      call stencil_starplus_pol_lapl(der%dim, der%order, polynomials)
      call get_rhs_lapl(rhs(:, 1))
      name = "Laplacian"
      call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(der%dim+1:der%dim+1), name)
      SAFE_DEALLOCATE_A(polynomials)
      SAFE_DEALLOCATE_A(rhs)

    case(DER_VARIATIONAL)
      ! we have the explicit coefficients
      call stencil_variational_coeff_lapl(der%dim, der%order, mesh%spacing, der%lapl, alpha = der%lapl_cutoff)

    end select

    ! Here the Laplacian is forced to be self-adjoint, and the gradient to be skew-self-adjoint
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

    POP_SUB(derivatives_build)

  contains

    ! ---------------------------------------------------------
    subroutine get_rhs_lapl(rhs)
      FLOAT, intent(out) :: rhs(:)

      integer :: i, j, k
      logical :: this_one

      PUSH_SUB(derivatives_build.get_rhs_lapl)

      ! find right-hand side for operator
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

      POP_SUB(derivatives_build.get_rhs_lapl)
    end subroutine get_rhs_lapl

    ! ---------------------------------------------------------
    subroutine get_rhs_grad(dir, rhs)
      integer, intent(in)  :: dir
      FLOAT,   intent(out) :: rhs(:)

      integer :: j, k
      logical :: this_one

      PUSH_SUB(derivatives_build.get_rhs_grad)

      ! find right-hand side for operator
      rhs(:) = M_ZERO
      do j = 1, der%grad(dir)%stencil%size
        this_one = .true.
        do k = 1, der%dim
          if(k == dir .and. polynomials(k, j).ne.1) this_one = .false.
          if(k.ne.dir .and. polynomials(k, j).ne.0) this_one = .false.
        end do
        if(this_one) rhs(j) = M_ONE
      end do

      POP_SUB(derivatives_build.get_rhs_grad)
    end subroutine get_rhs_grad

  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine derivatives_make_discretization(dim, mesh, masses, pol, rhs, n, op, name)
    integer,                intent(in)    :: dim
    type(mesh_t),           intent(in)    :: mesh
    FLOAT,                  intent(in)    :: masses(:)
    integer,                intent(in)    :: pol(:,:)
    integer,                intent(in)    :: n
    FLOAT,                  intent(inout) :: rhs(:,:)
    type(nl_operator_t),    intent(inout) :: op(:)
    character*32,           intent(in)    :: name

    integer :: p, p_max, i, j, k, pow_max
    FLOAT   :: x(MAX_DIM)
    FLOAT, allocatable :: mat(:,:), sol(:,:), powers(:,:)

    PUSH_SUB(derivatives_make_discretization)

    SAFE_ALLOCATE(mat(1:op(1)%stencil%size, 1:op(1)%stencil%size))
    SAFE_ALLOCATE(sol(1:op(1)%stencil%size, 1:n))

    message(1) = 'Info: Generating weights for finite-difference discretization of ' // trim(name)
    call messages_info(1)

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
          forall(j = 1:dim) x(j) = real(op(1)%stencil%points(j, i), REAL_PRECISION)*mesh%spacing(j)
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
    do i = 1, n
      call nl_operator_update_weights(op(i))
    end do

    SAFE_DEALLOCATE_A(mat)
    SAFE_DEALLOCATE_A(sol)
    SAFE_DEALLOCATE_A(powers)

    POP_SUB(derivatives_make_discretization)
  end subroutine derivatives_make_discretization

#ifdef HAVE_MPI    
  ! ---------------------------------------------------------
  logical function derivatives_overlap(this) result(overlap)
    type(derivatives_t), intent(in) :: this

    PUSH_SUB(derivatives_overlap)

    overlap = this%comm_method /= BLOCKING  

    POP_SUB(derivatives_overlap)
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
