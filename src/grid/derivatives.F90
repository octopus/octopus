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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module derivatives_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use global_oct_m
  use iso_c_binding
  use lalg_adv_oct_m
  use lattice_vectors_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use stencil_cube_oct_m
  use stencil_star_oct_m
  use stencil_starplus_oct_m
  use stencil_stargeneral_oct_m
  use stencil_variational_oct_m
  use transfer_table_oct_m
  use types_oct_m
  use utils_oct_m
  use varinfo_oct_m

!   debug purposes 
!   use io_binary_oct_m
!   use io_function_oct_m
!   use io_oct_m
!   use unit_oct_m
!   use unit_system_oct_m
  
  implicit none

  private
  public ::                             &
    derivatives_t,                      &
    derivatives_init,                   &
    derivatives_end,                    &
    derivatives_build,                  &
    derivatives_handle_batch_t,         &
    dderivatives_test,                  &
    zderivatives_test,                  &
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
    dderivatives_batch_grad,            &
    zderivatives_batch_grad,            &
    dderivatives_div,                   &
    zderivatives_div,                   &
    dderivatives_curl,                  &
    zderivatives_curl,                  &
    dderivatives_batch_curl,            &
    zderivatives_batch_curl,            &
    dderivatives_partial,               &
    zderivatives_partial,               &
    derivatives_get_lapl


  integer, parameter ::     &
    DER_BC_ZERO_F    = 0,   &  !< function is zero at the boundaries
    DER_BC_ZERO_DF   = 1,   &  !< first derivative of the function is zero
    DER_BC_PERIOD    = 2       !< boundary is periodic

  integer, parameter ::     &
    DER_STAR         = 1,   &
    DER_VARIATIONAL  = 2,   &
    DER_CUBE         = 3,   &
    DER_STARPLUS     = 4,   &
    DER_STARGENERAL  = 5

  integer, parameter ::     &
    BLOCKING = 1,           &
    NON_BLOCKING = 2 

  type derivatives_t
    ! Components are public by default
    type(boundaries_t)    :: boundaries
    type(mesh_t), pointer :: mesh => NULL()   !< pointer to the underlying mesh
    integer               :: dim = 0          !< dimensionality of the space (space%dim)
    integer               :: order = 0        !< order of the discretization (value depends on stencil)
    integer               :: stencil_type = 0 !< type of discretization

    FLOAT                 :: masses(MAX_DIM) = M_ZERO !< we can have different weights (masses) per space direction

    !> If the so-called variational discretization is used, this controls a
    !! possible filter on the Laplacian.
    FLOAT, private :: lapl_cutoff = M_ZERO

    type(nl_operator_t), allocatable, private :: op(:)  !< op(1:conf%dim) => gradient
                                                        !! op(conf%dim+1) => Laplacian
    type(nl_operator_t), pointer :: lapl => NULL()      !< these are just shortcuts for op
    type(nl_operator_t), pointer :: grad(:) => NULL()

    integer                      :: n_ghost(MAX_DIM) = 0   !< ghost points to add in each dimension
#if defined(HAVE_MPI)
    integer, private             :: comm_method = 0
#endif
    type(derivatives_t),    pointer :: finer  => NULL()
    type(derivatives_t),    pointer :: coarser => NULL()
    type(transfer_table_t), pointer :: to_finer => NULL()
    type(transfer_table_t), pointer :: to_coarser => NULL()
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
    FLOAT                        :: factor
  end type derivatives_handle_batch_t

  type(accel_kernel_t) :: kernel_uvw_xyz, kernel_dcurl, kernel_zcurl

contains

  ! ---------------------------------------------------------
  subroutine derivatives_init(der, namespace, space, latt, use_curvilinear, order)
    type(derivatives_t), target, intent(inout) :: der
    type(namespace_t),           intent(in)    :: namespace
    type(space_t),               intent(in)    :: space
    type(lattice_vectors_t),     intent(in)    :: latt
    logical,                     intent(in)    :: use_curvilinear
    integer, optional,           intent(in)    :: order

    integer :: idir
    integer :: default_stencil
    character(len=40) :: flags

    PUSH_SUB(derivatives_init)

    ! copy this value to my structure
    der%dim = space%dim

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
    !% since the cube typically needs far too much memory.
    !%Option stencil_star 1
    !% A star around each point (<i>i.e.</i>, only points on the axis).
    !%Option stencil_variational 2
    !% Same as the star, but with coefficients built in a different way.
    !%Option stencil_cube 3
    !% A cube of points around each point.
    !%Option stencil_starplus 4
    !% The star, plus a number of off-axis points.
    !%Option stencil_stargeneral 5
    !% The general star. Default for non-orthogonal grids.
    !%End
    default_stencil = DER_STAR
    if(use_curvilinear) default_stencil = DER_STARPLUS
    if (latt%nonorthogonal) default_stencil = DER_STARGENERAL

    call parse_variable(namespace, 'DerivativesStencil', default_stencil, der%stencil_type)
    
    if(.not.varinfo_valid_option('DerivativesStencil', der%stencil_type)) then
      call messages_input_error(namespace, 'DerivativesStencil')
    end if
    call messages_print_var_option(stdout, "DerivativesStencil", der%stencil_type)

    if(use_curvilinear  .and.  der%stencil_type < DER_CUBE) call messages_input_error(namespace, 'DerivativesStencil')
    if(der%stencil_type == DER_VARIATIONAL) then
      call parse_variable(namespace, 'DerivativesLaplacianFilter', M_ONE, der%lapl_cutoff)
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
    !% depends on the stencil used. Let <math>O</math> = <tt>DerivativesOrder</tt>, and <math>d</math> = <tt>Dimensions</tt>.
    !% <ul>
    !% <li> <tt>stencil_star</tt>: <math>2 O d + 1</math>
    !% <li> <tt>stencil_cube</tt>: <math>(2 O + 1)^d</math>
    !% <li> <tt>stencil_starplus</tt>: <math>2 O d + 1 + n</math> with <i>n</i> being 8
    !% in 2D and 24 in 3D.
    !% </ul>
    !%End
    call parse_variable(namespace, 'DerivativesOrder', 4, der%order)
    ! overwrite order if given as argument
    if(present(order)) then
      der%order = order
    end if

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
    !% Communication is based on non-blocking point-to-point communication.
    !%End
    
    call parse_variable(namespace, 'ParallelizationOfDerivatives', NON_BLOCKING, der%comm_method)
    
    if(.not. varinfo_valid_option('ParallelizationOfDerivatives', der%comm_method)) then
      call messages_input_error(namespace, 'ParallelizationOfDerivatives')
    end if

    call messages_obsolete_variable(namespace, 'OverlapDerivatives', 'ParallelizationOfDerivatives')
#endif

    ! if needed, der%masses should be initialized in modelmb_particles_init
    der%masses = M_ONE

    ! construct lapl and grad structures
    SAFE_ALLOCATE(der%op(1:der%dim + 1))
    der%grad => der%op
    der%lapl => der%op(der%dim + 1)

    call derivatives_get_stencil_lapl(der, space, latt)
    call derivatives_get_stencil_grad(der)

    ! find out how many ghost points we need in each dimension
    der%n_ghost(:) = 0
    do idir = 1, der%dim
      der%n_ghost(idir) = maxval(abs(der%lapl%stencil%points(idir, :)))
    end do

    nullify(der%coarser)
    nullify(der%finer)
    nullify(der%to_coarser)
    nullify(der%to_finer)

    if(accel_is_enabled()) then
      write(flags, '(A,I1.1)') ' -DDIMENSION=', der%dim
      call accel_kernel_build(kernel_uvw_xyz, 'uvw_to_xyz.cl', 'uvw_to_xyz', flags)
      call accel_kernel_build(kernel_dcurl, 'curl.cl', 'dcurl', flags = '-DRTYPE_DOUBLE')
      call accel_kernel_build(kernel_zcurl, 'curl.cl', 'zcurl', flags = '-DRTYPE_COMPLEX')
    end if

    POP_SUB(derivatives_init)
  end subroutine derivatives_init


  ! ---------------------------------------------------------
  subroutine derivatives_end(der)
    type(derivatives_t), intent(inout) :: der

    integer :: idim

    PUSH_SUB(derivatives_end)

    ASSERT(allocated(der%op))

    do idim = 1, der%dim+1
      call nl_operator_end(der%op(idim))
    end do

    SAFE_DEALLOCATE_A(der%op)
    nullify(der%lapl, der%grad)

    nullify(der%coarser)
    nullify(der%finer)
    nullify(der%to_coarser)
    nullify(der%to_finer)

    call boundaries_end(der%boundaries)

    POP_SUB(derivatives_end)
  end subroutine derivatives_end


  ! ---------------------------------------------------------
  subroutine derivatives_get_stencil_lapl(der, space, latt)
    type(derivatives_t),     intent(inout) :: der
    type(space_t),           intent(in)    :: space
    type(lattice_vectors_t), intent(in)    :: latt

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
    case(DER_STARGENERAL)
      call stencil_stargeneral_get_arms(der%lapl%stencil, space, latt)
      call stencil_stargeneral_get_lapl(der%lapl%stencil, der%dim, der%order)
    end select

    POP_SUB(derivatives_get_stencil_lapl)
  end subroutine derivatives_get_stencil_lapl


  ! ---------------------------------------------------------
  !> Returns the diagonal elements of the Laplacian, needed for preconditioning
  subroutine derivatives_lapl_diag(der, lapl)
    type(derivatives_t), intent(in)  :: der
    FLOAT,               intent(out) :: lapl(:)  !< lapl(mesh%np)

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
      case(DER_STARGENERAL)
      ! use the simple star stencil
        call stencil_star_get_grad(der%grad(ii)%stencil, ii, der%order)
      end select
    end do

    POP_SUB(derivatives_get_stencil_grad)

  end subroutine derivatives_get_stencil_grad

  ! ---------------------------------------------------------
  subroutine derivatives_build(der, namespace, space, mesh)
    type(derivatives_t),    intent(inout) :: der
    type(namespace_t),      intent(in)    :: namespace
    type(space_t),          intent(in)    :: space
    type(mesh_t),   target, intent(in)    :: mesh

    integer, allocatable :: polynomials(:,:)
    FLOAT,   allocatable :: rhs(:,:)
    integer :: i
    logical :: const_w_
    character(len=32) :: name
    type(nl_operator_t) :: auxop
    integer :: np_zero_bc

    PUSH_SUB(derivatives_build)

    call boundaries_init(der%boundaries, namespace, space, mesh)

    ASSERT(allocated(der%op))
    ASSERT(der%stencil_type>=DER_STAR .and. der%stencil_type<=DER_STARGENERAL)
    ASSERT(.not.(der%stencil_type==DER_VARIATIONAL .and. mesh%use_curvilinear))

    der%mesh => mesh    ! make a pointer to the underlying mesh

    const_w_  = .true.

    ! need non-constant weights for curvilinear and scattering meshes
    if(mesh%use_curvilinear) const_w_ = .false.

    np_zero_bc = 0

    ! build operators
    do i = 1, der%dim+1
      call nl_operator_build(space, mesh, der%op(i), der%mesh%np, const_w = const_w_)
      np_zero_bc = max(np_zero_bc, nl_operator_np_zero_bc(der%op(i)))
    end do

    ASSERT(np_zero_bc > mesh%np .and. np_zero_bc <= mesh%np_part)

    select case(der%stencil_type)

    case(DER_STAR) ! Laplacian and gradient have different stencils
      do i = 1, der%dim + 1
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%npoly))
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

      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(1)%stencil%npoly))
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
        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%npoly))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
        call stencil_starplus_pol_grad(der%dim, i, der%order, polynomials)
        call get_rhs_grad(i, rhs(:, 1))
        name = index2axis(i) // "-gradient"
        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(i:i), name)
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do
      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(der%dim+1)%stencil%npoly))
      SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
      call stencil_starplus_pol_lapl(der%dim, der%order, polynomials)
      call get_rhs_lapl(rhs(:, 1))
      name = "Laplacian"
      call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1, der%op(der%dim+1:der%dim+1), name)
      SAFE_DEALLOCATE_A(polynomials)
      SAFE_DEALLOCATE_A(rhs)

    case(DER_STARGENERAL)    
    
      do i = 1, der%dim        
        der%op(i)%stencil%npoly = der%op(i)%stencil%size ! for gradients we are fine

        SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(i)%stencil%npoly))
        SAFE_ALLOCATE(rhs(1:der%op(i)%stencil%size, 1:1))
        ! use simple star stencil polynomials
        call stencil_star_polynomials_grad(i, der%order, polynomials)
        call get_rhs_grad(i, rhs(:, 1))
        name = index2axis(i) // "-gradient"
        ! For directional derivatives the weights are the same as in the orthogonal case.
        ! Forcing the orthogonal case avoid incurring in ill-defined cases.
        ! NOTE: this is not so clean and also am not 100% sure is correct. It has to be tested. UDG
        call derivatives_make_discretization(der%dim, der%mesh, der%masses, polynomials, rhs, 1,&
                                              der%op(i:i), name, force_orthogonal = .true.)
        SAFE_DEALLOCATE_A(polynomials)
        SAFE_DEALLOCATE_A(rhs)
      end do

      der%op(der%dim+1)%stencil%npoly = der%op(der%dim+1)%stencil%size &
           + der%order*(2*der%order-1)*der%op(der%dim+1)%stencil%stargeneral%narms

      SAFE_ALLOCATE(polynomials(1:der%dim, 1:der%op(der%dim+1)%stencil%npoly))
      SAFE_ALLOCATE(rhs(1:der%op(der%dim+1)%stencil%size, 1:1))
      call stencil_stargeneral_pol_lapl(der%op(der%dim+1)%stencil, der%dim, der%order, polynomials)
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
    if(mesh%use_curvilinear) then
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

    ! useful for debug purposes
!     call dderivatives_test(der, 1, 1, 1)
!     call exit(1)


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
        do j = 1, der%lapl%stencil%npoly
          this_one = .true.
          do k = 1, der%dim
            if(k == i .and. polynomials(k, j) /= 2) this_one = .false.
            if(k /= i .and. polynomials(k, j) /= 0) this_one = .false.
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
      do j = 1, der%grad(dir)%stencil%npoly
        this_one = .true.
        do k = 1, der%dim
          if(k == dir .and. polynomials(k, j) /= 1) this_one = .false.
          if(k /= dir .and. polynomials(k, j) /= 0) this_one = .false.
        end do
        if(this_one) rhs(j) = M_ONE
      end do

      POP_SUB(derivatives_build.get_rhs_grad)
    end subroutine get_rhs_grad

  end subroutine derivatives_build


  ! ---------------------------------------------------------
  subroutine derivatives_make_discretization(dim, mesh, masses, pol, rhs, nderiv, op, name, force_orthogonal)
    integer,                intent(in)    :: dim
    type(mesh_t),           intent(in)    :: mesh
    FLOAT,                  intent(in)    :: masses(:)
    integer,                intent(in)    :: pol(:,:)
    integer,                intent(in)    :: nderiv
    FLOAT,                  intent(inout) :: rhs(:,:)
    type(nl_operator_t),    intent(inout) :: op(:)
    character(len=32),      intent(in)    :: name
    logical, optional,      intent(in)    :: force_orthogonal

    integer :: p, p_max, i, j, k, pow_max
    FLOAT   :: x(MAX_DIM)
    FLOAT, allocatable :: mat(:,:), sol(:,:), powers(:,:)

    PUSH_SUB(derivatives_make_discretization)

    SAFE_ALLOCATE(mat(1:op(1)%stencil%npoly, 1:op(1)%stencil%size))
    SAFE_ALLOCATE(sol(1:op(1)%stencil%size, 1:nderiv))

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
      ! first polynomial is just a constant
      mat(1,:) = M_ONE
      ! i indexes the point in the stencil
      do i = 1, op(1)%stencil%size
        if(mesh%use_curvilinear) then
          do j = 1, dim
            x(j) = mesh%x(p + op(1)%ri(i, op(1)%rimap(p)), j) - mesh%x(p, j)
          end do
        else
          do j = 1, dim
            x(j) = TOFLOAT(op(1)%stencil%points(j, i))*mesh%spacing(j)
          end do
          ! TODO : this internal if clause is inefficient - the condition is determined globally
          if (mesh%sb%latt%nonorthogonal .and. .not. optional_default(force_orthogonal, .false.))  & 
              x(1:dim) = matmul(mesh%sb%latt%rlattice_primitive(1:dim,1:dim), x(1:dim))
        end if
                         
! NB: these masses are applied on the cartesian directions. Should add a check for non-orthogonal axes
        do j = 1, dim
          x(j) = x(j)*sqrt(masses(j))
        end do

        ! calculate powers
        do j = 1, dim
          powers(j, 1) = x(j)
          do k = 2, pow_max
            powers(j, k) = x(j)*powers(j, k-1)
          end do
        end do

        ! generate the matrix
        ! j indexes the polynomial being used
        do j = 2, op(1)%stencil%npoly
          mat(j, i) = powers(1, pol(1, j))
          do k = 2, dim
            mat(j, i) = mat(j, i)*powers(k, pol(k, j))
          end do
        end do
      end do ! loop over i = point in stencil

      ! linear problem to solve for derivative weights:
      !   mat * sol = rhs

      if (op(1)%stencil%npoly==op(1)%stencil%size) then
        call lalg_linsyssolve(op(1)%stencil%size, nderiv, mat, rhs, sol)
      else
        call lalg_svd_inverse(op(1)%stencil%npoly, op(1)%stencil%size, mat, CNST(1.e-10))
        sol = matmul(transpose(mat(1:op(1)%stencil%size, 1:op(1)%stencil%size)), rhs)
      end if

      do i = 1, nderiv
!MJV 10 nov 2016: changed n to i below in sol(:,n): was only erroneous in cube case when several
!derivatives are calculated in 1 batch by this subroutine. All weights contained
!the last derivative`s weights (gradient z)
!op(1) is used systematically above to get dimensions, but here we have to save
!all operator stencil weights
! once we are happy and convinced, remove this comment
        op(i)%w(:, p) = sol(:, i)
      end do
              
    end do ! loop over points p

    do i = 1, nderiv
      call nl_operator_output_weights(op(i))
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

  ! ---------------------------------------------------------
  subroutine derivatives_get_lapl(this, op, space, name, order) 
    type(derivatives_t),         intent(in)    :: this
    type(nl_operator_t),         intent(inout) :: op(:) 
    type(space_t),               intent(in)    :: space
    character(len=32),           intent(in)    :: name
    integer,                     intent(in)    :: order

    integer, allocatable :: polynomials(:,:)
    FLOAT, allocatable :: rhs(:,:)
    integer :: i, j, k
    logical :: this_one
  
    PUSH_SUB(derivatives_get_lapl)

    call nl_operator_init(op(1), name)
    if(this%mesh%sb%latt%nonorthogonal) then
      call stencil_stargeneral_get_arms(op(1)%stencil, space, this%mesh%sb%latt)
      call stencil_stargeneral_get_lapl(op(1)%stencil, this%dim, order)
    else
      call stencil_star_get_lapl(op(1)%stencil, this%dim, order)
    end if
    call nl_operator_build(space, this%mesh, op(1), this%mesh%np, const_w = .not. this%mesh%use_curvilinear)

    !At the moment this code is almost copy-pasted from derivatives_build.
    if(this%mesh%sb%latt%nonorthogonal) then
      op(1)%stencil%npoly = op(1)%stencil%size &
         + order*(2*order-1)*op(1)%stencil%stargeneral%narms
    end if
    SAFE_ALLOCATE(polynomials(1:this%dim, 1:op(1)%stencil%npoly))
    SAFE_ALLOCATE(rhs(1:op(1)%stencil%size, 1:1))
    if(this%mesh%sb%latt%nonorthogonal) then
      call stencil_stargeneral_pol_lapl(op(1)%stencil, this%dim, order, polynomials)
    else
      call stencil_star_polynomials_lapl(this%dim, order, polynomials)
    end if
    ! find right-hand side for operator
    rhs(:,1) = M_ZERO
    do i = 1, this%dim
      do j = 1, op(1)%stencil%npoly
        this_one = .true.
        do k = 1, this%dim
          if(k == i .and. polynomials(k, j) /= 2) this_one = .false.
          if(k /= i .and. polynomials(k, j) /= 0) this_one = .false.
        end do
        if(this_one) rhs(j,1) = M_TWO
      end do
    end do
    call derivatives_make_discretization(this%dim, this%mesh, this%masses, &
             polynomials, rhs, 1, op(1:1), name, force_orthogonal = .true.)
    SAFE_DEALLOCATE_A(polynomials)
    SAFE_DEALLOCATE_A(rhs)

    POP_SUB(derivatives_get_lapl)
  end subroutine derivatives_get_lapl

  
#include "undef.F90"
#include "real.F90"
#include "derivatives_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "derivatives_inc.F90"

end module derivatives_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
