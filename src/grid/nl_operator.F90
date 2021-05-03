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

module nl_operator_oct_m
  use accel_oct_m
  use batch_oct_m
  use boundaries_oct_m
  use global_oct_m
  use index_oct_m
  use iso_c_binding
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use operate_f_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use stencil_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                     &
    nl_operator_t,              &
    nl_operator_index_t,        &
    nl_operator_global_init,    &
    nl_operator_global_end,     &
    nl_operator_init,           &
    nl_operator_copy,           &
    nl_operator_build,          &
    dnl_operator_operate,       &
    znl_operator_operate,       &
    dnl_operator_operate_batch, &
    znl_operator_operate_batch, &
    dnl_operator_operate_diag,  &
    znl_operator_operate_diag,  &
    nl_operator_end,            &
    nl_operator_skewadjoint,    &
    nl_operator_selfadjoint,    &
    nl_operator_get_index,      &
    nl_operator_output_weights, &
    nl_operator_np_zero_bc,         &
    nl_operator_compact_boundaries

  type nl_operator_index_t
    private
    integer              :: nri = 0
    integer, allocatable :: imin(:)
    integer, allocatable :: imax(:)
    integer, allocatable :: ri(:, :)
  end type nl_operator_index_t

  type nl_operator_t
    private
    type(stencil_t),   public :: stencil
    type(mesh_t), pointer     :: mesh => NULL()  !< pointer to the underlying mesh
    integer,      allocatable :: nn(:)     !< the size of the stencil at each point (for curvilinear coordinates)
    integer,          public :: np = 0     !< number of points in mesh
    !> When running in parallel mode, the next three arrays are unique on each node.
    integer, allocatable, public :: index(:,:)    !< index of the points. Unique on each parallel process.
    FLOAT,   allocatable, public :: w(:,:)        !< weights. Unique on each parallel process.

    logical,          public :: const_w = .true.   !< are the weights independent of index i

    character(len=40) :: label

    !> the compressed index of grid points
    integer, public :: nri = 0
    integer, allocatable, public :: ri(:,:)
    integer, allocatable, public :: rimap(:)
    integer, allocatable, public :: rimap_inv(:)

    integer                   :: ninner = 0
    integer                   :: nouter = 0

    type(nl_operator_index_t) :: inner
    type(nl_operator_index_t) :: outer

    type(accel_kernel_t) :: kernel
    type(accel_mem_t) :: buff_imin
    type(accel_mem_t) :: buff_imax
    type(accel_mem_t) :: buff_ri
    type(accel_mem_t) :: buff_map
    type(accel_mem_t) :: buff_all
    type(accel_mem_t) :: buff_inner
    type(accel_mem_t) :: buff_outer
    type(accel_mem_t) :: buff_stencil
    type(accel_mem_t) :: buff_ip_to_xyz
    type(accel_mem_t) :: buff_xyz_to_ip
  end type nl_operator_t

  integer, parameter :: &
       OP_FORTRAN = 0,  &
       OP_VEC     = 1,  &
       OP_MIN     = OP_FORTRAN, &
       OP_MAX     = OP_VEC

  integer, parameter ::     &
    OP_INVMAP    = 1,       &
    OP_MAP       = 2,       &
    OP_NOMAP     = 3

  integer, public, parameter :: OP_ALL = 3, OP_INNER = 1, OP_OUTER = 2

  logical :: compact_boundaries

  interface
    integer function op_is_available(opid, type)
      implicit none
      integer, intent(in) :: opid, type
    end function op_is_available
  end interface

  integer :: dfunction_global = -1
  integer :: zfunction_global = -1
  integer :: sfunction_global = -1
  integer :: cfunction_global = -1  
  integer :: function_opencl

contains
  
  ! ---------------------------------------------------------
  subroutine nl_operator_global_init(namespace)
    type(namespace_t),         intent(in)    :: namespace
    
    integer :: default

    PUSH_SUB(nl_operator_global_init)

    !%Variable OperateDouble
    !%Type integer
    !%Section Execution::Optimization
    !%Default optimized
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for real functions.
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End

    !%Variable OperateComplex
    !%Type integer
    !%Section Execution::Optimization
    !%Default optimized
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for complex functions. 
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End

    default = OP_VEC

    call parse_variable(namespace, 'OperateDouble', default, dfunction_global)
    if(.not.varinfo_valid_option('OperateDouble', dfunction_global)) call messages_input_error(namespace, 'OperateDouble')

    call parse_variable(namespace, 'OperateComplex', default, zfunction_global)
    if(.not.varinfo_valid_option('OperateComplex', zfunction_global)) call messages_input_error(namespace, 'OperateComplex')


    !%Variable OperateSingle
    !%Type integer
    !%Section Execution::Optimization
    !%Default optimized
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for single-precision real functions.
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End

    !%Variable OperateComplexSingle
    !%Type integer
    !%Section Execution::Optimization
    !%Default optimized
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for single-precision complex functions. 
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End
    
    call parse_variable(namespace, 'OperateSingle', OP_FORTRAN, sfunction_global)
    if(.not.varinfo_valid_option('OperateSingle', sfunction_global)) call messages_input_error(namespace, 'OperateSingle')
    
    call parse_variable(namespace, 'OperateComplexSingle', OP_FORTRAN, cfunction_global)
    if(.not.varinfo_valid_option('OperateComplexSingle', cfunction_global)) then
      call messages_input_error(namespace, 'OperateComplexSingle')
    end if

    if(accel_is_enabled()) then

      !%Variable OperateAccel
      !%Type integer
      !%Default map
      !%Section Execution::Optimization
      !%Description
      !% This variable selects the subroutine used to apply non-local
      !% operators over the grid when an accelerator device is used.
      !%Option invmap 1
      !% The standard implementation ported to OpenCL.
      !%Option map 2
      !% A different version, more suitable for GPUs.
      !%End
      call parse_variable(namespace, 'OperateAccel',  OP_MAP, function_opencl)

      call messages_obsolete_variable(namespace, 'OperateOpenCL', 'OperateAccel')

    end if

    !%Variable NLOperatorCompactBoundaries
    !%Type logical
    !%Default no
    !%Section Execution::Optimization
    !%Description
    !% (Experimental) When set to yes, for finite systems Octopus will
    !% map boundary points for finite-differences operators to a few
    !% memory locations. This increases performance, however it is
    !% experimental and has not been thoroughly tested.
    !%End

    call parse_variable(namespace, 'NLOperatorCompactBoundaries', .false., compact_boundaries)

    if(compact_boundaries) then
      call messages_experimental('NLOperatorCompactBoundaries')
    end if
      
    POP_SUB(nl_operator_global_init)
  end subroutine nl_operator_global_init

  ! ---------------------------------------------------------

  subroutine nl_operator_global_end()
    PUSH_SUB(nl_operator_global_end)

    POP_SUB(nl_operator_global_end)
  end subroutine nl_operator_global_end

  ! ---------------------------------------------------------
  subroutine nl_operator_init(op, label)
    type(nl_operator_t), intent(inout) :: op
    character(len=*),    intent(in)    :: label

    PUSH_SUB(nl_operator_init)

    op%label = label

    POP_SUB(nl_operator_init)
  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_copy(opo, opi)
    type(nl_operator_t),         intent(inout) :: opo
    type(nl_operator_t), target, intent(in)    :: opi

    PUSH_SUB(nl_operator_copy)

    ! We cannot currently copy the GPU kernel for the nl_operator
    ASSERT(.not. accel_is_enabled())

    call nl_operator_end(opo)
    call nl_operator_init(opo, opi%label)

    call stencil_copy(opi%stencil, opo%stencil)

    opo%np           =  opi%np
    opo%mesh         => opi%mesh

    SAFE_ALLOCATE_SOURCE_A(opo%nn, opi%nn)
    SAFE_ALLOCATE_SOURCE_A(opo%index, opi%index)
    SAFE_ALLOCATE_SOURCE_A(opo%w, opi%w)

    opo%const_w   = opi%const_w

    opo%nri       =  opi%nri
    ASSERT(allocated(opi%ri))

    SAFE_ALLOCATE_SOURCE_A(opo%ri, opi%ri)
    SAFE_ALLOCATE_SOURCE_A(opo%rimap, opi%rimap)
    SAFE_ALLOCATE_SOURCE_A(opo%rimap_inv, opi%rimap_inv)
    
    if(opi%mesh%parallel_in_domains) then
      opo%inner%nri = opi%inner%nri
      SAFE_ALLOCATE_SOURCE_A(opo%inner%imin, opi%inner%imin)
      SAFE_ALLOCATE_SOURCE_A(opo%inner%imax, opi%inner%imax)
      SAFE_ALLOCATE_SOURCE_A(opo%inner%ri,   opi%inner%ri)

      opo%outer%nri = opi%outer%nri
      SAFE_ALLOCATE_SOURCE_A(opo%outer%imin, opi%outer%imin)
      SAFE_ALLOCATE_SOURCE_A(opo%outer%imax, opi%outer%imax)
      SAFE_ALLOCATE_SOURCE_A(opo%outer%ri,   opi%outer%ri)
    end if


    POP_SUB(nl_operator_copy)
  end subroutine nl_operator_copy


  ! ---------------------------------------------------------
  subroutine nl_operator_build(space, mesh, op, np, const_w)
    type(space_t),        intent(in)    :: space
    type(mesh_t), target, intent(in)    :: mesh
    type(nl_operator_t),  intent(inout) :: op
    integer,              intent(in)    :: np       !< Number of (local) points.
    logical, optional,    intent(in)    :: const_w  !< are the weights constant (independent of the point)

    integer :: ii, jj, p1(MAX_DIM), time, current
    integer, allocatable :: st1(:), st2(:), st1r(:)
    integer :: nn
    integer :: ir, maxp, iinner, iouter
    logical :: change, force_change
    character(len=200) :: flags
    integer, allocatable :: inner_points(:), outer_points(:), all_points(:)    
    
    PUSH_SUB(nl_operator_build)

    if(mesh%parallel_in_domains .and. .not. const_w) then
      call messages_experimental('Domain parallelization with curvilinear coordinates')
    end if

    ASSERT(np > 0)

    ! store values in structure
    op%np       = np
    op%mesh     => mesh
    op%const_w  = .false.
    if(present(const_w )) op%const_w  = const_w

    ! allocate weights op%w
    if(op%const_w) then
      SAFE_ALLOCATE(op%w(1:op%stencil%size, 1:1))
      if(debug%info) then
        message(1) = 'Info: nl_operator_build: working with constant weights.'
        call messages_info(1)
      end if
    else
      SAFE_ALLOCATE(op%w(1:op%stencil%size, 1:op%np))
      if(debug%info) then
        message(1) = 'Info: nl_operator_build: working with non-constant weights.'
        call messages_info(1)
      end if
    end if

    ! set initially to zero
    op%w = M_ZERO

    ! Build lookup table
    SAFE_ALLOCATE(st1(1:op%stencil%size))
    SAFE_ALLOCATE(st1r(1:op%stencil%size))
    SAFE_ALLOCATE(st2(1:op%stencil%size))

    op%nri = 0
    do time = 1, 2
      st2 = 0
      do ii = 1, np
        p1 = 0
        call mesh_local_index_to_coords(mesh, ii, p1)

        do jj = 1, op%stencil%size
          ! Get local index of p1 plus current stencil point.
          st1(jj) = mesh_local_index_from_coords(mesh, p1(1:MAX_DIM) + op%stencil%points(1:MAX_DIM, jj))

          ! if boundary conditions are zero, we can remap boundary
          ! points to reduce memory accesses. We cannot do this for the
          ! first point, since it is used to build the weights, so it
          ! has to have the positions right
          if(ii > 1 .and. compact_boundaries .and. mesh_compact_boundaries(mesh) .and. .not. space%is_periodic()) then
            st1(jj) = min(st1(jj), mesh%np + 1)
          end if
          ASSERT(st1(jj) > 0)
        end do

        st1(1:op%stencil%size) = st1(1:op%stencil%size) - ii

        change = any(st1 /= st2) 

        !the next is to detect when we move from a point that does not
        !have boundary points as neighbours to one that has
        force_change = any(st1 + ii > mesh%np) .and. all(st2 + ii - 1 <= mesh%np) 

        if(change .and. compact_boundaries .and. mesh_compact_boundaries(mesh) .and. .not. space%is_periodic()) then 
          !try to repair it by changing the boundary points
          do jj = 1, op%stencil%size
            if(st1(jj) + ii > mesh%np .and. st2(jj) + ii - 1 > mesh%np .and. st2(jj) + ii <= mesh%np_part) then
              st1r(jj) = st2(jj)
            else
              st1r(jj) = st1(jj)
            end if
          end do

          change = any(st1r /= st2)

          if(.not. change) st1 = st1r
        end if

        ! if the stencil changes
        if (change .or. force_change) then 
          !store it
          st2 = st1

          !first time, just count
          if ( time == 1 ) op%nri = op%nri + 1

          !second time, store
          if( time == 2 ) then
            current = current + 1
            op%ri(1:op%stencil%size, current) = st1(1:op%stencil%size)
          end if
        end if

        if(time == 2) op%rimap(ii) = current

      end do

      !after counting, allocate
      if (time == 1 ) then 
        SAFE_ALLOCATE(op%ri(1:op%stencil%size, 1:op%nri))
        SAFE_ALLOCATE(op%rimap(1:op%np))
        SAFE_ALLOCATE(op%rimap_inv(1:op%nri + 1))
        op%ri        = 0
        op%rimap     = 0
        op%rimap_inv = 0
        current      = 0

        ! the sizes
        if(mesh%use_curvilinear) then
          SAFE_ALLOCATE(op%nn(1:op%nri))
          ! for the moment all the sizes are the same
          op%nn = op%stencil%size
        end if
      end if

    end do

    !the inverse mapping
    op%rimap_inv(1) = 0
    do jj = 1, op%np
      op%rimap_inv(op%rimap(jj) + 1) = jj
    end do
    op%rimap_inv(op%nri + 1) = op%np

    SAFE_DEALLOCATE_A(st1)
    SAFE_DEALLOCATE_A(st1r)
    SAFE_DEALLOCATE_A(st2)

    do jj = 1, op%nri
      nn = op%rimap_inv(jj + 1) - op%rimap_inv(jj)
    end do

    if(op%mesh%parallel_in_domains) then
      !now build the arrays required to apply the nl_operator by parts

      !count points
      op%inner%nri = 0
      op%outer%nri = 0
      do ir = 1, op%nri
        maxp = op%rimap_inv(ir + 1) + maxval(op%ri(1:op%stencil%size, ir))
        if (maxp <= np) then
          !inner point
          op%inner%nri = op%inner%nri + 1
          ASSERT(op%inner%nri <= op%nri)
        else
          !outer point
          op%outer%nri = op%outer%nri + 1
          ASSERT(op%outer%nri <= op%nri)
        end if
      end do
      
      ASSERT(op%inner%nri + op%outer%nri == op%nri)
      
      SAFE_ALLOCATE(op%inner%imin(1:op%inner%nri + 1))
      SAFE_ALLOCATE(op%inner%imax(1:op%inner%nri))
      SAFE_ALLOCATE(op%inner%ri(1:op%stencil%size, 1:op%inner%nri))

      SAFE_ALLOCATE(op%outer%imin(1:op%outer%nri + 1))
      SAFE_ALLOCATE(op%outer%imax(1:op%outer%nri))
      SAFE_ALLOCATE(op%outer%ri(1:op%stencil%size, 1:op%outer%nri))

      !now populate the arrays
      iinner = 0
      iouter = 0
      do ir = 1, op%nri
        maxp = op%rimap_inv(ir + 1) + maxval(op%ri(1:op%stencil%size, ir))
        if (maxp <= np) then
          !inner point
          iinner = iinner + 1
          op%inner%imin(iinner) = op%rimap_inv(ir)
          op%inner%imax(iinner) = op%rimap_inv(ir + 1)
          op%inner%ri(1:op%stencil%size, iinner) = op%ri(1:op%stencil%size, ir)
        else
          !outer point
          iouter = iouter + 1
          op%outer%imin(iouter) = op%rimap_inv(ir)
          op%outer%imax(iouter) = op%rimap_inv(ir + 1)
          op%outer%ri(1:op%stencil%size, iouter) = op%ri(1:op%stencil%size, ir)
        end if
      end do
      
      !verify that all points in the inner operator are actually inner
      do ir = 1, op%inner%nri
        do ii = op%inner%imin(ir) + 1, op%inner%imax(ir)
          ASSERT(all(ii + op%inner%ri(1:op%stencil%size, ir) <= mesh%np))
        end do
      end do
      
    end if

    if(accel_is_enabled() .and. op%const_w) then

      write(flags, '(i5)') op%stencil%size
      flags='-DNDIM=3 -DSTENCIL_SIZE='//trim(adjustl(flags))

      if(op%mesh%parallel_in_domains) flags = '-DINDIRECT '//trim(flags)
      
      select case(function_opencl)
      case(OP_INVMAP)
        call accel_kernel_build(op%kernel, 'operate.cl', 'operate', flags)
      case(OP_MAP)
        call accel_kernel_build(op%kernel, 'operate.cl', 'operate_map', flags)
      end select

      call accel_create_buffer(op%buff_ri, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, op%nri*op%stencil%size)
      call accel_write_buffer(op%buff_ri, op%nri*op%stencil%size, op%ri)

      select case(function_opencl)
      case(OP_INVMAP)
        call accel_create_buffer(op%buff_imin, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, op%nri)
        call accel_write_buffer(op%buff_imin, op%nri, op%rimap_inv(1:))
        call accel_create_buffer(op%buff_imax, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, op%nri)
        call accel_write_buffer(op%buff_imax, op%nri, op%rimap_inv(2:))

      case(OP_MAP)

        call accel_create_buffer(op%buff_map, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pad(op%mesh%np, accel_max_workgroup_size()))
        call accel_write_buffer(op%buff_map, op%mesh%np, (op%rimap - 1)*op%stencil%size)

        if(op%mesh%parallel_in_domains) then
          
          SAFE_ALLOCATE(inner_points(1:op%mesh%np))
          SAFE_ALLOCATE(outer_points(1:op%mesh%np))
          SAFE_ALLOCATE(all_points(1:op%mesh%np))
          
          op%ninner = 0
          op%nouter = 0
          
          do ii = 1, op%mesh%np
            all_points(ii) = ii - 1
            maxp = ii + maxval(op%ri(1:op%stencil%size, op%rimap(ii)))
            if(maxp <= op%mesh%np) then
              op%ninner = op%ninner + 1
              inner_points(op%ninner) = ii - 1
            else
              op%nouter = op%nouter + 1
              outer_points(op%nouter) = ii - 1
            end if
          end do
          
          call accel_create_buffer(op%buff_all, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pad(op%mesh%np, accel_max_workgroup_size()))
          call accel_write_buffer(op%buff_all, op%mesh%np, all_points)
        
          call accel_create_buffer(op%buff_inner, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pad(op%ninner, accel_max_workgroup_size()))
          call accel_write_buffer(op%buff_inner, op%ninner, inner_points)
          
          call accel_create_buffer(op%buff_outer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, pad(op%nouter, accel_max_workgroup_size()))
          call accel_write_buffer(op%buff_outer, op%nouter, outer_points)

          SAFE_DEALLOCATE_A(inner_points)
          SAFE_DEALLOCATE_A(outer_points)
          SAFE_DEALLOCATE_A(all_points)
  
        end if
      end select
    end if

    POP_SUB(nl_operator_build)

  end subroutine nl_operator_build

  ! ---------------------------------------------------------
  subroutine nl_operator_output_weights(this)
    type(nl_operator_t), intent(inout)  :: this

    integer :: istencil, idir

    PUSH_SUB(nl_operator_output_weights)

    if(debug%info) then

      write(message(1), '(3a)') 'Debug info: Finite difference weights for ', trim(this%label), '.'
      write(message(2), '(a)')  '            Spacing:'
      do idir = 1, this%mesh%sb%dim
        write(message(2), '(a,f16.8)') trim(message(2)), this%mesh%spacing(idir)
      end do
      call messages_info(2)
      
      do istencil = 1, this%stencil%size
        write(message(1), '(a,i3,3i4,f25.10)') '      ', istencil, this%stencil%points(1:3, istencil), this%w(istencil, 1)
        call messages_info(1)
      end do
      
    end if

    POP_SUB(nl_operator_output_weights)

  end subroutine nl_operator_output_weights

  ! ---------------------------------------------------------
  !> opt has to be initialised and built.
  subroutine nl_operator_skewadjoint(op, opt, mesh)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        target, intent(in)  :: mesh

    integer          :: ip, jp, kp, lp, index
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    PUSH_SUB(nl_operator_skewadjoint)

    call nl_operator_copy(opt, op)

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      call nl_operator_copy(opgt, opg)
      SAFE_ALLOCATE(vol_pp(1:mesh%np_global))
      call vec_allgather(mesh%vp, vol_pp, mesh%vol_pp)
    else
      opg  => op
      opgt => opt
      vol_pp => mesh%vol_pp
    end if

    opgt%w = M_ZERO
    do ip = 1, mesh%np_global
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(opg, jp, ip)
        if(index <= mesh%np_global) then
          do lp = 1, op%stencil%size
            kp = nl_operator_get_index(opg, lp, index)
            if( kp == ip ) then
              if(.not.op%const_w) then
                opgt%w(jp, ip) = M_HALF*opg%w(jp, ip) - M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w(lp, index)
              else
                opgt%w(jp, 1) = opg%w(lp, 1)
              end if
            end if
          end do
        end if
      end do
    end do

    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do ip = 1, mesh%np
        opt%w(:, ip) = opgt%w(:, mesh_local2global(mesh, ip))
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if

    POP_SUB(nl_operator_skewadjoint)
  end subroutine nl_operator_skewadjoint


  ! ---------------------------------------------------------
  subroutine nl_operator_selfadjoint(op, opt, mesh)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        target, intent(in)  :: mesh

    integer          :: ip, jp, kp, lp, index
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    PUSH_SUB(nl_operator_selfadjoint)

    call nl_operator_copy(opt, op)

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      opgt = opg
      SAFE_ALLOCATE(vol_pp(1:mesh%np_global))
      call vec_allgather(mesh%vp, vol_pp, mesh%vol_pp)
    else      
      opg  => op
      opgt => opt
      vol_pp => mesh%vol_pp
    end if

    opgt%w = M_ZERO
    do ip = 1, mesh%np_global
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(opg, jp, ip)

        if(index <= mesh%np_global) then
          do lp = 1, op%stencil%size
            kp = nl_operator_get_index(opg, lp, index)
            if( kp == ip ) then
              if(.not.op%const_w) then
                opgt%w(jp, ip) = M_HALF*opg%w(jp, ip) + M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w(lp, index)
              else
                opgt%w(jp, 1) = opg%w(lp, 1)
              end if
            end if
          end do
        end if

      end do
    end do

    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do ip = 1, mesh%np
        opt%w(:, ip) = opgt%w(:, mesh_local2global(mesh, ip))
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if

    POP_SUB(nl_operator_selfadjoint)
  end subroutine nl_operator_selfadjoint


  ! ---------------------------------------------------------
  !> Like nl_operator_gather but opg is present on all nodes
  !! (so do not forget to call nl_operator_end on all nodes
  !! afterwards).
  subroutine nl_operator_allgather(op, opg)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opg

    integer :: ip

    PUSH_SUB(nl_operator_allgather)

    ! Copy elements of op to opg that
    ! are independent from the partitions, i.e. everything
    ! except op%index and -- in the non-constant case -- op%w
    call nl_operator_common_copy(op, opg)

    ! Weights have to be collected only if they are not constant.
    if(.not.op%const_w) then
      do ip = 1, op%stencil%size
        call vec_allgather(op%mesh%vp, opg%w(ip, :), op%w(ip, :))
      end do
    end if

    POP_SUB(nl_operator_allgather)

  end subroutine nl_operator_allgather

  ! ---------------------------------------------------------
  ! The following are private routines.
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Copies all parts of op to opg that are independent of
  !! the partitions, i.e. everything except op%index and -- in the
  !! non-constant case -- op%w
  !! This can be considered as nl_operator_copy and
  !! reallocating w and i.
  !! \warning: this should be replaced by a normal copy with a flag.
  subroutine nl_operator_common_copy(op, opg)
    type(nl_operator_t), target, intent(in)    :: op
    type(nl_operator_t),         intent(inout) :: opg

    PUSH_SUB(nl_operator_common_copy)

    call nl_operator_end(opg)
    call nl_operator_init(opg, op%label)

    call stencil_copy(op%stencil, opg%stencil)

    if(op%const_w) then
      SAFE_ALLOCATE(opg%w(1:op%stencil%size, 1:1))
    else
      SAFE_ALLOCATE(opg%w(1:op%stencil%size, 1:op%mesh%np_global))
    end if
    opg%mesh     => op%mesh
    opg%np       =  op%mesh%np_global
    opg%const_w  =  op%const_w
    opg%nri      =  op%nri
    if(op%const_w) then
      opg%w = op%w
    end if

    POP_SUB(nl_operator_common_copy)

  end subroutine nl_operator_common_copy


  ! ---------------------------------------------------------
  ! End of private routines.
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_t), intent(inout) :: op

    PUSH_SUB(nl_operator_end)

    if(accel_is_enabled() .and. op%const_w) then

      call accel_release_buffer(op%buff_ri)
      select case(function_opencl)
      case(OP_INVMAP)
        call accel_release_buffer(op%buff_imin)
        call accel_release_buffer(op%buff_imax)

      case(OP_MAP)
        call accel_release_buffer(op%buff_map)
        if(op%mesh%parallel_in_domains) then
          call accel_release_buffer(op%buff_all)
          call accel_release_buffer(op%buff_inner)
          call accel_release_buffer(op%buff_outer)
        end if

      case(OP_NOMAP)
        call accel_release_buffer(op%buff_map)
        call accel_release_buffer(op%buff_stencil)
        call accel_release_buffer(op%buff_xyz_to_ip)
        call accel_release_buffer(op%buff_ip_to_xyz)
      end select
    end if

    SAFE_DEALLOCATE_A(op%inner%imin)
    SAFE_DEALLOCATE_A(op%inner%imax)
    SAFE_DEALLOCATE_A(op%inner%ri)
    SAFE_DEALLOCATE_A(op%outer%imin)
    SAFE_DEALLOCATE_A(op%outer%imax)
    SAFE_DEALLOCATE_A(op%outer%ri)

    SAFE_DEALLOCATE_A(op%index)
    SAFE_DEALLOCATE_A(op%w)

    SAFE_DEALLOCATE_A(op%ri)
    SAFE_DEALLOCATE_A(op%rimap)
    SAFE_DEALLOCATE_A(op%rimap_inv)
    SAFE_DEALLOCATE_A(op%nn)

    call stencil_end(op%stencil)

    POP_SUB(nl_operator_end)
  end subroutine nl_operator_end


  ! ---------------------------------------------------------
  integer pure function nl_operator_get_index(op, is, ip) result(res)
    type(nl_operator_t), intent(in)   :: op
    integer,             intent(in)   :: is
    integer,             intent(in)   :: ip
    
    res = ip + op%ri(is, op%rimap(ip))
  end function nl_operator_get_index

  ! ---------------------------------------------------------
  
  integer pure function nl_operator_np_zero_bc(op) result(np_bc)
    type(nl_operator_t), intent(in)   :: op

    integer :: jj, ii

    np_bc = 0
    do jj = 1, op%nri
      ii = op%rimap_inv(jj + 1) + maxval(op%ri(1:op%stencil%size, jj))
      np_bc = max(np_bc, ii)
    end do

  end function nl_operator_np_zero_bc

  ! ------------------------------------------------------

  logical pure function nl_operator_compact_boundaries(op)
    type(nl_operator_t), intent(in)   :: op

    nl_operator_compact_boundaries = compact_boundaries
  end function nl_operator_compact_boundaries


#include "undef.F90"
#include "real.F90"
#include "nl_operator_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "nl_operator_inc.F90"

end module nl_operator_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
