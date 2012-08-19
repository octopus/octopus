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

module nl_operator_m
  use batch_m
  use boundaries_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use c_pointer_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use math_m
  use index_m
  use mesh_m
  use messages_m
  use multicomm_m
  use mpi_m
  use octcl_kernel_m
  use opencl_m
  use par_vec_m
  use parser_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use types_m
  use varinfo_m

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
    nl_operator_transpose,      &
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
    nl_operator_write,          &
    nl_operator_update_weights, &
    nl_operator_op_to_matrix_cmplx

  type nl_operator_index_t
    private
    integer          :: nri
    integer, pointer :: imin(:)
    integer, pointer :: imax(:)
    integer, pointer :: ri(:, :)
  end type nl_operator_index_t

  type nl_operator_t
    type(stencil_t)       :: stencil
    type(mesh_t), pointer :: mesh      !< pointer to the underlying mesh
    integer, pointer      :: nn(:)     !< the size of the stencil at each point (for curvilinear coordinates)
    integer               :: np        !< number of points in mesh
    ! When running in parallel mode, the next three
    ! arrays are unique on each node.
    integer, pointer  :: i(:,:)    !< index of the points
    FLOAT,   pointer  :: w_re(:,:) !< weightsp, real part
    FLOAT,   pointer  :: w_im(:,:) !< weightsp, imaginary part

    logical               :: const_w   !< are the weights independent of index i
    logical               :: cmplx_op  !< .true. if we have also imaginary weights

    character(len=40) :: label

    !> the compressed index of grid points
    integer :: nri
    integer, pointer :: ri(:,:)
    integer, pointer :: rimap(:)
    integer, pointer :: rimap_inv(:)
    integer, pointer :: map_split(:, :)
    integer          :: n4
    integer          :: n1
    
    type(nl_operator_index_t) :: inner
    type(nl_operator_index_t) :: outer

#ifdef HAVE_OPENCL
    type(octcl_kernel_t) :: kernel
    type(opencl_mem_t) :: buff_imin
    type(opencl_mem_t) :: buff_imax
    type(opencl_mem_t) :: buff_ri
    type(opencl_mem_t) :: buff_map
    type(opencl_mem_t) :: buff_map_split
#endif
  end type nl_operator_t

  integer, parameter :: &
       OP_FORTRAN = 0,  &
       OP_VEC     = 1,  &
       OP_MIN     = OP_FORTRAN, &
       OP_MAX     = OP_VEC

#ifdef HAVE_OPENCL
  integer, parameter ::  &
    OP_INVMAP    = 1,       &
    OP_MAP       = 2,       &
    OP_MAP_SPLIT = 3
#endif

  integer, public, parameter :: OP_ALL = 3, OP_INNER = 1, OP_OUTER = 2

  logical :: initialized = .false.

  interface
    integer function op_is_available(opid, type)
      integer, intent(in) :: opid, type
    end function op_is_available
  end interface

  integer :: dfunction_global = -1
  integer :: zfunction_global = -1
#ifdef HAVE_OPENCL
  integer :: function_opencl
#endif

  type(profile_t), save :: nl_operate_profile
  type(profile_t), save :: operate_batch_prof


contains
  
  ! ---------------------------------------------------------
  subroutine nl_operator_global_init()
    integer :: default

    PUSH_SUB(nl_operator_global_init)

    !%Variable OperateDouble
    !%Type integer
    !%Default autodetect
    !%Section Execution::Optimization
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for real functions.
    !% By default the optimized version is used.
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End

    !%Variable OperateComplex
    !%Type integer
    !%Default autodetect
    !%Section Execution::Optimization
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for complex functions. 
    !% By default the optimized version is used.
    !%Option fortran 0
    !% The standard Fortran function.
    !%Option optimized 1
    !% This version is optimized using vector primitives (if available).
    !%End

#ifndef SINGLE_PRECISION
    default = OP_VEC
#else
    default = OP_FORTRAN
#endif

    call parse_integer(datasets_check('OperateDouble'),  default, dfunction_global)
    if(.not.varinfo_valid_option('OperateDouble', dfunction_global)) call input_error('OperateDouble')

    call parse_integer(datasets_check('OperateComplex'), default, zfunction_global)
    if(.not.varinfo_valid_option('OperateComplex', dfunction_global)) call input_error('OperateComplex')

#ifdef HAVE_OPENCL
    if(opencl_is_enabled()) then

      !%Variable OperateOpenCL
      !%Type integer
      !%Default split
      !%Section Execution::Optimization
      !%Description
      !% This variable selects the subroutine used to apply non-local
      !% operators over the grid when opencl is used. The default is
      !% map.
      !%Option invmap 1
      !% The standard implementation ported to OpenCL.
      !%Option map 2
      !% A different version, more suitable for GPUs.
      !%Option split 3
      !% (Experimental) This operator uses two different paths, one for points where
      !% the operator can be applied in blocks and other for single
      !% points.
      !%End
      call parse_integer(datasets_check('OperateOpenCL'),  OP_MAP, function_opencl)

      if(function_opencl == OP_MAP_SPLIT) then
        call messages_experimental('split non-local operator')
      end if
    end if
#endif

    POP_SUB(nl_operator_global_init)
  end subroutine nl_operator_global_init

  ! ---------------------------------------------------------

  subroutine nl_operator_global_end()
    PUSH_SUB(nl_operator_global_end)

    POP_SUB(nl_operator_global_end)
  end subroutine nl_operator_global_end

  ! ---------------------------------------------------------

  character(len=8) function op_function_name(id) result(str)
    integer, intent(in) :: id

    PUSH_SUB(op_function_name)
    
    str = 'unknown'
    if(id == OP_FORTRAN) str = 'Fortran'
    if(id == OP_VEC)     str = 'Vector'
    
    POP_SUB(op_function_name)
  end function op_function_name


  ! ---------------------------------------------------------
  subroutine nl_operator_init(op, label)
    type(nl_operator_t), intent(out) :: op
    character(len=*),    intent(in)  :: label

    PUSH_SUB(nl_operator_init)

    nullify(op%mesh, op%i, op%w_re, op%w_im, op%ri, op%rimap, op%rimap_inv)
    nullify(op%inner%imin, op%inner%imax, op%inner%ri)
    nullify(op%outer%imin, op%outer%imax, op%outer%ri)
    nullify(op%nn)

    op%label = label

    POP_SUB(nl_operator_init)
  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_copy(opo, opi)
    type(nl_operator_t), intent(out) :: opo
    type(nl_operator_t), intent(in)  :: opi

    PUSH_SUB(nl_operator_copy)

    call nl_operator_init(opo, opi%label)

    call stencil_copy(opi%stencil, opo%stencil)

    opo%np           =  opi%np
    opo%mesh         => opi%mesh

    call loct_pointer_copy(opo%nn, opi%nn)
    call loct_pointer_copy(opo%i,    opi%i)
    call loct_pointer_copy(opo%w_re, opi%w_re)
    call loct_pointer_copy(opo%w_im, opi%w_im)

    opo%const_w   = opi%const_w
    opo%cmplx_op  = opi%cmplx_op

    opo%nri       =  opi%nri
    ASSERT(associated(opi%ri))

    call loct_pointer_copy(opo%ri, opi%ri)
    call loct_pointer_copy(opo%rimap, opi%rimap)
    call loct_pointer_copy(opo%rimap_inv, opi%rimap_inv)
    
    if(opi%mesh%parallel_in_domains) then
      opo%inner%nri = opi%inner%nri
      call loct_pointer_copy(opo%inner%imin, opi%inner%imin)
      call loct_pointer_copy(opo%inner%imax, opi%inner%imax)
      call loct_pointer_copy(opo%inner%ri,   opi%inner%ri)      

      opo%outer%nri = opi%outer%nri
      call loct_pointer_copy(opo%outer%imin, opi%outer%imin)
      call loct_pointer_copy(opo%outer%imax, opi%outer%imax)
      call loct_pointer_copy(opo%outer%ri,   opi%outer%ri)
    end if

    POP_SUB(nl_operator_copy)
  end subroutine nl_operator_copy


  ! ---------------------------------------------------------
  subroutine nl_operator_build(mesh, op, np, const_w, cmplx_op)
    type(mesh_t), target, intent(in)    :: mesh
    type(nl_operator_t),  intent(inout) :: op
    integer,              intent(in)    :: np       ! Number of (local) points.
    logical, optional,    intent(in)    :: const_w  ! are the weights constant (independent of the point)
    logical, optional,    intent(in)    :: cmplx_op ! do we have complex weights?

    integer :: ii, jj, p1(MAX_DIM), time, current
    integer, allocatable :: st1(:), st2(:), st1r(:)
    integer :: nn, bl4, bl1
#ifdef HAVE_MPI
    integer :: ir, maxp, iinner, iouter
#endif
    logical :: change, force_change, iter_done
    character(len=20) :: flags

    PUSH_SUB(nl_operator_build)

#ifdef HAVE_MPI
    if(mesh%parallel_in_domains .and. .not. const_w) then
      call messages_experimental('Domain parallelization with curvilinear coordinates')
    end if
#endif

    ASSERT(np > 0)

    ! store values in structure
    op%np       = np
    op%mesh     => mesh
    op%const_w  = .false.
    op%cmplx_op = .false.
    if(present(const_w )) op%const_w  = const_w
    if(present(cmplx_op)) op%cmplx_op = cmplx_op

    ! allocate weights op%w
    if(op%const_w) then
      SAFE_ALLOCATE(op%w_re(1:op%stencil%size, 1:1))
      if (op%cmplx_op) then
        SAFE_ALLOCATE(op%w_im(1:op%stencil%size, 1:1))
      end if
      if(in_debug_mode) then
        message(1) = 'Info: nl_operator_build: working with constant weights.'
        call messages_info(1)
      end if
    else
      SAFE_ALLOCATE(op%w_re(1:op%stencil%size, 1:op%np))
      if (op%cmplx_op) then
        SAFE_ALLOCATE(op%w_im(1:op%stencil%size, 1:op%np))
      end if
      if(in_debug_mode) then
        message(1) = 'Info: nl_operator_build: working with non-constant weights.'
        call messages_info(1)
      end if
    end if

    ! set initially to zero
    op%w_re = M_ZERO
    if (op%cmplx_op) op%w_im = M_ZERO

    ! Build lookup table
    SAFE_ALLOCATE(st1(1:op%stencil%size))
    SAFE_ALLOCATE(st1r(1:op%stencil%size))
    SAFE_ALLOCATE(st2(1:op%stencil%size))

    op%nri = 0
    do time = 1, 2
      st2 = 0
      do ii = 1, np
        p1 = 0
        if(mesh%parallel_in_domains) then
          ! When running in parallel, get global number of
          ! point ii.
          call index_to_coords(mesh%idx, mesh%sb%dim, &
            mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ii - 1), p1)
        else
          call index_to_coords(mesh%idx, mesh%sb%dim, ii, p1)
        end if

        do jj = 1, op%stencil%size
          ! Get global index of p1 plus current stencil point.
          if(mesh%sb%mr_flag) then
            st1(jj) = index_from_coords(mesh%idx, mesh%sb%dim, &
                 p1(1:MAX_DIM) + mesh%resolution(p1(1), p1(2), p1(3))*op%stencil%points(1:MAX_DIM, jj))
          else
            st1(jj) = index_from_coords(mesh%idx, mesh%sb%dim, p1(1:MAX_DIM) + op%stencil%points(1:MAX_DIM, jj))
          end if
#ifdef HAVE_MPI
          if(mesh%parallel_in_domains) then
            ! When running parallel, translate this global
            ! number back to a local number.
            st1(jj) = vec_global2local(mesh%vp, st1(jj), mesh%vp%partno)
          end if
#endif
          ! if boundary conditions are zero, we can remap boundry
          ! points to reduce memory accesses. We cannot do this for the
          ! first point, since it is used to build the weights, so it
          ! has to have the positions right
          if(ii > 1 .and. mesh_compact_boundaries(mesh)) then
            st1(jj) = min(st1(jj), mesh%np + 1)
          end if
          ASSERT(st1(jj) > 0)
        end do

        st1(1:op%stencil%size) = st1(1:op%stencil%size) - ii

        change = any(st1 /= st2) 

        !the next is to detect when we move from a point that does not
        !have boundary points as neighbours to one that has
        force_change = any(st1 + ii > mesh%np) .and. all(st2 + ii - 1 <= mesh%np) 

        if(change .and. mesh_compact_boundaries(mesh)) then 
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

    op%n1 = 0
    op%n4 = 0
    do jj = 1, op%nri
      nn = op%rimap_inv(jj + 1) - op%rimap_inv(jj)
      op%n4 = op%n4 + nn/4
      op%n1 = op%n1 + mod(nn, 4)
    end do

    SAFE_ALLOCATE(op%map_split(1:2, 1:op%n1 + op%n4))

    bl1 = 1
    bl4 = 1
    ii = 1
    do
      iter_done = .false.

      if(ii + 3 <= op%np) then
        if(op%rimap(ii) == op%rimap(ii + 3)) then
          op%map_split(1, op%n1 + bl4) = ii - 1
          op%map_split(2, op%n1 + bl4) = (op%rimap(ii) - 1)*op%stencil%size

          bl4 = bl4 + 1
          ii = ii + 4
          iter_done = .true.
        end if
      end if

      if(.not. iter_done) then
        op%map_split(1, bl1) = ii - 1
        op%map_split(2, bl1) = (op%rimap(ii) - 1)*op%stencil%size

        bl1 = bl1 + 1
        ii = ii + 1
      end if

      if(ii > op%np) exit
    end do

    ASSERT(op%n4 == bl4 - 1)
    ASSERT(op%n1 == bl1 - 1)

#ifdef HAVE_MPI
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
#endif

#ifdef HAVE_OPENCL
    if(opencl_is_enabled() .and. op%const_w) then

      write(flags, '(i5)') op%stencil%size
      flags='-DSTENCIL_SIZE='//trim(adjustl(flags))

      select case(function_opencl)
      case(OP_MAP)
        call octcl_kernel_build(op%kernel, 'operate.cl', 'operate_map', flags)
      case(OP_MAP_SPLIT)
        call octcl_kernel_build(op%kernel, 'operate.cl', 'operate4', flags)
      case(OP_INVMAP)
        call octcl_kernel_build(op%kernel, 'operate.cl', 'operate', flags)
      end select

      call opencl_create_buffer(op%buff_ri, CL_MEM_READ_ONLY, TYPE_INTEGER, op%nri*op%stencil%size)
      call opencl_write_buffer(op%buff_ri, op%nri*op%stencil%size, op%ri)

      select case(function_opencl)
      case(OP_INVMAP)
        call opencl_create_buffer(op%buff_imin, CL_MEM_READ_ONLY, TYPE_INTEGER, op%nri)
        call opencl_write_buffer(op%buff_imin, op%nri, op%rimap_inv(1:))
        call opencl_create_buffer(op%buff_imax, CL_MEM_READ_ONLY, TYPE_INTEGER, op%nri)
        call opencl_write_buffer(op%buff_imax, op%nri, op%rimap_inv(2:))

      case(OP_MAP_SPLIT)
        nn = pad(op%n4 + op%n1, opencl_max_workgroup_size())*2
        call opencl_create_buffer(op%buff_map_split, CL_MEM_READ_ONLY, TYPE_INTEGER, nn)
        call opencl_write_buffer(op%buff_map_split, (op%n1 + op%n4)*2, op%map_split)

      case(OP_MAP)
        call opencl_create_buffer(op%buff_map, CL_MEM_READ_ONLY, TYPE_INTEGER, pad(op%mesh%np, opencl_max_workgroup_size()))
        call opencl_write_buffer(op%buff_map, op%mesh%np, (op%rimap - 1)*op%stencil%size)
      end select
    end if
#endif

    SAFE_DEALLOCATE_P(op%map_split)

    POP_SUB(nl_operator_build)

  end subroutine nl_operator_build

  ! ---------------------------------------------------------
  subroutine nl_operator_update_weights(this)
    type(nl_operator_t), intent(inout)  :: this

    integer :: istencil, idir

    PUSH_SUB(nl_operator_update_weights)

    if(in_debug_mode) then

      write(message(1), '(3a)') 'Debug info: Finite difference weights for ', trim(this%label), '.'
      write(message(2), '(a)')  '            Spacing:'
      do idir = 1, this%mesh%sb%dim
        write(message(2), '(a,f16.8)') trim(message(2)), this%mesh%spacing(idir)
      end do
      call messages_info(2)
      
      do istencil = 1, this%stencil%size
        write(message(1), '(a,i3,3i4,f16.10)') '      ', istencil, this%stencil%points(1:3, istencil), this%w_re(istencil, 1)
        call messages_info(1)
      end do
      
    end if

    POP_SUB(nl_operator_update_weights)

  end subroutine nl_operator_update_weights

  ! ---------------------------------------------------------
  subroutine nl_operator_transpose(op, opt)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opt

    integer :: ip, jp, kp, lp, index

    PUSH_SUB(nl_operator_transpose)

    call nl_operator_copy(opt, op)

    opt%label = trim(op%label)//' transposed'
    opt%w_re = M_ZERO
    if (op%cmplx_op) opt%w_im = M_ZERO
    do ip = 1, op%np
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(op, jp, ip)
        if(index <= op%np) then
          do lp = 1, op%stencil%size
            kp = nl_operator_get_index(op, lp, index)
            if( kp == ip ) then
              if(.not.op%const_w) then
                opt%w_re(jp, ip) = op%w_re(lp, index)
                if (op%cmplx_op) opt%w_im(jp, ip) = op%w_im(lp, index)
              else
                opt%w_re(jp, 1) = op%w_re(lp, 1)
                if (op%cmplx_op) opt%w_im(jp, 1) = op%w_im(lp, 1)
              end if
            end if
          end do
        end if
      end do
    end do

    POP_SUB(nl_operator_transpose)
  end subroutine nl_operator_transpose

  ! ---------------------------------------------------------
  ! opt has to be initialised and built.
  subroutine nl_operator_skewadjoint(op, opt, mesh)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        intent(in)  :: mesh

    integer          :: ip, jp, kp, lp, index
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    PUSH_SUB(nl_operator_skewadjoint)

    call nl_operator_copy(opt, op)

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      call nl_operator_copy(opgt, opg)
      SAFE_ALLOCATE(vol_pp(1:mesh%np_global))
      call dvec_allgather(mesh%vp, vol_pp, mesh%vol_pp)
    else
#endif
      opg  => op
      opgt => opt
      vol_pp => mesh%vol_pp
#if defined(HAVE_MPI)
    end if
#endif

    opgt%w_re = M_ZERO
    if (op%cmplx_op) opgt%w_im = M_ZERO
    do ip = 1, mesh%np_global
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(opg, jp, ip)
        if(index <= mesh%np_global) then
          do lp = 1, op%stencil%size
            kp = nl_operator_get_index(opg, lp, index)
            if( kp == ip ) then
              if(.not.op%const_w) then
                opgt%w_re(jp, ip) = M_HALF*opg%w_re(jp, ip) - M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w_re(lp, index)
                if (op%cmplx_op) &
                   opgt%w_im(jp, ip) = M_HALF*opg%w_im(jp, ip) - M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w_im(lp, index)
              else
                opgt%w_re(jp, 1) = opg%w_re(lp, 1)
                if (op%cmplx_op) opgt%w_im(jp, 1) = opg%w_im(lp, 1)
              end if
            end if
          end do
        end if
      end do
    end do

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do ip = 1, mesh%vp%np_local(mesh%vp%partno)
        opt%w_re(:, ip) = opgt%w_re(:, mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+ip-1))
        if(opt%cmplx_op) then
          opt%w_im(:, ip) = opgt%w_im(:, mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+ip-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if
#endif

    POP_SUB(nl_operator_skewadjoint)
  end subroutine nl_operator_skewadjoint


  ! ---------------------------------------------------------
  subroutine nl_operator_selfadjoint(op, opt, mesh)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        intent(in)  :: mesh

    integer          :: ip, jp, kp, lp, index
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    PUSH_SUB(nl_operator_selfadjoint)

    call nl_operator_copy(opt, op)

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      opgt = opg
      SAFE_ALLOCATE(vol_pp(1:mesh%np_global))
      call dvec_allgather(mesh%vp, vol_pp, mesh%vol_pp)
#else
      ASSERT(.false.)
#endif
    else      
      opg  => op
      opgt => opt
      vol_pp => mesh%vol_pp
    end if

    opgt%w_re = M_ZERO
    if (op%cmplx_op) opgt%w_im = M_ZERO
    do ip = 1, mesh%np_global
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(opg, jp, ip)

        if(index <= mesh%np_global) then
          do lp = 1, op%stencil%size
            kp = nl_operator_get_index(opg, lp, index)
            if( kp == ip ) then
              if(.not.op%const_w) then
                opgt%w_re(jp, ip) = M_HALF*opg%w_re(jp, ip) + M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w_re(lp, index)
                if (op%cmplx_op) &
                  opgt%w_im(jp, ip) = M_HALF*opg%w_im(jp, ip) + M_HALF*(vol_pp(index)/vol_pp(ip))*opg%w_im(lp, index)
              else
                opgt%w_re(jp, 1) = opg%w_re(lp, 1)
                if (op%cmplx_op) opgt%w_im(jp, 1) = opg%w_im(lp, 1)
              end if
            end if
          end do
        end if

      end do
    end do

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do ip = 1, mesh%vp%np_local(mesh%vp%partno)
        opt%w_re(:, ip) = opgt%w_re(:, mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+ip-1))
        if(opt%cmplx_op) then
          opt%w_im(:, ip) = opgt%w_im(:, mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+ip-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if
#endif

    POP_SUB(nl_operator_selfadjoint)
  end subroutine nl_operator_selfadjoint


#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  ! Collects a distributed non-local operator op into opg
  ! on the root node. nl_operator_end has to be called
  ! on opg when no longer needed.
  subroutine nl_operator_gather(op, opg)
    type(nl_operator_t), intent(in)  :: op  ! Local operator.
    type(nl_operator_t), intent(out) :: opg ! Global operator.

    integer :: ip

    PUSH_SUB(nl_operator_gather)

    ASSERT(associated(op%i))

    ! If root node, copy elements of op to opg that
    ! are independent from the partitions, i.e. everything
    ! except op%i and -- in the non-constant case -- op%w_re
    ! op%w_im.
    if(op%mesh%vp%rank.eq.op%mesh%vp%root) then
      call nl_operator_common_copy(op, opg)
    end if
    if(in_debug_mode) call messages_debug_newlines(4)

    ! Gather op%i and - if necessary - op%w_re and op%w_im.
    ! Collect for every point in the stencil in a single step.
    ! This permits to use ivec_gather.
    do ip = 1, op%stencil%size
      call ivec_gather(op%mesh%vp, op%mesh%vp%root, opg%i(ip, :), op%i(ip, :))
    end do
    if(op%mesh%vp%rank.eq.op%mesh%vp%root) then
      call nl_operator_translate_indices(opg)
    end if
    if(in_debug_mode) call messages_debug_newlines(2)

    ! Weights have to be collected only if they are not constant.
    if(.not.op%const_w) then
      do ip = 1, op%stencil%size
        call dvec_gather(op%mesh%vp, op%mesh%vp%root, opg%w_re(ip, :), op%w_re(ip, :))
        if(op%cmplx_op) call dvec_gather(op%mesh%vp, op%mesh%vp%root, opg%w_im(ip, :), op%w_im(ip, :))
      end do
    end if

    POP_SUB(nl_operator_gather)

  end subroutine nl_operator_gather


  ! ---------------------------------------------------------
  !> Reverse of nl_operator_gather. op is allocated, so
  !! it is necessary to call nl_operator_end on it some time.
  subroutine nl_operator_scatter(op, opg)
    type(nl_operator_t), intent(out) :: op
    type(nl_operator_t), intent(in)  :: opg

    integer :: ip

    PUSH_SUB(nl_operator_scatter)

    call nl_operator_init(op, op%label)
    call nl_operator_build(opg%mesh, op, opg%mesh%np, opg%const_w, opg%cmplx_op)

    do ip = 1, opg%stencil%size
      call dvec_scatter(opg%mesh%vp, op%mesh%vp%root, opg%w_re(ip, :), op%w_re(ip, :))
      if(opg%cmplx_op) then
        call dvec_scatter(opg%mesh%vp, op%mesh%vp%root, opg%w_im(ip, :), op%w_im(ip, :))
      end if
    end do

    POP_SUB(nl_operator_scatter)

  end subroutine nl_operator_scatter


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
    ! except op%i and -- in the non-constant case -- op%w_re
    ! op%w_im.
    call nl_operator_common_copy(op, opg)

    ! Gather op%i and -- if necessary -- op%w_re and op%w_im.
    ! Collect for every point in the stencil in a single step.
    ! This permits to use ivec_gather.
    do ip = 1, op%stencil%size
      call ivec_allgather(op%mesh%vp, opg%i(ip, :), op%i(ip, :))
    end do
    call nl_operator_translate_indices(opg)

    ! Weights have to be collected only if they are not constant.
    if(.not.op%const_w) then
      do ip = 1, op%stencil%size
        call dvec_allgather(op%mesh%vp, opg%w_re(ip, :), op%w_re(ip, :))
        if(op%cmplx_op) call dvec_allgather(op%mesh%vp, opg%w_im(ip, :), op%w_im(ip, :))
      end do
    end if

    POP_SUB(nl_operator_allgather)

  end subroutine nl_operator_allgather

  ! ---------------------------------------------------------
  ! The following are private routines.
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  !> Copies all parts of op to opg that are independent of
  !! the partitions, i.e. everything except op%i and -- in the
  !! non-constant case -- op%w_re op%w_im.
  !! This can be considered as nl_operator_copy and
  !! reallocating w_re, w_im and i.
  !! \warning: this should be replaced by a normal copy with a flag.
  subroutine nl_operator_common_copy(op, opg)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opg

    PUSH_SUB(nl_operator_common_copy)

    call nl_operator_init(opg, op%label)

    call stencil_copy(op%stencil, opg%stencil)

    SAFE_ALLOCATE(opg%i(1:op%stencil%size, 1:op%mesh%np_global))
    if(op%const_w) then
      SAFE_ALLOCATE(opg%w_re(1:op%stencil%size, 1:1))
      if(op%cmplx_op) then
        SAFE_ALLOCATE(opg%w_im(1:op%stencil%size, 1:1))
      end if
    else
      SAFE_ALLOCATE(opg%w_re(1:op%stencil%size, 1:op%mesh%np_global))
      if(op%cmplx_op) then
        SAFE_ALLOCATE(opg%w_im(1:op%stencil%size, 1:op%mesh%np_global))
      end if
    end if
    opg%mesh        => op%mesh
    opg%np       =  op%mesh%np_global
    opg%cmplx_op =  op%cmplx_op
    opg%const_w  =  op%const_w
    opg%nri      =  op%nri
    if(op%const_w) then
      opg%w_re = op%w_re
      if(op%cmplx_op) then
        opg%w_im = op%w_im
      end if
    end if

    POP_SUB(nl_operator_common_copy)

  end subroutine nl_operator_common_copy


  ! ---------------------------------------------------------
  !> Translates indices in i from local point numbers to
  !! global point numbers after gathering them.
  subroutine nl_operator_translate_indices(opg)
    type(nl_operator_t), intent(inout) :: opg

    integer :: ip, jp
    integer :: il, ig

    PUSH_SUB(nl_operator_translate_indices)

    ASSERT(associated(opg%i))

    do ip = 1, opg%stencil%size
      do jp = 1, opg%mesh%np_global
        il = opg%mesh%vp%np_local(opg%mesh%vp%part(jp))
        ig = il+opg%mesh%vp%np_ghost(opg%mesh%vp%part(jp))
        ! opg%i(ip, jp) is a local point number, i.e. it can be
        ! a real local point (i.e. the local point number
        ! is less or equal than the number of local points of
        ! the node which owns the point with global number jp):
        if(opg%i(ip, jp).le.il) then
          ! Write the global point number from the lookup
          ! table in op_(ip, jp).
          opg%i(ip, jp) = opg%mesh%vp%local(opg%mesh%vp%xlocal(opg%mesh%vp%part(jp)) &
            +opg%i(ip, jp)-1)
          ! Or a ghost point:
        else if(opg%i(ip, jp).gt.il.and.opg%i(ip, jp).le.ig) then
          opg%i(ip, jp) = opg%mesh%vp%ghost(opg%mesh%vp%xghost(opg%mesh%vp%part(jp)) &
            +opg%i(ip, jp)-1-il)
          ! Or a boundary point:
        else if(opg%i(ip, jp).gt.ig) then
          opg%i(ip, jp) = opg%mesh%vp%bndry(opg%mesh%vp%xbndry(opg%mesh%vp%part(jp)) &
            +opg%i(ip, jp)-1-ig)
        end if
      end do
    end do

    POP_SUB(nl_operator_translate_indices)

  end subroutine nl_operator_translate_indices

  ! ---------------------------------------------------------
  ! End of private routines.
  ! ---------------------------------------------------------
#endif


  ! ---------------------------------------------------------
  !> When running in parallel only the root node
  !! creates the matrix. But all nodes have to
  !! call this routine because the distributed operator has
  !! to be collected.
  subroutine nl_operator_op_to_matrix_cmplx(op, aa)
    type(nl_operator_t), target, intent(in) :: op
    CMPLX, intent(out)                      :: aa(:, :)

    integer          :: ip, jp, kp, index

    type(nl_operator_t), pointer :: opg

    PUSH_SUB(nl_operator_op_to_matrix_cmplx)

#if defined(HAVE_MPI)
    if(op%mesh%parallel_in_domains) then
      SAFE_ALLOCATE(opg)
      call nl_operator_gather(op, opg)
    else
#endif
      opg => op
#if defined(HAVE_MPI)
    end if
#endif

    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      kp = 1
      do ip = 1, op%mesh%np_global
        if(.not.op%const_w) kp = ip
        do jp = 1, op%stencil%size
          index = nl_operator_get_index(opg, jp, ip)
          if(index <= op%mesh%np_global) then
            aa(ip, index) = opg%w_re(jp, kp)
            if (op%cmplx_op) aa(ip, index) = aa(ip, index) + opg%w_im(jp, kp)
          end if
        end do
      end do
      
      if(op%mesh%parallel_in_domains) then 
        call nl_operator_end(opg)
        SAFE_DEALLOCATE_P(opg)
      end if
    end if
    if(in_debug_mode) call messages_debug_newlines(2)

    POP_SUB(nl_operator_op_to_matrix_cmplx)

  end subroutine nl_operator_op_to_matrix_cmplx

  ! ---------------------------------------------------------
  !> When running in parallel only the root node
  !! creates the matrix. But all nodes have to
  !! call this routine because the distributed operator has
  !! to be collected.
  subroutine nl_operator_op_to_matrix(op, aa, bb)
    type(nl_operator_t), target, intent(in) :: op
    FLOAT, intent(out)                      :: aa(:, :)
    FLOAT, optional, intent(out)            :: bb(:, :)

    integer          :: ip, jp, kp, index

    type(nl_operator_t), pointer :: opg

    PUSH_SUB(nl_operator_op_to_matrix)

    if(op%mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(opg)
      call nl_operator_gather(op, opg)
#else
      ASSERT(.false.)
#endif
    else
      opg => op
    end if

    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      kp = 1
      do ip = 1, op%mesh%np_global
        if(.not.op%const_w) kp = ip
        do jp = 1, op%stencil%size
          index = nl_operator_get_index(opg, jp, ip)
          if(index <= op%mesh%np_global) then
            aa(ip, index) = opg%w_re(jp, kp)
            if (op%cmplx_op) bb(ip, index) = opg%w_im(jp, kp)
          end if
        end do
      end do
      
      if(op%mesh%parallel_in_domains) then 
        call nl_operator_end(opg)
        SAFE_DEALLOCATE_P(opg)
      end if
    end if
    if(in_debug_mode) call messages_debug_newlines(2)

    POP_SUB(nl_operator_op_to_matrix)

  end subroutine nl_operator_op_to_matrix


  ! ---------------------------------------------------------
  subroutine nl_operator_matrix_to_op(op_ref, op, aa, bb)
    FLOAT, intent(in)                :: aa(:, :)
    FLOAT, optional, intent(in)      :: bb(:, :)
    type(nl_operator_t), intent(in)  :: op_ref
    type(nl_operator_t), intent(out) :: op

    integer :: ip, jp, index

    PUSH_SUB(nl_operator_matrix_to_op)

    ASSERT(associated(op_ref%i))

    call nl_operator_copy(op, op_ref)
    do ip = 1, op%np
      do jp = 1, op%stencil%size
        index = nl_operator_get_index(op, jp, ip)
        if(index <= op%np) &
          op%w_re(jp, ip) = aa(ip, index)
        if (op%cmplx_op) op%w_im(jp, ip) = bb(ip, index)
      end do
    end do

    POP_SUB(nl_operator_matrix_to_op)
  end subroutine nl_operator_matrix_to_op


  ! ---------------------------------------------------------
  !> Prints the matrix of the non-local operator op to
  !! a file (unit).
  !!
  !! When running in parallel, only the root node writes
  !! to the file but all nodes have to call this routine
  !! simultaneously because the distributed operator has
  !! to be collected.
  subroutine nl_operator_write(op, filename)
    character(len=*),    intent(in) :: filename
    type(nl_operator_t), intent(in) :: op

    integer            :: ip, jp
    integer            :: unit
    FLOAT, allocatable :: aa(:, :)

    PUSH_SUB(nl_operator_write)


    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      SAFE_ALLOCATE(aa(1:op%mesh%np_global, 1:op%mesh%np_global))
      aa = M_ZERO
    end if

    call nl_operator_op_to_matrix(op, aa)

    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      unit = io_open(filename, action='write')
      if(unit < 0) then
        message(1) = 'Could not open file '//filename//' to write operator.'
        call messages_fatal(1)
      end if
      do ip = 1, op%mesh%np_global
        do jp = 1, op%mesh%np_global - 1
          write(unit, fmt = '(f9.4)', advance ='no') aa(ip, jp)
        end do
        write(unit, fmt = '(f9.4)') aa(ip, op%mesh%np_global)
      end do
      call io_close(unit)

      SAFE_DEALLOCATE_A(aa)
    end if

    POP_SUB(nl_operator_write)

  end subroutine nl_operator_write


  ! ---------------------------------------------------------
  !> Same as nl_operator_write but transposed matrix.
  subroutine nl_operatorT_write(op, unit)
    type(nl_operator_t), intent(in) :: op
    integer, intent(in)             :: unit

    integer :: ip, jp
    FLOAT, allocatable :: aa(:, :)

    PUSH_SUB(nl_operatorT_write)


    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      SAFE_ALLOCATE(aa(1:op%mesh%np_global, 1:op%mesh%np_global))
      aa = M_ZERO
    end if

    call nl_operator_op_to_matrix(op, aa)

    if(mpi_grp_is_root(op%mesh%mpi_grp)) then
      do ip = 1, op%mesh%np_global
        do jp = 1, op%mesh%np_global - 1
          write(unit, fmt = '(f9.4)', advance ='no') aa(jp, ip)
        end do
        write(unit, fmt = '(f9.4)') aa(op%mesh%np_global, ip)
      end do

      SAFE_DEALLOCATE_A(aa)
    end if

    POP_SUB(nl_operatorT_write)

  end subroutine nl_operatorT_write


  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_t), intent(inout) :: op

    PUSH_SUB(nl_operator_end)

#ifdef HAVE_OPENCL
    if(opencl_is_enabled() .and. op%const_w) then

      call opencl_release_buffer(op%buff_ri)
      select case(function_opencl)
      case(OP_INVMAP)
        call opencl_release_buffer(op%buff_imin)
        call opencl_release_buffer(op%buff_imax)

      case(OP_MAP_SPLIT)
        call opencl_release_buffer(op%buff_map_split)

      case(OP_MAP)
        call opencl_release_buffer(op%buff_map)
      end select
    end if
#endif

    if(op%mesh%parallel_in_domains) then
      SAFE_DEALLOCATE_P(op%inner%imin)
      SAFE_DEALLOCATE_P(op%inner%imax)
      SAFE_DEALLOCATE_P(op%inner%ri)
      SAFE_DEALLOCATE_P(op%outer%imin)
      SAFE_DEALLOCATE_P(op%outer%imax)
      SAFE_DEALLOCATE_P(op%outer%ri)
    end if

    SAFE_DEALLOCATE_P(op%i)
    SAFE_DEALLOCATE_P(op%w_re)
    SAFE_DEALLOCATE_P(op%w_im)

    SAFE_DEALLOCATE_P(op%ri)
    SAFE_DEALLOCATE_P(op%rimap)
    SAFE_DEALLOCATE_P(op%rimap_inv)
    SAFE_DEALLOCATE_P(op%nn)

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

#include "undef.F90"
#include "real.F90"
#include "nl_operator_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "nl_operator_inc.F90"

end module nl_operator_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
