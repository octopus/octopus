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
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use loct_parser_m
  use math_m
  use index_m
  use mesh_m
  use messages_m
  use multicomm_m
  use mpi_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use stencil_m

  implicit none

  private
  public ::                     &
    nl_operator_t,              &
    nl_operator_index_t,        &
    nl_operator_global_init,    &
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
    type(mesh_t), pointer :: m         ! pointer to the underlying mesh
    integer               :: np        ! number of points in mesh
    integer               :: der_zero_max
    ! When running in parallel mode, the next three
    ! arrays are unique on each node.
    integer, pointer  :: i(:,:)    ! index of the points
    FLOAT,   pointer  :: w_re(:,:) ! weightsp, real part
    FLOAT,   pointer  :: w_im(:,:) ! weightsp, imaginary part

    logical               :: const_w   ! are the weights independent of index i
    logical               :: cmplx_op  ! .true. if we have also imaginary weights

    integer :: dfunction
    integer :: zfunction

    character(len=40) :: label

    !the compressed index of grid points
    integer :: nri
    integer, pointer :: ri(:,:)
    integer, pointer :: rimap(:)
    integer, pointer :: rimap_inv(:)
    integer(4), pointer :: ribit(:)
    
    type(nl_operator_index_t) :: inner
    type(nl_operator_index_t) :: outer
    
  end type nl_operator_t

  integer, parameter :: &
       OP_FORTRAN = 0,  &
       OP_C       = 1,  &
       OP_VEC     = 2,  &
       OP_AS      = 3,  &
       OP_BIT     = 4,  &
       OP_BG      = 5,  &
       OP_MIN     = OP_FORTRAN, &
       OP_MAX     = OP_BG
  
  integer, public, parameter :: OP_ALL = 3, OP_INNER = 1, OP_OUTER = 2

  logical :: initialized = .false.

  interface
    integer function op_is_available(opid, type)
      integer, intent(in) :: opid, type
    end function op_is_available
  end interface

  integer :: dfunction_global = -1
  integer :: zfunction_global = -1

  type(profile_t), save :: nl_operate_profile
  type(profile_t), save :: operate_batch_prof

contains
  
  ! ---------------------------------------------------------
  subroutine nl_operator_global_init()

    !%Variable OperateDouble
    !%Type integer
    !%Default autodetect
    !%Section Execution::Optimization
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for real functions. By default the best
    !% subroutine for your system is detected at runtime, but using
    !% this variable you might skip the detection.
    !%Option autodetect -1
    !% Automatically discover the fastest function
    !%Option fortran 0
    !% The standard fortran function.
    !%Option c 1
    !% The C version of the function unrolled by hand.
    !%Option vec 2
    !% This version has been optimized using vector primitives, it is
    !% available on x86 (with SSE2) and x86-64 systems (not supported
    !% by all compilers).
    !%Option as 3
    !% Hand-written assembler version, currently only available on Itanium systems.
    !%End

    !%Variable OperateComplex
    !%Type integer
    !%Default autodetect
    !%Section Execution::Optimization
    !%Description
    !% This variable selects the subroutine used to apply non-local
    !% operators over the grid for complex functions. By default the best
    !% subroutine for your system is detected at runtime, but using
    !% this variable you might skip the detection.
    !%Option autodetect -1
    !% Automatically discover the fastest function
    !%Option fortran 0
    !% The standard fortran function.
    !%Option c 1
    !% The C version of the function unrolled by hand.
    !%Option vec 2
    !% This version has been optimized using vector primitives, it is
    !% available on x86 (with SSE2) and x86-64 systems (not supported
    !% by all compilers).
    !%Option as 3
    !% Hand-written assembler version, currently only available on Itanium systems.
    !%End

    call loct_parse_int(datasets_check('OperateDouble'),  -1, dfunction_global)
    if(dfunction_global.ne.-1) then
      if(op_is_available(dfunction_global, M_REAL)  == 0) then
        message(1) = 'OperateDouble chosen in not available on this platform'
        call write_fatal(1)
      end if
    end if

    call loct_parse_int(datasets_check('OperateComplex'), -1, zfunction_global)
    if(zfunction_global.ne.-1) then
      if(op_is_available(zfunction_global, M_CMPLX) == 0) then
        message(1) = 'OperateComplex chosen in not available on this platform'
        call write_fatal(1)
      end if
    end if

  end subroutine nl_operator_global_init


  ! ---------------------------------------------------------
  character(len=8) function op_function_name(id) result(str)
    integer, intent(in) :: id
    
    str = 'unknown'
    if(id == OP_FORTRAN) str = 'Fortran'
    if(id == OP_C)       str = 'C'
    if(id == OP_VEC)     str = 'Vector'
    if(id == OP_AS)      str = 'AS'
    if(id == OP_BIT)     str = 'Bit'
    if(id == OP_BG)      str = 'BG'
    
  end function op_function_name


  ! ---------------------------------------------------------
  subroutine nl_operator_init(op, label)
    type(nl_operator_t), intent(out) :: op
    character(len=*),    intent(in)  :: label

    call push_sub('nl_operator.nl_operator_init')

    nullify(op%m, op%i, op%w_re, op%w_im, op%ri, op%rimap, op%rimap_inv)
    nullify(op%inner%imin, op%inner%imax, op%inner%ri)
    nullify(op%outer%imin, op%outer%imax, op%outer%ri)
    nullify(op%ribit)

    op%label = label

    call pop_sub()
  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_copy(opo, opi)
    type(nl_operator_t), intent(out) :: opo
    type(nl_operator_t), intent(in)  :: opi

    call push_sub('nl_operator.nl_operator_copy')

    call nl_operator_init(opo, opi%label)

    opo%np       =  opi%np
    opo%m        => opi%m
    opo%cmplx_op =  opi%cmplx_op
    opo%nri      =  opi%nri

    call stencil_copy(opi%stencil, opo%stencil)

    opo%const_w = opi%const_w

    call loct_pointer_copy(opo%w_re, opi%w_re)
    if (opi%cmplx_op) call loct_pointer_copy(opo%w_im, opi%w_im)

    opo%dfunction = opi%dfunction
    opo%zfunction = opi%zfunction

    ASSERT(associated(opi%ri))

    call loct_pointer_copy(opo%ri, opi%ri)
    call loct_pointer_copy(opo%ribit, opi%ribit)
    call loct_pointer_copy(opo%rimap, opi%rimap)
    call loct_pointer_copy(opo%rimap_inv, opi%rimap_inv)
    
    opo%dfunction = opi%dfunction
    opo%zfunction = opi%zfunction

    if(opi%m%parallel_in_domains) then
      opo%inner%nri = opi%inner%nri
      call loct_pointer_copy(opo%inner%imin, opi%inner%imin)
      call loct_pointer_copy(opo%inner%imax, opi%inner%imax)
      call loct_pointer_copy(opo%inner%ri,   opi%inner%ri)      

      opo%outer%nri = opi%outer%nri
      call loct_pointer_copy(opo%outer%imin, opi%outer%imin)
      call loct_pointer_copy(opo%outer%imax, opi%outer%imax)
      call loct_pointer_copy(opo%outer%ri,   opi%outer%ri)
    end if

    call pop_sub()
  end subroutine nl_operator_copy


  ! ---------------------------------------------------------
  subroutine nl_operator_build(m, op, np, const_w, cmplx_op)
    type(mesh_t), target, intent(in)    :: m
    type(nl_operator_t),  intent(inout) :: op
    integer,              intent(in)    :: np       ! Number of (local) points.
    logical, optional,    intent(in)    :: const_w  ! are the weights constant (independent of the point)
    logical, optional,    intent(in)    :: cmplx_op ! do we have complex weights?

    integer :: ii, jj, p1(MAX_DIM), time, current
    integer, allocatable :: st1(:), st2(:), st1r(:)
    FLOAT :: dbest, zbest
#ifdef HAVE_MPI
    integer :: ir, maxp, iinner, iouter
#endif
    logical :: change

    call push_sub('nl_operator.nl_operator_build')

#ifdef HAVE_MPI
    if(m%parallel_in_domains .and. .not. const_w) then
      call messages_devel_version('Domain parallelization with curvilinear coordinates')
    end if
#endif

    ASSERT(np > 0)

    ! store values in structure
    op%np       = np
    op%m        => m
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
      message(1) = 'Info: nl_operator_build: working with constant weights.'
      if(in_debug_mode) call write_info(1)
    else
      SAFE_ALLOCATE(op%w_re(1:op%stencil%size, 1:op%np))
      if (op%cmplx_op) then
        SAFE_ALLOCATE(op%w_im(1:op%stencil%size, 1:op%np))
      end if
      message(1) = 'Info: nl_operator_build: working with non-constant weights.'
      if(in_debug_mode) call write_info(1)
    end if

    ! set initially to zero
    op%w_re = M_ZERO
    if (op%cmplx_op) op%w_im = M_ZERO

    ! Build lookup table
    SAFE_ALLOCATE(st1(1:op%stencil%size))
    SAFE_ALLOCATE(st1r(1:op%stencil%size))
    SAFE_ALLOCATE(st2(1:op%stencil%size))

    op%nri = 0
    op%der_zero_max = m%np + 1
    do time = 1, 2
      st2 = 0
      do ii = 1, np
        p1 = 0
        if(m%parallel_in_domains) then
          ! When running in parallel, get global number of
          ! point ii.
          call index_to_coords(m%idx, m%sb%dim, m%vp%local(m%vp%xlocal(m%vp%partno) + ii - 1), p1)
        else
          call index_to_coords(m%idx, m%sb%dim, ii, p1)
        end if

        do jj = 1, op%stencil%size
          ! Get global index of p1 plus current stencil point.
          if(m%sb%mr_flag) then
            st1(jj) = index_from_coords(m%idx, m%sb%dim, &
                 p1(1:MAX_DIM) + m%resolution(p1(1), p1(2), p1(3))*op%stencil%points(1:MAX_DIM, jj))
          else
            st1(jj) = index_from_coords(m%idx, m%sb%dim, p1(1:MAX_DIM) + op%stencil%points(1:MAX_DIM, jj))
          end if
#ifdef HAVE_MPI
          if(m%parallel_in_domains) then
            ! When running parallel, translate this global
            ! number back to a local number.
            st1(jj) = vec_global2local(m%vp, st1(jj), m%vp%partno)
          end if
#endif
          if(mesh_compact_boundaries(m)) then
            st1(jj) = min(st1(jj), m%np + 1)
          end if
          ASSERT(st1(jj) > 0)
        end do

        st1(1:op%stencil%size) = st1(1:op%stencil%size) - ii

        change = any(st1 /= st2)
        
        if(change .and. mesh_compact_boundaries(m)) then 
          !try to repair it by changing the boundary points
          do jj = 1, op%stencil%size
            if(st1(jj) + ii > m%np .and. st2(jj) + ii - 1 > m%np .and. st2(jj) + ii <= m%np_part) then
              st1r(jj) = st2(jj)
            else
              st1r(jj) = st1(jj)
            end if
          end do

          change = any(st1r /= st2)

          if(.not. change) then
            st1 = st1r
            op%der_zero_max = max(op%der_zero_max, maxval(st1) + ii)
          end if
        end if

        ! if the stencil changes
        if (change) then 
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
      end if

    end do

    !the inverse mapping
    op%rimap_inv(1) = 0
    do jj = 1, op%np
      op%rimap_inv(op%rimap(jj) + 1) = jj
    end do
    op%rimap_inv(op%nri + 1) = op%np

    ! we allocate ribit with the maximum size it can take (the 1 +
    ! size/32 is to account for the descriptor that takes 1 bit per value)
    SAFE_ALLOCATE(op%ribit(1:(op%stencil%size + 1 + op%stencil%size/32)*op%nri))
    op%ribit = 0

    call generate_ribit(op%nri, op%stencil%size, op%ri, op%ribit)

    SAFE_DEALLOCATE_A(st1)
    SAFE_DEALLOCATE_A(st1r)
    SAFE_DEALLOCATE_A(st2)

#ifdef HAVE_MPI
    if(op%m%parallel_in_domains) then
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
          ASSERT(all(ii + op%inner%ri(1:op%stencil%size, ir) <= m%np))
        end do
      end do
      
    end if
#endif

    if(op%const_w .and. .not. op%cmplx_op) then
      message(1) = 'Info: '//trim(op%label)
      message(2) = '      Total throughput (MFlops)'
      if(dfunction_global == -1) then
        call dnl_operator_tune(op, dbest)
        write(message(2), '(2a,i7)') trim(message(2)), ' real = ', int(dbest)
      else
        op%dfunction = dfunction_global
      end if

      if(zfunction_global == -1) then
        call znl_operator_tune(op, zbest)
        write(message(2), '(2a,i7)') trim(message(2)), ' complex = ', int(zbest)
      else
        op%zfunction = zfunction_global
      end if

      if(dfunction_global == -1.or.zfunction_global == -1) call write_info(2)
    end if

    call pop_sub()

  end subroutine nl_operator_build

  ! ---------------------------------------------------------
  subroutine nl_operator_transpose(op, opt)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opt

    integer :: i, j, index, l, k

    call push_sub('nl_operator.nl_operator_transpose')

    call nl_operator_copy(opt, op)

    opt%label = trim(op%label)//' transposed'
    opt%w_re = M_ZERO
    if (op%cmplx_op) opt%w_im = M_ZERO
    do i = 1, op%np
      do j = 1, op%stencil%size
        index = nl_operator_get_index(op, j, i)
        if(index <= op%np) then
          do l = 1, op%stencil%size
            k = nl_operator_get_index(op, l, index)
            if( k == i ) then
              if(.not.op%const_w) then
                opt%w_re(j, i) = op%w_re(l, index)
                if (op%cmplx_op) opt%w_im(j, i) = op%w_im(l, index)
              else
                opt%w_re(j, 1) = op%w_re(l, 1)
                if (op%cmplx_op) opt%w_im(j, 1) = op%w_im(l, 1)
              end if
            end if
          end do
        end if
      end do
    end do

    if(opt%const_w .and. .not. opt%cmplx_op) then
      call dnl_operator_tune(opt)
      call znl_operator_tune(opt)
    end if

    call pop_sub()
  end subroutine nl_operator_transpose

  ! ---------------------------------------------------------
  ! opt has to be initialised and built.
  subroutine nl_operator_skewadjoint(op, opt, m)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        intent(in)  :: m

    integer          :: i, j, index, l, k
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    call push_sub('nl_operator.nl_operator_skewadjoint')

    call nl_operator_copy(opt, op)

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      call nl_operator_copy(opgt, opg)
      SAFE_ALLOCATE(vol_pp(1:m%np_global))
      call dvec_allgather(m%vp, vol_pp, m%vol_pp)
    else
#endif
      opg  => op
      opgt => opt
      vol_pp => m%vol_pp
#if defined(HAVE_MPI)
    end if
#endif

    opgt%w_re = M_ZERO
    if (op%cmplx_op) opgt%w_im = M_ZERO
    do i = 1, m%np_global
      do j = 1, op%stencil%size
        index = nl_operator_get_index(opg, j, i)
        if(index <= m%np_global) then
          do l = 1, op%stencil%size
            k = nl_operator_get_index(opg, l, index)
            if( k == i ) then
              if(.not.op%const_w) then
                opgt%w_re(j, i) = M_HALF*opg%w_re(j, i) - M_HALF*(vol_pp(index)/vol_pp(i))*opg%w_re(l, index)
                if (op%cmplx_op) &
                   opgt%w_im(j, i) = M_HALF*opg%w_im(j, i) - M_HALF*(vol_pp(index)/vol_pp(i))*opg%w_im(l, index)
              else
                opgt%w_re(j, 1) = opg%w_re(l, 1)
                if (op%cmplx_op) opgt%w_im(j, 1) = opg%w_im(l, 1)
              end if
            end if
          end do
        end if
      end do
    end do

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do i = 1, m%vp%np_local(m%vp%partno)
        opt%w_re(:, i) = opgt%w_re(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        if(opt%cmplx_op) then
          opt%w_im(:, i) = opgt%w_im(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if
#endif

    call pop_sub()
  end subroutine nl_operator_skewadjoint


  ! ---------------------------------------------------------
  subroutine nl_operator_selfadjoint(op, opt, m)
    type(nl_operator_t), target, intent(in)  :: op
    type(nl_operator_t), target, intent(out) :: opt
    type(mesh_t),        intent(in)  :: m

    integer          :: i, j, index, l, k
    FLOAT, pointer   :: vol_pp(:)

    type(nl_operator_t), pointer :: opg, opgt

    call push_sub('nl_operator.nl_operator_selfadjoint')

    call nl_operator_copy(opt, op)

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(opg)
      SAFE_ALLOCATE(opgt)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%label)
      opgt = opg
      SAFE_ALLOCATE(vol_pp(1:m%np_global))
      call dvec_allgather(m%vp, vol_pp, m%vol_pp)
#else
      ASSERT(.false.)
#endif
    else      
      opg  => op
      opgt => opt
      vol_pp => m%vol_pp
    end if

    opgt%w_re = M_ZERO
    if (op%cmplx_op) opgt%w_im = M_ZERO
    do i = 1, m%np_global
      do j = 1, op%stencil%size
        index = nl_operator_get_index(opg, j, i)

        if(index <= m%np_global) then
          do l = 1, op%stencil%size
            k = nl_operator_get_index(opg, l, index)
            if( k == i ) then
              if(.not.op%const_w) then
                opgt%w_re(j, i) = M_HALF*opg%w_re(j, i) + M_HALF*(vol_pp(index)/vol_pp(i))*opg%w_re(l, index)
                if (op%cmplx_op) &
                  opgt%w_im(j, i) = M_HALF*opg%w_im(j, i) + M_HALF*(vol_pp(index)/vol_pp(i))*opg%w_im(l, index)
              else
                opgt%w_re(j, 1) = opg%w_re(l, 1)
                if (op%cmplx_op) opgt%w_im(j, 1) = opg%w_im(l, 1)
              end if
            end if
          end do
        end if

      end do
    end do

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      SAFE_DEALLOCATE_P(vol_pp)
      do i = 1, m%vp%np_local(m%vp%partno)
        opt%w_re(:, i) = opgt%w_re(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        if(opt%cmplx_op) then
          opt%w_im(:, i) = opgt%w_im(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
      SAFE_DEALLOCATE_P(opg)
      SAFE_DEALLOCATE_P(opgt)
    end if
#endif

    call pop_sub()
  end subroutine nl_operator_selfadjoint


#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  ! Collects a distributed non-local operator op into opg
  ! on the root node. nl_operator_end has to be called
  ! on opg when no longer needed.
  subroutine nl_operator_gather(op, opg)
    type(nl_operator_t), intent(in)  :: op  ! Local operator.
    type(nl_operator_t), intent(out) :: opg ! Global operator.

    integer :: i

    call push_sub('nl_operator.nl_operator_gather')

    ASSERT(associated(op%i))

    ! If root node, copy elements of op to opg that
    ! are independent from the partitions, i. e. everything
    ! except op%i and - in the non constant case - op%w_re
    ! op%w_im.
    if(op%m%vp%rank.eq.op%m%vp%root) then
      call nl_operator_common_copy(op, opg)
    end if
    if(in_debug_mode) call write_debug_newlines(4)

    ! Gather op%i and - if necessary - op%w_re and op%w_im.
    ! Collect for every point in the stencil in a single step.
    ! This permits to use ivec_gather.
    do i = 1, op%stencil%size
      call ivec_gather(op%m%vp, opg%i(i, :), op%i(i, :))
    end do
    if(op%m%vp%rank.eq.op%m%vp%root) then
      call nl_operator_translate_indices(opg)
    end if
    if(in_debug_mode) call write_debug_newlines(2)

    ! Weights have to be collected only if they are
    ! not constant.
    if(.not.op%const_w) then
      do i = 1, op%stencil%size
        call dvec_gather(op%m%vp, opg%w_re(i, :), op%w_re(i, :))
        if(op%cmplx_op) call dvec_gather(op%m%vp, opg%w_im(i, :), op%w_im(i, :))
      end do
    end if

    call pop_sub()

  end subroutine nl_operator_gather


  ! ---------------------------------------------------------
  ! Reverse of nl_operator_gather. op is allocated, so
  ! it is necessary to call nl_operator_end on it some time.
  subroutine nl_operator_scatter(op, opg)
    type(nl_operator_t), intent(out) :: op
    type(nl_operator_t), intent(in)  :: opg

    integer :: i

    call push_sub('nl_operator.nl_operator_scatter')

    call nl_operator_init(op, op%label)
    call nl_operator_build(opg%m, op, opg%m%np, opg%const_w, opg%cmplx_op)

    do i = 1, opg%stencil%size
      call dvec_scatter(opg%m%vp, opg%w_re(i, :), op%w_re(i, :))
      if(opg%cmplx_op) then
        call dvec_scatter(opg%m%vp, opg%w_im(i, :), op%w_im(i, :))
      end if
    end do

    call pop_sub()

  end subroutine nl_operator_scatter


  ! ---------------------------------------------------------
  ! Like nl_operator_gather but opg is present on all nodes
  ! (so do not forget to call nl_operator_end on all nodes
  ! afterwards).
  subroutine nl_operator_allgather(op, opg)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opg

    integer :: i

    call push_sub('nl_operator.nl_operator_allgather')

    ! Copy elements of op to opg that
    ! are independent from the partitions, i. e. everything
    ! except op%i and - in the non constant case - op%w_re
    ! op%w_im.
    call nl_operator_common_copy(op, opg)

    ! Gather op%i and - if necessary - op%w_re and op%w_im.
    ! Collect for every point in the stencil in a single step.
    ! This permits to use ivec_gather.
    do i = 1, op%stencil%size
      call ivec_allgather(op%m%vp, opg%i(i, :), op%i(i, :))
    end do
    call nl_operator_translate_indices(opg)

    ! Weights have to be collected only if they are
    ! not constant.
    if(.not.op%const_w) then
      do i = 1, op%stencil%size
        call dvec_allgather(op%m%vp, opg%w_re(i, :), op%w_re(i, :))
        if(op%cmplx_op) call dvec_allgather(op%m%vp, opg%w_im(i, :), op%w_im(i, :))
      end do
    end if

    call pop_sub()

  end subroutine nl_operator_allgather

  ! ---------------------------------------------------------
  ! The following are private routines.
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  ! Copies all part of op to opg that are independent of
  ! the partitions,i. e. everything except op%i and - in the
  ! non constant case - op%w_re op%w_im.
  ! This can be considered as nl_operator_copy and
  ! reallocating w_re, w_im and i.
  subroutine nl_operator_common_copy(op, opg)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opg

    call push_sub('nl_operator.nl_operator_common_copy')

    call nl_operator_init(opg, op%label)

    call stencil_copy(op%stencil, opg%stencil)

    SAFE_ALLOCATE(opg%i(1:op%stencil%size, 1:op%m%np_global))
    if(op%const_w) then
      SAFE_ALLOCATE(opg%w_re(1:op%stencil%size, 1:1))
      if(op%cmplx_op) then
        SAFE_ALLOCATE(opg%w_im(1:op%stencil%size, 1:1))
      end if
    else
      SAFE_ALLOCATE(opg%w_re(1:op%stencil%size, 1:op%m%np_global))
      if(op%cmplx_op) then
        SAFE_ALLOCATE(opg%w_im(1:op%stencil%size, 1:op%m%np_global))
      end if
    end if
    opg%m        => op%m
    opg%np       =  op%m%np_global
    opg%cmplx_op =  op%cmplx_op
    opg%const_w  =  op%const_w
    opg%nri      =  op%nri
    if(op%const_w) then
      opg%w_re = op%w_re
      if(op%cmplx_op) then
        opg%w_im = op%w_im
      end if
    end if

    call pop_sub()

  end subroutine nl_operator_common_copy


  ! ---------------------------------------------------------
  ! Translates indices in i from local point numbers to
  ! global point numbers after gathering them.
  subroutine nl_operator_translate_indices(opg)
    type(nl_operator_t), intent(inout) :: opg

    integer :: i, j
    integer :: il, ig

    call push_sub('nl_operator.nl_operator_translate_indices')

    ASSERT(associated(opg%i))

    do i = 1, opg%stencil%size
      do j = 1, opg%m%np_global
        il = opg%m%vp%np_local(opg%m%vp%part(j))
        ig = il+opg%m%vp%np_ghost(opg%m%vp%part(j))
        ! opg%i(i, j) is a local point number, i. e. it can be
        ! a real local point (i. e. the local point number
        ! is less or equal than the number of local points of
        ! the node which owns the point with global number j):
        if(opg%i(i, j).le.il) then
          ! Write the global point number from the lookup
          ! table in op_(i, j).
          opg%i(i, j) = opg%m%vp%local(opg%m%vp%xlocal(opg%m%vp%part(j)) &
            +opg%i(i, j)-1)
          ! Or a ghost point:
        else if(opg%i(i, j).gt.il.and.opg%i(i, j).le.ig) then
          opg%i(i, j) = opg%m%vp%ghost(opg%m%vp%xghost(opg%m%vp%part(j)) &
            +opg%i(i, j)-1-il)
          ! Or a boundary point:
        else if(opg%i(i, j).gt.ig) then
          opg%i(i, j) = opg%m%vp%bndry(opg%m%vp%xbndry(opg%m%vp%part(j)) &
            +opg%i(i, j)-1-ig)
        end if
      end do
    end do

    call pop_sub()

  end subroutine nl_operator_translate_indices

  ! ---------------------------------------------------------
  ! End of private routines.
  ! ---------------------------------------------------------
#endif


  ! ---------------------------------------------------------
  ! When running in parallel only the root node
  ! creates the matrix. But all nodes have to
  ! call this routine because the distributed operator has
  ! to be collected.
  subroutine nl_operator_op_to_matrix_cmplx(op, a)
    type(nl_operator_t), target, intent(in) :: op
    CMPLX, intent(out)                      :: a(:, :)

    integer          :: i, j, k, index

    type(nl_operator_t), pointer :: opg

    call push_sub('nl_operator.nl_operator_op_to_matrix')

    if(op%m%parallel_in_domains) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(opg)
      call nl_operator_gather(op, opg)
#else
      ASSERT(.false.)
#endif
    else
      opg => op
    end if

    if(mpi_grp_is_root(op%m%mpi_grp)) then
      k = 1
      do i = 1, op%m%np_global
        if(.not.op%const_w) k = i
        do j = 1, op%stencil%size
          index = nl_operator_get_index(opg, j, i)
          if(index <= op%m%np_global) then
            a(i, index) = opg%w_re(j, k)
            if (op%cmplx_op) a(i, index) = a(i, index) + opg%w_im(j, k)
          end if
        end do
      end do
      
      if(op%m%parallel_in_domains) then 
        call nl_operator_end(opg)
        SAFE_DEALLOCATE_P(opg)
      end if
    end if
    if(in_debug_mode) call write_debug_newlines(2)

    call pop_sub()

  end subroutine nl_operator_op_to_matrix_cmplx

  ! ---------------------------------------------------------
  ! When running in parallel only the root node
  ! creates the matrix. But all nodes have to
  ! call this routine because the distributed operator has
  ! to be collected.
  subroutine nl_operator_op_to_matrix(op, a, b)
    type(nl_operator_t), target, intent(in) :: op
    FLOAT, intent(out)                 :: a(:, :)
    FLOAT, optional, intent(out)       :: b(:, :)

    integer          :: i, j, k, index

    type(nl_operator_t), pointer :: opg

    call push_sub('nl_operator.nl_operator_op_to_matrix')

    if(op%m%parallel_in_domains) then
#if defined(HAVE_MPI)
      SAFE_ALLOCATE(opg)
      call nl_operator_gather(op, opg)
#else
      ASSERT(.false.)
#endif
    else
      opg => op
    end if

    if(mpi_grp_is_root(op%m%mpi_grp)) then
      k = 1
      do i = 1, op%m%np_global
        if(.not.op%const_w) k = i
        do j = 1, op%stencil%size
          index = nl_operator_get_index(opg, j, i)
          if(index <= op%m%np_global) then
            a(i, index) = opg%w_re(j, k)
            if (op%cmplx_op) b(i, index) = opg%w_im(j, k)
          end if
        end do
      end do
      
      if(op%m%parallel_in_domains) then 
        call nl_operator_end(opg)
        SAFE_DEALLOCATE_P(opg)
      end if
    end if
    if(in_debug_mode) call write_debug_newlines(2)

    call pop_sub()

  end subroutine nl_operator_op_to_matrix


  ! ---------------------------------------------------------
  subroutine nl_operator_matrix_to_op(op_ref, op, a, b)
    FLOAT, intent(in)                   :: a(:, :)
    FLOAT, optional, intent(in)         :: b(:, :)
    type(nl_operator_t), intent(in)  :: op_ref
    type(nl_operator_t), intent(out) :: op

    integer :: i, j, index

    call push_sub('nl_operator.nl_operator_matrix_to_op')

    ASSERT(associated(op_ref%i))

    call nl_operator_copy(op, op_ref)
    do i = 1, op%np
      do j = 1, op%stencil%size
        index = nl_operator_get_index(op, j, i)
        if(index <= op%np) &
          op%w_re(j, i) = a(i, index)
        if (op%cmplx_op) op%w_im(j, i) = b(i, index)
      end do
    end do

    call pop_sub()
  end subroutine nl_operator_matrix_to_op


  ! ---------------------------------------------------------
  ! Prints the matrix of the non local operator op to
  ! a file (unit).
  ! When running in parallel, only the root node writes
  ! to the file but all nodes have to call this routine
  ! simultaneously because the distributed operator has
  ! to be collected.
  subroutine nl_operator_write(op, filename)
    character(len=*),       intent(in) :: filename
    type(nl_operator_t), intent(in) :: op

    integer            :: i, j
    integer            :: unit
    FLOAT, allocatable :: a(:, :)

    call push_sub('nl_operator.nl_operator_write')


    if(mpi_grp_is_root(op%m%mpi_grp)) then
      SAFE_ALLOCATE(a(1:op%m%np_global, 1:op%m%np_global))
      a = M_ZERO
    end if

    call nl_operator_op_to_matrix(op, a)

    if(mpi_grp_is_root(op%m%mpi_grp)) then
      unit = io_open(filename, action='write')
      if(unit < 0) then
        message(1) = 'Could not open file '//filename//' to write operator.'
        call write_fatal(1)
      end if
      do i = 1, op%m%np_global
        do j = 1, op%m%np_global - 1
          write(unit, fmt = '(f9.4)', advance ='no') a(i, j)
        end do
        write(unit, fmt = '(f9.4)') a(i, op%m%np_global)
      end do
      call io_close(unit)

      SAFE_DEALLOCATE_A(a)
    end if

    call pop_sub()

  end subroutine nl_operator_write


  ! ---------------------------------------------------------
  ! Same as nl_operator_write but transposed matrix.
  subroutine nl_operatorT_write(op, unit)
    type(nl_operator_t), intent(in) :: op
    integer, intent(in)                :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)

    call push_sub('nl_operator.nl_operator_write')


    if(mpi_grp_is_root(op%m%mpi_grp)) then
      SAFE_ALLOCATE(a(1:op%m%np_global, 1:op%m%np_global))
      a = M_ZERO
    end if

    call nl_operator_op_to_matrix(op, a)

    if(mpi_grp_is_root(op%m%mpi_grp)) then
      do i = 1, op%m%np_global
        do j = 1, op%m%np_global - 1
          write(unit, fmt = '(f9.4)', advance ='no') a(j, i)
        end do
        write(unit, fmt = '(f9.4)') a(op%m%np_global, i)
      end do

      SAFE_DEALLOCATE_A(a)
    end if

    call pop_sub()

  end subroutine nl_operatorT_write


  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_t), intent(inout) :: op

    call push_sub('nl_operator.nl_operator_end')

    if(op%m%parallel_in_domains) then
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
    SAFE_DEALLOCATE_P(op%ribit)
    SAFE_DEALLOCATE_P(op%rimap)
    SAFE_DEALLOCATE_P(op%rimap_inv)

    call stencil_end(op%stencil)

    call pop_sub()
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
