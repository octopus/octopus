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

module nl_operator_m
  use global_m
  use messages_m
  use mesh_m
  use mesh_lib_m
  use simul_box_m
  use io_m
  use profiling_m
  use par_vec_m
  use mpi_m

  implicit none

  private
  public ::                     &
    nl_operator_t,              &
    nl_operator_init,           &
    nl_operator_equal,          &
    nl_operator_build,          &
    nl_operator_transpose,      &
    dnl_operator_operate,       &
    znl_operator_operate,       &
    dnl_operator_operate_diag,  &
    znl_operator_operate_diag,  &
    nl_operator_end,            &
    nl_operator_skewadjoint,    &
    nl_operator_selfadjoint,    &
    nl_operator_write

  type nl_operator_t
    type(mesh_t), pointer :: m         ! pointer to the underlying mesh
    integer               :: n         ! number of points in discrete operator
    integer               :: np        ! number of points in mesh
    integer, pointer      :: stencil(:,:)
    integer               :: stencil_center ! center point of stencil

    ! When running in parallel mode, the next three
    ! arrays are unique on each node.
    integer, pointer      :: i(:,:)    ! index of the points
    FLOAT,   pointer      :: w_re(:,:) ! weightsp, real part
    FLOAT,   pointer      :: w_im(:,:) ! weightsp, imaginary part

    logical               :: const_w   ! are the weights independent of index i
    logical               :: cmplx_op  ! .true. if we have also imaginary weights
  end type nl_operator_t

  interface assignment (=)
    module procedure nl_operator_equal
  end interface

contains


  ! ---------------------------------------------------------
  subroutine nl_operator_init(op, n)
    type(nl_operator_t), intent(out) :: op
    integer,             intent(in)  :: n

    ASSERT(n  > 0)

    call push_sub('nl_operator.nl_operator_init')

    op%n  = n
    ALLOCATE(op%stencil(3, n), 3*n)

    call pop_sub()
  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_equal(opo, opi)
    type(nl_operator_t), intent(out) :: opo
    type(nl_operator_t), intent(in)  :: opi

    call push_sub('nl_operator.nl_operator_equal')

    call nl_operator_init(opo, opi%n)

    opo%np       =  opi%np
    opo%m        => opi%m
    opo%cmplx_op =  opi%cmplx_op
    opo%stencil(1:3, 1:opo%n) = opi%stencil(1:3, 1:opi%n)
    ALLOCATE(opo%i(opi%n, opi%np), opi%n*opi%np)

    if(opi%const_w) then
      ALLOCATE(opo%w_re(opi%n, 1), opi%n*1)
      if (opi%cmplx_op) then
        ALLOCATE(opo%w_im(opi%n, 1), opi%n*1)
      end if
    else
      ALLOCATE(opo%w_re(opi%n, opi%np), opi%n*opi%np)
      if (opi%cmplx_op) then
        ALLOCATE(opo%w_im(opi%n, opi%np), opi%n*opi%np)
      end if
    end if

    opo%const_w = opi%const_w
    opo%i       = opi%i
    opo%w_re    = opi%w_re
    if (opi%cmplx_op) then
      opo%w_im  = opi%w_im
    end if

    call pop_sub()
  end subroutine nl_operator_equal


  ! ---------------------------------------------------------
  subroutine nl_operator_build(m, op, np, const_w, cmplx_op)
    type(mesh_t), target, intent(in)    :: m
    type(nl_operator_t),  intent(inout) :: op
    integer,              intent(in)    :: np       ! Number of (local) points.
    logical, optional,    intent(in)    :: const_w  ! are the weights constant (independent of the point)
    logical, optional,    intent(in)    :: cmplx_op ! do we have complex weights?

    integer :: i, j, p1(MAX_DIM)

    call push_sub('nl_operator.nl_operator_build')

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
      ALLOCATE(op%w_re(op%n, 1), op%n*1)
      if (op%cmplx_op) then
        ALLOCATE(op%w_im(op%n, 1), op%n*1)
      end if
      message(1) = 'Info: nl_operator_build: working with constant weights.'
      if(in_debug_mode) call write_info(1)
    else
      ALLOCATE(op%w_re(op%n, op%np), op%n*op%np)
      if (op%cmplx_op) then
        ALLOCATE(op%w_im(op%n, op%np), op%n*op%np)
      end if
      message(1) = 'Info: nl_operator_build: working with non-constant weights.'
      if(in_debug_mode) call write_info(1)
    end if

    ! set initially to zero
    op%w_re = M_ZERO
    if (op%cmplx_op) op%w_im = M_ZERO

    ! store center point of the stencil
    do i = 1, op%n
      if(                              & 
        op%stencil(1, i) .eq. 0 .and.  &
        op%stencil(2, i) .eq. 0 .and.  &
        op%stencil(3, i) .eq. 0        &
        ) then
        op%stencil_center = i
        exit
      end if
    end do

    ! Build lookup table op%i from stencil.
    ALLOCATE(op%i(op%n, np), op%n*np)

    do i = 1, np
      if(m%parallel_in_domains) then
        ! When running in parallel, get global number of
        ! point i.
        p1(:) = m%Lxyz(m%vp%local(m%vp%xlocal(m%vp%partno)+i-1), :)
      else
        p1(:) = m%Lxyz(i, :)
      end if

      do j = 1, op%n
        ! Get global index of p1 plus current stencil point.
        op%i(j, i) = mesh_index(m%sb%dim, m%sb%periodic_dim, m%nr,    &
          m%Lxyz_inv, p1(:) + op%stencil(:, j))

        if(m%parallel_in_domains) then
          ! When running parallel, translate this global
          ! number back to a local number.

          op%i(j, i) = m%vp%global(op%i(j, i), m%vp%partno)
        end if

      end do
    end do

    call pop_sub()
  end subroutine nl_operator_build


  ! ---------------------------------------------------------
  subroutine nl_operator_transpose(op, opt)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opt

    integer :: np, i, j, index, l, k

    call push_sub('nl_operator.nl_operator_transpose')

    np = op%np
    opt = op
    opt%w_re = M_ZERO
    if (op%cmplx_op) opt%w_im = M_ZERO
    do i = 1, op%np
      do j = 1, op%n
        index = op%i(j, i)
        if(index <= op%np) then
          do l = 1, op%n
            k = op%i(l, index)
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

    call pop_sub()
  end subroutine nl_operator_transpose


  ! ---------------------------------------------------------
  ! opt has to be initialised and built.
  subroutine nl_operator_skewadjoint(op, opt, m)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opt
    type(mesh_t),        intent(in)  :: m

    integer          :: np, i, j, index, l, k
    integer, pointer :: op_i(:, :)
    FLOAT, pointer   :: vol_pp(:)
    FLOAT, pointer   :: w_re(:, :), w_re_t(:, :)
    FLOAT, pointer   :: w_im(:, :), w_im_t(:, :)

#if defined(HAVE_MPI)
    type(nl_operator_t) :: opg, opgt
#endif

    call push_sub('nl_operator.nl_operator_skewadjoint')

    opt = op

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%n)
      call nl_operator_equal(opgt, opg)
      ALLOCATE(vol_pp(m%np_global), m%np_global)
      call dvec_allgather(m%vp, vol_pp, m%vol_pp)

      op_i   => opg%i
      w_re   => opg%w_re
      w_re_t => opgt%w_re
      if(op%cmplx_op) then
        w_im   => opg%w_im
        w_im_t => opgt%w_im
      end if
#else
      ASSERT(.false.)
#endif
    else
      op_i   => op%i
      w_re   => op%w_re
      w_re_t => opt%w_re
      if(op%cmplx_op) then
        w_im   => op%w_im
        w_im_t => opt%w_im
      end if
      vol_pp => m%vol_pp
    end if

    np = m%np_global
    w_re_t = M_ZERO
    if (op%cmplx_op) w_im_t = M_ZERO
    do i = 1, m%np_global
      do j = 1, op%n
        index = op_i(j, i)
        if(index <= m%np_global) then
          do l = 1, op%n
            k = op_i(l, index)
            if( k == i ) then
              if(.not.op%const_w) then
                w_re_t(j, i) = M_HALF*w_re(j, i) - M_HALF*(vol_pp(index)/vol_pp(i))*w_re(l, index)
                if (op%cmplx_op) &
                   w_im_t(j, i) = M_HALF*w_im(j, i) - M_HALF*(vol_pp(index)/vol_pp(i))*w_im(l, index)
              else
                w_re_t(j, 1) = w_re(l, 1)
                if (op%cmplx_op) w_im_t(j, 1) = w_im(l, 1)
              end if
            end if
          end do
        end if
      end do
    end do

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      deallocate(vol_pp)
      do i = 1, m%vp%np_local(m%vp%partno)
        opt%w_re(:, i) = w_re_t(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        if(opt%cmplx_op) then
          opt%w_im(:, i) = w_im_t(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
#endif
    end if

    call pop_sub()
  end subroutine nl_operator_skewadjoint


  ! ---------------------------------------------------------
  subroutine nl_operator_selfadjoint(op, opt, m)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opt
    type(mesh_t),        intent(in)  :: m

    integer          :: np, i, j, index, l, k
    integer, pointer :: op_i(:, :)
    FLOAT, pointer   :: vol_pp(:)
    FLOAT, pointer   :: w_re(:, :), w_re_t(:, :)
    FLOAT, pointer   :: w_im(:, :), w_im_t(:, :)

#if defined(HAVE_MPI)
    type(nl_operator_t) :: opg, opgt
#endif

    call push_sub('nl_operator.nl_operator_selfadjoint')

    np = op%np
    opt = op

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      call nl_operator_allgather(op, opg)
      call nl_operator_init(opgt, op%n)
      opgt = opg
      ALLOCATE(vol_pp(m%np_global), m%np_global)
      call dvec_allgather(m%vp, vol_pp, m%vol_pp)

      op_i   => opg%i
      w_re   => opg%w_re
      w_re_t => opgt%w_re
      if(op%cmplx_op) then
        w_im   => opg%w_im
        w_im_t => opgt%w_im
      end if
#else
      ASSERT(.false.)
#endif
    else
      op_i   => op%i
      w_re   => op%w_re
      w_re_t => opt%w_re
      if(op%cmplx_op) then
        w_im   => op%w_im
        w_im_t => opt%w_im
      end if
      vol_pp => m%vol_pp
    end if

    np = m%np_global
    w_re_t = M_ZERO
    if (op%cmplx_op) w_im_t = M_ZERO
    do i = 1, m%np_global
      do j = 1, op%n
        index = op_i(j, i)
        if(index <= m%np_global) then
          do l = 1, op%n
            k = op_i(l, index)
            if( k == i ) then
              if(.not.op%const_w) then
                w_re_t(j, i) = M_HALF*w_re(j, i) + M_HALF*(vol_pp(index)/vol_pp(i))*w_re(l, index)
                if (op%cmplx_op) &
                  w_im_t(j, i) = M_HALF*w_im(j, i) + M_HALF*(vol_pp(index)/vol_pp(i))*w_im(l, index)
              else
                w_re_t(j, 1) = w_re(l, 1)
                if (op%cmplx_op) w_im_t(j, 1) = w_im(l, 1)
              end if
            end if
          end do
        end if
      end do
    end do

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      deallocate(vol_pp)
      do i = 1, m%vp%np_local(m%vp%partno)
        opt%w_re(:, i) = w_re_t(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        if(opt%cmplx_op) then
          opt%w_im(:, i) = w_im_t(:, m%vp%local(m%vp%xlocal(m%vp%partno)+i-1))
        end if
      end do
      call nl_operator_end(opg)
      call nl_operator_end(opgt)
    end if
#endif

    call pop_sub()
  end subroutine nl_operator_selfadjoint


#if defined(HAVE_MPI)
  ! ---------------------------------------------------------
  ! Collects a distributed non local operator op into opg
  ! on the root node. nl_operator_end has to be called
  ! on opg when no longer needed.
  subroutine nl_operator_gather(op, opg)
    type(nl_operator_t), intent(in)  :: op  ! Local operator.
    type(nl_operator_t), intent(out) :: opg ! Global operator.

    integer :: i

    call push_sub('nl_operator.nl_operator_gather')

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
    do i = 1, op%n
      call ivec_gather(op%m%vp, opg%i(i, :), op%i(i, :))
    end do
    if(op%m%vp%rank.eq.op%m%vp%root) then
      call nl_operator_translate_indices(opg)
    end if
    if(in_debug_mode) call write_debug_newlines(2)

    ! Weights have to be collected only if they are
    ! not constant.
    if(.not.op%const_w) then
      do i = 1, op%n
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

    call nl_operator_init(op, opg%n)
    call nl_operator_build(opg%m, op, opg%m%np, opg%const_w, opg%cmplx_op)

    do i = 1, opg%n
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
    do i = 1, op%n
      call ivec_allgather(op%m%vp, opg%i(i, :), op%i(i, :))
    end do
    call nl_operator_translate_indices(opg)

    ! Weights have to be collected only if they are
    ! not constant.
    if(.not.op%const_w) then
      do i = 1, op%n
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
  ! This can be considered as nl_operator_equal and
  ! reallocating w_re, w_im and i.
  subroutine nl_operator_common_copy(op, opg)
    type(nl_operator_t), intent(in)  :: op
    type(nl_operator_t), intent(out) :: opg

    call push_sub('nl_operator.nl_operator_common_copy')

    call nl_operator_init(opg, op%n)
    ALLOCATE(opg%i(op%n, op%m%np_global), op%n*op%m%np_global)
    if(op%const_w) then
      ALLOCATE(opg%w_re(op%n, 1), op%n*1)
      if(op%cmplx_op) then
        ALLOCATE(opg%w_im(op%n, 1), op%n*1)
      end if
    else
      ALLOCATE(opg%w_re(op%n, op%m%np_global), op%n*op%m%np_global)
      if(op%cmplx_op) then
        ALLOCATE(opg%w_im(op%n, op%m%np_global), op%n*op%m%np_global)
      end if
    end if
    opg%m        => op%m
    opg%np       =  op%m%np_global
    opg%stencil  =  op%stencil
    opg%cmplx_op =  op%cmplx_op
    opg%const_w  =  op%const_w
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

    do i = 1, opg%n
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
  subroutine nl_operator_op_to_matrix(op, a, b)
    type(nl_operator_t), intent(in) :: op
    FLOAT, intent(out)                 :: a(:, :)
    FLOAT, optional, intent(out)       :: b(:, :)

    integer          :: i, j, k, index
    integer, pointer :: op_i(:, :)
    FLOAT, pointer   :: w_re(:, :), w_im(:, :)
    type(nl_operator_t) :: opg

    call push_sub('nl_operator.nl_operator_op_to_matrix')

    if(op%m%parallel_in_domains) then
#if defined(HAVE_MPI)
      call nl_operator_gather(op, opg)
      ! Take these shortcuts.
      op_i => opg%i
      w_re => opg%w_re
      if(opg%cmplx_op) w_im => opg%w_im
#else
      ASSERT(.false.)
#endif
    else
      ! Take these shortcuts.
      op_i => op%i
      w_re => op%w_re
      if(op%cmplx_op) w_im => op%w_im
    end if

    if(mpi_grp_is_root(op%m%mpi_grp)) then
      k = 1
      do i = 1, op%m%np_global
        if(.not.op%const_w) k = i
        do j = 1, op%n
          index = op_i(j, i)
          if(index <= op%m%np_global) then
            a(i, index) = w_re(j, k)
            if (op%cmplx_op) b(i, index) = w_im(j, k)
          end if
        end do
      end do
      
      if(op%m%parallel_in_domains) call nl_operator_end(opg)
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

    op = op_ref
    do i = 1, op%np
      do j = 1, op%n
        index = op%i(j, i)
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
      ALLOCATE(a(op%m%np_global, op%m%np_global), op%m%np_global*op%m%np_global)
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

      deallocate(a)
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
      ALLOCATE(a(op%m%np_global, op%m%np_global), op%m%np_global*op%m%np_global)
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

      deallocate(a)
    end if

    call pop_sub()

  end subroutine nl_operatorT_write


  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_t), intent(inout) :: op

    call push_sub('nl_operator.nl_operator_end')

    ASSERT(associated(op%i))
    ASSERT(associated(op%w_re))
    ASSERT(associated(op%stencil))

    if (op%cmplx_op) then
      ASSERT(associated(op%w_im))
    end if

    deallocate(op%i, op%w_re, op%stencil)
    nullify   (op%i, op%w_re, op%stencil)

    if (op%cmplx_op) then
      deallocate(op%w_im)
      nullify   (op%w_im)
    end if

    call pop_sub()
  end subroutine nl_operator_end


  ! ---------------------------------------------------------
  ! calculates fo = op fi
  ! ---------------------------------------------------------
  subroutine dnl_operator_operate(op, fi, fo)
    FLOAT,               intent(inout) :: fi(:)  ! fi(op%np_part)
    type(nl_operator_t), intent(in)    :: op
    FLOAT,               intent(out)   :: fo(:)  ! fo(op%np_part)

    integer :: ii, nn

    call profiling_in(C_PROFILING_NL_OPERATOR)
    call push_sub('nl_operator.dnl_operator_operate')

#if defined(HAVE_MPI)
    if(op%m%parallel_in_domains) then
      call dvec_ghost_update(op%m%vp, fi)
    end if
#endif

    nn = op%n
    if(op%const_w) then
      do ii = 1, op%np
        fo(ii) = sum(op%w_re(1:nn, 1)  * fi(op%i(1:nn, ii)))
      end do
    else
      do ii = 1, op%np
        fo(ii) = sum(op%w_re(1:nn, ii) * fi(op%i(1:nn, ii)))
      end do
    end if
    do ii = op%np + 1, size(fo)
      fo(ii) = M_ZERO
    end do

    call pop_sub()
    call profiling_out(C_PROFILING_NL_OPERATOR)
  end subroutine dnl_operator_operate


  ! ---------------------------------------------------------
  subroutine znl_operator_operate(op, fi, fo)
    CMPLX,               intent(inout) :: fi(:)  ! fi(op%np)
    type(nl_operator_t), intent(in)    :: op
    CMPLX,               intent(out)   :: fo(:)  ! fo(op%np)

    integer :: ii, nn

    call profiling_in(C_PROFILING_NL_OPERATOR)
    call push_sub('nl_operator.znl_operator_operate')

#if defined(HAVE_MPI)
    if(op%m%parallel_in_domains) then
      call zvec_ghost_update(op%m%vp, fi)
    end if
#endif

    nn = op%n
    if(op%cmplx_op) then
      if(op%const_w) then
        do ii = 1, op%np
          fo(ii) = sum(cmplx(op%w_re(1:nn, 1),  op%w_im(1:nn, 1))  * fi(op%i(1:nn, ii)))
        end do
      else
        do ii = 1, op%np
          fo(ii) = sum(cmplx(op%w_re(1:nn, ii), op%w_im(1:nn, ii)) * fi(op%i(1:nn, ii)))
        end do
      end if
    else
      if(op%const_w) then
        do ii = 1, op%np
          fo(ii) = sum(op%w_re(1:nn, 1)  * fi(op%i(1:nn, ii)))
        end do
      else
        do ii = 1, op%np
          fo(ii) = sum(op%w_re(1:nn, ii) * fi(op%i(1:nn, ii)))
        end do
      end if
    end if
    do ii = op%np + 1, size(fo)
      fo(ii) = M_ZERO
    end do

    call pop_sub()
    call profiling_out(C_PROFILING_NL_OPERATOR)
  end subroutine znl_operator_operate


  ! ---------------------------------------------------------
  subroutine dnl_operator_operate_diag(op, fo)
    type(nl_operator_t), intent(in)    :: op
    FLOAT,               intent(out)   :: fo(:)  ! fo(op%np_part)

    integer :: ii, nn, jj

    call profiling_in(C_PROFILING_NL_OPERATOR)
    call push_sub('nl_operator.dnl_operator_operate_diag')

    nn = op%n
    if(op%const_w) then
      do ii = 1, nn
        if( 1 == op%i(ii,1) ) then
          fo(1:op%np) = op%w_re(ii, 1)
          exit
        end if
      end do
    else
      do ii = 1, op%np
        do jj = 1, nn
          if( ii == op%i(jj,ii) ) then
            fo(ii) = op%w_re(jj, ii)
            exit
          end if
        end do
      end do
    end if

    call pop_sub()
    call profiling_out(C_PROFILING_NL_OPERATOR)

  end subroutine dnl_operator_operate_diag


  ! ---------------------------------------------------------
  subroutine znl_operator_operate_diag(op, fo)
    type(nl_operator_t), intent(in)    :: op
    CMPLX,               intent(out)   :: fo(:)  ! fo(op%np)

    integer :: ii, nn, jj

    call profiling_in(C_PROFILING_NL_OPERATOR)
    call push_sub('nl_operator.znl_operator_operate_diag')

    nn = op%n

    if(op%cmplx_op) then
      if(op%const_w) then
        do ii = 1, nn
          if( 1 == op%i(ii,1) ) then
            fo(1:op%np) = cmplx(op%w_re(ii, 1), op%w_im(ii, 1))
            exit
          end if
        end do
      else
        do ii = 1, op%np
          do jj = 1, nn
            if( ii == op%i(jj,ii) ) then
              fo(ii) = cmplx(op%w_re(jj, ii), op%w_im(jj, ii))
              exit
            end if
          end do
        end do
      end if
    else
      if(op%const_w) then
        do ii = 1, nn
          if( 1 == op%i(ii,1) ) then
            fo(1:op%np) = op%w_re(ii, 1)
            exit
          end if
        end do
      else
        do ii = 1, op%np
          do jj = 1, nn
            if( ii == op%i(jj,ii) ) then
              fo(ii) = op%w_re(jj, ii)
              exit
            end if
          end do
        end do
      end if
    end if

    call pop_sub()
    call profiling_out(C_PROFILING_NL_OPERATOR)

  end subroutine znl_operator_operate_diag

end module nl_operator_m
