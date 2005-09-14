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

module nl_operator
  use global
  use messages
  use mesh
  use simul_box
  use io
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  use par_vec
#endif

  implicit none

  private
  public :: nl_operator_type,      &
       nl_operator_init,           &
       nl_operator_equal,          &
       nl_operator_build,          &
       nl_operator_transpose,      &
       dnl_operator_operate,       &
       znl_operator_operate,       &
       znl_operator_operate_cmplx, &
       nl_operator_end,            &
       nl_operator_skewadjoint,    &
       nl_operator_selfadjoint, &
nl_operator_write

  type nl_operator_type
     type(mesh_type), pointer :: m         ! pointer to the underlying mesh
     integer                  :: n         ! number of points in discrete operator
     integer                  :: np        ! number of points in mesh
     integer, pointer         :: stencil(:,:)

     integer, pointer         :: i(:,:)    ! index of the points
     FLOAT,   pointer         :: w_re(:,:) ! weightsp, real part
     FLOAT,   pointer         :: w_im(:,:) ! weightsp, imaginary part

     logical                  :: const_w   ! are the weights independent of i
     logical                  :: cmplx_op  ! .true. if we have also imaginary weights
  end type nl_operator_type

  interface assignment (=)
     module procedure nl_operator_equal
  end interface

contains


  ! ---------------------------------------------------------
  subroutine nl_operator_init(op, n)
    type(nl_operator_type), intent(out) :: op
    integer,                intent(in)  :: n

    ASSERT(n  > 0)

    call push_sub('nl_operator.nl_operator_init')

    op%n  = n
    allocate(op%stencil(3, n))

    call pop_sub()
  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_equal(opo, opi)
    type(nl_operator_type), intent(out) :: opo
    type(nl_operator_type), intent(in)  :: opi

    call push_sub('nl_operator.nl_operator_equal')

    call nl_operator_init(opo, opi%n)

    opo%np = opi%np
    opo%m  => opi%m
    opo%stencil(1:3, 1:opo%n) = opi%stencil(1:3, 1:opi%n)
    allocate(opo%i(opi%n, opi%np))
    if(opi%const_w) then
       allocate(opo%w_re(opi%n, 1))
       if (opi%cmplx_op) then
          allocate(opo%w_im(opi%n, 1))
       endif
    else
       allocate(opo%w_re(opi%n, opi%np))
       if (opi%cmplx_op) then
          allocate(opo%w_im(opi%n, opi%np))
       endif
    endif

    opo%const_w = opi%const_w
    opo%i       = opi%i
    opo%w_re    = opi%w_re
    if (opi%cmplx_op) then
       opo%w_im  = opi%w_im
    endif

    call pop_sub()
  end subroutine nl_operator_equal


  ! ---------------------------------------------------------
  subroutine nl_operator_build(m, op, np, const_w, cmplx_op)
    type(mesh_type), target, intent(in)    :: m
    type(nl_operator_type),  intent(inout) :: op       
    integer,                 intent(in)    :: np       ! Number of (local)
                                                      ! points.
    logical, optional,      intent(in)    :: const_w  ! are the weights constant (independent of the point)
    logical, optional,      intent(in)    :: cmplx_op ! do we have complex weights?

    integer :: i, j, p1(3)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    integer :: ierr ! MPI errorcode.
    integer :: rank ! Current nodes rank.
#endif

    call push_sub('nl_operator.nl_operator_build')

    ASSERT(np > 0)

    ! store values in structure
    op%np       = np
    op%m        => m
    op%const_w  = .false.
    op%cmplx_op = .false.
    if(present(const_w )) op%const_w = const_w
    if(present(cmplx_op)) op%const_w = cmplx_op

    ! allocate weights op%w
    if(op%const_w) then
       allocate(op%w_re(op%n, 1))
       if (op%cmplx_op) then
          allocate(op%w_im(op%n, 1))
       endif
       message(1) = 'Info: nl_operator_build: working with constant weights.'
       if(conf%verbose > VERBOSE_DEBUG) call write_info(1)
    else
       allocate(op%w_re(op%n, op%np))
       if (op%cmplx_op) then
          allocate(op%w_im(op%n, op%np))
       endif
       message(1) = 'Info: nl_operator_build: working with non-constant weights.'
       if(conf%verbose > VERBOSE_DEBUG) call write_info(1)
    end if

    ! set initially to zero
    op%w_re = M_ZERO
    if (op%cmplx_op) op%w_im = M_ZERO

    ! Build lookup table op%i from stencil.
    allocate(op%i(op%n, np))

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    call MPI_Comm_rank(m%vp%comm, rank, ierr)
#endif

    do i = 1, np
#if defined(HAVE_MPI) && defined(HAVE_METIS)
      ! When running in parallel, get global number of
      ! point i.
      p1(:) = m%Lxyz(m%vp%local(m%vp%xlocal(rank+1)+i-1), :)
#else
      p1(:) = m%Lxyz(i, :)
#endif

      do j = 1, op%n
        ! Get global index of p1 plus current stencil point.
        op%i(j, i) = mesh_index(m, p1(:) + op%stencil(:, j), 1)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
        ! When running parallel, translate this global
        ! number back to a local number.
        op%i(j, i) = m%vp%global(op%i(j, i), rank+1)
#endif
      end do
    end do

    call pop_sub()

  end subroutine nl_operator_build


  ! ---------------------------------------------------------
  subroutine nl_operator_transpose(op, opt)
    type(nl_operator_type), intent(in)  :: op
    type(nl_operator_type), intent(out) :: opt

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
                   endif
                endif
             enddo
          endif
       enddo
    enddo

    call pop_sub()
  end subroutine nl_operator_transpose

  ! ---------------------------------------------------------
  subroutine nl_operator_skewadjoint(op, opt, m)
    type(nl_operator_type), intent(in)  :: op
    type(nl_operator_type), intent(out) :: opt
    type(mesh_type),        intent(in)  :: m

    integer :: np, i, j, index, l, k

    call push_sub('nl_operator.nl_operator_skewadjoint')

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
                      opt%w_re(j, i) = M_HALF*op%w_re(j, i) - M_HALF*(m%vol_pp(index)/m%vol_pp(i))*op%w_re(l, index)
                      if (op%cmplx_op) &
                           opt%w_im(j, i) = M_HALF*op%w_im(j, i) - M_HALF*(m%vol_pp(index)/m%vol_pp(i))*op%w_im(l, index)
                   else
                      opt%w_re(j, 1) = op%w_re(l, 1)
                      if (op%cmplx_op) opt%w_im(j, 1) = op%w_im(l, 1)
                   endif
                endif
             enddo
          endif
       enddo
    enddo

    call pop_sub()
  end subroutine nl_operator_skewadjoint

  ! ---------------------------------------------------------
  subroutine nl_operator_selfadjoint(op, opt, m)
    type(nl_operator_type), intent(in)  :: op
    type(nl_operator_type), intent(out) :: opt
    type(mesh_type),        intent(in)  :: m

    integer :: np, i, j, index, l, k

    call push_sub('nl_operator.nl_operator_selfadjoint')

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
                      opt%w_re(j, i) = M_HALF*op%w_re(j, i) + M_HALF*(m%vol_pp(index)/m%vol_pp(i))*op%w_re(l, index)
                      if (op%cmplx_op) &
                           opt%w_im(j, i) = M_HALF*op%w_im(j, i) + M_HALF*(m%vol_pp(index)/m%vol_pp(i))*op%w_im(l, index)
                   else
                      opt%w_re(j, 1) = op%w_re(l, 1)
                      if (op%cmplx_op) opt%w_im(j, 1) = op%w_im(l, 1)
                   endif
                endif
             enddo
          endif
       enddo
    enddo

    call pop_sub()
  end subroutine nl_operator_selfadjoint


  ! ---------------------------------------------------------
  ! When running in parallel only the root node
  ! creates the matrix. But all nodes have to
  ! call this routine because the distributed operator has
  ! to be collected.
  subroutine nl_operator_op_to_matrix(op, a, b)
    type(nl_operator_type), intent(in) :: op
    FLOAT, intent(out)                 :: a(:, :)
    FLOAT, optional, intent(out)       :: b(:, :)

    integer          :: i, j, k, index
    integer, pointer :: op_i(:, :)
    FLOAT, pointer   :: w_re(:, :), w_im(:, :)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    integer :: il, ig
    integer :: rank, ierr
#endif

    call push_sub('nl_operator.nl_operator_op_to_matrix')

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    ! When running parallel, we have to collect the operator first.
    call MPI_Comm_rank(op%m%vp%comm, rank, ierr)

    ! Collect the distributed index table op%i into op_i.
    ! Only root collects the index table.
    if(rank.eq.op%m%vp%root) then
      allocate(op_i(op%n, op%m%np_glob))
    end if
    ! Collect for every stencil point.
    do i = 1, op%n
      call ivec_gather(op%m%vp, op_i(i, :), op%i(i, :))
      ! All point numbers gathered for stencil point i
      ! have to be translated to global points numbers
      ! This, again, is work for root.
      if(rank.eq.op%m%vp%root) then
        do j = 1, op%m%np_glob
          il = op%m%vp%np_local(op%m%part(j))
          ig = il+op%m%vp%np_ghost(op%m%part(j))
          ! op_i(i, j) is a local point number, i. e. it can be
          ! a real local point (i. e. the local point number
          ! is less or equal than the number of local points of
          ! the node which owns the point with global number j):
          if(op_i(i, j).le.il) then
            ! Write the global point number from the lookup
            ! table in op_(i, j).
            op_i(i, j) = op%m%vp%local(op%m%vp%xlocal(op%m%part(j)) &
                         +op_i(i, j)-1)
          ! Or a ghost point:
          else if(op_i(i, j).gt.il.and.op_i(i, j).le.ig) then
            op_i(i, j) = op%m%vp%ghost(op%m%vp%xghost(op%m%part(j)) &
                         +op_i(i, j)-1-il)
          ! Or a boundary point:
          else if(op_i(i, j).gt.ig) then
            op_i(i, j) = op%m%vp%bndry(op%m%vp%xbndry(op%m%part(j)) &
                         +op_i(i, j)-1-ig)
            !op_i(i, j) = 0
          end if
        end do
      end if
    end do


    ! Weights have to be collected only if they are
    ! non constant.
    if(.not.op%const_w) then
      if(rank.eq.op%m%vp%root) then
        allocate(w_re(op%n, op%m%np_glob))
        if(op%cmplx_op) allocate(w_im(op%n, op%m%np_glob))
      end if
      do i = 1, op%n
        call dvec_gather(op%m%vp, w_re(i, :), op%w_re(i, :))
        if(op%cmplx_op) call dvec_gather(op%m%vp, w_im(i, :), op%w_im(i, :))
      end do
    ! Otherwise, those of root are just fine.
    else
      w_re => op%w_re
      if(op%cmplx_op) w_im => op%w_im
    end if
#else
    ! Take these shortcuts.
    op_i => op%i
    w_re => op%w_re
    if(op%cmplx_op) w_im => op%w_im
#endif

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(rank.eq.op%m%vp%root) then
#endif
      k = 1
      do i = 1, op%m%np_glob
         if(.not.op%const_w) k = i
         do j = 1, op%n
            index = op_i(j, i)
            if(index <= op%m%np_glob) then
               a(i, index) = w_re(j, k)
               if (op%cmplx_op) b(i, index) = w_im(j, k)
            endif
         enddo
      enddo
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    endif
#endif

    call pop_sub()

  end subroutine nl_operator_op_to_matrix


  ! ---------------------------------------------------------
  subroutine nl_operator_matrix_to_op(op_ref, op, a, b)
    FLOAT, intent(in)                   :: a(:, :)
    FLOAT, optional, intent(in)         :: b(:, :)
    type(nl_operator_type), intent(in)  :: op_ref
    type(nl_operator_type), intent(out) :: op

    integer :: i, j, index

    call push_sub('nl_operator.nl_operator_matrix_to_op')

    op = op_ref
    do i = 1, op%np
       do j = 1, op%n
          index = op%i(j, i)
          if(index <= op%np) &
               op%w_re(j, i) = a(i, index)
          if (op%cmplx_op) op%w_im(j, i) = b(i, index)
       enddo
    enddo

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
    type(nl_operator_type), intent(in) :: op

    integer            :: i, j
    integer            :: unit
    FLOAT, allocatable :: a(:, :)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    integer :: rank, ierr
#endif

    call push_sub('nl_operator.nl_operator_write')
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    call MPI_Comm_rank(op%m%vp%comm, rank, ierr)
#endif

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(rank.eq.op%m%vp%root) then
#endif
      allocate(a(op%m%np_glob, op%m%np_glob))
      a = M_ZERO
#if defined(HAVE_MPI) && defined(HAVE_METIS) 
    end if
#endif

      call nl_operator_op_to_matrix(op, a)

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(rank.eq.op%m%vp%root) then
#endif
      unit = io_open(filename, action='write')
      if(unit < 0) then
        message(1) = 'Could not open file '//filename//' to write operator.'
        call write_fatal(1)
      end if
      do i = 1, op%m%np_glob
         do j = 1, op%m%np_glob - 1
            write(unit, fmt = '(f9.4)', advance ='no') a(i, j)
         enddo
         write(unit, fmt = '(f9.4)') a(i, op%m%np_glob)
      enddo
      call io_close(unit)

      deallocate(a)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    end if
#endif

    call pop_sub()

  end subroutine nl_operator_write


  ! ---------------------------------------------------------
  ! Same as nl_operator_write but transposed matrix.
  subroutine nl_operatorT_write(op, unit)

    type(nl_operator_type), intent(in) :: op
    integer, intent(in)                :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    integer :: rank, ierr
#endif

    call push_sub('nl_operator.nl_operator_write')
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    call MPI_Comm_rank(op%m%vp%comm, rank, ierr)
#endif

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(rank.eq.op%m%vp%root) then
#endif
      allocate(a(op%m%np_glob, op%m%np_glob))
      a = M_ZERO
#if defined(HAVE_MPI) && defined(HAVE_METIS) 
    end if
#endif

      call nl_operator_op_to_matrix(op, a)

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(rank.eq.op%m%vp%root) then
#endif
      do i = 1, op%m%np_glob
         do j = 1, op%m%np_glob - 1
            write(unit, fmt = '(f9.4)', advance ='no') a(j, i)
         enddo
         write(unit, fmt = '(f9.4)') a(op%m%np_glob, i)
      enddo

      deallocate(a)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
    end if
#endif

    call pop_sub()

  end subroutine nl_operatorT_write


  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_type), intent(inout) :: op

    call push_sub('nl_operator.nl_operator_end')

    ASSERT(associated(op%i))
    ASSERT(associated(op%w_re))
    ASSERT(associated(op%stencil))

    if (op%cmplx_op) then
       ASSERT(associated(op%w_im))
    endif

    deallocate(op%i, op%w_re, op%stencil)
    nullify   (op%i, op%w_re, op%stencil)

    if (op%cmplx_op) then
       deallocate(op%w_im)
       nullify   (op%w_im)
    endif

    call pop_sub()
  end subroutine nl_operator_end


  ! ---------------------------------------------------------
  ! calculates fo = op fi
  ! ---------------------------------------------------------
  subroutine dnl_operator_operate(op, fi, fo)
    FLOAT,                  intent(inout) :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)    :: op
    FLOAT,                  intent(out)   :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:)

    call push_sub('nl_operator.dnl_operator_operate')

#if defined(HAVE_MPI) && defined(HAVE_METIS) 
    call dvec_ghost_update(op%m%vp, fi)
#endif

    n = op%n
    if(op%const_w) then
       allocate(w_re(n))
       w_re(1:n) = op%w_re(1:n, 1)
       do i = 1, op%np
          fo(i) = sum(w_re(1:n)*fi(op%i(1:n,i)))
       end do
       deallocate(w_re)
    else
       do i = 1, op%np
          fo(i) = sum(op%w_re(1:n,i)*fi(op%i(1:n,i)))
       end do
    end if

    call pop_sub()
  end subroutine dnl_operator_operate


  ! ---------------------------------------------------------
  subroutine znl_operator_operate(op, fi, fo)
    CMPLX,                  intent(inout) :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)    :: op
    CMPLX,                  intent(out)   :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:)

    call push_sub('nl_operator.znl_operator_operate')

#if defined(HAVE_MPI) && defined(HAVE_METIS) 
    call zvec_ghost_update(op%m%vp, fi)
#endif

    n = op%n
    if(op%const_w) then
       allocate(w_re(n))
       w_re(1:n) = op%w_re(1:n, 1)
       do i = 1, op%np
          fo(i) = sum(  cmplx(  w_re(1:n)*real(fi(op%i(1:n,i))),  w_re(1:n)*aimag(fi(op%i(1:n, i))), PRECISION  )   )
       end do
       deallocate(w_re)
    else
       do i = 1, op%np
          fo(i) = sum(  cmplx(  op%w_re(1:n, i)*real(fi(op%i(1:n,i))),  op%w_re(1:n, i)*aimag(fi(op%i(1:n, i))), PRECISION  )   )
       end do
    end if

    call pop_sub()
  end subroutine znl_operator_operate


  ! ---------------------------------------------------------
  ! allow for complex operators
  ! ---------------------------------------------------------
  subroutine znl_operator_operate_cmplx(op, fi, fo)
    CMPLX,                  intent(inout) :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)    :: op
    CMPLX,                  intent(out)   :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:), w_im(:)

    call push_sub('nl_operator.znl_operator_operate_complex')

#if defined(HAVE_MPI) && defined(HAVE_METIS) 
    call zvec_ghost_update(op%m%vp, fi)
#endif

    n = op%n
    if(op%const_w) then
       allocate(w_re(n),w_im(n))
       w_re(1:n) = op%w_re(1:n, 1)
       w_im(1:n) = op%w_im(1:n, 1)
       do i = 1, op%np
          fo(i) = sum(  cmplx(  w_re(1:n)* real(fi(op%i(1:n,i))),  w_re(1:n)*aimag(fi(op%i(1:n, i))), PRECISION  )  ) &
               + sum(   cmplx( -w_im(1:n)*aimag(fi(op%i(1:n,i))),  w_im(1:n)* real(fi(op%i(1:n, i))), PRECISION  )  )
       end do
       deallocate(w_re, w_im)
    else
       do i = 1, op%np
          fo(i) = sum(  cmplx(  op%w_re(1:n, i)* real(fi(op%i(1:n,i))),  op%w_re(1:n, i)*aimag(fi(op%i(1:n, i))), PRECISION  )  )  &
               + sum(   cmplx( -op%w_im(1:n, i)*aimag(fi(op%i(1:n,i))),  op%w_im(1:n, i)* real(fi(op%i(1:n, i))), PRECISION  )  )
       end do
    end if

    call pop_sub()
  end subroutine znl_operator_operate_cmplx

end module nl_operator
