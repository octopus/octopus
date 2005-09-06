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
       nl_operator_selfadjoint

  type nl_operator_type
     integer          :: n          ! number of points in discrete operator
     integer          :: np         ! number of points in mesh
     integer, pointer :: stencil(:,:)

     integer, pointer :: i(:,:)     ! index of the points
     FLOAT,   pointer :: w_re(:,:)  ! weightsp, real part
     FLOAT,   pointer :: w_im(:,:)  ! weightsp, imaginary part

     logical          :: const_w    ! are the weights independent of i
     logical          :: cmplx_op   ! .true. if we have also imaginary weights
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
    type(mesh_type),        intent(in)    :: m
    type(nl_operator_type), intent(inout) :: op       
    integer,                intent(in)    :: np       ! Number of (local)
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
  subroutine nl_operator_op_to_matrix(op, a, b)
    type(nl_operator_type), intent(in) :: op
    FLOAT, intent(out)                 :: a(:, :)
    FLOAT, optional, intent(out)       :: b(:, :)

    integer :: i, j, k, index

    call push_sub('nl_operator.nl_operator_op_to_matrix')

    k = 1
    do i = 1, op%np
       if(.not.op%const_w) k = i
       do j = 1, op%n
          index = op%i(j, i)
          if(index <= op%np) then
             a(i, index) = op%w_re(j, k)
             if (op%cmplx_op) b(i, index) = op%w_im(j, k)
          endif
       enddo
    enddo
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
  subroutine nl_operator_write(op, unit)
    type(nl_operator_type), intent(in) :: op
    integer, intent(in)                :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)

    call push_sub('nl_operator.nl_operator_write')

    allocate(a(op%np, op%np))
    a = M_ZERO

    call nl_operator_op_to_matrix(op, a)

    do i = 1, op%np
       do j = 1, op%np - 1
          write(unit, fmt = '(f9.4)', advance ='no') a(i, j)
       enddo
       write(unit, fmt = '(f9.4)') a(i, op%np)
    enddo

    deallocate(a)
    call pop_sub()
  end subroutine nl_operator_write


  ! ---------------------------------------------------------
  subroutine nl_operatorT_write(op, unit)
    type(nl_operator_type), intent(in) :: op
    integer, intent(in)                :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)

    call push_sub('nl_operator.nl_operator_Twrite')

    allocate(a(op%np, op%np))
    a = M_ZERO

    call nl_operator_op_to_matrix(op, a)

    do i = 1, op%np
       do j = 1, op%np - 1
          write(unit, fmt = '(f9.4)', advance ='no') a(j, i)
       enddo
       write(unit, fmt = '(f9.4)') a(op%np, i)
    enddo

    deallocate(a)
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
    FLOAT,                  intent(in)  :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)  :: op
    FLOAT,                  intent(out) :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:)

    call push_sub('nl_operator.dnl_operator_operate')

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
    CMPLX,                  intent(in)  :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)  :: op
    CMPLX,                  intent(out) :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:)

    call push_sub('nl_operator.znl_operator_operate')

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
    CMPLX,                  intent(in)  :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)  :: op
    CMPLX,                  intent(out) :: fo(:)  ! fo(op%np)

    integer :: i, n
    FLOAT, allocatable :: w_re(:), w_im(:)

    call push_sub('nl_operator.znl_operator_operate_complex')

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
