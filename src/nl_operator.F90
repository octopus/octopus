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

#include "global.h"

module nl_operator
  use mesh

  implicit none

  type nl_operator_type
    integer          :: n         ! number of points in discrete operator
    integer          :: np        ! number of points in mesh
    integer, pointer :: stencil(:,:)

    integer, pointer :: i(:,:)    ! index of the points
    FLOAT,   pointer :: w(:,:)    ! weightsp

    logical          :: const_w   ! are the weights independent of i
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

    op%n  = n
    allocate(op%stencil(3, n))

  end subroutine nl_operator_init


  ! ---------------------------------------------------------
  subroutine nl_operator_equal(opo, opi)
    type(nl_operator_type), intent(out) :: opo
    type(nl_operator_type), intent(in)  :: opi
    call nl_operator_init(opo, opi%n)
    opo%np = opi%np
    opo%stencil(1:3, 1:opo%n) = opi%stencil(1:3, 1:opi%n)
!!$    allocate(opo%i(size(opi%i, 1), size(opi%i, 2)))
!!$    allocate(opo%w(size(opi%w, 1), size(opi%w, 2)))
    allocate(opo%i(opi%n, opi%np))
    if(opi%const_w) then
      allocate(opo%w(opi%n, 1))
    else
      allocate(opo%w(opi%n, opi%np))
    endif
    opo%i = opi%i
    opo%w = opi%w
    opo%const_w = opi%const_w
  end subroutine nl_operator_equal


  ! ---------------------------------------------------------
  subroutine nl_operator_build(m, op, np, const_w)
    type(mesh_type),        intent(in)    :: m
    type(nl_operator_type), intent(inout) :: op
    integer,                intent(in)    :: np
    logical, optional,      intent(in)    :: const_w  ! are the weights constant (independent of the point)

    integer :: i, j, ix(conf%dim)

    ASSERT(np > 0)

    ! store values in structure
    op%np = np
    op%const_w = .false.
    if(present(const_w)) op%const_w = const_w

    ! allocate weights op%w
    if(op%const_w) then
      allocate(op%w(op%n, 1))
    else
      allocate(op%w(op%n, op%np))
    end if

    ! build lookup table op%i from stencil
    allocate(op%i(op%n, op%np))

    do i = 1, m%np     ! for all points in mesh
      ix(:) = m%Lxyz(i,:)

      do j = 1, op%n   ! for all points in stencil
        op%i(j, i) = mesh_index(m, ix(1:conf%dim)+op%stencil(1:conf%dim,j), 1)
      end do
    end do
  end subroutine nl_operator_build


  ! ---------------------------------------------------------
  subroutine nl_operator_transpose(op, opt)
    type(nl_operator_type), intent(in)  :: op
    type(nl_operator_type), intent(out) :: opt

    integer :: np, i, j, index, l, k

    np = op%np
    opt = op
    !call nl_operator_equal(opt, op)
    opt%w = M_ZERO
    do i = 1, op%np
       do j = 1, op%n
          index = op%i(j, i)
          if(index<=op%np) then
             do l = 1, op%n
                k = op%i(l, index)
                if(k==i) then
                   if(.not.op%const_w) then
                      opt%w(j, i) = op%w(l, index)
                   else
                      opt%w(j, 1) = op%w(l, 1)
                   endif
                endif
             enddo
          endif
       enddo
    enddo

  end subroutine nl_operator_transpose


  ! ---------------------------------------------------------
  subroutine nl_operator_op_to_matrix(op, a)
    type(nl_operator_type), intent(in) :: op
    FLOAT, intent(out)                 :: a(:, :)

    integer :: i, j, k, index

    k = 1
    do i = 1, op%np
       if(.not.op%const_w) k = i
       do j = 1, op%n
          index = op%i(j, i)
          if(index<=op%np) then
            a(i, index) = op%w(j, k)
          endif
       enddo
    enddo
  end subroutine nl_operator_op_to_matrix


  ! ---------------------------------------------------------
  subroutine nl_operator_matrix_to_op(a, op_ref, op)
    FLOAT, intent(out)                  :: a(:, :)
    type(nl_operator_type), intent(in)  :: op_ref
    type(nl_operator_type), intent(out) :: op

    integer :: i, j, index

    op = op_ref
    do i = 1, op%np
       do j = 1, op%n
          index = op%i(j, i)
          if(index<=op%np) &
          op%w(j, i) = a(i, index)
       enddo
    enddo

  end subroutine nl_operator_matrix_to_op


  ! ---------------------------------------------------------
  subroutine nl_operator_write(op, unit)
    type(nl_operator_type), intent(in) :: op
    integer, intent(in)                  :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)
    character(len=50) :: fmt

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
  end subroutine nl_operator_write


  subroutine nl_operatorT_write(op, unit)
    type(nl_operator_type), intent(in) :: op
    integer, intent(in)                  :: unit

    integer :: i, j
    FLOAT, allocatable :: a(:, :)
    character(len=50) :: fmt

    allocate(a(op%np, op%np))
    a = M_ZERO

    call nl_operator_op_to_matrix(op, a)

    do i = 1, op%np
       do j = 1, op%np - 1
!!$          write(unit, fmt = '(f9.4)', advance ='no') a(i, j)
          write(unit, fmt = '(f9.4)', advance ='no') a(j, i)
       enddo
!!$       write(unit, fmt = '(f9.4)') a(i, op%np)
       write(unit, fmt = '(f9.4)') a(op%np, i)
    enddo

    deallocate(a)
  end subroutine nl_operatorT_write


  ! ---------------------------------------------------------
  subroutine nl_operator_end(op)
    type(nl_operator_type), intent(inout) :: op

    ASSERT(associated(op%i))
    ASSERT(associated(op%w))
    ASSERT(associated(op%stencil))

    deallocate(op%i, op%w, op%stencil)
    nullify   (op%i, op%w, op%stencil)
  end subroutine nl_operator_end


  ! ---------------------------------------------------------
  ! calculates fo = op fi
  ! ---------------------------------------------------------
  subroutine dnl_operator_operate(op, fi, fo)
    FLOAT,                  intent(in)  :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)  :: op
    FLOAT,                  intent(out) :: fo(:)  ! fo(op%np)

    integer :: i, n

    n = op%n
    if(op%const_w) then
      do i = 1, op%np
        fo(i) = sum(op%w(1:n,1)*fi(op%i(1:n,i)))
      end do
    else
      do i = 1, op%np
        fo(i) = sum(op%w(1:n,i)*fi(op%i(1:n,i)))
      end do
    end if

  end subroutine dnl_operator_operate


  subroutine znl_operator_operate(op, fi, fo)
    CMPLX,                  intent(in)  :: fi(:)  ! fi(op%np)
    type(nl_operator_type), intent(in)  :: op
    CMPLX,                  intent(out) :: fo(:)  ! fo(op%np)

    integer :: i, n

    n = op%n
    if(op%const_w) then
      do i = 1, op%np
        fo(i) = sum(op%w(1:n,1)*fi(op%i(1:n,i)))
      end do
    else
      do i = 1, op%np
        fo(i) = sum(op%w(1:n,i)*fi(op%i(1:n,i)))
      end do
    end if

  end subroutine znl_operator_operate

end module nl_operator
