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
    FLOAT,   pointer :: w(:,:)    ! weights

    logical          :: const_w   ! are the weights independent of i
  end type nl_operator_type

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
