!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module mesh_function
  use global
  use lib_basic_alg
  use mesh
  use cube_function

  implicit none

contains

subroutine mf_add_der(m)
  type(mesh_type) :: m

  integer :: unit, j, i
  call push_sub('mf_add_der')

  allocate(m%grad(conf%dim))

  ! Sets the order
  m%lapl%norder             = m%der_order
  m%grad(1:conf%dim)%norder = m%der_order

  ! Sets the number of points of the mesh
  m%lapl%np             = m%np
  m%grad(1:conf%dim)%np = m%np

  ! Sets the number of points per point for the operators
  m%grad(1:conf%dim)%n = (2*m%der_order + 1)
  m%lapl%n             = (2*m%der_order + 1)*conf%dim

  ! Now we can allocate
  allocate(m%lapl%i(m%lapl%n, m%lapl%np)); m%lapl%i = 1
  allocate(m%lapl%w(m%lapl%n, m%lapl%np)); m%lapl%w = M_ZERO
  do j = 1, conf%dim
     allocate(m%grad(j)%i(m%grad(j)%n, m%grad(j)%np)); m%grad(j)%i = 1
     allocate(m%grad(j)%w(m%grad(j)%n, m%grad(j)%np)); m%grad(j)%w = M_ZERO
  enddo

  call derivatives_init_diff(m, m%der_order)

  call pop_sub()
end subroutine mf_add_der

subroutine mf_end_der(m)
  type(mesh_type) :: m
  integer :: ierr, j

  call push_sub('mf_end')

  deallocate(m%lapl%i, m%lapl%w); nullify(m%lapl%i); nullify(m%lapl%w)
  do j = 1, conf%dim
     deallocate(m%grad(j)%i, m%grad(j)%w); nullify(m%grad(j)%i); nullify(m%grad(j)%w)
  enddo
  deallocate(m%grad); nullify(m%grad)

  call pop_sub()
end subroutine mf_end_der

subroutine derivatives_init_diff(m, order)
  type(mesh_type)     :: m
  integer, intent(in) :: order

  FLOAT, allocatable :: dgidfj(:) ! the coefficients for the gradient
  FLOAT, allocatable :: dlidfj(:) ! the coefficients for the laplacian
  integer :: i, j, ix(3), ik, in, nl, ng(3)

  call push_sub('derivatives_init_diff')

  call derivatives_coeff()

  do i = 1, m%np
    ix(:) = m%Lxyz(:, i)
    
    ! fill in the table
    nl = 0; ng = 0
    do j = 1, conf%dim
      do ik = -order, order
        ix(j) = ix(j) + ik
        
        in = mesh_index(m, ix, sign(j, ik))
        if(in > 0) then
          nl    = nl    + 1
          ng(j) = ng(j) + 1
          
          m%lapl%i(nl, i) = in
          m%lapl%w(nl, i) = dlidfj(abs(ik))/m%h(j)**2

          m%grad(j)%i(ng(j), i) = in
          m%grad(j)%w(ng(j), i) = dgidfj(ik)/m%h(j)
        end if
        
        ix(j) = ix(j) - ik
      end do
    end do
    
  end do        

  deallocate(dlidfj,dgidfj)
  call pop_sub()
contains

  subroutine derivatives_coeff
    use math, only: weights

    integer :: k, i, j, morder
    FLOAT, allocatable :: cc(:,:,:)

    call push_sub('derivatives_coeff')
    k = order
    if (k < 1) then
      write(message(1), '(a,i4,a)') "Input: '", k, "' is not a valid OrderDerivatives"
      message(2) = '(1 <= OrderDerivatives)'
      call write_fatal(2)
    end if
    allocate(dlidfj(-k:k), dgidfj(-k:k))
    morder = 2*k
    allocate(cc(0:morder, 0:morder, 0:2))
    call weights(2, morder, cc)
    dgidfj(0) = cc(0, morder, 1)
    dlidfj(0) = cc(0, morder, 2)
  
    j = 1
    do i = 1, k
      dgidfj(-i) = cc(j, morder, 1)
      dlidfj(-i) = cc(j, morder, 2)
      j = j + 1
      dgidfj( i) = cc(j, morder, 1)
      dlidfj( i) = cc(j, morder, 2)
      j = j + 1
    end do
    deallocate(cc)

    call pop_sub()
  end subroutine derivatives_coeff

end subroutine derivatives_init_diff

#include "undef.F90"
#include "real.F90"
#include "mf_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mf_inc.F90"

end module mesh_function
