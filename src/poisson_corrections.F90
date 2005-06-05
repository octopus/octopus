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

module poisson_corrections
  use global
  use lib_oct
  use derivatives
  use mesh

  implicit none

  private
  public :: phi, &
            aux, &
            maxl, &
            der_pointer, &
            mesh_pointer, &
            build_phi, &
            build_aux, &
            correct_rho, &
            get_multipoles, &
            op, opt

  type(der_discr_type), pointer :: der_pointer
  type(mesh_type),      pointer :: mesh_pointer

  integer :: maxl
  FLOAT, allocatable :: phi(:, :)
  FLOAT, allocatable :: aux(:, :)

  FLOAT, parameter :: alpha_ = M_FIVE

  contains

  subroutine correct_rho(m, ml, rho, rho_corrected, vh_correction)
    implicit none
    type(mesh_type), intent(in)  :: m
    integer,         intent(in)  :: ml
    FLOAT,           intent(in)  :: rho(:)
    FLOAT,           intent(out) :: rho_corrected(:)
    FLOAT,           intent(out) :: vh_correction(:)

    integer :: i, add_lm, l, mm, lldfac, j
    FLOAT   :: alpha, r2
    FLOAT, allocatable :: mult(:)
    FLOAT, allocatable :: betal(:)


    allocate(mult((ml+1)**2))
    call get_multipoles(m, rho, ml, mult)

    alpha = alpha_*m%h(1)

    allocate(betal((ml+1)**2))
    add_lm = 1
    do l = 0, ml
       do mm = -l, l
          lldfac = 1; do j = 1, 2*l+1, 2; lldfac = lldfac * j; enddo
          betal(add_lm) = (2**(l+2))/( alpha**(2*l+3) * sqrt(M_PI) * lldfac )
          add_lm = add_lm + 1
       enddo
    enddo

    rho_corrected = rho
    vh_correction = M_ZERO
    do i = 1, m%np
       r2 = dot_product(m%x(i, 1:3), m%x(i, 1:3)) ! mesh_r could be used, but it wastes time.
       r2 = exp(-r2/alpha**2)
       add_lm = 1
       do l = 0, ml
          do mm = -l, l
             rho_corrected(i) = rho_corrected(i) - mult(add_lm) * betal(add_lm) * aux(add_lm, i) * r2
             vh_correction(i) = vh_correction(i) + mult(add_lm) * phi(add_lm, i)   
             add_lm = add_lm + 1
          enddo
       enddo
    enddo

    deallocate(mult, betal)
  end subroutine correct_rho

  subroutine get_multipoles(m, rho, ml, mult)
    implicit none
    type(mesh_type), intent(in)  :: m
    FLOAT,           intent(in)  :: rho(:)  ! rho(m%np)
    integer,         intent(in)  :: ml
    FLOAT,           intent(out) :: mult((ml+1)**2)

    integer :: i, add_lm, l, mm

    mult(:) = M_ZERO
    do i = 1, m%np
      add_lm = 1
      do l = 0, ml
         do mm = -l, l
            mult(add_lm) = mult(add_lm) + rho(i)*aux(add_lm, i)*m%vol_pp(i)
            add_lm = add_lm + 1
         enddo
      enddo
    enddo

  end subroutine get_multipoles

  subroutine build_phi(m)
    type(mesh_type), intent(in) :: m

    FLOAT :: alpha, beta, gamma, ylm, r, x(3)
    integer :: i, l, add_lm, lldfac, j, mm

    allocate(phi((maxl+1)**2, m%np))
    alpha = alpha_*m%h(1)
    do i = 1, m%np
       call mesh_r(m, i, r, x = x)
       add_lm = 1
       do l = 0, maxl
          lldfac = 1; do j = 1, 2*l+1, 2; lldfac = lldfac * j; enddo
          beta = (2**(l+2))/( alpha**(2*l+3) * sqrt(M_PI) * lldfac )
          gamma = ( sqrt(M_PI)*2**(l+3) ) / lldfac
          do mm = -l, l
             ylm  = loct_ylm(x(1), x(2), x(3), l, mm)
             if(r .ne. M_ZERO) then
                phi(add_lm, i) = gamma * isubl( l, r/alpha) * ylm / r**(l+1)
             else
                phi(add_lm, i) = gamma * ylm / alpha
             endif
             add_lm = add_lm + 1
          enddo
       enddo
    enddo

    contains

     FLOAT function isubl(l, x)
       integer, intent(in) :: l
       FLOAT, intent(in) :: x
       isubl = M_HALF*loct_gamma(l + M_HALF)*(M_ONE - loct_incomplete_gamma(l+M_HALF, x**2) )
     end function isubl

  end subroutine build_phi

  subroutine build_aux(m)
    type(mesh_type), intent(in) :: m

    FLOAT :: ylm, r, x(3)
    integer :: i, l, add_lm, mm

    allocate(aux((maxl+1)**2, m%np))
    do i = 1, m%np
       call mesh_r(m, i, r, x = x)
       add_lm = 1
       do l = 0, maxl
          do mm = -l, l
             ylm  = loct_ylm(x(1), x(2), x(3), l, mm)
             if(r .ne. M_ZERO) then
                aux(add_lm, i) = r**l*ylm
             else
                if(l==0) then
                   aux(add_lm, i) = ylm
                else
                   aux(add_lm, i) = M_ZERO
                endif
             endif
             add_lm = add_lm + 1
          enddo
       enddo
    enddo

  end subroutine build_aux

  subroutine op(x, y)
    FLOAT, intent(in)  :: x(:)
    FLOAT, intent(out) :: y(:)
    integer :: i
    FLOAT, allocatable :: u(:)
    allocate(u(mesh_pointer%np))
    do i = 1, mesh_pointer%np
       u(i) = x(i)/sqrt(mesh_pointer%vol_pp(i))
    enddo
    call dderivatives_lapl(der_pointer, u, y)
    do i = 1, mesh_pointer%np
       y(i) = y(i)*sqrt(mesh_pointer%vol_pp(i))
    enddo
    deallocate(u)
  end subroutine op
  
  subroutine opt(x, y)
    FLOAT, intent(in)  :: x(:)
    FLOAT, intent(out) :: y(:)
    integer :: i
    FLOAT, allocatable :: u(:)
    allocate(u(mesh_pointer%np))
    do i = 1, mesh_pointer%np
       u(i) = x(i)*sqrt(mesh_pointer%vol_pp(i))
    enddo
    call dderivatives_laplt(der_pointer, u, y)
    do i = 1, mesh_pointer%np
       y(i) = y(i)/sqrt(mesh_pointer%vol_pp(i))
    enddo
    deallocate(u)
  end subroutine opt


end module poisson_corrections
