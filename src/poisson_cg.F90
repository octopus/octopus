!! Copyright (C) 2004 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module poisson_cg
  use mesh
  use functions

  implicit none

  public :: poisson_cg1, &
            poisson_cg2

contains
  

  subroutine poisson_cg1(m, der, pot, rho)
    type(mesh_type),      intent(in)    :: m
    type(der_discr_type), intent(in)    :: der
    FLOAT,                intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                intent(in)    :: rho(:) ! rho(m%np)
  
    integer, parameter :: ml = 4 ! maximum multipole moment to include

    integer :: i, iter
    FLOAT :: s1, s2, s3, ak, ck
    FLOAT, allocatable :: zk(:), tk(:), pk(:), mult(:), wk(:), lwk(:)
    
    call push_sub('poisson_cg')

    allocate(zk(m%np), tk(m%np), pk(m%np)) ! for the cg
    allocate(wk(m%np_tot), lwk(m%np_tot))
    
    ! build initial guess for the potential
    wk(1:m%np) = pot(1:m%np)
    call boundary_conditions(m, rho, ml, wk)
    
    call dderivatives_lapl(der, wk, lwk, .true.)
    
    zk(:) = -M_FOUR*M_PI*rho(:) ! this will be needed later
    zk(:) = zk(:) - lwk(1:m%np)
    deallocate(wk, lwk) ! they are no longer needed
    
    pk = zk
    s1 = dmf_nrm2(m, zk)**2
    
    ! now we start the conjugate gradient cycles
    iter = 0
    do 
      iter = iter + 1
      call dderivatives_lapl(der, pk, tk)
      
      s2 = dmf_dotp(m, zk, tk)
      ak = s1/s2
      pot = pot + ak*pk
      zk = zk - ak*tk
      s3 = dmf_nrm2(m, zk)**2
      
      if(iter >= 400 .or. abs(s3) < CNST(1.0e-6)) exit
      
      ck = s3/s1
      s1 = s3
      pk = zk + ck*pk
    end do
    
    if(iter >= 400) then
      message(1) = "Poisson using conjugated gradients: Not converged!"
      write(message(2),*) "iter = ", iter, " s3 = ", s3
      call write_warning(2)
    endif
    
    deallocate(zk, tk, pk)
    
    call pop_sub()
  end subroutine poisson_cg1

  subroutine poisson_cg2(m, der, pot, rho)
    implicit none
    type(mesh_type),      intent(in)    :: m
    type(der_discr_type), intent(in)    :: der
    FLOAT,                intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                intent(in)    :: rho(:) ! rho(m%np)
  
    integer, parameter :: ml = 4 ! maximum multipole moment to include

    integer :: i, iter
    FLOAT :: s1, s2, s3, ak, ck
    FLOAT, allocatable :: zk(:), tk(:), pk(:), mult(:), wk(:), lwk(:)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)
    
    call push_sub('poisson_cg')

    allocate(zk(m%np), tk(m%np), pk(m%np)) ! for the cg
    allocate(wk(m%np), lwk(m%np))
    allocate(rho_corrected(m%np), vh_correction(m%np))

    call correct_rho(m, ml, rho, rho_corrected, vh_correction)

    ! build initial guess for the potential
    wk(1:m%np) = pot(1:m%np)
    
    call dderivatives_lapl(der, wk, lwk)
    
    zk(:) = -M_FOUR*M_PI*rho_corrected(:) ! this will be needed later
    zk(:) = zk(:) - lwk(1:m%np)
    deallocate(wk, lwk) ! they are no longer needed
    
    pk = zk
    s1 = dmf_nrm2(m, zk)**2
    
    ! now we start the conjugate gradient cycles
    iter = 0
    do 
      iter = iter + 1
      call dderivatives_lapl(der, pk, tk)
      
      s2 = dmf_dotp(m, zk, tk)
      ak = s1/s2
      pot = pot + ak*pk
      zk = zk - ak*tk
      s3 = dmf_nrm2(m, zk)**2
      
      if(iter >= 400 .or. abs(s3) < CNST(1.0e-6)) exit
      
      ck = s3/s1
      s1 = s3
      pk = zk + ck*pk
    end do
    
    if(iter >= 400) then
      message(1) = "Poisson using conjugated gradients: Not converged!"
      write(message(2),*) "iter = ", iter, " s3 = ", s3
      call write_warning(2)
    endif

    pot = pot + vh_correction
    
    deallocate(zk, tk, pk)
    deallocate(rho_corrected, vh_correction)
    
    call pop_sub()
  end subroutine poisson_cg2

  subroutine get_multipoles(m, rho, ml, mult)
    implicit none
    type(mesh_type), intent(in)  :: m
    FLOAT,           intent(in)  :: rho(:)  ! rho(m%np)
    integer,         intent(in)  :: ml
    FLOAT,           intent(out) :: mult((ml+1)**2)

    integer :: i, add_lm, l, mm
    FLOAT   :: x(3), r, s1, sa, rr

    mult(:) = M_ZERO
    do i = 1, m%np
      call mesh_r(m, i, r, x=x)
      add_lm = 1
      rr = M_ONE
      do l = 0, ml
         s1 = rho(i)*rr
         do mm = -l, l
            sa = loct_ylm(x(1), x(2), x(3), l, mm)
            mult(add_lm) = mult(add_lm) + s1*sa*m%vol_pp(i)
            add_lm = add_lm + 1
         enddo
         rr = rr*r
      enddo
    enddo

  end subroutine get_multipoles

  ! calculate multipole moments until ml
  subroutine boundary_conditions(m, rho, ml, pot)
    implicit none
    type(mesh_type), intent(in)  :: m
    FLOAT,           intent(in)  :: rho(:)  ! rho(m%np)
    integer,         intent(in)  :: ml
    FLOAT,           intent(out) :: pot(:)  ! pot(m%np_tot)

    integer :: i, add_lm, l, mm
    FLOAT   :: x(3), r, s1, sa
    FLOAT, allocatable :: mult(:)
    
    allocate(mult((ml+1)**2))
    call get_multipoles(m, rho, ml, mult)
    do i = m%np+1, m%np_tot ! boundary conditions
      call mesh_r(m, i, r, x=x)
      add_lm = 1
      do l = 0, ml
        s1 = M_FOUR*M_PI/((M_TWO*l + M_ONE)*r**(l + 1))
        do mm = -l, l
          sa = loct_ylm(x(1), x(2), x(3), l, mm)
          pot(i) = pot(i) + sa * mult(add_lm) * s1
          add_lm = add_lm+1
        end do
      end do
    end do

    deallocate(mult)
  end subroutine boundary_conditions

  subroutine correct_rho(m, ml, rho, rho_corrected, vh_correction)
    implicit none
    type(mesh_type), intent(in)  :: m
    integer,         intent(in)  :: ml
    FLOAT,           intent(in)  :: rho(:)
    FLOAT,           intent(out) :: rho_corrected(:)
    FLOAT,           intent(out) :: vh_correction(:)

    integer :: i, add_lm, l, mm, lldfac, j
    FLOAT   :: x(3), r, s1, sa, ylm, alpha, beta, gamma
    FLOAT, allocatable :: mult(:)

    allocate(mult((ml+1)**2))
    call get_multipoles(m, rho, ml, mult)

    rho_corrected = rho
    vh_correction = M_ZERO
    alpha = M_TWO
    do i = 1, m%np
       call mesh_r(m, i, r, x = x)
       add_lm = 1
       do l = 0, ml
          lldfac = 1; do j = 1, 2*l+1, 2; lldfac = lldfac * j; enddo
          beta = (2**(l+2))/( alpha**(2*l+3) * sqrt(M_PI) * lldfac )
          gamma = ( sqrt(M_PI)*2**(l+3) ) / lldfac
          do mm = -l, l
             ylm  = loct_ylm(x(1), x(2), x(3), l, mm)
             rho_corrected(i) = rho_corrected(i) - mult(add_lm)*beta * r**l * exp(-(r/alpha)**2)*ylm
             if(r .ne. M_ZERO) then
                vh_correction(i) = vh_correction(i) + mult(add_lm)*gamma * isubl( l, r/alpha) * ylm / r**(l+1)
             else
                vh_correction(i) = vh_correction(i) + mult(add_lm)*gamma * ylm / alpha
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

  end subroutine correct_rho


end module poisson_cg
