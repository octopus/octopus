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

  public :: poisson_cg1

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


  ! calculate multipole moments until ml
  subroutine boundary_conditions(m, rho, ml, pot)
    type(mesh_type), intent(in)  :: m
    FLOAT,           intent(in)  :: rho(:)  ! rho(m%np)
    integer,         intent(in)  :: ml
    FLOAT,           intent(out) :: pot(:)  ! pot(m%np_tot)

    integer :: i, add_lm, l, mm
    FLOAT   :: x(3), r, s1, sa
    FLOAT, allocatable :: mult(:)
    
    allocate(mult(2*ml+1))

    mult(:) = M_ZERO
    do i = 1, m%np
      call mesh_r(m, i, r, x=x)
      
      add_lm = 1
      do l = 0, ml
        if(l.ne.0) then 
          s1 = rho(i)*m%vol_pp(i)*r**l
        else
          s1 = rho(i)*m%vol_pp(i)
        end if
        
        do mm = -l, l
          sa = loct_ylm(x(1), x(2), x(3), l, mm)
          mult(add_lm) = mult(add_lm) + s1*sa
          add_lm = add_lm+1
        end do
      end do
    end do

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

end module poisson_cg
