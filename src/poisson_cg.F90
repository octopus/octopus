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
!!
!! $Id$

#include "global.h"

module poisson_cg
  use global
  use math
  use mesh
  use derivatives
  use functions
  use poisson_corrections

  implicit none

  public :: poisson_cg1_init, &
            poisson_cg1,      &
            poisson_cg1_end,  &
            poisson_cg2_init, &
            poisson_cg2,      &
            poisson_cg2_end

  FLOAT :: threshold

contains

  subroutine poisson_cg1_init(m, ml, thr)
    type(mesh_type), intent(in) :: m
    integer, intent(in) :: ml
    FLOAT, intent(in)   :: thr
    call push_sub('poisson_cg1_init')
    maxl = ml
    threshold = thr
    call build_aux(m)
    call pop_sub()
  end subroutine poisson_cg1_init

  subroutine poisson_cg1_end
    deallocate(aux)
  end subroutine poisson_cg1_end

  subroutine poisson_cg2_init(m, ml, thr)
    type(mesh_type), intent(in) :: m
    integer, intent(in) :: ml
    FLOAT, intent(in) :: thr
    call push_sub('poisson_cg2_init')
    maxl = ml
    threshold = thr
    call build_aux(m)
    call build_phi(m)
    call pop_sub()
  end subroutine poisson_cg2_init

  subroutine poisson_cg2_end
    deallocate(aux, phi)
  end subroutine poisson_cg2_end

  subroutine poisson_cg1(m, der, pot, rho)
    type(mesh_type), target, intent(in)         :: m
    type(der_discr_type), target, intent(in)    :: der
    FLOAT,                intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                intent(in)    :: rho(:) ! rho(m%np)
  
    integer :: iter
    FLOAT :: res
    FLOAT, allocatable :: wk(:), lwk(:), zk(:), pk(:)

    call push_sub('poisson_cg')

    allocate(wk(m%np_tot), lwk(m%np_tot), zk(m%np), pk(m%np))
    
    ! build initial guess for the potential
    wk(1:m%np) = pot(1:m%np)
    call boundary_conditions(m, rho, maxl, wk)
    call dderivatives_lapl(der, wk, lwk, .true.)
    
    zk(1:m%np) = -M_FOUR*M_PI*rho(1:m%np) - lwk(1:m%np)
    deallocate(wk, lwk) ! they are no longer needed

    der_pointer  => der
    mesh_pointer => m
    pk = zk
    iter = 400
    call dconjugate_gradients(m%np, pk, zk, op, iter, res, threshold)
    if(res >= threshold) then
       message(1) = 'Conjugate gradients Poisson solver did not converge.'
       write(message(2), '(a,i8)')    '  Iter = ',iter
       write(message(3), '(a,e14.6)') '  Res = ', res
       call write_warning(3)
    endif 
    nullify(der_pointer, mesh_pointer)
    pot = pot + pk
    
    deallocate(zk, pk)
    call pop_sub()
  end subroutine poisson_cg1

  subroutine poisson_cg2(m, der, pot, rho)
    implicit none
    type(mesh_type), target, intent(in)    :: m
    type(der_discr_type), target, intent(in)    :: der
    FLOAT,                intent(inout) :: pot(:) ! pot(m%np)
    FLOAT,                intent(in)    :: rho(:) ! rho(m%np)

    integer :: i, iter
    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)
    FLOAT :: res

    call push_sub('poisson_cg')

    der_pointer  => der
    mesh_pointer => m

    allocate(rho_corrected(m%np), vh_correction(m%np))
    call correct_rho(m, maxl, rho, rho_corrected, vh_correction)
    rho_corrected = - M_FOUR*M_PI*rho_corrected
    pot = pot - vh_correction
    ! This is the change of base, which for the moment is not done.
    do i = 1, m%np
       pot(i) = pot(i)*sqrt(m%vol_pp(i))
       rho_corrected(i) = rho_corrected(i)*sqrt(m%vol_pp(i))
    enddo
    iter = 400
    call dconjugate_gradients(m%np, pot, rho_corrected, op, opt, &
                             iter, res, threshold)
    if(res >= threshold) then
       message(1) = 'Conjugate gradients Poisson solver did not converge.'
       write(message(2), '(a,i8)')    '  Iter = ',iter
       write(message(3), '(a,e14.6)') '  Res = ', res
       call write_warning(3)
    endif 
    do i = 1, m%np
       pot(i) = pot(i)/sqrt(m%vol_pp(i))
    enddo
    pot = pot + vh_correction

    nullify(der_pointer, mesh_pointer)
    deallocate(rho_corrected, vh_correction)
    call pop_sub()
  end subroutine poisson_cg2

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

end module poisson_cg
