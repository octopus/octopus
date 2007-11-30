!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! ---------------------------------------------------------
subroutine poisson2D_init(gr)
  type(grid_t), intent(inout) :: gr

  call push_sub('poisson2D.poisson2D_init')


  ASSERT(poisson_solver == FFT_SPH .or. poisson_solver == DIRECT_SUM_2D)
  if (poisson_solver == FFT_SPH) call poisson_fft_build_2d(gr, poisson_solver)

  call pop_sub()

end subroutine poisson2D_init


! ---------------------------------------------------------
subroutine poisson2D_solve(m, pot, rho)
  type(mesh_t), intent(in) :: m
  FLOAT, intent(out)       :: pot(m%np)
  FLOAT, intent(in)        :: rho(m%np)

  integer  :: i, j
  FLOAT    :: x(2), y(2), tmp
  FLOAT, allocatable :: pvec(:)

  ASSERT(poisson_solver == -2)

  call push_sub('poisson2D.poisson2D_solve')

  if(m%parallel_in_domains) then
    ALLOCATE(pvec(m%np), m%np)

    pot = M_ZERO
    do i = 1, m%np_global
      x(:) = m%x_global(i,1:2)
      do j = 1, m%np
        if(m%vp%global(i, m%vp%partno) == j) then
          pvec(j) = M_TWO*sqrt(M_PI)*rho(j)/m%h(1)
        else
          y(:) = m%x(j,1:2)
          pvec(j) = rho(j)/sqrt(sum((x-y)**2))
        end if
      end do
      tmp = dmf_integrate(m, pvec)
      if (m%vp%part(i).eq.m%vp%partno) then
        pot(m%vp%global(i, m%vp%partno)) = tmp
      end if
    end do

    deallocate(pvec)

  else ! serial mode

    pot = M_ZERO
    do i = 1, m%np
      x(:) = m%x(i,1:2)
      do j = 1, m%np
        if(i == j) then
          pot(i) = pot(i) + M_TWO*sqrt(M_PI)*rho(i)/m%h(1)*m%vol_pp(j)
        else
          y(:) = m%x(j,1:2)
          pot(i) = pot(i) + rho(j)/sqrt(sum((x-y)**2))*m%vol_pp(j)
        end if
      end do
    end do

  end if

  call pop_sub()
end subroutine poisson2D_solve

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
