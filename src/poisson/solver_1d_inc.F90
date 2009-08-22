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


!-----------------------------------------------------------------
subroutine poisson1d_init(gr)
  type(grid_t), intent(inout) :: gr
  call push_sub('solver_1D.poisson1d_init')

  call loct_parse_float(datasets_check('Poisson1DSoftCoulombParam'), &
    M_ONE, poisson_soft_coulomb_param)

  select case(poisson_solver)
  case(POISSON_FFT_SPH)
    call poisson_fft_build_1d_0d(gr, poisson_soft_coulomb_param)
  case(POISSON_FFT_NOCUT)
    call poisson_fft_build_1d_1d(gr, poisson_soft_coulomb_param)
  end select

  call pop_sub()
end subroutine poisson1d_init
!-----------------------------------------------------------------


!-----------------------------------------------------------------
subroutine poisson1D_solve(m, pot, rho)
  type(mesh_t), intent(in) :: m
  FLOAT, intent(out) :: pot(m%np)
  FLOAT, intent(in)  :: rho(m%np)

  integer  :: i, j
  FLOAT    :: x, y
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(poisson_solver == -1)

  call push_sub('poisson1D.poisson1D_solve')

#ifdef HAVE_MPI
  if(m%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:m%np))

    pot = M_ZERO
    do i = 1, m%np_global
      xg = mesh_x_global(m, i)
      x = xg(1)
      do j = 1, m%np
        y = m%x(j, 1)
        pvec(j) = rho(j)/sqrt(poisson_soft_coulomb_param**2 + (x-y)**2)
      end do
      tmp = dmf_integrate(m, pvec)
      if (m%vp%part(i).eq.m%vp%partno) then
        pot(vec_global2local(m%vp, i, m%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else  ! running in serial
#endif
    pot = M_ZERO
    do i = 1, m%np
      x = m%x(i, 1)
      do j = 1, m%np
        y = m%x(j, 1)
        pot(i) = pot(i) + rho(j)/sqrt(poisson_soft_coulomb_param**2 + (x-y)**2) * m%vol_pp(j)
      end do
    end do
#ifdef HAVE_MPI
  end if
#endif

  call pop_sub()
end subroutine poisson1D_solve
!-----------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
