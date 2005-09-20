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

subroutine poisson1D_solve(m, pot, rho)
  type(mesh_type), intent(in) :: m
  FLOAT, intent(out) :: pot(m%np)
  FLOAT, intent(in)  :: rho(m%np)

  integer  :: i, j
  FLOAT    :: x, y
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  FLOAT    :: tmp
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(poisson_solver == -1)

  call push_sub('poisson1D.poisson1D_solve')

#if defined(HAVE_MPI) && defined(HAVE_METIS)

  allocate(pvec(1:m%np))

  pot = M_ZERO
  do i = 1, m%np_global
     x = m%x_global(i, 1)
     do j = 1, m%np
        y = m%x(j, 1)
        pvec(j) = rho(j)/sqrt(M_ONE + (x-y)**2)
     end do
     tmp = dmf_integrate(m, pvec)
     if (m%part(i).eq.m%vp%partno) then
        pot(m%vp%global(i, m%vp%partno)) = tmp
     endif
  end do

  deallocate(pvec)

#else

  pot = M_ZERO
  do i = 1, m%np
     x = m%x(i, 1)
     do j=1, m%np
        y = m%x(j, 1)
        pot(i) = pot(i) + rho(j)/sqrt(M_ONE + (x-y)**2) * m%vol_pp(j)
     end do
     pot(i) = pot(i)
  end do

#endif

  call pop_sub()
end subroutine poisson1D_solve
