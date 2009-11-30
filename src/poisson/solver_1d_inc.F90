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
subroutine poisson1d_init(mesh)
  type(mesh_t), intent(inout) :: mesh

  call push_sub('solver_1D.poisson1d_init')

  call parse_float(datasets_check('Poisson1DSoftCoulombParam'), &
    M_ONE, poisson_soft_coulomb_param)

  select case(poisson_solver)
  case(POISSON_FFT_SPH)
    call poisson_fft_build_1d_0d(mesh, poisson_soft_coulomb_param)
  case(POISSON_FFT_NOCUT)
    call poisson_fft_build_1d_1d(mesh, poisson_soft_coulomb_param)
  end select

  call pop_sub()
end subroutine poisson1d_init

!-----------------------------------------------------------------

subroutine poisson1D_solve(mesh, pot, rho)
  type(mesh_t), intent(in)  :: mesh
  FLOAT,        intent(out) :: pot(:)
  FLOAT,        intent(in)  :: rho(:)

  integer  :: ip, jp
  FLOAT    :: xx, yy
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(poisson_solver == -1)

  call push_sub('poisson1D.poisson1D_solve')

#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:mesh%np))

    pot = M_ZERO
    do ip = 1, mesh%np_global
      xg = mesh_x_global(mesh, ip)
      xx = xg(1)
      do jp = 1, mesh%np
        yy = mesh%x(jp, 1)
        pvec(jp) = rho(jp)/sqrt(poisson_soft_coulomb_param**2 + (xx-yy)**2)
      end do
      tmp = dmf_integrate(mesh, pvec)
      if (mesh%vp%part(ip).eq.mesh%vp%partno) then
        pot(vec_global2local(mesh%vp, ip, mesh%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else  ! running in serial
#endif
    pot = M_ZERO
    do ip = 1, mesh%np
      xx = mesh%x(ip, 1)
      do jp = 1, mesh%np
        yy = mesh%x(jp, 1)
        pot(ip) = pot(ip) + rho(jp)/sqrt(poisson_soft_coulomb_param**2 + (xx-yy)**2)*mesh%vol_pp(jp)
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
