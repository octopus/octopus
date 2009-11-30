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
subroutine poisson2D_init(mesh)
  type(mesh_t), intent(inout) :: mesh

  call push_sub('poisson2D.poisson2D_init')

  select case(poisson_solver)
  case(POISSON_FFT_SPH)
    call poisson_fft_build_2d_0d(mesh)
  case(POISSON_FFT_CYL)
    call poisson_fft_build_2d_1d(mesh)
  case(POISSON_FFT_NOCUT)
    call poisson_fft_build_2d_2d(mesh)
  end select

  call pop_sub()

end subroutine poisson2D_init

! ---------------------------------------------------------

subroutine poisson2D_solve(mesh, pot, rho)
  type(mesh_t), intent(in)  :: mesh
  FLOAT,        intent(out) :: pot(:)
  FLOAT,        intent(in)  :: rho(:)

  integer  :: i, j
  FLOAT    :: x(2), y(2)
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(poisson_solver == -2)

  call push_sub('poisson2D.poisson2D_solve')
#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:mesh%np))

    pot = M_ZERO
    do i = 1, mesh%np_global
      xg = mesh_x_global(m, i)
      x(1:2) = xg(1:2)
      do j = 1, mesh%np
        if(vec_global2local(mesh%vp, i, mesh%vp%partno) == j) then
          pvec(j) = M_TWO*sqrt(M_PI)*rho(j)/mesh%h(1)
        else
          y(:) = mesh%x(j,1:2)
          pvec(j) = rho(j)/sqrt(sum((x-y)**2))
        end if
      end do
      tmp = dmf_integrate(m, pvec)
      if (mesh%vp%part(i).eq.mesh%vp%partno) then
        pot(vec_global2local(mesh%vp, i, mesh%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else ! serial mode
#endif
    pot = M_ZERO
    do i = 1, mesh%np
      x(:) = mesh%x(i,1:2)
      do j = 1, mesh%np
        if(i == j) then
          pot(i) = pot(i) + M_TWO*sqrt(M_PI)*rho(i)/mesh%h(1)*mesh%vol_pp(j)
        else
          y(:) = mesh%x(j,1:2)
          pot(i) = pot(i) + rho(j)/sqrt(sum((x-y)**2))*mesh%vol_pp(j)
        end if
      end do
    end do
#ifdef HAVE_MPI
  end if
#endif

  call pop_sub()
end subroutine poisson2D_solve

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
