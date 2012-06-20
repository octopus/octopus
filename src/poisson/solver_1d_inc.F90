!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
subroutine poisson1d_init(this)
  type(poisson_t), intent(inout) :: this

  PUSH_SUB(poisson1d_init)

  call parse_float(datasets_check('Poisson1DSoftCoulombParam'), &
    M_ONE, this%poisson_soft_coulomb_param)

  if(this%method == POISSON_FFT) then
    call poisson_fft_init(this%fft_solver, this%der%mesh, this%cube, this%kernel, &
      soft_coulb_param = this%poisson_soft_coulomb_param)
  end if

  POP_SUB(poisson1d_init)
end subroutine poisson1d_init

!-----------------------------------------------------------------

subroutine poisson1D_solve(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  integer  :: ip, jp
  FLOAT    :: xx, yy
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(this%method == -1)

  PUSH_SUB(poisson1D_solve)

#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))

    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx = xg(1)
      do jp = 1, this%der%mesh%np
        yy = this%der%mesh%x(jp, 1)
        pvec(jp) = rho(jp)/sqrt(this%poisson_soft_coulomb_param**2 + (xx-yy)**2)
      end do
      tmp = dmf_integrate(this%der%mesh, pvec)
      if (this%der%mesh%vp%part(ip).eq.this%der%mesh%vp%partno) then
        pot(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else  ! running in serial
#endif
    pot = M_ZERO
    do ip = 1, this%der%mesh%np
      xx = this%der%mesh%x(ip, 1)
      do jp = 1, this%der%mesh%np
        yy = this%der%mesh%x(jp, 1)
        pot(ip) = pot(ip) + rho(jp)/sqrt(this%poisson_soft_coulomb_param**2 + &
	         (xx-yy)**2)*this%der%mesh%vol_pp(1)
      end do
    end do
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(poisson1D_solve)
end subroutine poisson1D_solve
!-----------------------------------------------------------------

!
! Complex scaled soft Coulomb Hartree Solver
!
subroutine zpoisson1D_solve(this, pot, rho, theta)
  type(poisson_t), intent(in)  :: this
  CMPLX,           intent(out) :: pot(:)
  CMPLX,           intent(in)  :: rho(:)
  FLOAT,           intent(in)  :: theta !< complex scaling angle

  integer  :: ip, jp
  CMPLX    :: xx, yy
#ifdef HAVE_MPI
  CMPLX    :: tmp, xg(1:MAX_DIM)
  CMPLX, allocatable :: pvec(:)
#endif

  ASSERT(this%method == -1)

  PUSH_SUB(zpoisson1D_solve)

#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))

    pot = M_z0
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx = xg(1)
      do jp = 1, this%der%mesh%np
        yy = this%der%mesh%x(jp, 1)
        pvec(jp) = rho(jp)/sqrt(this%poisson_soft_coulomb_param**2 +&
         (xx-yy)**2 * exp(M_zI*M_TWO*theta))
      end do
      tmp = dmf_integrate(this%der%mesh, pvec)
      if (this%der%mesh%vp%part(ip).eq.this%der%mesh%vp%partno) then
        pot(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else  ! running in serial
#endif
    pot = M_z0
    do ip = 1, this%der%mesh%np
      xx = this%der%mesh%x(ip, 1)
      do jp = 1, this%der%mesh%np
        yy = this%der%mesh%x(jp, 1)
        pot(ip) = pot(ip) + rho(jp)/sqrt(this%poisson_soft_coulomb_param**2 + &
	         (xx-yy)**2 * exp(M_zI*M_TWO*theta))*this%der%mesh%vol_pp(1)
      end do
    end do
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(zpoisson1D_solve)
end subroutine zpoisson1D_solve
!-----------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
