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

! ---------------------------------------------------------
subroutine poisson2D_init(this)
  type(poisson_t), intent(inout) :: this

  PUSH_SUB(poisson2D_init)

  if(this%method == POISSON_FFT) then
    select case(this%kernel)
    case(POISSON_FFT_KERNEL_SPH)
      call poisson_fft_build_2d_0d(this%der%mesh, this%cube)
    case(POISSON_FFT_KERNEL_CYL)
      call poisson_fft_build_2d_1d(this%der%mesh, this%cube)
    case(POISSON_FFT_KERNEL_NOCUT)
      call poisson_fft_build_2d_2d(this%der%mesh, this%cube)
    end select
  end if

  POP_SUB(poisson2D_init)

end subroutine poisson2D_init

! ---------------------------------------------------------

subroutine poisson2D_solve(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  integer  :: ip, jp
  FLOAT    :: xx(2), yy(2)
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:)
#endif

  ASSERT(this%method == -2)

  PUSH_SUB(poisson2D_solve)
#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))

    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx(1:2) = xg(1:2)
      do jp = 1, this%der%mesh%np
        if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp) then
          pvec(jp) = M_TWO*sqrt(M_PI)*rho(jp)/this%der%mesh%spacing(1)
        else
          yy(:) = this%der%mesh%x(jp,1:2)
          pvec(jp) = rho(jp)/sqrt(sum((xx-yy)**2))
        end if
      end do
      tmp = dmf_integrate(this%der%mesh, pvec)
      if (this%der%mesh%vp%part(ip).eq.this%der%mesh%vp%partno) then
        pot(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno)) = tmp
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)

  else ! serial mode
#endif
    pot = M_ZERO
    if(this%der%mesh%use_curvilinear) then
      do ip = 1, this%der%mesh%np
        xx(:) = this%der%mesh%x(ip, 1:2)
        do jp = 1, this%der%mesh%np
          if(ip == jp) then
            pot(ip) = pot(ip) + M_TWO*sqrt(M_PI)*rho(ip)/this%der%mesh%spacing(1)*this%der%mesh%vol_pp(jp)
          else
            yy(:) = this%der%mesh%x(jp, 1:2)
            pot(ip) = pot(ip) + rho(jp)/sqrt(sum((xx-yy)**2))*this%der%mesh%vol_pp(jp)
          end if
        end do
      end do
    else
      do ip = 1, this%der%mesh%np
        xx(:) = this%der%mesh%x(ip, 1:2)
        do jp = 1, this%der%mesh%np
          if(ip == jp) then
            pot(ip) = pot(ip) + M_TWO*sqrt(M_PI)*rho(ip)/this%der%mesh%spacing(1)*this%der%mesh%vol_pp(1)
          else
            yy(:) = this%der%mesh%x(jp, 1:2)
            pot(ip) = pot(ip) + rho(jp)/sqrt(sum((xx-yy)**2))*this%der%mesh%vol_pp(1)
          end if
        end do
      end do
    end if
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(poisson2D_solve)
end subroutine poisson2D_solve

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
