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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: solver_1d_inc.F90 9854 2013-01-19 23:28:12Z dstrubbe $

subroutine X(poisson1D_solve_direct)(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  R_TYPE,          intent(out) :: pot(:)
  R_TYPE,          intent(in)  :: rho(:)
  !! Using theta is not needed for a normal Poisson solver due to linearity,
  !! which makes it possible to solve separately for real/imaginary parts.
  !! But the soft-Coulomb kernel makes it necessary.
  !! Note that we don`t divide by e^(i theta) here, we do it "outside" as with the other Poisson solvers.

  integer             :: ip, jp
  integer, allocatable :: ip_v(:), part_v(:)
  FLOAT               :: xx, yy
  R_TYPE              :: soft_coulomb_param_squared
#ifdef HAVE_MPI
  R_TYPE              :: tmp
  FLOAT               :: xg(1:MAX_DIM)
  R_TYPE, allocatable :: pvec(:)
#endif

  PUSH_SUB(X(poisson1D_solve_direct))

  ASSERT(this%method == POISSON_DIRECT_SUM)

  soft_coulomb_param_squared = this%poisson_soft_coulomb_param**2 * exp(-M_TWO * M_zI * this%theta)
  ! This will discard imaginary part for R_TYPE real.
  ! But theta won`t be there unless we already use complex scaling, so only cmplx is relevant.

#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))
    SAFE_ALLOCATE(part_v(1:this%der%mesh%np_global))
    SAFE_ALLOCATE(ip_v(1:this%der%mesh%np_global))
    do ip = 1, this%der%mesh%np_global
      ip_v(ip) = ip
    end do
    call partition_get_partition_number(this%der%mesh%inner_partition, this%der%mesh%np_global, ip_v, part_v)

    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx = xg(1)
      do jp = 1, this%der%mesh%np
        yy = this%der%mesh%x(jp, 1)
        pvec(jp) = rho(jp)/sqrt(soft_coulomb_param_squared + (xx - yy)**2)
      end do
      tmp = X(mf_integrate)(this%der%mesh, pvec) 
      ip_v(1) = ip
      call partition_get_partition_number(this%der%mesh%inner_partition, 1, ip_v, part_v)

      if (part_v(1) == this%der%mesh%vp%partno) then
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
        if(this%der%mesh%use_curvilinear) then
          pot(ip) = pot(ip) + rho(jp)/sqrt(soft_coulomb_param_squared + (xx - yy)**2)*this%der%mesh%vol_pp(jp)
        else
          pot(ip) = pot(ip) + rho(jp)/sqrt(soft_coulomb_param_squared + (xx - yy)**2)
        endif
      end do
    end do
    if(.not. this%der%mesh%use_curvilinear) pot(:) = pot(:) * this%der%mesh%volume_element
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(X(poisson1D_solve_direct))
end subroutine X(poisson1D_solve_direct)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
