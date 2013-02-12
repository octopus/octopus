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
!> This routine calculates the SIC exchange functional.
subroutine X(oep_sic) (xcs, gr, st, is, oep, ex, ec)
  type(xc_t),     intent(in)    :: xcs
  type(grid_t),   intent(inout) :: gr
  type(states_t), intent(inout) :: st
  integer,        intent(in)    :: is
  type(xc_oep_t), intent(inout) :: oep
  FLOAT,          intent(inout) :: ex, ec

  integer  :: ist
  FLOAT :: ex2, ec2, ex_, ec_, edummy
  FLOAT, allocatable :: vxc(:, :), rho(:,:)
  R_TYPE, allocatable :: psi(:, :)

  call profiling_in(C_PROFILING_XC_SIC)
  PUSH_SUB(X(oep_sic))

  ASSERT(st%d%dim == 1)

  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(rho(1:gr%mesh%np, 1:2))
  SAFE_ALLOCATE(Vxc(1:gr%mesh%np, 1:2))
  rho(1:gr%mesh%np, 2) = M_ZERO

  ! loop over states
  ex_ = M_ZERO
  ec_ = M_ZERO
  do ist = st%st_start, st%st_end
    if(st%occ(ist, is) .gt. small) then ! we only need the occupied states

      call states_get_state(st, gr%mesh, ist, is, psi)

      ! get orbital density
      rho(1:gr%mesh%np, 1) = oep%socc*st%occ(ist, is)*R_ABS(psi(1:gr%mesh%np, 1))**2

      ! initialize before calling get_vxc
      vxc = M_ZERO
      ex2  = M_ZERO
      ec2  = M_ZERO

      ! calculate LDA/GGA contribution to the SIC (does not work for LB94)
      edummy = M_ZERO
      call xc_get_vxc(gr%fine%der, xcs, st, rho, SPIN_POLARIZED, edummy, edummy, vxc, ex=ex2, ec=ec2)

      ex_ = ex_ - oep%sfact*ex2
      ec_ = ec_ - oep%sfact*ec2

      oep%X(lxc)(1:gr%mesh%np, ist, is) = oep%X(lxc)(1:gr%mesh%np, ist, is) - vxc(1:gr%mesh%np, 1)*R_CONJ(psi(1:gr%mesh%np, 1))

      ! calculate the Hartree contribution using Poisson equation
      vxc(1:gr%mesh%np, 1) = M_ZERO
      call dpoisson_solve(psolver, vxc(:, 1), rho(:, 1), all_nodes=.false.)

      ! The exchange energy.
      ex_ = ex_ - M_HALF*oep%sfact*oep%socc*st%occ(ist, is)* &
        dmf_dotp(gr%mesh, vxc(1:gr%mesh%np, 1), R_ABS(psi(1:gr%mesh%np, 1))**2)

      oep%X(lxc)(1:gr%mesh%np, ist, is) = oep%X(lxc)(1:gr%mesh%np, ist, is) - &
        vxc(1:gr%mesh%np, 1)*R_CONJ(psi(1:gr%mesh%np, 1))
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(ec_, edummy, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    ec_ = edummy
    call MPI_Allreduce(ex_, edummy, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    ex_ = edummy
  end if
#endif

  ec = ec + ec_
  ex = ex + ex_

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(Vxc)
  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(oep_sic))
  call profiling_out(C_PROFILING_XC_SIC)
end subroutine X(oep_sic)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
