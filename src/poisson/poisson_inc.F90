!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio,
!! G. Bertsch, M. Oliveira
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
!! $Id$

!------------------------------------------------------------------

subroutine X(poisson_solve_start)(this, rho)
  type(poisson_t),      intent(in) :: this
  R_TYPE,               intent(in) :: rho(:)

#ifdef HAVE_MPI2    
  integer :: islave
  type(profile_t), save :: prof

  PUSH_SUB(X(poisson_solve_start))
  call profiling_in(prof, "POISSON_START")

  ! we assume all nodes have a copy of the density
  do islave = this%local_grp%rank, this%nslaves - 1, this%local_grp%size !all nodes are used for communication
    call MPI_Send(rho(1), this%der%mesh%np, R_MPITYPE, islave, CMD_POISSON_SOLVE, this%intercomm, mpi_err)
  end do
  
  call profiling_out(prof)
  POP_SUB(X(poisson_solve_start))
#endif
    
end subroutine X(poisson_solve_start)
  
!----------------------------------------------------------------

subroutine X(poisson_solve_finish)(this, pot)
  type(poisson_t),  intent(in)    :: this
  R_TYPE,           intent(inout) :: pot(:)

#ifdef HAVE_MPI2
  type(profile_t), save :: prof

  PUSH_SUB(X(poisson_solve_finish))
  call profiling_in(prof, "POISSON_FINISH")

  call MPI_Bcast(pot(1), this%der%mesh%np, R_MPITYPE, 0, this%intercomm, mpi_err)

  call profiling_out(prof)
  POP_SUB(X(poisson_solve_finish))
#endif

end subroutine X(poisson_solve_finish)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
