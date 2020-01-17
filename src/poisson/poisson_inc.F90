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


 !> Calculates the Poisson equation.
 !! Given the density returns the corresponding potential.
 !!
 !! Different solvers are available that can be chosen in the input file
 !! with the "PoissonSolver" parameter
subroutine X(poisson_solve_sm)(this, sm, pot, rho, all_nodes)
  type(poisson_t),      intent(in)    :: this
  type(submesh_t),      intent(in)    :: sm
  R_TYPE,               intent(inout) :: pot(:) !< Local size of the \b potential vector.
  R_TYPE,               intent(inout) :: rho(:) !< Local size of the \b density (rho) vector.
  !> Is the Poisson solver allowed to utilise
  !! all nodes or only the domain nodes for
  !! its calculations? (Defaults to .true.)
  logical, optional,    intent(in)    :: all_nodes
  type(derivatives_t), pointer :: der

  logical               :: all_nodes_value
  type(profile_t), save :: prof

  call profiling_in(prof, 'POISSON_SOLVE_SM')
  PUSH_SUB(X(poisson_solve_sm))

  der => this%der

  ASSERT(ubound(pot, dim = 1) == sm%np_part .or. ubound(pot, dim = 1) == sm%np)
  ASSERT(ubound(rho, dim = 1) == sm%np_part .or. ubound(rho, dim = 1) == sm%np)

  ! Check optional argument and set to default if necessary.
  if(present(all_nodes)) then
    all_nodes_value = all_nodes
  else
    all_nodes_value = this%all_nodes_default
  end if

  ASSERT(this%method /= POISSON_NULL)
  ASSERT(this%der%mesh%sb%dim /= 1)
  call X(poisson_solve_direct_sm)(this, sm, pot, rho)

  POP_SUB(X(poisson_solve_sm))
  call profiling_out(prof)
end subroutine X(poisson_solve_sm)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
