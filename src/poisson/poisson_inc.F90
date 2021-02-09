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
  call profiling_in(prof, TOSTRING(X(POISSON_START)))

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
  call profiling_in(prof, TOSTRING(X(POISSON_FINISH)))

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
#ifdef R_TCOMPLEX
  FLOAT, allocatable :: aux1(:), aux2(:)
#endif

  call profiling_in(prof, TOSTRING(X(POISSON_SOLVE_SM)))
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
  ASSERT(this%der%dim /= 1)

  select case(this%method)
  case(POISSON_DIRECT_SUM)
    call X(poisson_solve_direct_sm)(this, sm, pot, rho)
  case(POISSON_ISF)
#ifdef R_TREAL
    call poisson_isf_solve(this%isf_solver, sm%mesh, this%cube, pot, rho, all_nodes_value, sm) 
#else
    SAFE_ALLOCATE(aux1(1:sm%np))
    SAFE_ALLOCATE(aux2(1:sm%np))
    ! first the real part
    aux1(1:sm%np) = real(rho(1:sm%np))
    aux2(1:sm%np) = real(pot(1:sm%np))
    call poisson_isf_solve(this%isf_solver, der%mesh, this%cube, aux2, aux1, all_nodes_value, sm)
    pot(1:sm%np)  = aux2(1:sm%np)

    ! now the imaginary part
    aux1(1:sm%np) = aimag(rho(1:sm%np))
    aux2(1:sm%np) = aimag(pot(1:sm%np))
    call poisson_isf_solve(this%isf_solver, der%mesh, this%cube, aux2, aux1, all_nodes_value, sm)
    pot(1:sm%np) = pot(1:sm%np) + M_zI*aux2(1:sm%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)
#endif
  case(POISSON_PSOLVER)
#ifdef R_TREAL
    if (this%psolver_solver%datacode == "G") then
      ! Global version
      call poisson_psolver_global_solve(this%psolver_solver, sm%mesh, this%cube, pot, rho, sm)
    else ! "D" Distributed version
      ASSERT(.false.)
      call poisson_psolver_parallel_solve(this%psolver_solver, sm%mesh, this%cube, pot, rho, this%mesh_cube_map)
    end if
#else
    SAFE_ALLOCATE(aux1(1:sm%np))
    SAFE_ALLOCATE(aux2(1:sm%np))
    ! first the real part
    aux1(1:sm%np) = real(rho(1:sm%np))
    aux2(1:sm%np) = real(pot(1:sm%np))
    if (this%psolver_solver%datacode == "G") then
      ! Global version
      call poisson_psolver_global_solve(this%psolver_solver, sm%mesh, this%cube, aux2, aux1, sm)
    else ! "D" Distributed version
      ASSERT(.false.)
      call poisson_psolver_parallel_solve(this%psolver_solver, sm%mesh, this%cube, aux2, aux1, this%mesh_cube_map)
    end if

    pot(1:sm%np)  = aux2(1:sm%np)

    ! now the imaginary part
    aux1(1:sm%np) = aimag(rho(1:sm%np))
    aux2(1:sm%np) = aimag(pot(1:sm%np))
    if (this%psolver_solver%datacode == "G") then
      ! Global version
      call poisson_psolver_global_solve(this%psolver_solver, sm%mesh, this%cube, aux2, aux1, sm)
    else ! "D" Distributed version
      ASSERT(.false.)
      call poisson_psolver_parallel_solve(this%psolver_solver, sm%mesh, this%cube, aux2, aux1, this%mesh_cube_map)
    end if
    pot(1:sm%np) = pot(1:sm%np) + M_zI*aux2(1:sm%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)
#endif
  case(POISSON_FFT)
    call X(poisson_fft_solve)(this%fft_solver, sm%mesh, this%cube, pot, rho, this%mesh_cube_map, sm=sm)
  end select

  POP_SUB(X(poisson_solve_sm))
  call profiling_out(prof)
end subroutine X(poisson_solve_sm)
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
