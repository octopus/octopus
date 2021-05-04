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

! ---------------------------------------------------------
subroutine poisson_kernel_init(this, namespace, space, all_nodes_comm)
  type(poisson_t),   intent(inout) :: this
  type(namespace_t), intent(in)    :: namespace
  type(space_t),     intent(in)    :: space
  integer,           intent(in)    :: all_nodes_comm

  integer :: maxl, iter, dim_electronic
  logical :: valid_solver
  FLOAT :: qq(1:MAX_DIM)
  
  PUSH_SUB(poisson_kernel_init)

  select case(this%method)
  case(POISSON_DIRECT_SUM, POISSON_FMM, POISSON_FFT, POISSON_CG, POISSON_CG_CORRECTED)
    valid_solver = .true.
  case(POISSON_MULTIGRID, POISSON_ISF, POISSON_PSOLVER, POISSON_POKE)
    valid_solver = .true.
  case(POISSON_NO)
    write(message(1),'(a)')'Info: you have elected to not use a Poisson solver.'
    write(message(2),'(a)')' Hartree potential and energy will be 0'
    call messages_info(2)
    valid_solver = .true.
  case default
    valid_solver = .false.
  end select

  ASSERT(valid_solver)

  !%Variable PoissonSolverMaxMultipole
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% Order of the multipolar expansion for boundary corrections. 
  !%
  !% The Poisson solvers <tt>multigrid</tt>, <tt>cg</tt>, and <tt>cg_corrected</tt>
  !% (and <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>)
  !% do a multipolar expansion of the given
  !% charge density, such that <math>\rho = \rho_{multip.expansion}+\Delta
  !% \rho</math>. The Hartree potential due to the <math>\rho_{multip.expansion}</math> is
  !% calculated analytically, while the Hartree potential due to <math>\Delta \rho</math>
  !% is calculated with either a multigrid or cg solver.
  !% The order of the multipolar expansion is set by this variable.
  !%
  !% Default is 4 for <tt>PoissonSolver = cg_corrected</tt> and <tt>multigrid</tt>, and 2
  !% for <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>.
  !%End

  !%Variable PoissonSolverMaxIter
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Default 500
  !%Description
  !% The maximum number of iterations for conjugate-gradient
  !% Poisson solvers.
  !%End

  !%Variable PoissonSolverThreshold
  !%Type float
  !%Section Hamiltonian::Poisson
  !%Default 1e-6
  !%Description
  !% The tolerance for the Poisson solution, used by the <tt>cg</tt>,
  !% <tt>cg_corrected</tt>, and <tt>multigrid</tt> solvers.
  !%End

  !! This variable is disabled for the moment
  !!
  !!Variable PoissonSolverIncreaseBox
  !!Type logical
  !!Section Hamiltonian::Poisson
  !!Description
  !! (experimental) If the selected Poisson solver is
  !! <tt>cg_corrected</tt> the boundary conditions have to be
  !! calculated by means of performing a multipole
  !! expansion. Unfortunately, if the charge distribution is not
  !! contained in a simulation box of approximately spherical shape,
  !! the error can be quite large. Good cases are the spherical box,
  !! the parallelepiped when all dimensions are of similar magnitude,
  !! or the cylinder when the height is not too different to the
  !! diameter of the base. Bad cases are the rest, including the
  !! <tt>minimum</tt> box, when the geometry of the molecule is not
  !! compact enough.
  !!
  !! In order to cure this problem, the Hartree problem may be solved
  !! in an auxiliary simulation box, which will contain the original
  !! one, but which will be a sphere.  This implies some extra
  !! computational effort -- since the density and potential have to
  !! be transferred between boxes -- and extra memory consumption --
  !! since a new grid has to be stored, which may need quite a lot of
  !! memory if you use curvilinear coordinates.
  !!End
  
  if(this%is_dressed) then
    dim_electronic = space%dim -1
  else
    dim_electronic = space%dim
  end if
  
  if(dim_electronic == 1) then
    !%Variable Poisson1DSoftCoulombParam
    !%Type float
    !%Default 1.0 bohr
    !%Section Hamiltonian::Poisson
    !%Description
    !% When <tt>Dimensions = 1</tt>, to prevent divergence, the Coulomb interaction treated by the Poisson
    !% solver is not <math>1/r</math> but <math>1/\sqrt{a^2 + r^2}</math>, where this variable sets the value of <math>a</math>.
    !%End
    call parse_variable(namespace, 'Poisson1DSoftCoulombParam', M_ONE, this%poisson_soft_coulomb_param, units_inp%length)
  else
    this%poisson_soft_coulomb_param = M_ZERO
  end if

  select case(this%method)
  case(POISSON_FMM)
    call poisson_fmm_init(this%params_fmm, space, this%der, all_nodes_comm)

  case(POISSON_CG)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
    call messages_info(1)
    call parse_variable(namespace, 'PoissonSolverMaxIter', 500, iter)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    call poisson_corrections_init(this%corrector, namespace, space, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_CG_CORRECTED)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    call parse_variable(namespace, 'PoissonSolverMaxIter', 500, iter)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)
    call poisson_corrections_init(this%corrector, namespace, space, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_MULTIGRID)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)

    call poisson_multigrid_init(this%mg, namespace, space, this%der%mesh, maxl, threshold)
     
  case(POISSON_ISF)
    call poisson_isf_init(this%isf_solver, namespace, this%der%mesh, this%cube, all_nodes_comm, init_world = this%all_nodes_default)
    
  case(POISSON_PSOLVER)
    !! We`ll use the MPI_WORLD_COMM, to use all the available processes for the
    !! Poisson solver
    if(this%all_nodes_default) then
      this%cube%mpi_grp = mpi_world
    else
      this%cube%mpi_grp = this%der%mesh%mpi_grp
    end if
    qq(1:MAX_DIM) = M_ZERO
    call poisson_psolver_init(this%psolver_solver, namespace, space, this%der%mesh, this%cube, M_ZERO, qq)
    call poisson_psolver_get_dims(this%psolver_solver, this%cube)
    this%cube%parallel_in_domains = this%psolver_solver%datacode == "D" .and. mpi_world%size > 1
    if (this%cube%parallel_in_domains) then
      call mesh_cube_parallel_map_init(this%mesh_cube_map, this%der%mesh, this%cube)
    end if
    
  case(POISSON_FFT)

    call poisson_fft_init(this%fft_solver, namespace, space, this%der%mesh, this%cube, this%kernel, &
      soft_coulb_param = this%poisson_soft_coulomb_param)
    ! soft parameter has no effect unless in 1D

    if (this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
      call parse_variable(namespace, 'PoissonSolverMaxMultipole', 2, maxl)
      write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
      call messages_info(1)
      call poisson_corrections_init(this%corrector, namespace, space, maxl, this%der%mesh)
    end if

  case(POISSON_NO)
    call poisson_no_init(this%no_solver, this%der%mesh, this%cube)
  end select

  POP_SUB(poisson_kernel_init)
end subroutine poisson_kernel_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
