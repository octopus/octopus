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
subroutine poisson3D_init(gr, geo)
  type(grid_t), intent(inout) :: gr
  type(geometry_t), intent(in) :: geo

  integer :: maxl, iter
  integer :: nx, ny, nz
  FLOAT :: xl, yl, zl

  call push_sub('poisson3D.poisson3D_init')

  ASSERT(poisson_solver >= POISSON_FFT_SPH .and. poisson_solver <= POISSON_SETE)

  !%Variable PoissonSolverMaxMultipole
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% Order of the multipolar expansion for boundary
  !% corrections. Default is 4 for <tt>PoissonSolver = cg_corrected</tt> and <tt>multigrid</tt> and 2
  !% for <tt>fft_corrected</tt>.
  !%End

  !%Variable PoissonSolverMaxIter
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% The maximum number of iterations for conjugate gradient
  !% Poisson solvers. Default is 400.
  !%End

  !%Variable PoissonSolverThreshold
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% The tolerance for the Poisson solution, used by the <tt>cg</tt> and
  !% <tt>multigrid</tt> solvers. Default is <math>10^{-5}</math>.
  !%End

  !%Variable PoissonSolverIncreaseBox
  !%Type logical
  !%Section Hamiltonian::Poisson
  !%Description
  !% If the selected Poisson solver is <tt>cg</tt> or <tt>cg_corrected</tt>, the boundary conditions
  !% have to be calculated by means of performing a multipole expansion. Unfortunately,
  !% if the charge distribution is not contained in a simulation box of approximately
  !% spherical shape, the error can be quite large. Good cases are the spherical box,
  !% the parallelepiped when all dimensions are of similar magnitude, or the cylinder
  !% when the height is not too different to the diameter of the base. Bad cases
  !% are the rest, including the <tt>minimum</tt> box, when the geometry of the molecule is
  !% not compact enough.
  !%
  !% In order to cure this problem, the Hartree problem may be solved in an auxiliary
  !% simulation box, which will contain the original one, but which will be a sphere.
  !% This implies some extra computational effort -- since the density and potential have
  !% to be transferred between boxes -- and extra memory consumption -- since a new
  !% grid has to be stored, which may need quite a lot of memory if you use curvilinear
  !% coordinates.
  !%End

  select case(poisson_solver)
  case(POISSON_CG)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
     call write_info(1)
     call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_corrections_init(corrector, maxl, gr%mesh)
     call poisson_cg_init(gr%mesh, maxl, threshold, iter)

  case(POISSON_CG_CORRECTED)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)
     call parse_logical(datasets_check('PoissonSolverIncreaseBox'), .false., hartree_integrator%increase_box)
     if(gr%mesh%sb%box_shape .eq. SPHERE) hartree_integrator%increase_box = .false.
     if(hartree_integrator%increase_box) then
       write(message(1),'(a)') "Info: Poisson equation will be solved in a larger grid."
       call write_info(1)
       SAFE_ALLOCATE(hartree_integrator%grid)
       call grid_create_largergrid(gr, geo, hartree_integrator%grid)
     end if
     call poisson_corrections_init(corrector, maxl, gr%mesh)
     call poisson_cg_init(gr%mesh, maxl, threshold, iter)


  case(POISSON_MULTIGRID)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)

     call poisson_multigrid_init(mg, gr%mesh, maxl, threshold)

     call grid_create_multigrid(gr, geo)
     
  case(POISSON_ISF)
    call poisson_isf_init(gr%mesh, init_world = all_nodes_default)

  case(POISSON_FFT_SPH)
    call poisson_fft_build_3d_0d(gr, poisson_solver)

  case(POISSON_FFT_CYL)
    call poisson_fft_build_3d_1d(gr)

  case(POISSON_FFT_PLA)
    call poisson_fft_build_3d_2d(gr)

  case(POISSON_FFT_NOCUT)
    call poisson_fft_build_3d_3d(gr)

  case(POISSON_FFT_CORRECTED)
    call poisson_fft_build_3d_0d(gr, poisson_solver)
    call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 2, maxl)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call write_info(1)
    call poisson_corrections_init(corrector, maxl, gr%mesh)

  case(POISSON_SETE)
    nx = gr%mesh%idx%nr(2, 1) - gr%mesh%idx%nr(1, 1) + 1 - 2*gr%mesh%idx%enlarge(1)
    ny = gr%mesh%idx%nr(2, 2) - gr%mesh%idx%nr(1, 2) + 1 - 2*gr%mesh%idx%enlarge(2)
    nz = gr%mesh%idx%nr(2, 3) - gr%mesh%idx%nr(1, 3) + 1 - 2*gr%mesh%idx%enlarge(3)

    xl = 2*gr%mesh%sb%lsize(1)
    yl = 2*gr%mesh%sb%lsize(2)
    zl = 2*gr%mesh%sb%lsize(3)

    call poisson_sete_init(sete_solver, nx, ny, nz, xl, yl, zl, geo%natoms)
     
  end select

  call pop_sub()
end subroutine poisson3D_init




!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
