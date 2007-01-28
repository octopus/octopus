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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine poisson3D_init(gr, geo)
  type(grid_t), intent(inout) :: gr
  type(geometry_t), intent(in) :: geo

  integer :: maxl

  call push_sub('poisson3D.poisson3D_init')

#if defined(HAVE_FFT)
  ASSERT(poisson_solver >= FFT_SPH .and. poisson_solver <= ISF)
#else
  ASSERT(poisson_solver >= CG .and. poisson_solver <= ISF)
#endif

  !%Variable PoissonSolverMaxMultipole
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% Order of the multipolar expansion for boundary
  !% corrections. Default is 4 for cg_corrected and multigrid and 2
  !% for fft_corrected.
  !%End


  !%Variable PoissonSolverThreshold
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% The tolerance for the poisson solution, used by the cg and
  !% multigrid solvers. Default is <math>10^{-5}</math>.
  !%End

  !%Variable PoissonSolverIncreaseBox
  !%Type logical
  !%Section Hamiltonian::Poisson
  !%Description
  !% In case the selected Poisson solver is cg or cg_corrected, the boundary conditions
  !% have to be calculated by means of performing a multipole expansion. Unfortunately,
  !% if the charge distribution is not contained in a simulation box of approximately
  !% spherical shape, the error can be quite large. Good cases are the spherical box,
  !% the parallelpiped when all dimension are of similar magnitude, or the cylinder
  !% in case the height is not too different to the diameter of the base. Bad cases
  !% are the rest, including the "minimum" box, when the geometry of the molecule is
  !% not compact enough.
  !%
  !% In order to cure this problem, the Hartree problem may me solved in an auxiliary
  !% simulation box, which will contain the original one, but which will be a sphere.
  !% This implies some extra computational effort -- since the density and potential have
  !% to be transeferred between boxes --, and extra memory consumption -- since a new
  !% grid has to be stored, which may need quite a lot of memory if you use curvilinear
  !% coordinates.
  !%End

  select case(poisson_solver)
  case(CG)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
     call write_info(1)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_corrections_init(corrector, maxl, gr%m)
     call poisson_cg_init(gr%m, maxl, threshold)

  case(CG_CORRECTED)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)
     call loct_parse_logical(check_inp('PoissonSolverIncreaseBox'), .false., hartree_integrator%increase_box)
     if(gr%m%sb%box_shape .eq. SPHERE) hartree_integrator%increase_box = .false.
     if(hartree_integrator%increase_box) then
       write(message(1),'(a)') "Info: Poisson' equation will be solved in a larger grid."
       call write_info(1)
       call grid_create_largergrid(gr, geo, hartree_integrator%grid)
     end if
     call poisson_corrections_init(corrector, maxl, gr%m)
     call poisson_cg_init(gr%m, maxl, threshold)


  case(MULTIGRILLA)
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 4, maxl)
     call loct_parse_float(check_inp('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)

     call poisson_multigrid_init(mg, gr%m, maxl, threshold)

     call grid_create_multigrid(gr, geo)
     
  case(ISF)
    call poisson_isf_init(gr%m)
     
  end select

#ifdef HAVE_FFT
  if (poisson_solver <= FFT_CORRECTED) call poisson_fft_build_3d(gr, poisson_solver)

  if (poisson_solver == FFT_CORRECTED) then
     call loct_parse_int(check_inp('PoissonSolverMaxMultipole'), 2, maxl)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call write_info(1)
     call poisson_corrections_init(corrector, maxl, gr%m)
  end if
#endif

  call pop_sub()

end subroutine poisson3D_init



