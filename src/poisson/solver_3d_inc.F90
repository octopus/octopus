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
subroutine poisson3D_init(this, geo, all_nodes_comm)
  type(poisson_t),  intent(inout) :: this
  type(geometry_t), intent(in)    :: geo
  integer,          intent(in)    :: all_nodes_comm

  integer :: maxl, iter
  integer :: nx, ny, nz
  FLOAT   :: xl, yl, zl

  logical :: valid_solver

  PUSH_SUB(poisson3D_init)

  select case(this%method)
  case(POISSON_DIRECT_SUM, POISSON_FMM, POISSON_FFT, POISSON_CG, POISSON_CG_CORRECTED)
    valid_solver = .true.
  case(POISSON_MULTIGRID, POISSON_ISF, POISSON_SETE)
    valid_solver = .true.
  case default
    valid_solver = .false.
  end select

  ASSERT(valid_solver)

  !%Variable PoissonSolverMaxMultipole
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Description
  !% Order of the multipolar expansion for boundary
  !% corrections. Default is 4 for <tt>PoissonSolver = cg_corrected</tt> and <tt>multigrid</tt>, and 2
  !% for <tt>fft</tt> with <tt>PoissonFFTKernel = multipole_correction</tt>.
  !%End

  !%Variable PoissonSolverMaxIter
  !%Type integer
  !%Section Hamiltonian::Poisson
  !%Default 400
  !%Description
  !% The maximum number of iterations for conjugate-gradient
  !% Poisson solvers.
  !%End

  !%Variable PoissonSolverThreshold
  !%Type float
  !%Section Hamiltonian::Poisson
  !%Default 1e-5
  !%Description
  !% The tolerance for the Poisson solution, used by the <tt>cg</tt> and
  !% <tt>multigrid</tt> solvers.
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

  select case(this%method)
  case(POISSON_FMM)
    call poisson_fmm_init(this%params_fmm, this%der, all_nodes_comm)

  case(POISSON_CG)
    call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
    write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
    call messages_info(1)
    call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
    call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
    call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_CG_CORRECTED)
    call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
    call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
    call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)
    call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_MULTIGRID)
    call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
    call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)

    call poisson_multigrid_init(this%mg, this%der%mesh, maxl, threshold)
     
  case(POISSON_ISF)
    call poisson_isf_init(this%isf_solver, this%der%mesh, this%cube, all_nodes_comm, init_world = this%all_nodes_default)

  case(POISSON_FFT)

    call poisson_fft_init(this%fft_solver, this%der%mesh, this%cube, this%kernel)

    if (this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
      call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 2, maxl)
      write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
      call messages_info(1)
      call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
    end if

  case(POISSON_SETE)
    nx = this%der%mesh%idx%ll(1) 
    ny = this%der%mesh%idx%ll(2)
    nz = this%der%mesh%idx%ll(3)
    xl = 2*this%der%mesh%sb%lsize(1)
    yl = 2*this%der%mesh%sb%lsize(2)
    zl = 2*this%der%mesh%sb%lsize(3)

    call poisson_sete_init(this%sete_solver, nx, ny, nz, xl, yl, zl, geo%natoms)
     
  end select

  POP_SUB(poisson3D_init)
end subroutine poisson3D_init


subroutine poisson_solve_direct(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  FLOAT :: prefactor
  integer  :: ip, jp
  FLOAT    :: xx(this%der%mesh%sb%dim), yy(this%der%mesh%sb%dim)
#ifdef HAVE_MPI
  FLOAT    :: tmp
  FLOAT, allocatable :: pvec(:) 
#endif

  PUSH_SUB(poisson_solve_direct)

  ASSERT(this%method == POISSON_DIRECT_SUM)

  select case(this%der%mesh%sb%dim)
  case(3)
    prefactor = M_TWO*M_PI*(M_THREE/(M_PI*M_FOUR))**(M_TWOTHIRD)
  case(2)
    prefactor = M_TWO*sqrt(M_PI)
  case default
    message(1) = "Internal error: poisson_solve_direct can only be called for 2D or 3D."
    ! why not? all that is needed is the appropriate prefactors to be defined above, actually. then 1D, 4D etc. can be done
    call messages_fatal(1)
  end select

#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))

    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xx(:) = mesh_x_global(this%der%mesh, ip)
      do jp = 1, this%der%mesh%np
        if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp) then
          if(this%der%mesh%use_curvilinear) then
            pvec(jp) = rho(jp)*prefactor / (this%der%mesh%vol_pp(jp)**M_THIRD)
          else
            pvec(jp) = rho(jp)*prefactor / (this%der%mesh%volume_element**M_THIRD)
          endif
       else
          yy(:) = this%der%mesh%x(jp,:)
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
    do ip = 1, this%der%mesh%np
      xx(:) = this%der%mesh%x(ip,:)
      do jp = 1, this%der%mesh%np
        if(this%der%mesh%use_curvilinear) then
          if(ip == jp) then
            pot(ip) = pot(ip) + prefactor*rho(ip)/(this%der%mesh%vol_pp(jp)**M_THIRD)*this%der%mesh%vol_pp(jp)
          else
            yy(:) = this%der%mesh%x(jp,:)
            pot(ip) = pot(ip) + rho(jp)/sqrt(sum((xx-yy)**2))*this%der%mesh%vol_pp(jp)
          endif
        else
          if(ip == jp) then
            pot(ip) = pot(ip) + prefactor*rho(ip)/(this%der%mesh%volume_element**M_THIRD)
          else
            yy(:) = this%der%mesh%x(jp,1:3)
            pot(ip) = pot(ip) + rho(jp)/sqrt(sum((xx-yy)**2))
          endif
        end if
      end do
    end do
    if(.not. this%der%mesh%use_curvilinear) then
      pot(:) = pot(:) * this%der%mesh%volume_element
    endif
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(poisson_solve_direct) 
end subroutine poisson_solve_direct



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
