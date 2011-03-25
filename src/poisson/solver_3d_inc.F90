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
subroutine poisson3D_init(this, geo, all_nodes_comm)
  type(poisson_t),  intent(inout) :: this
  type(geometry_t), intent(in)    :: geo
  integer,          intent(in)    :: all_nodes_comm

  integer :: maxl, iter
  integer :: nx, ny, nz
  FLOAT   :: xl, yl, zl

  type(poisson_fmm_t) :: fmm_params
  logical :: valid_solver

  PUSH_SUB(poisson3D_init)

  valid_solver = (this%method >= POISSON_FFT_SPH .and. this%method<=POISSON_SETE) &
    .or. this%method == POISSON_DIRECT_SUM_3D .or. this%method == POISSON_FMM

  ASSERT(valid_solver)

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
#ifdef HAVE_LIBFM
    call poisson_fmm_init(this%params_fmm, this%der%mesh, all_nodes_comm)
#endif

  case(POISSON_CG)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
     call messages_info(1)
     call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
     call poisson_cg_init(this%der%mesh, maxl, threshold, iter)

  case(POISSON_CG_CORRECTED)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     call parse_integer(datasets_check('PoissonSolverMaxIter'), 400, iter)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call messages_info(1)
     call poisson_corrections_init(this%corrector, maxl, this%der%mesh)
     call poisson_cg_init(this%der%mesh, maxl, threshold, iter)

  case(POISSON_MULTIGRID)
     call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 4, maxl)
     call parse_float(datasets_check('PoissonSolverThreshold'), CNST(1.0e-6), threshold)
     write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
     call messages_info(1)

     call poisson_multigrid_init(this%mg, this%der%mesh, maxl, threshold)
     
  case(POISSON_ISF)
     call poisson_isf_init(this%isf_solver, this%der%mesh, all_nodes_comm, init_world = this%all_nodes_default)

  case(POISSON_FFT_SPH)
    call poisson_fft_build_3d_0d(this%der%mesh, this%method)

  case(POISSON_FFT_CYL)
    call poisson_fft_build_3d_1d(this%der%mesh)

  case(POISSON_FFT_PLA)
    call poisson_fft_build_3d_2d(this%der%mesh)

  case(POISSON_FFT_NOCUT)
    call poisson_fft_build_3d_3d(this%der%mesh)

  case(POISSON_FFT_CORRECTED)
    call poisson_fft_build_3d_0d(this%der%mesh, this%method)
    call parse_integer(datasets_check('PoissonSolverMaxMultipole'), 2, maxl)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)
    call poisson_corrections_init(this%corrector, maxl, this%der%mesh)

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


subroutine poisson3D_solve_direct(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  integer  :: ip, jp
  FLOAT    :: xx(3), yy(3)
#ifdef HAVE_MPI
  FLOAT    :: tmp, xg(1:MAX_DIM)
  FLOAT, allocatable :: pvec(:) 
#endif

  ASSERT(this%method == POISSON_DIRECT_SUM_3D)

  PUSH_SUB(poisson3D_solve_direct)
#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))

    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx(1:3) = xg(1:3) 
      do jp = 1, this%der%mesh%np
        if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp) then
          if((this%der%mesh%use_curvilinear .eqv. .false.) .and. (this%der%mesh%spacing(1)==this%der%mesh%spacing(2)) .and. &
               (this%der%mesh%spacing(2)==this%der%mesh%spacing(3)) .and.&
               (this%der%mesh%spacing(1)==this%der%mesh%spacing(3))) then
             pvec(jp) = rho(jp)*2.380077363979553356918/(this%der%mesh%spacing(1)) 
          else
             pvec(jp) = (1.0/(this%der%mesh%spacing(1)*this%der%mesh%spacing(2)*&
                  this%der%mesh%spacing(3))**(1./3.))*rho(jp)*M_TWO*M_PI*(3./(M_PI*4.))**(2./3.) 
          end if 
       else
          yy(:) = this%der%mesh%x(jp,1:3)
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
      xx(:) = this%der%mesh%x(ip,1:3)
      do jp = 1, this%der%mesh%np
        if(ip == jp) then
          pot(ip) = pot(ip) + M_TWO*sqrt(M_PI)*rho(ip)/this%der%mesh%spacing(1)*this%der%mesh%vol_pp(jp)
        else
          yy(:) = this%der%mesh%x(jp,1:3)
          pot(ip) = pot(ip) + rho(jp)/sqrt(sum((xx-yy)**2))*this%der%mesh%vol_pp(jp)
        end if
      end do
    end do
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(poisson3D_solve_direct) 
end subroutine poisson3D_solve_direct



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
