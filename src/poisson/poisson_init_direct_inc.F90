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
subroutine poisson_kernel_init(this, namespace, all_nodes_comm)
  type(poisson_t),   intent(inout) :: this
  type(namespace_t), intent(in)    :: namespace
  integer,           intent(in)    :: all_nodes_comm

  integer :: maxl, iter
  logical :: valid_solver

  PUSH_SUB(poisson_kernel_init)

  select case(this%method)
  case(POISSON_DIRECT_SUM, POISSON_FMM, POISSON_FFT, POISSON_CG, POISSON_CG_CORRECTED)
    valid_solver = .true.
  case(POISSON_MULTIGRID, POISSON_ISF, POISSON_LIBISF, POISSON_POKE)
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

  if(this%der%mesh%sb%dim == 1) then
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
    call poisson_fmm_init(this%params_fmm, this%der, all_nodes_comm)

  case(POISSON_CG)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    write(message(1),'(a,i2)')'Info: Boundary conditions fixed up to L =',  maxl
    call messages_info(1)
    call parse_variable(namespace, 'PoissonSolverMaxIter', 400, iter)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    call poisson_corrections_init(this%corrector, namespace, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_CG_CORRECTED)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    call parse_variable(namespace, 'PoissonSolverMaxIter', 400, iter)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)
    call poisson_corrections_init(this%corrector, namespace, maxl, this%der%mesh)
    call poisson_cg_init(threshold, iter)

  case(POISSON_MULTIGRID)
    call parse_variable(namespace, 'PoissonSolverMaxMultipole', 4, maxl)
    call parse_variable(namespace, 'PoissonSolverThreshold', CNST(1.0e-6), threshold)
    write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
    call messages_info(1)

    call poisson_multigrid_init(this%mg, namespace, this%der%mesh, maxl, threshold)
     
  case(POISSON_ISF)
    call poisson_isf_init(this%isf_solver, namespace, this%der%mesh, this%cube, all_nodes_comm, init_world = this%all_nodes_default)
    
  case(POISSON_LIBISF)
    !! We`ll use the MPI_WORLD_COMM, to use all the available processes for the
    !! Poisson solver
    if(this%all_nodes_default) then
      this%cube%mpi_grp = mpi_world
    else
      this%cube%mpi_grp = this%der%mesh%mpi_grp
    end if
    call poisson_libisf_init(this%libisf_solver, namespace, this%der%mesh, this%cube)
    call poisson_libisf_get_dims(this%libisf_solver, this%cube)
    this%cube%parallel_in_domains = this%libisf_solver%datacode == "D" .and. mpi_world%size > 1
    if (this%cube%parallel_in_domains) then
      call mesh_cube_parallel_map_init(this%mesh_cube_map, this%der%mesh, this%cube)
    end if
    
  case(POISSON_FFT)
    this%fft_solver%precalculated_g = .false.

    call poisson_fft_init(this%fft_solver, namespace, this%der%mesh, this%cube, this%kernel, &
      soft_coulb_param = this%poisson_soft_coulomb_param, qq = this%qq, singul = M_ZERO)
    ! soft parameter has no effect unless in 1D

    if (this%kernel == POISSON_FFT_KERNEL_CORRECTED) then
      call parse_variable(namespace, 'PoissonSolverMaxMultipole', 2, maxl)
      write(message(1),'(a,i2)')'Info: Multipoles corrected up to L =',  maxl
      call messages_info(1)
      call poisson_corrections_init(this%corrector, namespace, maxl, this%der%mesh)
    end if

  case(POISSON_NO)
    call poisson_no_init(this%no_solver, this%der%mesh, this%cube)
  end select

  POP_SUB(poisson_kernel_init)
end subroutine poisson_kernel_init


!-----------------------------------------------------------------
subroutine poisson_kernel_reinit(this, namespace, qq, singul)
  type(poisson_t),   intent(inout) :: this
  type(namespace_t), intent(in)    :: namespace
  FLOAT,             intent(in)    :: qq(:)
  FLOAT,             intent(in)    :: singul

  type(profile_t), save :: prof

  PUSH_SUB(poisson_kernel_reinit)

  call profiling_in(prof, 'POISSON_REINIT')

  !We only reinitialize the poisson sover if needed
  if(any(abs(this%qq(1:this%der%mesh%sb%periodic_dim) - qq(1:this%der%mesh%sb%periodic_dim)) > M_EPSILON)) then
    this%qq(1:this%der%mesh%sb%periodic_dim) = qq(1:this%der%mesh%sb%periodic_dim)
   
    select case(this%method)
    case(POISSON_FFT)
      if(.not. this%fft_solver%precalculated_g) then
        call poisson_fft_precalculate_g(this%fft_solver, this%der%mesh, this%cube)
      end if
      call poisson_fft_end(this%fft_solver)
      call poisson_fft_init(this%fft_solver, namespace, this%der%mesh, this%cube, this%kernel, &
        soft_coulb_param = this%poisson_soft_coulomb_param, qq = this%qq, singul = singul)
    case default
      call messages_not_implemented("poisson_kernel_reinit with other methods than FFT")
    end select

  end if

  call profiling_out(prof)

  POP_SUB(poisson_kernel_reinit)
end subroutine poisson_kernel_reinit


!-----------------------------------------------------------------
subroutine poisson_solve_direct(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  FLOAT                :: prefactor, aa1, aa2, aa3, aa4
  integer              :: ip, jp
  integer              :: dim           !< physical dimensions
  integer              :: dim_effective !< effective dimensions (= dim,  if no soft coulomb)
                                        !<                      (= dim+1 if soft coulomb)

  FLOAT                :: xx1(1:MAX_DIM+1), xx2(1:MAX_DIM+1), xx3(1:MAX_DIM+1), xx4(1:MAX_DIM+1)
  FLOAT                :: xx(1:MAX_DIM+1), yy(1:MAX_DIM+1) 

  logical              :: include_diag

#ifdef HAVE_MPI
  FLOAT                :: xg(MAX_DIM)
  integer, allocatable :: ip_v(:), part_v(:)
  FLOAT, allocatable   :: pvec(:), tmp(:)
#endif

  PUSH_SUB(poisson_solve_direct)

  dim = this%der%mesh%sb%dim
  ASSERT(this%method == POISSON_DIRECT_SUM)


  if (this%poisson_soft_coulomb_param**2 > M_ZERO) then
    dim_effective = dim+1
    xx(dim_effective) = this%poisson_soft_coulomb_param
    xx1(dim_effective) = this%poisson_soft_coulomb_param
    xx2(dim_effective) = this%poisson_soft_coulomb_param
    xx3(dim_effective) = this%poisson_soft_coulomb_param
    xx4(dim_effective) = this%poisson_soft_coulomb_param
    yy(dim_effective) = M_ZERO
    include_diag = .true.
  else
    dim_effective = dim
    include_diag = .false.
  endif


  select case(dim)
  case(3)
    prefactor = M_TWO*M_PI*(M_THREE/(M_PI*M_FOUR))**(M_TWOTHIRD)
  case(2)
    prefactor = M_TWO*sqrt(M_PI)
  case(1)
    prefactor = M_ONE 
  case default
    message(1) = "Internal error: poisson_solve_direct can only be called for 1D, 2D or 3D."
    ! why not? all that is needed is the appropriate prefactors to be defined above, actually. then 1D, 4D etc. can be done
    call messages_fatal(1)
  end select

  if(.not. this%der%mesh%use_curvilinear) then
    prefactor = prefactor / (this%der%mesh%volume_element**(M_ONE/this%der%mesh%sb%dim))
  end if

#ifdef HAVE_MPI
  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))
    SAFE_ALLOCATE(part_v(1:this%der%mesh%np_global))
    SAFE_ALLOCATE(ip_v(1:this%der%mesh%np_global))
    SAFE_ALLOCATE(tmp(1:this%der%mesh%np_global))
    do ip = 1, this%der%mesh%np_global
      ip_v(ip) = ip
    end do
    call partition_get_partition_number(this%der%mesh%inner_partition, this%der%mesh%np_global, ip_v, part_v)
    
    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx(1:dim) = xg(1:dim)
      if(this%der%mesh%use_curvilinear) then
        do jp = 1, this%der%mesh%np
          if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp .and. .not. include_diag) then
            pvec(jp) = rho(jp)*prefactor**(M_ONE - M_ONE/this%der%mesh%sb%dim)
          else
            yy(1:dim) = this%der%mesh%x(jp, 1:dim)
            pvec(jp) = rho(jp)/sqrt(sum((xx(1:dim_effective) - yy(1:dim_effective))**2))
          end if
        end do
      else
        do jp = 1, this%der%mesh%np
          if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp .and. .not. include_diag) then
            pvec(jp) = rho(jp)*prefactor
          else
            yy(1:dim) = this%der%mesh%x(jp, 1:dim)
            pvec(jp) = rho(jp)/sqrt(sum((xx(1:dim_effective) - yy(1:dim_effective))**2))
          end if
        end do
      end if
      tmp(ip) = dmf_integrate(this%der%mesh, pvec, reduce = .false.)
    end do

    call comm_allreduce(this%der%mesh%mpi_grp%comm, tmp)

    do ip = 1, this%der%mesh%np_global
      if (part_v(ip) == this%der%mesh%vp%partno) then
        pot(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno)) = tmp(ip)
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)
    SAFE_DEALLOCATE_A(tmp)

  else ! serial mode
#endif
    
    do ip = 1, this%der%mesh%np - 4 + 1, 4

      xx1(1:dim) = this%der%mesh%x(ip    , 1:dim)
      xx2(1:dim) = this%der%mesh%x(ip + 1, 1:dim)
      xx3(1:dim) = this%der%mesh%x(ip + 2, 1:dim)
      xx4(1:dim) = this%der%mesh%x(ip + 3, 1:dim)
      
      if(this%der%mesh%use_curvilinear) then

        if(.not. include_diag) then
          aa1 = prefactor*rho(ip    )*this%der%mesh%vol_pp(ip    )**(M_ONE - M_ONE/this%der%mesh%sb%dim)
          aa2 = prefactor*rho(ip + 1)*this%der%mesh%vol_pp(ip + 1)**(M_ONE - M_ONE/this%der%mesh%sb%dim)
          aa3 = prefactor*rho(ip + 2)*this%der%mesh%vol_pp(ip + 2)**(M_ONE - M_ONE/this%der%mesh%sb%dim)
          aa4 = prefactor*rho(ip + 3)*this%der%mesh%vol_pp(ip + 3)**(M_ONE - M_ONE/this%der%mesh%sb%dim)
        end if
        
        !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)

        do jp = 1, this%der%mesh%np
          yy(1:dim) = this%der%mesh%x(jp, 1:dim)
          if(ip     /= jp .or. include_diag) aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_effective) - yy(1:dim_effective))**2)) & 
                                                         *this%der%mesh%vol_pp(jp)
          if(ip + 1 /= jp .or. include_diag) aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim_effective) - yy(1:dim_effective))**2)) & 
                                                         *this%der%mesh%vol_pp(jp)
          if(ip + 2 /= jp .or. include_diag) aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim_effective) - yy(1:dim_effective))**2)) & 
                                                         *this%der%mesh%vol_pp(jp)
          if(ip + 3 /= jp .or. include_diag) aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim_effective) - yy(1:dim_effective))**2)) & 
                                                         *this%der%mesh%vol_pp(jp)
        end do

      else

        if(.not. include_diag) then
          aa1 = prefactor*rho(ip    )
          aa2 = prefactor*rho(ip + 1)
          aa3 = prefactor*rho(ip + 2)
          aa4 = prefactor*rho(ip + 3)
        else
          aa1 = M_ZERO
          aa2 = M_ZERO
          aa3 = M_ZERO
          aa4 = M_ZERO
        end if

        !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)

        do jp = 1, this%der%mesh%np
          yy(1:dim) = this%der%mesh%x(jp, 1:dim)
          if(ip     /= jp .or. include_diag) aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_effective) - yy(1:dim_effective))**2))
          if(ip + 1 /= jp .or. include_diag) aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim_effective) - yy(1:dim_effective))**2))
          if(ip + 2 /= jp .or. include_diag) aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim_effective) - yy(1:dim_effective))**2))
          if(ip + 3 /= jp .or. include_diag) aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim_effective) - yy(1:dim_effective))**2))
        end do
        
      end if

      pot(ip    ) = this%der%mesh%volume_element*aa1
      pot(ip + 1) = this%der%mesh%volume_element*aa2
      pot(ip + 2) = this%der%mesh%volume_element*aa3
      pot(ip + 3) = this%der%mesh%volume_element*aa4
      
    end do
    
    do ip = ip, this%der%mesh%np

      aa1 = CNST(0.0)

      xx1(1:dim) = this%der%mesh%x(ip,1:dim)
      if(this%der%mesh%use_curvilinear) then
        do jp = 1, this%der%mesh%np
          yy(1:dim) = this%der%mesh%x(jp, 1:dim)
          if(ip == jp .and. .not. include_diag) then
            aa1 = aa1 + prefactor*rho(ip)*this%der%mesh%vol_pp(jp)**(M_ONE - M_ONE/this%der%mesh%sb%dim)
          else
            aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_effective) - yy(1:dim_effective))**2))*this%der%mesh%vol_pp(jp)
          end if
        end do
      else
        do jp = 1, this%der%mesh%np
          yy(1:dim) = this%der%mesh%x(jp, 1:dim)
          if(ip == jp) then
            aa1 = aa1 + prefactor*rho(ip)
          else
            aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_effective) - yy(1:dim_effective))**2))
          end if
        end do
      end if

      pot(ip) = this%der%mesh%volume_element*aa1
      
    end do
    
#ifdef HAVE_MPI
  end if
#endif

  POP_SUB(poisson_solve_direct) 
end subroutine poisson_solve_direct

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
