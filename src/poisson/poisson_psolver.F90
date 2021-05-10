!! Copyright (C) 2013 J. Alberdi-Rodriguez
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

#include "global.h"

module poisson_psolver_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use kpoints_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_oct_m
#ifdef HAVE_PSOLVER
  use mesh_function_oct_m
#endif
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use space_oct_m
  use submesh_oct_m

  !! Support for PSolver from BigDFT

#ifdef HAVE_LIBISF
  use poisson_solver
  use dynamic_memory
#endif

#ifdef HAVE_PSOLVER
  use poisson_solver
  use dictionaries, dict_set => set
  use yaml_output, only: yaml_map
  use wrapper_MPI, only: mpi_environment, mpi_environment_set
#endif


  implicit none

  private
  public ::                         &
    poisson_psolver_t,              &
    poisson_psolver_init,           &
    poisson_psolver_end,            &
    poisson_psolver_reinit,         &
    poisson_psolver_global_solve,   &
    poisson_psolver_parallel_solve, &
    poisson_psolver_get_dims

  type poisson_psolver_t
    private
    type(fourier_space_op_t) :: coulb  !< object for Fourier space operations
#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)
    type(coulomb_operator) :: kernel !< choice of kernel, one of options above
#endif
    !> Indicates the boundary conditions (BC) of the problem:
    !!            'F' free BC, isolated systems.
    !!                The program calculates the solution as if the given density is
    !!                "alone" in R^3 space.
    !!            'S' surface BC, isolated in y direction, periodic in xz plane                
    !!                The given density is supposed to be periodic in the xz plane,
    !!                so the dimensions in these direction mus be compatible with the FFT
    !!                Beware of the fact that the isolated direction is y!
    !!            'P' periodic BC.
    !!                The density is supposed to be periodic in all the three directions,
    !!                then all the dimensions must be compatible with the FFT.
    !!                No need for setting up the kernel.
    character(len = 1) :: geocode  = "F" !< 'F' free boundary contition
    !> Indicates the distribution of the data of the input/output array:
    !!            'G' global data. Each process has the whole array of the density 
    !!                which will be overwritten with the whole array of the potential
    !!            'D' distributed data. Each process has only the needed part of the density
    !!                and of the potential. The data distribution is such that each processor
    !!                has the xy planes needed for the calculation AND for the evaluation of the 
    !!                gradient, needed for XC part, and for the White-Bird correction, which
    !!                may lead up to 8 planes more on each side. Due to this fact, the information
    !!                between the processors may overlap.
    character(len = 1), public :: datacode = "G" 

    integer :: isf_order           !< order of the interpolating scaling functions used in the decomposition 
    integer :: localnscatterarr(5)
    integer :: rs_n_global(3)      !< total size of the fft in each direction in real space
    integer :: rs_istart(1:3)      !< where does the local portion of the function start in real space
#ifdef HAVE_PSOLVER
    type(dictionary), pointer :: inputs !input parameters
#endif
    double precision          :: offset
  end type poisson_psolver_t

#ifdef HAVE_PSOLVER
  logical, save :: flib_initialized = .false.
#endif

contains

  !-----------------------------------------------------------------
  subroutine poisson_psolver_init(this, namespace, space, mesh, cube, mu, qq, force_isolated)
    type(poisson_psolver_t), intent(out)   :: this
    type(namespace_t),       intent(in)    :: namespace
    type(space_t),           intent(in)    :: space
    type(mesh_t),            intent(inout) :: mesh
    type(cube_t),            intent(inout) :: cube
    FLOAT,                   intent(in)    :: mu
    FLOAT,                   intent(in)    :: qq(1:MAX_DIM)
    logical, optional,       intent(in)    :: force_isolated

    logical data_is_parallel
#ifdef HAVE_PSOLVER
    FLOAT :: alpha, beta, gamma
    FLOAT :: modq2
    type(mpi_environment) mpi_env
#endif
 
    PUSH_SUB(poisson_psolver_init)
    
#ifdef HAVE_PSOLVER
    if(.not.flib_initialized) then
      call f_lib_initialize()
      flib_initialized = .true.
    end if

    call dict_init(this%inputs)
#elif HAVE_LIBISF
    call f_lib_initialize()
#endif

    if(optional_default(force_isolated, .false.)) then
      this%geocode = "F"
    else
      select case(space%periodic_dim)
      case(0)
        ! Free BC
        this%geocode = "F"
      case(1)
        ! Wire BC
        this%geocode = "W"
        call messages_not_implemented("PSolver support for 1D periodic boundary conditions.")
      case(2)
        ! Surface BC
        this%geocode = "S"
        call messages_not_implemented("PSolver support for 2D periodic boundary conditions.")
      case(3)
        ! Periodic BC
        this%geocode = "P"
        call messages_experimental("PSolver support for 3D periodic boundary conditions.")
      end select
    end if


#ifdef HAVE_PSOLVER
    !Verbosity switch
    call dict_set(this%inputs//'setup'//'verbose', .false.)
    !Order of the Interpolating Scaling Function family
    call dict_set(this%inputs//'kernel'//'isf_order', 16)
    !Mu screening parameter
    call dict_set(this%inputs//'kernel'//'screening', mu)
    !Calculation of the stress tensor
    call dict_set(this%inputs//'kernel'//'stress_tensor', .false.)
#else
    this%isf_order = 16
#endif

    !%Variable PoissonSolverPSolverParallelData
    !%Type logical
    !%Section Hamiltonian::Poisson::PSolver
    !%Default yes
    !%Description
    !% Indicates whether data is partitioned within the PSolver library.
    !% If data is distributed among processes, Octopus uses parallel data-structures 
    !% and, thus, less memory.
    !% If "yes", data is parallelized. The <i>z</i>-axis of the input vector
    !% is split among the MPI processes.
    !% If "no", entire input and output vector is saved in all the MPI processes.
    !% If k-points parallelization is used, "no" must be selected.
    !%End
    call parse_variable(namespace, 'PoissonSolverPSolverParallelData', .true., data_is_parallel)

    call messages_obsolete_variable(namespace, 'PoissonSolverISFParallelData', 'PoissonSolverPSolverParallelData')

    if (data_is_parallel) then
#ifdef HAVE_PSOLVER
      call dict_set(this%inputs//'setup'//'global_data', .false.)
#endif
      this%datacode = "D"
    else 
#ifdef HAVE_PSOLVER
      call dict_set(this%inputs//'setup'//'global_data', .true.)
#endif
      this%datacode = "G"
    end if

#ifdef HAVE_PSOLVER
    call dict_set(this%inputs//'setup'//'verbose', debug%info)

    alpha = mesh%sb%latt%alpha*M_PI/(CNST(180.0))
    beta  = mesh%sb%latt%beta*M_PI/(CNST(180.0))
    gamma = mesh%sb%latt%gamma*M_PI/(CNST(180.0))

    ! Previously, pkernel_init set the communicator used within PSolver to comm_world.
    ! This can be overwritten by passing an optional argument of type(mpi_environment)
    ! to pkernel_init(). This data type is defined within the wrapper_MPI module of 
    ! the Futile library. Futile is a prerequisit for PSolver.

    ! TODO: check that cube%mpi_grp corresponds to the correct mpi group, when parallelizing, e.g., over systems !!
    
    call mpi_environment_set(mpi_env, cube%mpi_grp%rank, cube%mpi_grp%size, cube%mpi_grp%comm, cube%mpi_grp%size )

    this%kernel = pkernel_init(cube%mpi_grp%rank, cube%mpi_grp%size, this%inputs, this%geocode, cube%rs_n_global, &
      mesh%spacing, alpha_bc = alpha, beta_ac = beta, gamma_ab = gamma, mpi_env = mpi_env)
    call pkernel_set(this%kernel, verbose=debug%info)

    !G=0 component
    modq2 = sum(qq(1:space%periodic_dim)**2)
    if (modq2 > M_EPSILON) then
      this%offset = M_ONE/modq2
    else
      this%offset = M_ZERO
    end if

    !Screened coulomb potential (erfc function)
    if (mu > M_EPSILON) then
      if(modq2 > M_EPSILON) then
        this%offset = this%offset*(M_ONE - exp(-modq2/((M_TWO*mu)**2)))
      else
        !Analytical limit of 1/|q|^2*(1-exp(-|q|^2/4mu^2))
        this%offset = M_ONE/((M_TWO*mu)**2)
      end if
    end if
    this%offset = this%offset*M_FOUR*M_PI

#elif HAVE_LIBISF
    this%kernel = pkernel_init(.false., cube%mpi_grp%rank, cube%mpi_grp%size, 0, this%geocode, cube%rs_n_global, &
      mesh%spacing, this%isf_order)
    call pkernel_set(this%kernel, .false.)
#endif

    POP_SUB(poisson_psolver_init)
  end subroutine poisson_psolver_init
  
  !-----------------------------------------------------------------
  subroutine poisson_psolver_end(this)
    type(poisson_psolver_t), intent(inout) :: this

    PUSH_SUB(poisson_psolver_end)

#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)
    call pkernel_free(this%kernel)
    call f_lib_finalize()
#endif
#ifdef HAVE_PSOLVER
    call dict_free(this%inputs)
#endif
    
    POP_SUB(poisson_psolver_end)
  end subroutine poisson_psolver_end

  !-----------------------------------------------------------------
  subroutine poisson_psolver_reinit(this, space, mesh, cube, mu, qq_in)
    type(poisson_psolver_t), intent(inout) :: this
    type(space_t),           intent(in)    :: space
    type(cube_t),            intent(inout) :: cube
    type(mesh_t),            intent(inout) :: mesh
    FLOAT,                   intent(in)    :: mu
    FLOAT,                   intent(in)    :: qq_in(1:MAX_DIM)

    FLOAT :: alpha, beta, gamma
    FLOAT :: qq_abs(1:space%dim), qq(1:space%dim)
    FLOAT :: modq2
    integer :: idim

    PUSH_SUB(poisson_psolver_reinit)

#ifdef HAVE_PSOLVER
    call pkernel_free(this%kernel)
#endif

    !We might change the cell angles
    alpha = mesh%sb%latt%alpha*M_PI/(CNST(180.0))
    beta  = mesh%sb%latt%beta*M_PI/(CNST(180.0))
    gamma = mesh%sb%latt%gamma*M_PI/(CNST(180.0))

#ifdef HAVE_PSOLVER
    call dict_set(this%inputs//'kernel'//'screening',mu)

    this%kernel = pkernel_init(cube%mpi_grp%rank, cube%mpi_grp%size, this%inputs,&
        this%geocode,cube%rs_n_global,mesh%spacing, &
        alpha_bc = alpha, beta_ac = beta, gamma_ab = gamma)
    call pkernel_set(this%kernel,verbose=debug%info)
#endif

    !G=0 component
    !We remove potential umklapp
    do idim = 1, space%periodic_dim
      qq(idim) = qq_in(idim) - anint(qq_in(idim) + M_HALF*CNST(1e-8))
    end do
    qq(space%periodic_dim + 1:space%dim) = M_ZERO
    call kpoints_to_absolute(mesh%sb%latt, qq, qq_abs)
    modq2 = norm2(qq_abs)
    if (modq2 > M_EPSILON) then
      this%offset = M_ONE/modq2
    else
      this%offset = M_ZERO
    end if

    !Screened coulomb potential (erfc function)
    if (mu > M_EPSILON) then
      if (modq2 > M_EPSILON) then
        this%offset = this%offset*(M_ONE - exp(-modq2/((M_TWO*mu)**2)))
      else
        !Analytical limit of 1/|q|^2*(1-exp(-|q|^2/4mu^2))
        this%offset = M_ONE/((M_TWO*mu)**2)
      end if
    end if
    this%offset = this%offset*M_FOUR*M_PI

    POP_SUB(poisson_psolver_reinit)
  end subroutine poisson_psolver_reinit

  !-----------------------------------------------------------------
  subroutine poisson_psolver_parallel_solve(this, mesh, cube, pot, rho,  mesh_cube_map)
    type(poisson_psolver_t),        intent(in), target :: this
    type(mesh_t),                   intent(in)         :: mesh
    type(cube_t),                   intent(in)         :: cube
    FLOAT,                          intent(out)        :: pot(:)
    FLOAT,                          intent(in)         :: rho(:)
    type(mesh_cube_parallel_map_t), intent(in)         :: mesh_cube_map

    type(profile_t), save :: prof
    type(cube_function_t) :: cf   
#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)    
    double precision :: hartree_energy !<  Hartree energy
#endif 
    !> offset:  Total integral on the supercell of the final potential on output.
    !! To be used only in the periodic case, ignored for other boundary conditions.
#ifdef HAVE_PSOLVER
    FLOAT :: offset
#elif HAVE_LIBISF
    double precision :: offset
#endif

    !> pot_ion:  additional external potential that is added to the output
    !! when the XC parameter ixc/=0 and sumpion=.true.
    !! When sumpion=.true., it is always provided in the distributed form,
    !! clearly without the overlapping terms which are needed only for the XC part
    double precision, allocatable :: pot_ion(:,:,:) 

    character(len=3) :: quiet

#ifdef HAVE_PSOLVER
    type(coulomb_operator), pointer :: kernel_pointer
#elif HAVE_LIBISF
    double precision :: strten(6)
#endif

    PUSH_SUB(poisson_psolver_parallel_solve)

    call dcube_function_alloc_RS(cube, cf)

    call dmesh_to_cube_parallel(mesh, rho, cube, cf, mesh_cube_map)
  
    SAFE_ALLOCATE(pot_ion(1:cube%rs_n(1),1:cube%rs_n(2),1:cube%rs_n(3)))

    if(.not.debug%info) then
      quiet = "YES"
    else
      quiet = "NO "
    end if

#ifdef HAVE_PSOLVER
    !The offset is the integral over space of the potential
    !this%offset is the G=0 component of the (screened) Coulomb potential
    !The G=0 component of the Hartree therefore needs to be
    ! multiplied by the G=0 component of the density
    if(this%offset > M_EPSILON) then
      offset = this%offset*dmf_integrate(mesh,rho)
    end if
#endif

    call profiling_in(prof,"PSOLVER_LIBRARY")
#ifdef HAVE_PSOLVER
    kernel_pointer => this%kernel
    call H_potential(this%datacode, kernel_pointer, cf%dRS, pot_ion, hartree_energy, offset, .false., quiet = quiet)
#elif HAVE_LIBISF
    call H_potential(this%datacode, this%kernel, cf%dRS, pot_ion, hartree_energy, offset, .false., quiet = quiet, &
      stress_tensor = strten)
#endif
    call profiling_out(prof)
    SAFE_DEALLOCATE_A(pot_ion)

    call dcube_to_mesh_parallel(cube, cf, mesh, pot, mesh_cube_map)

    call dcube_function_free_RS(cube, cf)

    POP_SUB(poisson_psolver_parallel_solve)
  end subroutine poisson_psolver_parallel_solve

  !-----------------------------------------------------------------
  subroutine poisson_psolver_global_solve(this, mesh, cube, pot, rho, sm)
    type(poisson_psolver_t), intent(in), target :: this
    type(mesh_t),            intent(in)         :: mesh
    type(cube_t),            intent(in)         :: cube
    FLOAT,                   intent(out)        :: pot(:)
    FLOAT,                   intent(in)         :: rho(:)
    type(submesh_t),     optional,  intent(in)    :: sm  !< If present pot and rho are assumed to come from it

    character(len=3) :: quiet
    type(profile_t), save :: prof
    type(cube_function_t) :: cf
#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)
    double precision :: hartree_energy !<  Hartree energy
#endif
    !> offset:  Total integral on the supercell of the final potential on output
    !! To be used only in the periodic case, ignored for other boundary conditions.
#ifdef HAVE_PSOLVER
    FLOAT :: offset
#elif HAVE_LIBISF
    double precision :: offset
#endif

    !> pot_ion:  additional external potential
    !! that is added to the output when the XC parameter ixc/=0 and sumpion=.true.
    !! When sumpion=.true., it is always provided in the distributed form,
    !! clearly without the overlapping terms which are needed only for the XC part
    double precision, allocatable ::  pot_ion(:,:,:) 

#ifdef HAVE_PSOLVER
    type(coulomb_operator), pointer :: kernel_pointer
#elif HAVE_LIBISF
    double precision :: strten(6)
#endif

    PUSH_SUB(poisson_psolver_global_solve)
    
    call dcube_function_alloc_RS(cube, cf)

    if(present(sm)) then
      call dsubmesh_to_cube(sm, rho, cube, cf)
    else
      if(mesh%parallel_in_domains) then
        call dmesh_to_cube(mesh, rho, cube, cf, local=.true.)
      else
        call dmesh_to_cube(mesh, rho, cube, cf)
      end if
    end if

    SAFE_ALLOCATE(pot_ion(1:cube%rs_n(1),1:cube%rs_n(2),1:cube%rs_n(3)))

    if (.not. debug%info) then
      quiet = "YES"
    else
      quiet = "NO "
    end if

#ifdef HAVE_PSOLVER
    !The offset is the integral over space of the potential
    !this%offset is the G=0 component of the (screened) Coulomb potential
    !The G=0 component of the Hartree therefore needs to be
    ! multiplied by the G=0 component of the density
    if(this%offset > M_EPSILON) then
      offset = this%offset*dmf_integrate(mesh, rho)
    end if
#endif

    call profiling_in(prof,"PSOLVER_LIBRARY")

#ifdef HAVE_PSOLVER
    kernel_pointer => this%kernel
    call H_potential(this%datacode, kernel_pointer, &
         cf%dRS,  pot_ion, hartree_energy, offset, .false., &
         quiet = quiet) !optional argument
#elif HAVE_LIBISF
    call H_potential(this%datacode, this%kernel, &
         cf%dRS,  pot_ion, hartree_energy, offset, .false., &
         quiet = quiet, stress_tensor = strten) !optional argument
#endif
    SAFE_DEALLOCATE_A(pot_ion)

    if(present(sm)) then
      call dcube_to_submesh(cube, cf, sm, pot)
    else
      if(mesh%parallel_in_domains) then
        call dcube_to_mesh(cube, cf, mesh, pot, local=.true.)
      else
        call dcube_to_mesh(cube, cf, mesh, pot)
      end if
    end if

    call dcube_function_free_RS(cube, cf)

    POP_SUB(poisson_psolver_global_solve)
  end subroutine poisson_psolver_global_solve

  subroutine poisson_psolver_get_dims(this, cube) 
    type(poisson_psolver_t), intent(inout) :: this
    type(cube_t),            intent(inout) :: cube

#if (defined HAVE_LIBISF) || (defined HAVE_PSOLVER)
    !>    ixc         eXchange-Correlation code. Indicates the XC functional to be used 
    !!                for calculating XC energies and potential. 
    !!                ixc=0 indicates that no XC terms are computed. The XC functional codes follow
    !!                the ABINIT convention.   
    !>    n3d         third dimension of the density. For distributed data, it takes into account 
    !!                the enlarging needed for calculating the XC functionals.
    !!                For global data it is simply equal to n03. 
    !!                When there are too many processes and there is no room for the density n3d=0
    integer :: n3d 
    !>    n3p         third dimension for the potential. The same as n3d, but without 
    !!                taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
    integer :: n3p
    !>    n3pi        Dimension of the pot_ion array, always with distributed data. 
    !!                For distributed data n3pi=n3p
    integer :: n3pi
    !>     i3xcsh     Shift of the density that must be performed to enter in the 
    !!                non-overlapping region. Useful for recovering the values of the potential
    !!                when using GGA XC functionals. If the density starts from rhopot(1,1,1),
    !!                the potential starts from rhopot(1,1,i3xcsh+1). 
    !!                For non-GGA XCs and for global distribution data i3xcsh=0
    integer :: i3xcsh
    !>    i3s         Starting point of the density effectively treated by each processor 
    !!                in the third direction.
    !!                It takes into account also the XC enlarging. The array rhopot will correspond
    !!                To the planes of third coordinate from i3s to i3s+n3d-1. 
    !!                The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
    !!                The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
    !!                For global disposition i3s is equal to distributed case with i3xcsh=0.
    integer :: i3s

    !> use_gradient:  .true. if functional is using the gradient.
    logical :: use_gradient = .false.
    !> use_wb_corr:  .true. if functional is using WB corrections.
    logical :: use_wb_corr = .false.
#endif

    PUSH_SUB(poisson_psolver_get_dims)

    !! Get the dimensions of the cube

#ifdef HAVE_PSOLVER
    call PS_dim4allocation(this%geocode, this%datacode, cube%mpi_grp%rank, cube%mpi_grp%size, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         use_gradient, use_wb_corr, &
         0, n3d, n3p, n3pi, i3xcsh, i3s)
    this%localnscatterarr(:) = (/ n3d, n3p, n3pi, i3xcsh, i3s /)
#elif HAVE_LIBISF
    call PS_dim4allocation(this%geocode, this%datacode, cube%mpi_grp%rank, cube%mpi_grp%size, &
         cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3), &
         use_gradient, use_wb_corr, &
         n3d, n3p, n3pi, i3xcsh, i3s)
    this%localnscatterarr(:) = (/ n3d, n3p, n3pi, i3xcsh, i3s /)
#endif
    cube%rs_n(1:2)      = cube%rs_n_global(1:2)
    cube%rs_n(3)        = this%localnscatterarr(1)
    cube%rs_istart(1:2) = 1
    cube%rs_istart(3)   = this%localnscatterarr(5)
    
    !! With PSolver we don`t care about the Fourier space and its dimensions
    !! We`ll put as in RS
    cube%fs_n_global(1) = cube%rs_n_global(1)
    cube%fs_n_global(2) = cube%rs_n_global(2)
    cube%fs_n_global(3) = cube%rs_n_global(3)
    cube%fs_n(1:2)      = cube%rs_n_global(1:2)
    cube%fs_n(3)        = this%localnscatterarr(1)
    cube%fs_istart(1:2) = 1
    cube%fs_istart(3)   = this%localnscatterarr(5)

    POP_SUB(poisson_psolver_get_dims)
  end subroutine poisson_psolver_get_dims
  
end module poisson_psolver_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
