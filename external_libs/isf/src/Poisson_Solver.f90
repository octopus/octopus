!> @file
!!    Define the module for Poisson Solver
!!
!! @author
!!    Luigi Genovese (February 2007)
!!    PSolverNC added by Anders Bergman, March 2008
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module used by the Poisson Solver library.
!! It must be used in the parent routine. 
!! @details
!!    In the main routine in which the Poisson Solver is called
!!    -# The Poisson kernel must be declared as a pointer, then the 
!!       routine createKernel can be called. On exit, the kernel will be allocated and
!!       ready to use. See the documentation of the createKernel routine for more details
!!    -# The correct sizes for allocating the density/potential and the pot_ion arrays
!!       are given from the routine PS_dim4allocation (see routine documentation for details).
!!       Its argument MUST be in agreement with the arguments of the PSolver routine. 
!!       WARNING: No cross-check of the arguments is performed!
!!    -# The PSolver routine can then be called. On exit, the Hartree potential is computed
!!       and summed (following ixc value) to XC and external potential. 
!!       The input array is overwritten. Again, see routine documentation for details.
!!    -# QUICK INSTRUCTION FOR THE IMPATIENT:If you want to use the Poisson Solver in the 
!!       "ordinary" way, for a grid of dimensions nx,ny,nz and grid spacings hx,hy,hz, 
!!       just create the Kernel with
!!           call createKernel(geocode,nx,ny,nz,hx,hy,hz,14,0,1,kernel)
!!       where kernel is a pointer as described above; 
!!       geocode is 'F','S' or 'P' for Free, Surfaces of Periodic BC respectively.
!!       (Beware that for Surfaces BC the isolated direction is y!)
!!       After that you can calculate the potential with
!!           call PSolver(geocode,'G',0,1,nx,ny,nz,0,hx,hy,hz,&
!!                rhopot,kernel,fake_arr,energy,fake_exc,fake_vxc,0.d0,.false.,1)
!!       where:
!!          rhopot    is the density on input, and the electrostatic potential on output
!!                    (an array of dimension(nx,ny,nz))
!!          energy    is the result of @f$ 1/2 \int dx rho(x) potA @f$(x)
!!          fake_arr  is an array of dimension(1), untouched
!!          fake_*xc  values of the XC energies, automatically zero in that case
!!
!!       Any other changment of the arguments require reading of the documentation.
!!       See documentations of the Public routines
!! @warning
!!    This module REQUIRES the module of XC functional from ABINIT, defs_xc, which
!!    require defs_basis and defs_datatypes. 
!!    Such routines are provided inside the abinit directory of this bundle.
!!    They are based on XC functionals present in ABINIT 5.x
!!    If you want to use this Poisson Solver without the XC functionals, you can comment out
!!    the XC part in the PSolver routine
!!    Search for
!!
!! @ingroup PSOLVER
module Poisson_Solver
   use wrapper_linalg
   use wrapper_MPI
   use dynamic_memory
   use time_profiling, only: TIMING_UNINITIALIZED
   !use m_profiling
   ! TO BE REMOVED with f_malloc
   
   implicit none
   
   private
   
   ! General precision, density and the potential types
   integer, parameter, public :: gp=kind(1.0d0)  !< general-type precision
   integer, parameter, public :: dp=kind(1.0d0)  !< density-type precision
   integer, parameter, public :: wp=kind(1.0d0)  !< potential-type precision
   ! Associated MPI precisions.
   integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
   integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
   integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION

   !timing categories
   integer, public, save :: TCAT_PSOLV_COMPUT=TIMING_UNINITIALIZED
   integer, public, save :: TCAT_PSOLV_COMMUN=TIMING_UNINITIALIZED
   integer, public, save :: TCAT_PSOLV_KERNEL=TIMING_UNINITIALIZED
   
   include 'configure.inc'
   
   !> Defines the fundamental structure for the kernel
   type, public :: coulomb_operator
      !variables with physical meaning
      integer :: itype_scf             !< Order of the ISF family to be used
      real(gp) :: mu                   !< Inverse screening length for the Helmholtz Eq. (Poisson Eq. -> mu=0)
       !> geocode is used in all the code to specify the boundary conditions (BC) the problem:
       !!          - 'F' free BC, isolated systems.
       !!                The program calculates the solution as if the given density is
       !!                "alone" in R^3 space.
       !!          - 'S' surface BC, isolated in y direction, periodic in xz plane                
       !!                The given density is supposed to be periodic in the xz plane,
       !!                so the dimensions in these direction mus be compatible with the FFT
       !!                Beware of the fact that the isolated direction is y!
       !!          - 'P' periodic BC.
       !!                The density is supposed to be periodic in all the three directions,
       !!                then all the dimensions must be compatible with the FFT.
       !!                No need for setting up the kernel (in principle for Plane Waves)
       !!          - 'W' Wires BC.
       !!                The density is supposed to be periodic in z direction, 
       !!                which has to be compatible with the FFT.
       !!          - 'H' Helmholtz Equation Solver
      character(len=1) :: geocode
      integer, dimension(3) :: ndims   !< dimension of the box of the density
      real(gp), dimension(3) :: hgrids !<grid spacings in each direction
      real(gp), dimension(3) :: angrad !< angles in radiants between each of the axis
      real(dp), dimension(:), pointer :: kernel !< kernel of the Poisson Solver
      real(dp) :: work1_GPU,work2_GPU,k_GPU
      integer, dimension(5) :: plan
      integer, dimension(3) :: geo
      !variables with computational meaning
      type(mpi_environment) :: mpi_env !< complete environment for the POisson Solver
      type(mpi_environment) :: inplane_mpi,part_mpi !<mpi_environment for internal ini-plane parallelization
      integer :: igpu !< control the usage of the GPU
      integer :: initCufftPlan
      integer :: keepGPUmemory
   end type coulomb_operator

   !intialization of the timings
   public :: PS_initialize_timing_categories
   ! Calculate the allocation dimensions
   public :: PS_dim4allocation, PS_getVersion
   ! Routine that creates the kernel
   public :: pkernel_init, pkernel_set, pkernel_free
   ! Calculate the poisson solver
   public :: H_potential 
   ! Calculate the allocation dimensions
   public :: P_FFT_dimensions, S_FFT_dimensions, F_FFT_dimensions, W_FFT_dimensions, xc_dimensions

   !> This structure is used to indicate the arguments of the routine which are used commonly
   !! Doxygen will duplicate the documentation for the arguments
   type doc
      character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode 
                                  !! @ingroup RESERVED
      !> Indicates the distribution of the data of the input/output array:
      !!    - 'G' global data. Each process has the whole array of the density 
      !!          which will be overwritten with the whole array of the potential.
      !!    - 'D' distributed data. Each process has only the needed part of the density
      !!          and of the potential. The data distribution is such that each processor
      !!          has the xy planes needed for the calculation AND for the evaluation of the 
      !!          gradient, needed for XC part, and for the White-Bird correction, which
      !!          may lead up to 8 planes more on each side. Due to this fact, the information
      !!          between the processors may overlap.
      !!          @ingroup RESERVED
      character(len=1) :: datacode
   end type doc

contains

  pure function pkernel_null() result(k)
    type(coulomb_operator) :: k
    k%itype_scf=0
    k%geocode='F'
    k%mu=0.0_gp
    k%ndims=(/0,0,0/)
    k%hgrids=(/0.0_gp,0.0_gp,0.0_gp/)
    k%angrad=(/0.0_gp,0.0_gp,0.0_gp/)
    nullify(k%kernel)
    k%work1_GPU=0.d0
    k%work2_GPU=0.d0
    k%k_GPU=0.d0
    k%plan=(/0,0,0,0,0/)
    k%geo=(/0,0,0/)
    k%mpi_env=mpi_environment_null()
    k%inplane_mpi=mpi_environment_null()
    k%part_mpi=mpi_environment_null()
    k%igpu=0
    k%initCufftPlan=0
    k%keepGPUmemory=1
  end function pkernel_null

  !> switch on the timing categories for the Poisson Solver
  !! shuold be called if the time_profiling module has to be used for profiling the routines
  subroutine PS_initialize_timing_categories()
    use time_profiling, only: f_timing_category_group,f_timing_category
    use wrapper_mpi, only: comm => tgrp_mpi_name, mpi_initialize_timing_categories
    use wrapper_linalg, only: linalg_initialize_timing_categories
    implicit none
    character(len=*), parameter :: pscpt='PS Computation'

    call mpi_initialize_timing_categories()

    call linalg_initialize_timing_categories()
    !group of Poisson Solver operations, separate category
    call f_timing_category_group(pscpt,&
         'Computations of Poisson Solver operations')

  !define the timing categories
  call f_timing_category('PSolver Computation',pscpt,&
       '3D SG_FFT and related operations',&
       TCAT_PSOLV_COMPUT)
  call f_timing_category('PSolver Kernel Creation',pscpt,&
       'ISF operations and creation of the kernel',&
       TCAT_PSOLV_KERNEL)
  call f_timing_category('PSolver Communication',comm,&
       'MPI_ALLTOALL and MPI_ALLGATHERV',&
       TCAT_PSOLV_COMMUN)

  end subroutine PS_initialize_timing_categories

  function PS_getVersion() result(str)
    character(len = 128) :: str

    write(str, "(A)") package_version
  end function PS_getVersion

  include 'PSolver_Main.f90'
  include 'createKernel.f90'
end module Poisson_Solver
