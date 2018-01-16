!> @file
!!    Define the module for Poisson Solver
!!
!! @author
!!    Luigi Genovese (February 2007)<br/>
!!    PSolverNC added by Anders Bergman, March 2008<br/>
!!    Copyright (C) 2002-2017 BigDFT group<br/>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module used by the Poisson Solver library.
!! It must be used in the parent routine.
!! @details
!!    In the main routine in which the Poisson Solver is called:
!!
!!    1. The Poisson kernel must be declared as a pointer, then the
!!       routine createKernel can be called. On exit, the kernel will be allocated and
!!       ready to use. See the documentation of the createKernel routine for more details
!!    2. The correct sizes for allocating the density/potential and the pot_ion arrays
!!       are given from the routine PS_dim4allocation (see routine documentation for details).
!!       Its argument MUST be in agreement with the arguments of the PSolver routine.
!!       WARNING: No cross-check of the arguments is performed!
!!    3. The PSolver routine can then be called. On exit, the Hartree potential is computed
!!       and summed (following ixc value) to XC and external potential.
!!       The input array is overwritten. Again, see routine documentation for details.
!!    4. QUICK INSTRUCTION FOR THE IMPATIENT:If you want to use the Poisson Solver in the
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
!!       Any other changement of the arguments require reading of the documentation.
!!       See documentations of the Public routines
!!
!! @warning
!!    This module REQUIRES the module of XC functional from ABINIT, defs_xc, which
!!    require defs_basis and defs_datatypes.
!!    Such routines are provided inside the abinit directory of this bundle.
!!    They are based on XC functionals present in ABINIT 5.x
!!    If you want to use this Poisson Solver without the XC functionals, you can comment out
!!    the XC part in the PSolver routine
!!    Search for
!!
module Poisson_Solver
   use dictionaries, only: f_err_throw
   use f_utils
   use f_enums
   use PSbase
   use wrapper_linalg
   use wrapper_MPI
   use dynamic_memory
   use time_profiling, only: TIMING_UNINITIALIZED, f_timing
   use yaml_output
   use yaml_strings
   use psolver_environment
   use PStypes
   use PSbox
   !use m_profiling
   ! TO BE REMOVED with f_malloc

   implicit none

   private

   ! Associated MPI precisions.
   integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
   integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
   integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION

   ! Timing categories
   integer, public, save :: TCAT_PSOLV_COMPUT=TIMING_UNINITIALIZED
   integer, public, save :: TCAT_PSOLV_COMMUN=TIMING_UNINITIALIZED
   integer, public, save :: TCAT_PSOLV_KERNEL=TIMING_UNINITIALIZED

   include 'configure.inc'

   ! Intialization of the timings
   public :: PS_initialize_timing_categories,coulomb_operator,PSolver_energies
   ! Calculate the allocation dimensions
   public :: PS_dim4allocation,PSolver_logo,ps_soft_PCM_forces
   ! Routine that creates the kernel
   public :: pkernel_init, pkernel_set, pkernel_free, pkernel_set_epsilon, pkernel_allocate_cavity,pkernel_get_radius
   ! Calculate the poisson solver
   public :: H_potential,Electrostatic_Solver,PS_set_options
   ! Calculate the allocation dimensions
   public :: P_FFT_dimensions, S_FFT_dimensions, F_FFT_dimensions, W_FFT_dimensions
   public :: dp,gp,PS_dump_coulomb_operator, xc_dimensions

   !> This structure is used to indicate the arguments of the routine which are used commonly
   !! Doxygen will duplicate the documentation for the arguments
   type doc
      character(len=1) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
      !> Indicates the distribution of the data of the input/output array:
      !!    - 'G' global data. Each process has the whole array of the density
      !!          which will be overwritten with the whole array of the potential.
      !!    - 'D' distributed data. Each process has only the needed part of the density
      !!          and of the potential. The data distribution is such that each processor
      !!          has the xy planes needed for the calculation AND for the evaluation of the
      !!          gradient, needed for XC part, and for the White-Bird correction, which
      !!          may lead up to 8 planes more on each side. Due to this fact, the information
      !!          between the processors may overlap.
      character(len=1) :: datacode
   end type doc

contains

  !> Switch on the timing categories for the Poisson Solver
  !! should be called if the time_profiling module has to be used for profiling the routines
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

  subroutine PSolver_logo()
    implicit none
    call yaml_map('Reference Paper','The Journal of Chemical Physics 137, 134108 (2012)')
    call yaml_map('Version Number', "PSolver " // trim(PS_getVersion()))
    call yaml_map('Timestamp of this run',yaml_date_and_time_toa())
      call yaml_map('Root process Hostname',mpihostname())
  end subroutine PSolver_logo
 

  include 'PSolver_Main.f90'

  include 'createKernel.f90'

end module Poisson_Solver
