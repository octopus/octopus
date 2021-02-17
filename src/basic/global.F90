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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module global_oct_m
  use hardware_oct_m
  use loct_oct_m
  use mpi_oct_m
  use varinfo_oct_m
#ifdef HAVE_OPENMP
  use omp_lib
#endif

  implicit none

  private

  ! ---------------------------------------------------------
  !> Public types, variables and procedures.
  public ::           &
    conf_t,           &
    global_init,      &
    global_end,       & 
    optional_default, &
    assert_die,       &
    not_in_openmp,    &
    operator(+),      &
    bitand


  integer, public, parameter :: MAX_PATH_LEN=256
  integer, public, parameter :: MAX_OUTPUT_TYPES=40

  type conf_t
    ! Components are public by default
    logical :: devel_version !< If true then allow unstable parts of the code
    logical :: report_memory
    character(len=256) :: share       !< Name of the share dir
    character(len=256) :: git_commit  !< hash of latest git commit
    character(len=50)  :: build_time  !< time octopus was compiled
    character(len=20)  :: version     !< version number
    character(len=256) :: cc
    character(len=256) :: cflags
    character(len=256) :: cxx
    character(len=256) :: cxxflags
    character(len=256) :: fc
    character(len=256) :: fcflags
    integer            :: target_states_block_size
  end type conf_t

  type(conf_t),      public :: conf
   
  FLOAT, public, parameter :: R_SMALL = CNST(0.0001)

  !> some mathematical constants
  FLOAT, public, parameter :: M_Pi        = CNST(3.1415926535897932384626433832795029)
  FLOAT, public, parameter :: M_E         = CNST(2.7182818284590452353602874713526625)  
  FLOAT, public, parameter :: M_ZERO      = CNST(0.0)
  FLOAT, public, parameter :: M_ONE       = CNST(1.0)
  FLOAT, public, parameter :: M_TWO       = CNST(2.0)
  FLOAT, public, parameter :: M_THREE     = CNST(3.0)
  FLOAT, public, parameter :: M_FOUR      = CNST(4.0)
  FLOAT, public, parameter :: M_FIVE      = CNST(5.0)
  FLOAT, public, parameter :: M_HALF      = CNST(0.5)
  FLOAT, public, parameter :: M_THIRD     = M_ONE/M_THREE
  FLOAT, public, parameter :: M_TWOTHIRD  = M_TWO/M_THREE
  FLOAT, public, parameter :: M_FOURTH    = M_ONE/M_FOUR
  CMPLX, public, parameter :: M_z0        = (CNST(0.0), CNST(0.0))
  CMPLX, public, parameter :: M_z1        = (CNST(1.0), CNST(0.0))
  CMPLX, public, parameter :: M_z2        = (CNST(2.0), CNST(0.0))
  CMPLX, public, parameter :: M_z2I       = (CNST(0.0), CNST(2.0))
  CMPLX, public, parameter :: M_zI        = (CNST(0.0), CNST(1.0))

  FLOAT, public, parameter :: M_EPSILON   =  epsilon(M_ONE)
  FLOAT, public, parameter :: M_TINY      =  tiny(M_ONE)
  FLOAT, public, parameter :: M_HUGE      =  huge(M_ONE)

  !> some physical constants
  FLOAT, public, parameter :: P_a_B =  CNST(0.52917720859)
  FLOAT, public, parameter :: P_Ang =  M_ONE / P_a_B
  FLOAT, public, parameter :: P_Ry  =  CNST(13.60569193)
  FLOAT, public, parameter :: P_eV  =  M_ONE / P_Ry
  FLOAT, public, parameter :: P_Kb  =  CNST(8.617343e-5)/(M_TWO*P_Ry)  !< Boltzmann constant in Ha/K
  FLOAT, public, parameter :: P_c   =  CNST(137.035999679)
  FLOAT, public, parameter :: P_g   =  CNST(2.0023193043768)   !< Electron gyromagnetic ratio
  FLOAT, public, parameter :: P_PROTON_CHARGE = CNST(-1.0)
  FLOAT, public, parameter :: P_ep  =  M_ONE/(M_FOUR*M_Pi)
  FLOAT, public, parameter :: P_mu  =  M_FOUR*M_PI/(P_c**2)

  !> the standard input and output
  integer, public :: stderr, stdin, stdout

  !> global epoch time (time at startup)
  integer, public :: s_epoch_sec, s_epoch_usec

  !> The stack.
  character(len=80), public          :: sub_stack(50)
  FLOAT, public                      :: time_stack(50)
  integer, public                    :: no_sub_stack = 0

  !> Same for profiling mode.
  logical, public :: in_profiling_mode = .false.

  integer,    public :: global_alloc_err
  integer(8), public :: global_sizeof
  
  ! The code directories should be defined here, and not hard coded in the Fortran files.
  character(len=*), public, parameter :: GS_DIR = "gs/"
  character(len=*), public, parameter :: TD_DIR = "td/"
  character(len=*), public, parameter :: STATIC_DIR = "static/"
  character(len=*), public, parameter :: EM_RESP_DIR = "em_resp/"
  character(len=*), public, parameter :: EM_RESP_FD_DIR = "em_resp_fd/"
  character(len=*), public, parameter :: KDOTP_DIR = "kdotp/"
  character(len=*), public, parameter :: VIB_MODES_DIR = "vib_modes/"
  character(len=*), public, parameter :: VDW_DIR = "vdw/"
  character(len=*), public, parameter :: CASIDA_DIR = "casida/"
  character(len=*), public, parameter :: OCT_DIR = "opt-control/"
  character(len=*), public, parameter :: PCM_DIR = "pcm/"
  character(len=*), public, parameter :: PARTITION_DIR = "partition/"

  ! End of declaration of public objects.
  ! ---------------------------------------------------------

  interface optional_default
    module procedure doptional_default, zoptional_default, ioptional_default, loptional_default, soptional_default
  end interface optional_default


  !> This function is defined in messages.F90
  interface 
    subroutine assert_die(s, f, l)
      implicit none
      character(len=*), intent(in) :: s, f
      integer, intent(in) :: l
    end subroutine assert_die
  end interface

  interface operator (+)
    module procedure cat
  end interface operator (+)

  interface bitand
    module procedure bitand48
    module procedure bitand84
    module procedure bitand88
    module procedure bitand44
  end interface bitand
  
contains

  ! ---------------------------------------------------------
  subroutine global_init(is_serial)
    logical, optional, intent(in) :: is_serial !< if .true., do not call MPI_Init

    character(len=256) :: share

    ! initialize mpi
    call mpi_mod_init(optional_default(is_serial, .false.))

    ! Get epoch time at node startup, just after the barrier to synchronize nodes first.
    call loct_gettimeofday(s_epoch_sec, s_epoch_usec)

    call hardware_init()

    ! Get the environment variable OCTOPUS_SHARE that overrides SHARE_DIR/share/octopus.
    call loct_getenv("OCTOPUS_SHARE", share)

    if(share /= "") then
      conf%share = trim(share)
    else
      conf%share = SHARE_DIR
    end if
    conf%git_commit = GIT_COMMIT
    conf%build_time = BUILD_TIME
    conf%version    = PACKAGE_VERSION
    conf%cc         = CC
    ! not indented to have the whole line in case it is long
    conf%cflags     = &
CFLAGS
    conf%cxx        = CXX
    conf%cxxflags   = &
CXXFLAGS
    conf%fc         = FC
    conf%fcflags    = &
FCFLAGS

    ! initialize info for the input variables
    call varinfo_init(trim(conf%share)//'/varinfo')

    conf%target_states_block_size = -1

  end subroutine global_init


  ! ---------------------------------------------------------
  subroutine global_end()

    call hardware_end()
    call mpi_mod_end()
    call varinfo_end()

  end subroutine global_end

  !----------------------------------------------------------

  FLOAT pure function doptional_default(opt, def) result(val)
    FLOAT, optional, intent(in) :: opt
    FLOAT,           intent(in) :: def
    
    val = def
    if(present(opt)) val = opt
  end function doptional_default

  !----------------------------------------------------------

  CMPLX pure function zoptional_default(opt, def) result(val)
    CMPLX, optional, intent(in) :: opt
    CMPLX,           intent(in) :: def
    
    val = def
    if(present(opt)) val = opt
  end function zoptional_default

  !----------------------------------------------------------

  integer pure function ioptional_default(opt, def) result(val)
    integer, optional, intent(in) :: opt
    integer,           intent(in) :: def
    
    val = def
    if(present(opt)) val = opt
  end function ioptional_default
  
  !----------------------------------------------------------

  logical pure function loptional_default(opt, def) result(val)
    logical, optional, intent(in) :: opt
    logical,           intent(in) :: def
    
    val = def
    if(present(opt)) val = opt
  end function loptional_default

  !----------------------------------------------------------

  character(len=80) pure function soptional_default(opt, def) result(val)
    character(len=*), optional, intent(in) :: opt
    character(len=*),           intent(in) :: def

    val = def
    if(present(opt)) val = opt
  end function soptional_default

  !-----------------------------------------------------------

  logical &
#ifndef HAVE_OPENMP
    pure & 
#endif
    function not_in_openmp()
    
#ifdef HAVE_OPENMP
    not_in_openmp = .not. omp_in_parallel()
#else
    not_in_openmp = .true.
#endif

  end function not_in_openmp

  !-----------------------------------------------------------

  function cat(str1, str2)
    character(len=*), intent(in) :: str1
    character(len=*), intent(in) :: str2

    character(len=len(str1) + len(str2)) :: cat
    cat = str1//str2
    
  end function cat

  ! -----------------------------------------------------------

  integer(8) pure function bitand48(val1, val2)
    integer(4), intent(in) :: val1
    integer(8), intent(in) :: val2

    bitand48 = iand(int(val1, 8), val2)
    
  end function bitand48

  ! -----------------------------------------------------------

  integer(8) pure function bitand84(val1, val2)
    integer(8), intent(in) :: val1
    integer(4), intent(in) :: val2

    bitand84 = iand(val1, int(val2, 8))
        
  end function bitand84

  ! -----------------------------------------------------------

  integer(8) pure function bitand88(val1, val2)
    integer(8), intent(in) :: val1
    integer(8), intent(in) :: val2
    
    bitand88 = iand(val1, val2)
    
  end function bitand88

  ! -----------------------------------------------------------

  integer(4) pure function bitand44(val1, val2)
    integer(4), intent(in) :: val1
    integer(4), intent(in) :: val2
    
    bitand44 = iand(val1, val2)
    
  end function bitand44

end module global_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
