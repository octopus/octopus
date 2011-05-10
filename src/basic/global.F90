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

#include "global.h"

module global_m
  use hardware_m
  use loct_m
  use mpi_m
  use varinfo_m

  implicit none

  private

  ! ---------------------------------------------------------
  ! Public types, variables and procedures.
  public ::          &
    conf_t,          &
    global_init,     &
    global_end,      & 
    assert_die,      &
    optional_default

  type conf_t
    integer :: debug_level   !< How much debug should print
    logical :: devel_version !< If true then allow unstable parts of the code
    logical :: report_memory
    character(len=256) :: share       !< Name of the share dir
    character(len=256) :: latest_svn  !< rcs info of latest svn commit
    character(len=50)  :: build_time  !< time octopus was compiled
    character(len=20)  :: version     !< version number
    character(len=256) :: cc
    character(len=256) :: cflags
    character(len=256) :: fc
    character(len=256) :: fcflags
  end type conf_t

  type(conf_t),      public :: conf

  ! the kinds used in the program
  integer, public, parameter  ::  r8 = selected_real_kind(12,256)
  integer, public, parameter  ::  r4 = selected_real_kind(6,37)
  integer, public, parameter  ::  i4 = selected_int_kind(9)
  integer, public, parameter  ::  i2 = selected_int_kind(3)

  FLOAT, public, parameter :: r_small = CNST(0.0001)

  ! some mathematical constants
  FLOAT, public, parameter :: M_Pi        = CNST(3.1415926535897932384626433832795029)
  FLOAT, public, parameter :: M_E         = CNST(2.7182818284590452353602874713526625)  
  FLOAT, public, parameter :: M_TWO       = CNST(2.0)
  FLOAT, public, parameter :: M_THREE     = CNST(3.0)
  FLOAT, public, parameter :: M_FOUR      = CNST(4.0)
  FLOAT, public, parameter :: M_FIVE      = CNST(5.0)
  FLOAT, public, parameter :: M_SIX       = CNST(6.0)
  FLOAT, public, parameter :: M_SEVEN     = CNST(7.0)
  FLOAT, public, parameter :: M_EIGHT     = CNST(8.0)
  FLOAT, public, parameter :: M_NINE      = CNST(9.0)
  FLOAT, public, parameter :: M_TEN       = CNST(10.0)
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
  FLOAT, public, parameter :: M_HUGE      =  huge(M_ONE)

  ! some physical constants
  FLOAT, public, parameter :: P_a_B =  CNST(0.52917720859)
  FLOAT, public, parameter :: P_Ang =  M_ONE / P_a_B
  FLOAT, public, parameter :: P_Ry  =  CNST(13.60569193)
  FLOAT, public, parameter :: P_eV  =  M_ONE / P_Ry
  FLOAT, public, parameter :: P_Kb  =  CNST(8.617343e-5)/(M_TWO*P_Ry)  ! Boltzmann constant in Ha/K
  FLOAT, public, parameter :: P_c   =  CNST(137.035999679)
  FLOAT, public, parameter :: P_g   =  CNST(2.0023193043768)   ! Electron gyromagnetic ratio
  FLOAT, public, parameter :: P_PROTON_CHARGE = CNST(-1.0)

  integer, public  :: ifinal_sete = 0

  ! the standard input and output
  integer, public :: stderr, stdin, stdout
  
  ! used to store return values of mpi calls
  integer, public :: mpi_err  

  ! global epoch time (time at startup)
  integer, public :: s_epoch_sec, s_epoch_usec

  ! The stack.
  character(len=80), public          :: sub_stack(50)
  real(8), public                    :: time_stack(50)
  integer, public                    :: no_sub_stack = 0

  ! should we run in debug mode
  logical, public :: in_debug_mode = .false.
  ! Same for profiling mode.
  logical, public :: in_profiling_mode = .false.

  ! End of declaration of public objects.
  ! ---------------------------------------------------------

  interface optional_default
    module procedure doptional_default, zoptional_default, ioptional_default, loptional_default, soptional_default
  end interface optional_default

contains

  ! ---------------------------------------------------------
  subroutine global_init()

    character(len=256) :: share

    ! initialize mpi
    call mpi_mod_init()

    ! Get epoch time at node startup, just after the barrier to synchronize nodes first.
    call loct_gettimeofday(s_epoch_sec, s_epoch_usec)

    call hardware_init()

    ! Get the environment variable OCTOPUS_SHARE that overrides SHARE_OCTOPUS/share/octopus.
    call loct_getenv("OCTOPUS_SHARE", share)

    if(share.ne."") then
      conf%share = trim(share)
    else
      conf%share = SHARE_OCTOPUS
    end if
    conf%latest_svn = LATEST_SVN
    conf%build_time = BUILD_TIME
    conf%version    = PACKAGE_VERSION
    conf%cc         = CC
#if defined (LONG_LINES)
    conf%cflags     = CFLAGS
#else
    conf%cflags     = "No flags information available."
#endif
    conf%fc         = FC
#if defined (LONG_LINES)
    conf%fcflags    = FCFLAGS
#else
    conf%fcflags    = "No flags information available."
#endif

    ! initialize info for the input variables
    call varinfo_init(trim(conf%share)//'/varinfo');

  end subroutine global_init


  ! ---------------------------------------------------------
  subroutine global_end()

    call hardware_end()
    call mpi_mod_end()
    call varinfo_end()

  end subroutine global_end


  ! ---------------------------------------------------------
  ! This subroutine is called by the assert macro
  ! ---------------------------------------------------------
  subroutine assert_die(s, f, l)
    character(len=*), intent(in) :: s, f
    integer, intent(in) :: l

    write(stderr, '(a,i5,3a,i5,3a)') 'Node ', mpi_world%rank, ': Assertion "', trim(s), '" failed in line ', l, ' in file "', f, '"'
    stop
  end subroutine assert_die

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

  pure function soptional_default(opt, def) result(val)
    character(len=*), optional, intent(in) :: opt
    character(len=*),           intent(in) :: def
    character(len=max(len(opt),len(def)))  :: val

    val = def
    if(present(opt)) val = opt
  end function soptional_default

end module global_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
