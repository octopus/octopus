!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module global
  use varinfo
  use string
  use lib_oct
  use mpi_mod

  implicit none

  private

  ! ---------------------------------------------------------
  ! Public types, variables and procedures.
  public ::         &
    conf_type,      &
    global_init,    &
    global_end,     &
    assert_die

  type conf_type
    integer :: debug_level   ! How much debug should print
    logical :: devel_version ! If true then allow unstable parts of the code
    character(len=256) :: share       ! Name of the share dir
    character(len=256) :: latest_cvs  ! rcs info of latest cvs commit
    character(len=50)  :: build_time  ! time octopus was compiled
    character(len=10)  :: version     ! version number
    character(len=256) :: compiler
    character(len=256) :: fcflags
  end type conf_type

  type(conf_type),      public :: conf

  ! the kinds used in the program
  integer, public, parameter  ::  r8 = selected_real_kind(12,256)
  integer, public, parameter  ::  r4 = selected_real_kind(6,37)
  integer, public, parameter  ::  i4 = selected_int_kind(9)
  integer, public, parameter  ::  i2 = selected_int_kind(3)

  FLOAT, public, parameter :: r_small = CNST(0.0001)

  ! some mathematical constants
  FLOAT, public, parameter :: M_Pi        = CNST(3.141592653589793)
  FLOAT, public, parameter :: M_ZERO      = CNST(0.0)
  FLOAT, public, parameter :: M_ONE       = CNST(1.0)
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
  CMPLX, public, parameter :: M_z2I       = (CNST(0.0), CNST(2.0))
  CMPLX, public, parameter :: M_zI        = (CNST(0.0), CNST(1.0))


  ! some physical constants
  FLOAT, public, parameter :: P_Ang =  M_ONE / CNST(0.529177)
  FLOAT, public, parameter :: P_eV  =  M_ONE / CNST(13.60580)
  FLOAT, public, parameter :: P_E   =  CNST(3.79470065)        ! (electron charge)
  FLOAT, public, parameter :: P_M   =  CNST(0.13123)           ! (electrons mass)
  FLOAT, public, parameter :: P_E2  =  CNST(14.399753)         ! = 2/(Ang*eV)
  FLOAT, public, parameter :: P_a_B =  CNST(0.529177)
  FLOAT, public, parameter :: P_Ry  =  CNST(13.6058)
  FLOAT, public, parameter :: P_Kb  =  CNST(3.166815104e-6)    ! Boltzmann constant in Ha/K
  FLOAT, public, parameter :: P_c   =  CNST(137.036)
  FLOAT, public, parameter :: P_g   =  CNST(2.0023193043768)   ! Electron gyromagnetic ratio

  integer, public  :: calc_mode
  integer, public  :: calc_dim

  ! the standard input and output
  integer, public :: stderr, stdin, stdout

  ! global epoch time (time at startup)
  integer, public :: s_epoch_sec, s_epoch_usec

  ! The stack.
  character(len=40), public          :: sub_stack(50)
  FLOAT, public                      :: time_stack(50)
  integer, public                    :: no_sub_stack = 0

  ! should we run in debug mode
  logical, public :: in_debug_mode = .false.
  ! Same for profiling mode.
  logical, public :: in_profiling_mode = .false.

  ! End of declaration of public objects.
  ! ---------------------------------------------------------

contains

  ! ---------------------------------------------------------
  subroutine global_init()

    ! initialize mpi
    call mpi_mod_init()

    ! Get epoch time at node startup, just after the barrier to synchronize nodes first.
    call loct_gettimeofday(s_epoch_sec, s_epoch_usec)

    conf%share      = SHARE_OCTOPUS
    conf%latest_cvs = LATEST_CVS
    conf%build_time = BUILD_TIME
    conf%version    = OCTOPUS_VERSION
    conf%compiler   = FC
#if defined (LONG_LINES)
    conf%fcflags    = FCFLAGS
#else
    conf%fcflags = "No flags information available."
#endif

    ! initialize info for the input variables
    call varinfo_init(trim(conf%share)//'/varinfo');

  end subroutine global_init


  ! ---------------------------------------------------------
  subroutine global_end()

    call mpi_mod_end()
    call varinfo_end()

  end subroutine global_end


  ! ---------------------------------------------------------
  ! This subroutine is called by the assert macro
  ! ---------------------------------------------------------
  subroutine assert_die(s, f, l)
    character(len=*), intent(in) :: s, f
    integer, intent(in) :: l

    write(stderr, '(3a,i5,3a)') 'Assertion "', trim(s), '" failed in line ', l, ' in file "', f, '"'
    stop
  end subroutine assert_die

end module global
