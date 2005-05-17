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
  use lib_oct_parser
  use lib_oct
  use io

#if defined(HAVE_MPI) && !defined(MPI_H)
  use mpi
#endif

  implicit none

#if defined(HAVE_MPI) && defined(MPI_H)
# include "mpif.h"
#endif

type conf_type
  integer :: verbose     ! <= 0  -> silent, no output except fatal errors
                         ! > 0   -> warning only
                         ! > 20  -> normal program info
                         ! > 999 -> debug
  integer :: debug_level ! How much debug should print

  character(len=256) :: share       ! Name of the share dir
  character(len=50)  :: build_time  ! time octopus was compiled
  character(len=10)  :: version     ! version number

  integer :: periodic_dim
  integer :: dim
  logical :: boundary_zero_derivative
  logical :: only_user_def
end type conf_type

type mpi_type
  integer :: numprocs ! how many are we
  integer :: node ! who am I
end type mpi_type

type(mpi_type) :: mpiv
type(conf_type) :: conf

! the kinds used in the program
integer, parameter  ::  r8 = selected_real_kind(12,256)
integer, parameter  ::  r4 = selected_real_kind(6,37)
integer, parameter  ::  i4 = selected_int_kind(9)
integer, parameter  ::  i2 = selected_int_kind(3)

FLOAT, parameter :: r_small = CNST(0.0001)

! some mathematical constants
FLOAT, parameter :: M_Pi        = CNST(3.141592653589793)
FLOAT, parameter :: M_ZERO      = CNST(0.0)
FLOAT, parameter :: M_ONE       = CNST(1.0)
FLOAT, parameter :: M_TWO       = CNST(2.0)
FLOAT, parameter :: M_THREE     = CNST(3.0)
FLOAT, parameter :: M_FOUR      = CNST(4.0)
FLOAT, parameter :: M_FIVE      = CNST(5.0)
FLOAT, parameter :: M_SIX       = CNST(6.0)
FLOAT, parameter :: M_SEVEN     = CNST(7.0)
FLOAT, parameter :: M_EIGHT     = CNST(8.0)
FLOAT, parameter :: M_NINE      = CNST(9.0)
FLOAT, parameter :: M_TEN       = CNST(10.0)
FLOAT, parameter :: M_HALF      = CNST(0.5)
FLOAT, parameter :: M_THIRD     = M_ONE/M_THREE
FLOAT, parameter :: M_TWOTHIRD  = M_TWO/M_THREE
FLOAT, parameter :: M_FOURTH    = M_ONE/M_FOUR
CMPLX, parameter :: M_z0        = (CNST(0.0), CNST(0.0))
CMPLX, parameter :: M_z1        = (CNST(1.0), CNST(0.0))
CMPLX, parameter :: M_z2I       = (CNST(0.0), CNST(2.0))
CMPLX, parameter :: M_zI        = (CNST(0.0), CNST(1.0))


! some physical constants
FLOAT, parameter :: P_Ang =  M_ONE / CNST(0.529177)
FLOAT, parameter :: P_eV  =  M_ONE / CNST(13.60580)
FLOAT, parameter :: P_E   =  CNST(3.79470065)        ! (electron charge)
FLOAT, parameter :: P_M   =  CNST(0.13123)           ! (electrons mass)
FLOAT, parameter :: P_E2  =  CNST(14.399753)         ! = 2/(Ang*eV)
FLOAT, parameter :: P_a_B =  CNST(0.529177)
FLOAT, parameter :: P_Ry  =  CNST(13.6058)
FLOAT, parameter :: P_Kb  =  CNST(3.166815104e-6)    ! Boltzmann constant in Ha/K
FLOAT, parameter :: P_c   =  CNST(137.036)

integer, parameter :: &
   VERBOSE_QUIET   = -9999, &
   VERBOSE_WARNING = 0,     &
   VERBOSE_NORMAL  = 30,    &
   VERBOSE_DEBUG   = 999

integer :: calc_mode

! some variables to be used everywhere
character(len=256), dimension(20) :: message ! to be output by fatal, warning

! variables to treat multi subsytems
character(len=32), allocatable :: subsys_label(:)
integer,           allocatable :: subsys_runmode(:), subsys_run_order(:)
integer, parameter             :: multi_subsys_mode = 1000
integer                        :: current_subsystem = 1
integer                        :: no_syslabels, no_subsys_runmodes
character(len=32)              :: current_label, tmpdir

! some private variables to this module
#ifdef DEBUG
character(len=40), private :: sub_stack(50)
FLOAT, private             :: time_stack(50)
integer, private           :: no_sub_stack = 0
#endif

character(len=70), parameter :: stars =  &
    '**********************************************************************'
character(len=68), parameter :: hyphens = &
    '--------------------------------------------------------------------'

contains

subroutine global_init()
  integer :: ierr

#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiv%node, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiv%numprocs, ierr)
  write(stdout,'(a,i4,a,i4,a)') 'Process ', mpiv%node, ' of ', mpiv%numprocs, ' is alive'  
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  mpiv%node = 0
  mpiv%numprocs = 1
#endif

  conf%share      = SHARE_OCTOPUS
  conf%build_time = BUILD_TIME
  conf%version    = OCTOPUS_VERSION

  ! initialize the parser
  ierr = loct_parse_init('out.oct')
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    message(2) = "Do you have write permissions in this directory?"
    call write_fatal(2)
  end if

  ! read in default variables
  ierr = loct_parse_input(trim(conf%share)//'/variables')

  ! setup standard input
  ierr = loct_parse_input('inp')
  if(ierr .ne. 0) then
    ierr = loct_parse_input("-")
    if(ierr .ne. 0) then
      message(1) = "Error initializing liboct"
      message(2) = "Can not open input file or standard input!"
      call write_fatal(2)
    end if
  end if

  ! initialize input/output system
  call io_init()

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems
  call loct_parse_int('CalculationMode', 1, calc_mode)
  if(calc_mode == multi_subsys_mode) then
     call read_system_labels()
  else
     call init_default_system_labels()
  endif

  ! verbosity level
  call loct_parse_int(check_inp('Verbose'), VERBOSE_NORMAL, conf%verbose)
  if(conf%verbose > VERBOSE_DEBUG .and. mpiv%node == 0) then
    call loct_parse_int(check_inp('DebugLevel'), 3, conf%debug_level)
    message(1) = 'Entering DEBUG mode'
    call write_warning(1)
  end if

  ! Sets the dimensionaliy of the problem.
  call loct_parse_int(check_inp('Dimensions'), 3, conf%dim)
  if(conf%dim<1 .or. conf%dim>3) then
    message(1) = 'Dimensions must be either 1, 2, or 3'
    call write_fatal(1)
  end if

  ! handle periodic directions
  call loct_parse_int(check_inp('PeriodicDimensions'), 0, conf%periodic_dim)
  if ((conf%periodic_dim < 0) .or. (conf%periodic_dim > 3)) then
    message(1) = 'Periodic dimensions must be either 0, 1, 2, or 3'
    call write_fatal(1)
  endif
  if(conf%periodic_dim > conf%dim) then
    message(1) = 'PeriodicDimensions must be <= Dimensions'
    call write_fatal(1)
  end if

  call loct_parse_logical(check_inp('BoundaryZeroDerivative'), .false., conf%boundary_zero_derivative)

end subroutine global_init

subroutine global_end()
#ifdef HAVE_MPI
  integer :: ierr
  call MPI_FINALIZE(ierr)
#endif

  call io_end()
  call loct_parse_end()
  
end subroutine global_end

! This subroutine is called by the assert macro
#ifdef DEBUG
subroutine assert_die(s, f, l)
  character(len=*), intent(in) :: s, f
  integer, intent(in) :: l
  
  write(stderr, '(3a,i5,3a)') 'Assertion "', trim(s), '" failed in line ', l, ' in file "', f, '"'
  stop
end subroutine assert_die
#endif

subroutine write_fatal(no_lines)
  integer, intent(in) :: no_lines
  integer :: i

  write(stderr, '(/,a,/,a)') stars, '*** Fatal Error (description follows)'
#ifdef HAVE_MPI
  write(stderr, '(a,a)') '*', hyphens
  write(stderr, '(a,i4)') "* From node = ", mpiv%node
#endif
  write(stderr, '(a,a)') '*', hyphens
  do i=1,no_lines
    write(stderr, '(a,1x,a)') '*', trim(message(i))
  end do
  write(stderr, '(a,a)') '*', hyphens

#ifdef DEBUG
  write(stderr, '(a)', advance='no') '* Stack: '
  do i=1,no_sub_stack
    write(stderr, '(a,a)', advance='no') ' > ', trim(sub_stack(i))
  end do
  write(stderr, '(/,a,/)') stars
  call io_status(stderr)
#endif

#ifdef HAVE_MPI
  call MPI_FINALIZE(i)
#endif

  stop
end subroutine write_fatal

subroutine write_warning(no_lines)
  integer, intent(in) :: no_lines
  integer :: i

  ! this always writes from ALL nodes

  if(conf%verbose >= VERBOSE_WARNING) then
    write(stderr, '(/,a)') '** Warning:'
#ifdef HAVE_MPI
    write(stderr, '(a,i4)') "** From node = ", mpiv%node
#endif
    do i=1,no_lines
      write(stderr, '(a,3x,a)') '**', trim(message(i))
    end do
  end if
#ifdef HAVE_FLUSH
  call flush(stderr)
#endif 
 
  return
end subroutine write_warning

subroutine write_info(no_lines, iunit, verbose_limit, stress)
  integer, intent(in) :: no_lines
  integer, intent(in), optional :: iunit
  integer, intent(in), optional :: verbose_limit
  logical, optional, intent(in) :: stress

  integer :: i, iu

#ifdef HAVE_MPI
  if(mpiv%node .ne. 0) return
#endif

  if(present(iunit)) then
    iu = iunit
  else
    iu = stdout
  end if

  if(conf%verbose >= VERBOSE_NORMAL) then
    if(present(stress)) write(iu, '(a)') stars
    do i = 1, no_lines
      if(.not.present(verbose_limit)) then
        write(iu, '(a)') trim(message(i))
      else if(conf%verbose>verbose_limit) then
        write(iu, '(a)') trim(message(i))
      endif
    enddo
    if(present(stress)) write(iu, '(a,/)') stars
  end if

#ifdef HAVE_FLUSH
  call flush(iu)
#endif
end subroutine write_info

#ifdef DEBUG
subroutine push_sub(sub_name)
  character(len=*), intent(in) :: sub_name
  integer i

  no_sub_stack = no_sub_stack + 1
  if(no_sub_stack > 49) then
    sub_stack(50) = 'push_sub'
    message(1) = 'Too many recursion levels (max=50)'
    call write_fatal(1)
  else
    sub_stack(no_sub_stack) = trim(sub_name)
    time_stack(no_sub_stack) = loct_clock()

    if(conf%verbose >= VERBOSE_DEBUG .and. no_sub_stack <= conf%debug_level .and. mpiv%node == 0) then
      write(stderr,'(a,f10.3,i10, a)', advance='no') "* I ", loct_clock()/CNST(1e6), &
           loct_getmem(), " | "
      do i = no_sub_stack-1, 1, -1
        write(stderr,'(a)', advance='no') "..|"
      end do
      write(stderr,'(a)') trim(sub_name)
    end if
  end if

  return
end subroutine push_sub

subroutine pop_sub()
  integer i

  if(no_sub_stack > 0) then
    if(conf%verbose > VERBOSE_DEBUG .and. no_sub_stack <= conf%debug_level .and. mpiv%node == 0) then
      
      ! It seems in std C libraries the number of clock ticks per second is 1e6...
      write(stderr,'(a,f10.3,i10, a)', advance='no') "* O ", &
          (loct_clock()-time_stack(no_sub_stack))/CNST(1e6), &
          loct_getmem(), " | "
      do i = no_sub_stack-1, 1, -1
        write(stderr,'(a)', advance='no') "..|"
      end do
      write(stderr,'(a)') trim(sub_stack(no_sub_stack))
    end if
    no_sub_stack = no_sub_stack - 1
  else
    no_sub_stack = 1
    sub_stack(1) = 'pop_sub'
    message(1) = 'Too few recursion levels'
    call write_fatal(1)    
  end if

end subroutine pop_sub
#endif


subroutine read_system_labels()
  integer               :: i
  integer(POINTER_SIZE) :: blk

  ! first we read the required information from the input file
  ! and prompt the user for possible errors in the input

  ! find out how many subsystem we want to treat ...
  if(loct_parse_block('SystemLabels', blk) == 0) then
     no_syslabels = loct_parse_block_cols(blk,0)
  else
     message(1) = "Could not find block SystemLabels in the input file."
     message(2) = "This block is mandatory for run mode MultiSubsystem."
     call write_fatal(2)
  endif
  
  allocate(subsys_label(no_syslabels), subsys_runmode(no_syslabels))
  allocate(subsys_run_order(no_syslabels))
  
  ! ... and how the user would like to call them.
  do i = 1, no_syslabels
     call loct_parse_block_string(blk, 0, i-1, subsys_label(i))
  enddo
  call loct_parse_block_end(blk)
  
  !  now we check what we have to run in the respective subsystem
  if(loct_parse_block('SystemRunModes', blk) == 0) then
     no_subsys_runmodes = loct_parse_block_cols(blk,0)
  else
     message(1) = "Could not find block SystemRunModes in the input file."
     message(2) = "This block is mandatory for run mode MultiSubsystem."
     call write_fatal(2)
  endif
  
  if(no_subsys_runmodes/=no_syslabels) then
     message(1) = "The blocks SystemLabels and SystemRunModes do not have"
     message(2) = "the same size. Please correct your input file."
     call write_fatal(1)
  endif
  
  do i = 1, no_subsys_runmodes
     call loct_parse_block_int(blk, 0, i-1, subsys_runmode(i))
  enddo
  call loct_parse_block_end(blk)
  
  
  ! ... and in which order 
  if(loct_parse_block('SystemRunOrder', blk) == 0) then
     no_subsys_runmodes = loct_parse_block_cols(blk,0)
  else
     message(1) = "Could not find block SystemRunOrder in the input file."
     message(2) = "This block is mandatory for run mode MultiSubsystem."
     call write_fatal(2)
  endif
  
  if(no_subsys_runmodes/=no_syslabels) then
     message(1) = "The blocks SystemLabels and SystemRunOrder do not have"
     message(2) = "the same size. Please correct your input file."
     call write_fatal(1)
  endif
  
  do i = 1, no_subsys_runmodes
     call loct_parse_block_int(blk, 0, i-1, subsys_run_order(i))
  enddo
  call loct_parse_block_end(blk)
  
end subroutine read_system_labels

subroutine init_default_system_labels()
  
  no_syslabels = 1
  allocate(subsys_label(no_syslabels), subsys_runmode(no_syslabels))
  allocate(subsys_run_order(no_syslabels))
  current_subsystem = 1
  subsys_label(current_subsystem) = ""
  subsys_runmode(current_subsystem) = calc_mode
  subsys_run_order(current_subsystem) = 1

end subroutine init_default_system_labels

! returns true if a file named stop exists
function clean_stop()
  logical clean_stop, file_exists
  
  clean_stop = .false.
  inquire(file='stop', exist=file_exists)
  if(file_exists) then
    message(1) = 'Clean STOP'
    message(2) = "(don't forget to remove the file 'stop' ;)"
    call write_warning(2)
    clean_stop = .true.
  end if

  return
end function clean_stop

character(len=64) function check_inp(variable) result(var_name)
  character(len = * ),           intent(in)  :: variable
!  character(len = * ), optional, intent(in)  :: prefix
  character(len = 64)                        :: composite_name

  composite_name = trim(subsys_label(current_subsystem))//trim(variable)

  if(loct_parse_isdef(composite_name).ne.0) then
!    composite name has been defined in the input file
     var_name = composite_name
  else
!    could not find composite name in the input; 
!    will use bare variable name
     var_name = variable
  endif

end function check_inp

! Given a path, it returns the extension (if it exists) of the file
! (that is, the part of the name that comes after its last point)
! If the filename does not have an extension, it returns the empty string.
character(len=8) function get_extension(path) result(ext)
  character(len = * ), intent(in)  :: path
  integer :: i, j

  i = index(path, ".", back = .true.)
  j = index(path(i+1:), "/")
  if(i.eq.0 .or. j.ne.0) then
    ext = ""
  else
    ext = path(i+1:)
  endif
end function get_extension


end module global
