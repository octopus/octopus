#include "config_F90.h"

module global

#if defined(HAVE_MPI) && !defined(MPI_H)
use mpi
#endif

implicit none

#if defined(HAVE_MPI) && defined(MPI_H)
#include "mpif.h"
#endif

type conf_type
  integer :: verbose ! <= 0  -> silent, no output except fatal errors
                     ! > 0   -> warning only
                     ! > 20  -> normal program info
                     ! >= 999 -> debug
  integer :: dim
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

real(r8), parameter :: eps_r8 = epsilon(1._r8)
real(r8), parameter :: inf_r8 = 10.0_r8 ** (12)
real(r8), parameter :: r_small = 0.0001_r8

! some mathematical constants
real(r8), parameter    :: M_Pi=3.141592653589793_r8
complex(r8), parameter :: M_zI=(0.0_r8,1.0_r8)
complex(r8), parameter :: M_z0=(0.0_r8,0.0_r8)
complex(r8), parameter :: M_z1=(1.0_r8,0.0_r8)

! some physical constants
real(r8), parameter :: P_Ang =  1._r8 / 0.529177_r8
real(r8), parameter :: P_eV  =  1._r8 / 13.60580_r8
real(r8), parameter :: P_E   =  3.79470065_r8! (electrons charge)
real(r8), parameter :: P_M   =  0.13123_r8   ! (electrons mass)
real(r8), parameter :: P_E2  =  14.399753_r8 ! ( = 2.0_r8/(Ang*eV) )
real(r8), parameter :: P_H2M = 3.81_r8
real(r8), parameter :: P_a_B = 0.529177_r8
real(r8), parameter :: P_Ry = 13.6058_r8

! the standard input and output
integer :: stderr = 0
integer :: stdin  = 5
integer :: stdout = 6

! some variables to be used everywhere
character(len=75), dimension(20) :: message ! to be output by fatal, warning
character(len=40) :: sub_name ! the sub name to be used by push/pop_sub

! some private variables to this module
character(len=40), private :: sub_stack(50)
real(r8), private          :: time_stack(50)
integer, private           :: no_sub_stack = 0
character(len=70), parameter, private :: stars =  &
    '**********************************************************************'
character(len=68), parameter, private :: hyphens = &
    '--------------------------------------------------------------------'

contains

subroutine write_fatal(no_lines)
  integer, intent(in) :: no_lines
  integer :: i

  write(stdout, '(/,a,/,a)') stars, '*** Fatal Error (description follows)'
#ifdef HAVE_MPI
  write(stdout, '(a,a)') '*', hyphens
  write('(a,i4)') "* From node = ", mpiv%node
#endif
  write(stdout, '(a,a)') '*', hyphens
  do i=1,no_lines
    write(stdout, '(a,1x,a)') '*', trim(message(i))
  end do
  write(stdout, '(a,a)') '*', hyphens
  write(stdout, '(a)', advance='no') '* Stack: '
  do i=1,no_sub_stack
    write(stdout, '(a,a)', advance='no') ' > ', trim(sub_stack(i))
  end do
  write(stdout, '(/,a,/)') stars

#ifdef HAVE_MPI
  call MPI_FINALIZE(i)
#endif

  stop
end subroutine write_fatal

subroutine write_warning(no_lines)
  integer, intent(in) :: no_lines
  integer :: i

  ! this always writes from ALL nodes

  if(conf%verbose>0) then
    write(stdout, '(/,a)') '** Warning:'
#ifdef HAVE_MPI
    write('(a,i4)') "** From node = ", mpiv%node
#endif
    do i=1,no_lines
      write(stdout, '(a,3x,a)') '**', trim(message(i))
    end do
  end if
  
  return
end subroutine write_warning

subroutine write_info(no_lines, iunit)
  integer, intent(in) :: no_lines
  integer, intent(in), optional :: iunit

  integer :: i, iu

#ifdef HAVE_MPI
  if(mpiv%node .ne. 0) return
#endif

  if(present(iunit)) then
    iu = iunit
  else
    iu = stdout
  end if

  if(conf%verbose>20) then
    do i=1,no_lines
      write(iu, '(a)') trim(message(i))
    end do
  end if

  return
end subroutine write_info

subroutine push_sub()
  integer i

  character(len=8)  :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer, dimension(8) :: values
  real(r8) :: t

  no_sub_stack = no_sub_stack + 1
  if(no_sub_stack > 49) then
    sub_stack(50) = 'push_sub'
    message(1) = 'Too many recursion levels (max=50)'
    call write_fatal(1)
  else
    sub_stack(no_sub_stack) = trim(sub_name)

    call date_and_time(date, time, zone, values)
    t = values(8)/1000._r8 + (values(7) + &
      60*(values(6) + 60*(values(5))))

    time_stack(no_sub_stack) = t

    if(conf%verbose > 999 .and. mpiv%node == 0) then
      write(stdout,'(a)', advance='no') "* Debug: In: "
      do i = no_sub_stack-1, 1, -1
        write(stdout,'(a)', advance='no') "  "
      end do
      write(stdout,'(a)') trim(sub_name)
    end if
  end if

  return
end subroutine push_sub

subroutine pop_sub()
  integer i

  character(len=8)  :: date
  character(len=10) :: time
  character(len=5) :: zone
  integer, dimension(8) :: values
  real(r8) :: t

  if(no_sub_stack > 0) then
    if(conf%verbose > 999 .and. mpiv%node == 0) then
      call date_and_time(date, time, zone, values)
      t = values(8)/1000._r8 + (values(7) + &
          60*(values(6) + 60*(values(5))))
 
      write(stdout,'(a)', advance='no') "* Debug: Out:"
      do i = no_sub_stack-1, 1, -1
        write(stdout,'(a)', advance='no') "  "
      end do
      write(stdout,'(a)', advance='no') trim(sub_stack(no_sub_stack))
      write(stdout, '(a,f10.3,a)') ' (', t-time_stack(no_sub_stack), 's)'
    end if
    no_sub_stack = no_sub_stack - 1
  else
    no_sub_stack = 1
    sub_stack(1) = 'pop_sub'
    message(1) = 'Too few recursion levels'
    call write_fatal(1)    
  end if

end subroutine pop_sub

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

! Upcases a string
!  15-OCT-2000: First version, Fernando Nogueira
SUBROUTINE upcase(str)
  
  CHARACTER (LEN=*), INTENT(INOUT) :: str
  INTEGER :: i, s
  
  DO i = 1, LEN(str)
    s = IACHAR(str(i:i))
    IF ( (s<=122) .AND. (s>=97) ) s = s-32
    str(i:i) = ACHAR(s)
  ENDDO
  
END SUBROUTINE upcase

! Lowcases a string
!  15-OCT-2000: First version, Fernando Nogueira
SUBROUTINE lowcase(str)

  CHARACTER (LEN=*), INTENT(INOUT) :: str
  INTEGER :: i, s
  
  DO i = 1, LEN(str)
    s = IACHAR(str(i:i))
    IF ( (s<=90) .AND. (s>=65) ) s = s+32
    str(i:i) = ACHAR(s)
  ENDDO
  
END SUBROUTINE lowcase

! Removes all spaces from a string
!  15-OCT-2000: First version, Fernando Nogueira
SUBROUTINE compact(str)

  CHARACTER (LEN=*), INTENT(INOUT) :: str
  INTEGER :: i, j
  
  DO i = 1, LEN(str)
    IF ( str(i:i) == ' ' ) THEN
      DO j = i, LEN(str)-1
        str(j:j) = str(j+1:j+1)
      ENDDO
      str(LEN(str):LEN(str)) = ' '
    ENDIF
  ENDDO
  
END SUBROUTINE compact

subroutine str_trim(str)
  character (len=*), intent(inout) :: str
  integer :: i, j, k, l
  
  l = len(str)
  do i = 1, l
    if(str(i:i) .ne. ' ') exit
  end do

  do j = 1, l - i + 1
    str(j:j) = str(i:i)
    i = i + 1
  end do

  do i = j, l
    str(j:j) = ' '
  end do

end subroutine str_trim

! puts space around string, so that it is centered
character(len=100) function str_center(s_in, l) result(s_out)
  character(len=*), intent(IN) :: s_in
  integer, intent(in) :: l

  integer :: pad, i, li

  li = len(s_in)
  if(l < li) then
    s_out(1:l) = s_in(1:l)
    return
  end if
  
  pad = (l - li)/2

  s_out = ""
  do i = 1, pad
    s_out(i:i) = " ";
  end do
  s_out(pad + 1:pad + li + 1) = s_in(1:li)
  do i = pad + li + 1, l
    s_out(i:i) = " ";
  end do
  
end function str_center

end module global
