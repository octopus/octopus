#include "config_F90.h"

program octopus
  use global
  use liboct
  use run_prog

  implicit none

  integer :: ierr, val(8)

#ifdef HAVE_MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiv%node, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiv%numprocs, ierr)
  write(stdout,'(a,i4,a,i4,a)') 'Process ', mpiv%node, ' of ', mpiv%numprocs, ' is alive'  
#else
  mpiv%node = 0
  mpiv%numprocs = 1
#endif

  ! init some of the stuff
  ierr = oct_parse_init(C_string('inp'), C_string('out.oct'))
  if(ierr .ne. 0) then
    message(1) = "Error initializing liboct"
    call write_fatal(1)
  end if
  
  call oct_parse_int(C_string('verbose'), 30, conf%verbose)
  
  if(conf%verbose >= 999 .and. mpiv%node == 0) then
    message(1) = 'Entering DEBUG mode'
    call write_warning(1)
  end if
  
  ! Sets the dimensionaliy of the problem.
#ifdef ONE_D
  conf%dim=1
#else
  conf%dim=3
#endif

  ! Let us print our logo
  ! it is damn hard to print ascii-art in FORTRAN ;((
  if(conf%verbose > 20 .and. mpiv%node == 0) then
    write(stdout, *)
    write(stdout,'(4x,a)')"<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
    write(stdout,'(9x,a)')"                       ___"
    write(stdout,'(9x,a)')"                    .-'   `'."
    write(stdout,'(9x,a)')"                   /         \ "
    write(stdout,'(9x,a)')"                   |         ;"
    write(stdout,'(9x,a)')"                   |         |           ___.--,"
    write(stdout,'(9x,a)')"          _.._     |0) ~ (0) |    _.---'`__.-( (_."
    write(stdout,'(9x,2a)')"   __.--'`_.. '.__.\    '--. \_.-' ,.--'`     `",'""`'
    write(stdout,'(9x,a)')"  ( ,.--'`   ',__ /./;   ;, '.__.'`    __"
    write(stdout,'(9x,2a)')"  _`) )  .---.__.' / |   |\   \__..--",'""  """--.,_'
    write(stdout,'(9x,a)')" `---' .'.''-._.-'`_./  /\ '.  \ _.-~~~````~~~-._`-.__.'"
    write(stdout,'(9x,a)')"       | |  .' _.-' |  |  \  \  '.               `~---`"
    write(stdout,'(9x,a)')"        \ \/ .'     \  \   '. '-._)"
    write(stdout,'(9x,a)')"         \/ /        \  \    `=.__`~-."
    write(stdout,'(9x,2a)')"    jgs  / /\         `) )    / / `",'"".`\ '
    write(stdout,'(9x,a)')"   , _.-'.'\ \        / /    ( (     / /"
    write(stdout,'(9x,a)')"    `--~`   ) )    .-'.'      '.'.  | ("
    write(stdout,'(9x,a)')"           (/`    ( (`          ) )  '-;"
    write(stdout,'(9x,a)')"            `      '-;         (-'"
    write(stdout,'(4x,a)')"<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
    write(stdout, *)
  end if

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation started on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

  ! now we really start
  !call test()
  !stop
  call run()

  ! print date
  call date_and_time(values=val)
  write(message(1),'(a,i4,a1,i2.2,a1,i2.2,a,i2.2,a1,i2.2,a1,i2.2)') &
       "Info: Calculation ended on ", val(1), "/", val(2), "/", val(3), &
       " at ", val(5), ":", val(6), ":", val(7)
  call write_info(1)

#ifdef HAVE_MPI
  call MPI_FINALIZE(ierr)
#endif
  
  call oct_parse_end()
  stop
end program octopus

subroutine test()
  use global

  character(len=100) :: c1, c2

  c1 = 'ola'
  write(*,'(a2,a7,a2)') 'aa', str_center(trim(c1), 7), 'aa'

end subroutine test
