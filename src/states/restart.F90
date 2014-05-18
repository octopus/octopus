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
!! $Id$

#include "global.h"

module restart_m
  use batch_m
  use datasets_m
  use global_m
  use io_binary_m
  use io_function_m
  use loct_m
  use mesh_m
  use mesh_batch_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use unit_system_m

  implicit none

  private
  public ::                  &
    restart_init,            &
    restart_dump_start,      &
    restart_dump_finish,     &
    clean_stop,              &
    restart_dir,             &
    restart_format,          &
    drestart_write_function, &
    zrestart_write_function, &
    drestart_read_function,  &
    zrestart_read_function,  &
    write_restart

  logical           :: restart_write_files
  integer           :: restart_format
  character(len=32) :: restart_dir

 !> from signals.c
  interface
    subroutine block_signals()
      implicit none
    end subroutine block_signals

    subroutine unblock_signals()
      implicit none
    end subroutine unblock_signals
  end interface

contains

  !> returns true if a file named stop exists
  function clean_stop(comm)
    integer, intent(in) :: comm !< communicator spanning all nodes that will call this function, i.e. not any slaves

    logical :: clean_stop, file_exists

    PUSH_SUB(clean_stop)

    clean_stop = .false.

    if(mpi_grp_is_root(mpi_world)) then
      inquire(file='stop', exist=file_exists)
      if(file_exists) then
        call loct_rm('stop')
        clean_stop = .true.
      endif
    end if

#ifdef HAVE_MPI
    ! make sure all nodes agree on whether this condition occurred
    call MPI_Bcast(clean_stop, 1, MPI_LOGICAL, 0, comm, mpi_err)
#endif

    if(clean_stop) then
      message(1) = 'Clean STOP'
      call messages_warning(1)
    endif

    POP_SUB(clean_stop)
  end function clean_stop


  ! ---------------------------------------------------------
  !> read restart format information
  subroutine restart_init()

    PUSH_SUB(restart_init)

    call messages_obsolete_variable('RestartFileFormat', 'RestartWrite')

    !%Variable RestartWrite
    !%Type logical
    !%Default true
    !%Section Execution::IO
    !%Description
    !% If this variable is set to no, restart information is not
    !% written. The default is yes.
    !%End

    call parse_logical(datasets_check('RestartWrite'), .true., restart_write_files)

    if(restart_write_files) then
      restart_format = io_function_fill_how("Binary")
    else
      message(1) = 'Restart information will not be written.'
      call messages_warning(1)
    end if

    !%Variable RestartDir
    !%Type string
    !%Default ''
    !%Section Execution::IO
    !%Description
    !% When <tt>Octopus</tt> reads restart files, e.g. when running a time-propagation
    !% after a ground-state calculation, these files will be read from
    !% <tt>&lt;RestartDir&gt/</tt>. Usually, <tt>RestartDir</tt> is
    !% <tt>TmpDir</tt> but in a transport calculation, the output of
    !% a periodic dataset is required to calculate the extended ground state.
    !%End
    call parse_string(datasets_check('RestartDir'), tmpdir, restart_dir)
    ! Append "/" if necessary.
    if(scan(restart_dir, '/', .true.) < 1) then
      restart_dir = trim(restart_dir)//'/'
    end if

    POP_SUB(restart_init)
  end subroutine restart_init


  subroutine restart_dump_start()
    PUSH_SUB(restart_dump_start)

    if(restart_write_files) call block_signals()

    POP_SUB(restart_dump_start)
  end subroutine restart_dump_start


  subroutine restart_dump_finish()
    PUSH_SUB(restart_dump_finish)

    if(restart_write_files) call unblock_signals()

    POP_SUB(restart_dump_finish)
  end subroutine restart_dump_finish


  ! --------------------------------------------
  logical pure function write_restart()
    
    write_restart = restart_write_files

  end function write_restart


#include "undef.F90"
#include "real.F90"
#include "restart_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "restart_inc.F90"

end module restart_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
