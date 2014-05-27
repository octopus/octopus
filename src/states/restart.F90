!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2014 M. Oliveira
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
  use index_m
  use io_m
  use io_binary_m
  use io_function_m
  use loct_m
  use mesh_m
  use mesh_batch_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use simul_box_m
  use unit_system_m

  implicit none

  private
  public ::                  &
    restart_t,               &
    clean_stop,              &
    restart_init,            &
    restart_end,             &
    restart_dir,             &
    restart_cd,              &
    restart_mkdir,           &
    restart_rm,              &
    restart_open,            &
    restart_close,           &
    restart_block_signals,   &
    restart_unblock_signals, &
    restart_skip,            &
    RESTART_TYPE_DUMP,       &
    RESTART_TYPE_LOAD,       &
    drestart_write_function, &
    zrestart_write_function, &
    drestart_read_function,  &
    zrestart_read_function

  type restart_t
    private
    integer           :: type    !< Restart type: RESTART_TYPE_DUMP or RESTART_TYPE_LOAD
    logical           :: skip    !< If set to .true., no restart information should be loaded or dumped.
    integer           :: format  !< Format used to store the restart information.
    character(len=80) :: dir     !< Directory where the restart information is stored.
    character(len=80) :: pwd     !< The current directory where the restart information is being loaded from or dumped to.
                                 !! It can be either dir or a subdirectory of dir.
    type(mpi_grp_t)   :: mpi_grp !< Some operations require an mpi group to be used.

    integer, pointer  :: map(:)  !< Map between the points of the stored mesh and the mesh used in the current calculations.
  end type restart_t

  integer, parameter :: RESTART_TYPE_DUMP = 1, &
                        RESTART_TYPE_LOAD = 2

 !> from signals.c
  interface restart_block_signals
    subroutine block_signals()
      implicit none
    end subroutine block_signals
  end interface restart_block_signals

  interface restart_unblock_signals
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
  !> Initializes a restart object.
  subroutine restart_init(restart, type, dirname, mpi_grp, mesh, sb, basedir, exact)
    type(restart_t),             intent(out) :: restart !< Restart information.
    integer,                     intent(in)  :: type    !< Is this restart used for dumping (type = RESTART_TYPE_DUMP)
                                                        !! or for loading (type = RESTART_TYPE_LOAD)?
    character(len=*),            intent(in)  :: dirname !< Directory where restart information is going to be loaded 
                                                        !! from or dumped to.
    type(mpi_grp_t),             intent(in)  :: mpi_grp !> The mpi group in charge of handling this restart.
    type(mesh_t),      optional, intent(in)  :: mesh    !< If present, depending on the type of restart, the mesh 
                                                        !! information is either dumped or the mesh compatibility is checked.
    type(simul_box_t), optional, intent(in)  :: sb      !< If present and type = RESTART_TYPE_DUMP, the simulation box 
                                                        !! information will be dumped.
    character(len=*),  optional, intent(in)  :: basedir !< Parent directory of dirname. If not present, it is obtained from 
                                                        !! the input file.
    logical,           optional, intent(in)  :: exact   !< If loading the restart information, should the mesh be 
                                                        !! exactly the same or not?

    logical :: grid_changed, grid_reordered, restart_write, dir_exists
    integer :: ierr, iunit
    character(len=80) :: basedir_, dirname_

    PUSH_SUB(restart_init)

    if (present(exact) .and. .not. present(mesh)) then
      message(1) = "Error in restart_init: the 'exact' optional argument requires a mesh."
      call messages_fatal(1)
    end if

    ! Some initializations
    restart%type = type
    nullify(restart%map)
    restart%mpi_grp = mpi_grp
    restart%format = io_function_fill_how("Binary")

    select case (restart%type)
    case (RESTART_TYPE_DUMP)

      if (.not. present(basedir)) then
        !%Variable TmpDir
        !%Default "restart/"
        !%Type string
        !%Section Execution::IO
        !%Description
        !% The name of the directory where <tt>Octopus</tt> stores binary information
        !% such as the wavefunctions.
        !%End
        call parse_string('TmpDir', trim(current_label)//'restart', basedir_)

        call messages_obsolete_variable('RestartFileFormat', 'RestartWrite')

      else
        basedir_ = basedir
      end if

      !%Variable RestartWrite
      !%Type logical
      !%Default true
      !%Section Execution::IO
      !%Description
      !% If this variable is set to no, restart information is not
      !% written. The default is yes.
      !%End

      call parse_logical(datasets_check('RestartWrite'), .true., restart_write)
      restart%skip = .not. restart_write

      if(restart%skip) then
        message(1) = 'Restart information will not be written.'
        call messages_warning(1)
      end if

    case (RESTART_TYPE_LOAD)
      ! We should never skip anything when loading the restart information
      restart%skip = .false.

      if (.not. present(basedir)) then
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
        call parse_string(datasets_check('RestartDir'),  trim(current_label)//'restart', basedir_)
      else
        basedir_ = basedir
      end if

    case default
      message(1) = "Unknown restart type in restart_init"
      call messages_fatal(1)
    end select

    ! Remove any trailing "/" from the paths (all the routines from this module should add the trailing "/" when needed)
    if (index(basedir_, '/', .true.) == len_trim(basedir_)) then
      basedir_ = basedir_(1:len_trim(basedir_)-1)
    end if
    if (index(dirname, '/', .true.) == len_trim(dirname)) then
      dirname_ = dirname(1:len_trim(dirname)-1)
    else
      dirname_ = dirname
    end if

    ! Set final paths
    restart%dir = trim(basedir_)//"/"//trim(dirname_)
    restart%pwd = restart%dir

    ! Check if the directory already exists and create it if necessary
    dir_exists = loct_dir_exists(trim(restart%pwd))
    if (restart%type == RESTART_TYPE_DUMP .and. .not. dir_exists) then
      call io_mkdir(trim(restart%pwd), is_tmp=.true., parents=.true.)
    end if

    select case (restart%type)
    case (RESTART_TYPE_DUMP)
      if (.not. restart%skip) then
        message(1) = "Info: Restart information will be written to '"//trim(restart%pwd)//"'."
        call messages_info(1)

        if (present(mesh)) then
          ! Maybe the restart directory is also being used to load restart information. In that case, the
          ! mesh must be exactly the same. If not, stop the code, as this might result in inconsistent
          ! restart data.
          iunit = io_open(trim(restart%pwd)//'/loading', action='read', status='old', die=.false., is_tmp=.true.)
          if (iunit > 0) then
            call mesh_check_dump_compatibility(restart%pwd, mesh, grid_changed, grid_reordered, restart%map, ierr)
            if (ierr >= 0 .and. (grid_changed .or. grid_reordered)) then
              message(1) = "Internal error: trying to write to a directory with previous"
              message(2) = "restart information written with a different mesh."
              call messages_fatal(2)
            end if
          end if
          close(iunit)
        end if

        ! Dump the grid information. The main parameters of the grid should not change
        ! during the calculation, so we should only need to dump it once.
        if (mpi_grp_is_root(restart%mpi_grp)) then
          if (present(mesh)) then
            iunit = io_open(trim(restart%pwd)//'/mesh', action='write', is_tmp=.true.)
            write(iunit,'(a)') '# This file contains the necessary information to generate the'
            write(iunit,'(a)') '# grid with which the functions in this directory were calculated,'
            write(iunit,'(a)') '# except for the geometry of the system.'
            call io_close(iunit)

            call mesh_dump(mesh, restart%pwd, "mesh")
            call index_dump_lxyz(mesh%idx, mesh%np_part_global, restart%pwd, ierr)

            call mesh_write_fingerprint(mesh, trim(restart%pwd)//"/"//"grid")
          end if

          if (present(sb)) then
            call simul_box_dump(sb, restart%pwd, "mesh")
          end if
        end if

        ! Mark the directory as been used for dumping. Note that only one restart instance can 
        ! dump to a given directory at the same time.
        if (mpi_grp_is_root(restart%mpi_grp)) then
          iunit = io_open(trim(restart%pwd)//'/dumping', action='write', status='new', die=.false., is_tmp=.true.)
          if (iunit < 0) then
            message(1) = "Internal error: directory '"//trim(restart%pwd)//"' already been used"
            message(2) = "for restart dumping."
            call messages_fatal(2)
          end if
          call io_close(iunit)
        end if

      end if

    case (RESTART_TYPE_LOAD)
      if (.not. dir_exists) then
        restart%skip = .true.

        message(1) = "Info: Could not find '"//trim(restart%pwd)//"' directory for restart."
        message(2) = "Info: No restart information will be read."
        call messages_warning(2)

      else
        message(1) = "Info: Restart information will be read from '"//trim(restart%pwd)//"'."
        call messages_info(1)

        ! Mark the directory as been used for loading. Note that only one restart instance can 
        ! load from a given directory at the same time.
        if(mpi_grp_is_root(restart%mpi_grp)) then
          iunit = io_open(trim(restart%pwd)//'/loading', action='write', status='new', die=.false., is_tmp=.true.)
          if (iunit < 0) then
            message(1) = "Internal error: directory '"//trim(restart%pwd)//"' already been used for"
            message(2) = "restart loading."
            call messages_fatal(2)
          end if
          call io_close(iunit)
        end if

        if (present(mesh)) then
          call mesh_check_dump_compatibility(restart%pwd, mesh, grid_changed, grid_reordered, restart%map, ierr)

          ! Check if an error occurred. If so, stop the calculation, because at the moment we really need a compatible mesh.
          if (ierr /= 0) then
            if (ierr == -1) then
              message(1) = "Unable to check mesh compatibility: unable to open mesh fingerprint"
              message(2) = "in '"//trim(restart%pwd)//"'."
            else if (ierr > 0) then
              message(1) = "Mesh from current calculation is not compatible with mesh found in"
              message(2) = "'"//trim(restart%pwd)//"'."
            end if
            call messages_fatal(2)
          end if

          ! Print some warnings in case the mesh is compatible, but changed.
          if (grid_changed) then
            if (grid_reordered) then
              message(1) = "Octopus is attempting to restart from a mesh with a different order of points."
            else
              message(1) = "Octopus is attempting to restart from a different mesh."
            end if
            call messages_warning(1)
          end if

          if (present(exact)) then
            restart%skip = grid_changed .and. .not. grid_reordered .and. exact
            if (restart%skip) then
              message(1) = "This calculation requires the exact same mesh to restart."
              message(2) = "No restart information will be read from '"//trim(restart%pwd)//"'."
              call messages_warning(2)
            end if
          else
            restart%skip = .false.
          end if
        end if
      end if

    end select


    ! Make sure all the processes have finished reading/writing all the grid information,
    ! as there might be some subsequent calls to this function where that information will
    ! be written/read to/from the same directory.
#if defined(HAVE_MPI)
        call MPI_Barrier(restart%mpi_grp, mpi_err)
#endif

    POP_SUB(restart_init)
  end subroutine restart_init


  ! ---------------------------------------------------------
  subroutine restart_end(restart)
    type(restart_t),  intent(inout) :: restart

    PUSH_SUB(restart_end)

    if(mpi_grp_is_root(restart%mpi_grp) .and. .not. restart%skip) then
      select case (restart%type)
      case (RESTART_TYPE_LOAD)
        message(1) = "Info: Finished reading information from '"//trim(restart%dir)//"'."
        call loct_rm(trim(restart%pwd)//"/loading")
      case (RESTART_TYPE_DUMP)
        call loct_rm(trim(restart%pwd)//"/dumping")
        message(1) = "Info: Finished writing information to '"//trim(restart%dir)//"'."
      end select
      call messages_info(1)
    end if

    restart%type = 0
    restart%skip = .true.
    SAFE_DEALLOCATE_P(restart%map)

    POP_SUB(restart_end)
  end subroutine restart_end


  ! ---------------------------------------------------------
  !> Returns the name of the directory containing the restart 
  !! information. The use of this function should be avoided, 
  !! as the access to the restart data should always be done 
  !! through this module, and it is only provided to allow 
  !! some older parts of the code to keep functioning until 
  !! someone fixes them.
  function restart_dir(restart)
    type(restart_t), intent(in) :: restart
    character(len=80) :: restart_dir

    PUSH_SUB(restart_dir)

    restart_dir = restart%pwd

    POP_SUB(restart_dir)
  end function restart_dir


  ! ---------------------------------------------------------
  !> If "dirname" is present, change the restart directory to
  !! dirname, where "dirname" is a subdirectory of the current 
  !! restart directory. If "dirname" is not present, change back
  !! to the initial directory.
  subroutine restart_cd(restart, dirname)
    type(restart_t),            intent(inout) :: restart
    character(len=*), optional, intent(in)    :: dirname

    PUSH_SUB(restart_cd)

    ASSERT(.not. restart%skip)

    if (present(dirname)) then
      if (restart%type == RESTART_TYPE_DUMP) then
        call restart_mkdir(restart, dirname)
      end if

      if (index(dirname, '/', .true.) == len_trim(dirname)) then
        restart%pwd = trim(restart%dir)//"/"//dirname(1:len_trim(dirname)-1)
      else
        restart%pwd = trim(restart%dir)//"/"//trim(dirname)
      end if

    else
      restart%pwd = restart%dir
    end if

    POP_SUB(restart_cd)
  end subroutine restart_cd


  ! ---------------------------------------------------------
  !> Make directory "dirname" inside the current restart directory.
  subroutine restart_mkdir(restart, dirname)
    type(restart_t),  intent(in) :: restart
    character(len=*), intent(in) :: dirname

    PUSH_SUB(restart_mkdir)

    ASSERT(.not. restart%skip)

    ASSERT (restart%type == RESTART_TYPE_DUMP)

    call io_mkdir(trim(restart%pwd)//"/"//trim(dirname), is_tmp=.true., parents=.true.)

    POP_SUB(restart_mkdir)
  end subroutine restart_mkdir


  ! ---------------------------------------------------------  
  !> Remove directory "dirname" that is located inside the current restart directory.
  subroutine restart_rm(restart, name)
    type(restart_t),  intent(in) :: restart
    character(len=*), intent(in) :: name

    ASSERT(.not. restart%skip)
    ASSERT(restart%type == RESTART_TYPE_DUMP)

    PUSH_SUB(restart_rm)

    call loct_rm(trim(restart%pwd)//"/"//trim(name))

    POP_SUB(restart_rm)
  end subroutine restart_rm


  ! ---------------------------------------------------------
  !> Open file "filename" found inside the current restart
  !! directory. Depending on the type of restart, the file
  !! will be used for reading  (RESTART_TYPE_LOAD), or for
  !! writing (RESTART_TYPE_DUMP).
  function restart_open(restart, filename, status, position, form)
    type(restart_t),            intent(in) :: restart
    character(len=*),           intent(in) :: filename
    character(len=*), optional, intent(in) :: status
    character(len=*), optional, intent(in) :: position
    character(len=*), optional, intent(in) :: form
    integer :: restart_open

    logical :: die
    character(len=20) :: action, status_

    PUSH_SUB(restart_open)

    ASSERT(.not. restart%skip)

    select case (restart%type)
    case (RESTART_TYPE_DUMP)
      status_ = 'unknown'
      action = 'write'
      die = .true.

    case (RESTART_TYPE_LOAD)
      status_ = 'old'
      action = 'read'
      die = .false.

    case default
      message(1) = "Error in restart_open: illegal restart type"
    end select

    if (present(status)) status_ = status

    restart_open = io_open(trim(restart%pwd)//"/"//trim(filename), action=trim(action), status=trim(status_), &
                           die=die, is_tmp=.true., position=position, form=form, grp=restart%mpi_grp)

    POP_SUB(restart_open)
  end function restart_open


  ! ---------------------------------------------------------
  !> Close a file previously opened with restart_open.
  subroutine restart_close(restart, iunit)
    type(restart_t), intent(in)    :: restart
    integer,         intent(inout) :: iunit

    PUSH_SUB(restart_close)

    call io_close(iunit, restart%mpi_grp)

    POP_SUB(restart_close)
  end subroutine restart_close


  ! ---------------------------------------------------------
  !> Returns true if the restart information should not be 
  !! read nor written. This might happen because the user
  !! chose not to write any restart information, or because
  !! the restart information is not available for reading.
  logical pure function restart_skip(restart)
    type(restart_t), intent(in) :: restart

    restart_skip = restart%skip

  end function restart_skip


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
