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

#include "global.h"

module restart_oct_m
  use batch_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use space_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                       &
    restart_t,                    &
    clean_stop,                   &
    restart_module_init,          &
    restart_init,                 &
    restart_end,                  &
    restart_dir,                  &
    restart_open_dir,             &
    restart_close_dir,            &
    restart_mkdir,                &
    restart_rm,                   &
    restart_open,                 &
    restart_write,                &
    restart_read,                 &
    restart_close,                &
    restart_block_signals,        &
    restart_unblock_signals,      &
    restart_skip,                 &
    restart_has_flag,             &
    restart_has_map,              &
    drestart_write_mesh_function, &
    zrestart_write_mesh_function, &
    drestart_read_mesh_function,  &
    zrestart_read_mesh_function,  &
    drestart_write_binary,        &
    zrestart_write_binary,        &
    drestart_read_binary,         &
    zrestart_read_binary,         &
    restart_are_basedirs_equal

  interface drestart_write_binary
    module procedure drestart_write_binary1, drestart_write_binary2, drestart_write_binary3, drestart_write_binary5
  end interface drestart_write_binary

  interface zrestart_write_binary
    module procedure zrestart_write_binary1, zrestart_write_binary2, zrestart_write_binary3, zrestart_write_binary5
  end interface zrestart_write_binary

  interface drestart_read_binary
    module procedure drestart_read_binary1, drestart_read_binary2, drestart_read_binary3, drestart_read_binary5
  end interface drestart_read_binary

  interface zrestart_read_binary
    module procedure zrestart_read_binary1, zrestart_read_binary2, zrestart_read_binary3, zrestart_read_binary5
  end interface zrestart_read_binary


  type restart_t
    private
    integer           :: data_type !< Type of information that the restart is supposed to read/write (GS, TD, etc)
    integer           :: type      !< Restart type: RESTART_TYPE_DUMP or RESTART_TYPE_LOAD
    logical           :: skip      !< If set to .true., no restart information should be loaded or dumped.
    integer(8)        :: format    !< Format used to store the restart information.
    character(len=MAX_PATH_LEN) :: dir !< Directory where the restart information is stored.
    character(len=MAX_PATH_LEN) :: pwd !< The current directory where the restart information is being loaded from or dumped to.
                                       !! It can be either dir or a subdirectory of dir.
    type(namespace_t), pointer :: namespace !< namespace depending on system to modify path
    type(mpi_grp_t)   :: mpi_grp   !< Some operations require an mpi group to be used.
    type(multicomm_t), pointer :: mc
    logical           :: has_mesh  !< If no, mesh info is not written or read, and mesh functions cannot be written or read.
    integer, allocatable :: map(:) !< Map between the points of the stored mesh and the mesh used in the current calculations.
  end type restart_t


  type restart_data_t
    private
    character(len=20) :: tag
    character(len=MAX_PATH_LEN) :: basedir
    character(len=MAX_PATH_LEN) :: dir
    integer :: flags
  end type restart_data_t


  integer, parameter, public :: RESTART_TYPE_DUMP = 1, &
                                RESTART_TYPE_LOAD = 2

  integer, parameter, public :: RESTART_UNDEFINED  = -1,  &
                                RESTART_ALL        =  0,  &                                
                                RESTART_GS         =  1,  &
                                RESTART_UNOCC      =  2,  &
                                RESTART_TD         =  3,  &
                                RESTART_EM_RESP    =  4,  &
                                RESTART_EM_RESP_FD =  5,  &
                                RESTART_KDOTP      =  6,  &
                                RESTART_VIB_MODES  =  7,  &
                                RESTART_VDW        =  8,  &
                                RESTART_CASIDA     =  9,  &
                                RESTART_OCT        =  10, &
                                RESTART_PARTITION  =  11, &
                                RESTART_PROJ       =  12, &
                                RESTART_MAXWELL    =  13, &
                                RESTART_TD_MAXWELL =  14

  integer, parameter :: RESTART_N_DATA_TYPES = 12

  integer, parameter, public :: RESTART_FLAG_STATES = 1,  &
                                RESTART_FLAG_RHO    = 2,  &
                                RESTART_FLAG_VHXC   = 4,  &
                                RESTART_FLAG_MIX    = 8,  &
                                RESTART_FLAG_SKIP   = 16

  type(restart_data_t) :: info(RESTART_N_DATA_TYPES)


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
      end if
    end if

#ifdef HAVE_MPI
    ! make sure all nodes agree on whether this condition occurred
    call MPI_Bcast(clean_stop, 1, MPI_LOGICAL, 0, comm, mpi_err)
#endif

    if(clean_stop) then
      message(1) = 'Clean STOP'
      call messages_warning(1)
    end if

    POP_SUB(clean_stop)
  end function clean_stop


  ! ---------------------------------------------------------
  subroutine restart_module_init(namespace)
    type(namespace_t),         intent(in)    :: namespace

    logical :: set(RESTART_N_DATA_TYPES)
    integer :: iline, n_cols, data_type
    character(len=MAX_PATH_LEN) :: default_basedir
    type(block_t) :: blk

    PUSH_SUB(restart_module_init)

    ! Each data type should have a tag
    info(RESTART_GS)%tag = "Ground-state"
    info(RESTART_UNOCC)%tag = "Unoccupied states"
    info(RESTART_TD)%tag = "Time-dependent"
    info(RESTART_EM_RESP)%tag = "EM Resp."
    info(RESTART_EM_RESP_FD)%tag = "EM Resp. FD"
    info(RESTART_KDOTP)%tag = "KdotP"
    info(RESTART_VIB_MODES)%tag = "Vib. Modes"
    info(RESTART_VDW)%tag = "VdW"
    info(RESTART_CASIDA)%tag = "Casida"
    info(RESTART_OCT)%tag = "Optimal Control"
    info(RESTART_PROJ)%tag = "GS for TDOutput"
    info(RESTART_PARTITION)%tag = "Mesh Partition"

    ! Default flags and directories (flags not yet used)
    info(:)%basedir = 'restart'
    info(:)%flags = 0
    info(RESTART_PARTITION)%flags = RESTART_FLAG_SKIP
    
    info(RESTART_GS)%dir = GS_DIR
    info(RESTART_UNOCC)%dir = GS_DIR
    info(RESTART_TD)%dir = TD_DIR
    info(RESTART_EM_RESP)%dir = EM_RESP_DIR
    info(RESTART_EM_RESP_FD)%dir = EM_RESP_FD_DIR
    info(RESTART_KDOTP)%dir = KDOTP_DIR
    info(RESTART_VIB_MODES)%dir = VIB_MODES_DIR
    info(RESTART_VDW)%dir = VDW_DIR
    info(RESTART_CASIDA)%dir = CASIDA_DIR
    info(RESTART_OCT)%dir = OCT_DIR
    info(RESTART_PROJ)%dir = GS_DIR
    info(RESTART_PARTITION)%dir = PARTITION_DIR

    ! Read input
    call messages_obsolete_variable(namespace, 'RestartFileFormat', 'RestartOptions')
    call messages_obsolete_variable(namespace, 'TmpDir', 'RestartOptions')
    call messages_obsolete_variable(namespace, 'RestartDir', 'RestartOptions')
    call messages_obsolete_variable(namespace, 'MeshPartitionRead', 'RestartOptions')
    call messages_obsolete_variable(namespace, 'MeshPartitionWrite', 'RestartOptions')
    call messages_obsolete_variable(namespace, 'MeshPartitionDir', 'RestartOptions')

    !%Variable RestartOptions
    !%Type block
    !%Section Execution::IO
    !%Description
    !% <tt>Octopus</tt> usually stores binary information, such as the wavefunctions, to be used
    !% in subsequent calculations. The most common example is the ground-state states
    !% that are used to start a time-dependent calculation. This variable allows to control
    !% where this information is written to or read from. The format of this block is the following:
    !% for each line, the first column indicates the type of data, the second column indicates
    !% the path to the directory that should be used to read and write that restart information, and the
    !% third column, which is optional, allows one to set some flags to modify the way how the data
    !% is read or written. For example, if you are running a time-dependent calculation, you can 
    !% indicate where <tt>Octopus</tt> can find the ground-state information in the following way:
    !%
    !% <tt>%RestartOptions
    !% <br>&nbsp;&nbsp;restart_gs | "gs_restart"
    !% <br>&nbsp;&nbsp;restart_td | "td_restart"
    !% <br>%</tt>
    !%
    !% The second line of the above example also tells <tt>Octopus</tt> that the time-dependent restart data
    !% should be read from and written to the "td_restart" directory.
    !%
    !% In case you want to change the path of all the restart directories, you can use the <tt>restart_all</tt> option.
    !% When using the <tt>restart_all</tt> option, it is still possible to have a different restart directory for specific
    !% data types. For example, when including the following block in your input file:
    !%
    !% <tt>%RestartOptions
    !% <br>&nbsp;&nbsp;restart_all | "my_restart"
    !% <br>&nbsp;&nbsp;restart_td&nbsp;  | "td_restart"
    !% <br>%</tt>
    !%
    !% the time-dependent restart information will be stored in the "td_restart" directory, while all the remaining 
    !% restart information will be stored in the "my_restart" directory.
    !%
    !% By default, the name of the "restart_all" directory is set to "restart".
    !% 
    !% Some <tt>CalculationMode</tt>s also take into account specific flags set in the third column of the <tt>RestartOptions</tt>
    !% block. These are used to determine if some specific part of the restart data is to be taken into account
    !% or not when reading the restart information. For example, when restarting a ground-state calculation, one can
    !% set the <tt>restart_rho</tt> flags, so that the density used is not built from the saved wavefunctions, but is
    !% instead read from the restart directory. In this case, the block should look like this:
    !% 
    !% <tt>%RestartOptions
    !% <br>&nbsp;&nbsp;restart_gs | "restart" | restart_rho
    !% <br>%</tt>
    !%
    !% A list of available flags is given below, but note that the code might ignore some of them, which will happen if they
    !% are not available for that particular calculation, or might assume some of them always present, which will happen
    !% in case they are mandatory.
    !% 
    !% Finally, note that all the restart information of a given data type is always stored in a subdirectory of the
    !% specified path. The name of this subdirectory is fixed and cannot be changed. For example, ground-state information 
    !% will always be stored in a subdirectory named "gs". This makes it safe in most situations to use the same path for
    !% all the data types. The name of these subdirectories is indicated in the description of the data types below.
    !%
    !% Currently, the available restart data types and flags are the following:
    !%Option restart_all 0
    !% (data type)
    !% Option to globally change the path of all the restart information.
    !%Option restart_gs  1
    !% (data type) 
    !% The data resulting from a ground-state calculation.
    !% This information is stored under the "gs" subdirectory.
    !%Option restart_unocc 2
    !% (data type) 
    !% The data resulting from an unoccupied states calculation. This information also corresponds to a ground state and 
    !% can be used as such, so it is stored under the same subdirectory as the one of restart_gs.
    !%Option restart_td 3
    !% (data type) 
    !% The data resulting from a real-time time-dependent calculation. 
    !% This information is stored under the "td" subdirectory.
    !%Option restart_em_resp 4
    !% (data type) 
    !% The data resulting from the calculation of the electromagnetic response using the Sternheimer approach. 
    !% This information is stored under the "em_resp" subdirectory.
    !%Option restart_em_resp_fd 5
    !% (data type) 
    !% The data resulting from the calculation of the electromagnetic response using finite-differences. 
    !% This information is stored under the "em_resp_fd" subdirectory.
    !%Option restart_kdotp 6
    !% (data type) 
    !% The data resulting from the calculation of effective masses by k.p perturbation theory.
    !% This information is stored under the "kdotp" subdirectory.
    !%Option restart_vib_modes 7
    !% (data type) 
    !% The data resulting from the calculation of vibrational modes.
    !% This information is stored under the "vib_modes" subdirectory.
    !%Option restart_vdw 8
    !% (data type) 
    !% The data resulting from the calculation of van der Waals coefficients.
    !% This information is stored under the "vdw" subdirectory.
    !%Option restart_casida 9
    !% (data type) 
    !% The data resulting from a Casida calculation.
    !% This information is stored under the "casida" subdirectory.
    !%Option restart_oct 10
    !% (data type) 
    !% The data for optimal control calculations.
    !% This information is stored under the "opt-control" subdirectory.
    !%Option restart_partition 11
    !% (data type) 
    !% The data for the mesh partitioning.
    !% This information is stored under the "partition" subdirectory.
    !%Option restart_proj 12
    !% (data type)
    !% The ground-state to be used with the td_occup and populations options of <tt>TDOutput</tt>.
    !% This information should be a ground state, so the "gs" subdirectory is used.
    !%Option restart_states 1
    !% (flag)
    !% Read the electronic states. (not yet implemented)
    !%Option restart_rho 2
    !% (flag)
    !% Read the electronic density.
    !%Option restart_vhxc 4
    !% (flag)
    !% Read the Hartree and XC potentials.
    !%Option restart_mix 8
    !% (flag)
    !% Read the SCF mixing information.   
    !%Option restart_skip 16
    !% (flag)
    !% This flag allows to selectively skip the reading and writting of specific restart information.
    !%End
    set = .false.
    if(parse_block(namespace, 'RestartOptions', blk) == 0) then

      default_basedir = 'restart'

      do iline = 1, parse_block_n(blk)
        n_cols = parse_block_cols(blk,iline-1)

        call parse_block_integer(blk, iline-1, 0, data_type)
        if (data_type < 0 .or. data_type > RESTART_N_DATA_TYPES) then
          call messages_input_error(namespace, 'RestartOptions', "Invalid data type", row=iline-1, column=0)
        end if
        if (data_type == 0) then
          call parse_block_string(blk, iline-1, 1, default_basedir)
        else
          set(data_type) = .true.
          call parse_block_string(blk, iline-1, 1, info(data_type)%basedir)

          if (n_cols > 2) call parse_block_integer(blk, iline-1, 2, info(data_type)%flags)
        end if

      end do
      call parse_block_end(blk)

      where (.not. set)
        info(:)%basedir = default_basedir
      end where
    end if

    POP_SUB(restart_module_init)
  end subroutine restart_module_init


  ! ---------------------------------------------------------
  !> Initializes a restart object.
  subroutine restart_init(restart, namespace, data_type, type, mc, ierr, mesh, dir, exact)
    type(restart_t),             intent(out) :: restart   !< Restart information
    type(namespace_t), target,   intent(in)  :: namespace
    integer,                     intent(in)  :: data_type !< Restart data type (RESTART_GS, RESTART_TD, etc)
    integer,                     intent(in)  :: type      !< Is this restart used for dumping (type = RESTART_TYPE_DUMP)
                                                          !! or for loading (type = RESTART_TYPE_LOAD)?
    type(multicomm_t), target,   intent(in)  :: mc        !< The multicommunicator in charge of handling this restart.
    integer,                     intent(out) :: ierr      !< Error code, if any. Required for LOAD, should not be present for DUMP.
    type(mesh_t),      optional, intent(in)  :: mesh      !< If present, depending on the type of restart, the mesh 
                                                          !! information is either dumped or the mesh compatibility is checked.
    character(len=*),  optional, intent(in)  :: dir       !< Directory where to find the restart data. It is mandatory if 
                                                          !! data_type=RESTART_UNDEFINED and is ignored in all the other cases.
    logical,           optional, intent(in)  :: exact     !< If loading the restart information, should the mesh be 
                                                          !! exactly the same or not?

    logical :: grid_changed, grid_reordered, restart_write, dir_exists
    character(len=20) :: tag
    character(len=MAX_PATH_LEN) :: basedir, dirname

    PUSH_SUB(restart_init)

    ierr = 0

    ! Sanity checks
    if (present(exact) .and. .not. present(mesh)) then
      message(1) = "Error in restart_init: the 'exact' optional argument requires a mesh."
      call messages_fatal(1)
    end if

    restart%has_mesh = present(mesh)

    ! Some initializations
    restart%type = type
    restart%mc => mc
    call mpi_grp_init(restart%mpi_grp, mc%master_comm)
    restart%format = io_function_fill_how("Binary")
    if (data_type < RESTART_UNDEFINED .and. data_type > RESTART_N_DATA_TYPES) then
      message(1) = "Illegal data_type in restart_init"
      call messages_fatal(1)
    end if
    restart%data_type = data_type
    restart%namespace => namespace

    select case (restart%type)
    case (RESTART_TYPE_DUMP)
      !%Variable RestartWrite
      !%Type logical
      !%Default true
      !%Section Execution::IO
      !%Description
      !% If this variable is set to no, restart information is not
      !% written. Note that some run modes will ignore this
      !% option and write some restart information anyway.
      !%End

      call parse_variable(namespace, 'RestartWrite', .true., restart_write)
      restart%skip = .not. restart_write

      if(restart%skip) then
        message(1) = 'Restart information will not be written.'
        call messages_warning(1)
      end if
        
    case (RESTART_TYPE_LOAD)
      ! This is set to true as an error condition, checked by assertions in some routines.
      restart%skip = .false.
      
    case default
      message(1) = "Unknown restart type in restart_init"
      call messages_fatal(1)
    end select


    ! If the restart data type is not defined, the directories should be set explicitly
    if (restart%data_type == RESTART_UNDEFINED) then
      ASSERT(present(dir))
      basedir = dir
      dirname = ""
    else
      basedir = info(restart%data_type)%basedir
      if (index(basedir, '/', .true.) /= len_trim(basedir)) then
        basedir = trim(basedir)//"/"
      end if
      dirname = info(restart%data_type)%dir
    end if

    ! Set final path
    restart%dir = trim(basedir)//trim(dirname)

    ! Remove any trailing "/" from the path (all the routines from this module should add the trailing "/" when needed)
    if (index(restart%dir, '/', .true.) == len_trim(restart%dir)) then
      restart%dir = restart%dir(1:len_trim(restart%dir)-1)
    end if

    ! Set initial path to the working directory
    restart%pwd = restart%dir

    ! Check if the directory already exists and create it if necessary
    dir_exists = io_dir_exists(trim(restart%pwd), namespace)
    if (restart%type == RESTART_TYPE_DUMP .and. .not. dir_exists) then
      call io_mkdir(trim(restart%pwd), namespace, parents=.true.)
    end if

    if (restart%data_type == RESTART_UNDEFINED) then
      tag = "Some "
    else
      tag = info(data_type)%tag
    end if

    select case (restart%type)
    case (RESTART_TYPE_DUMP)
      if (.not. restart%skip) then
        message(1) = "Info: "//trim(tag)//" restart information will be written to '"//trim(restart%pwd)//"'."
        call messages_info(1)

        ! Dump the grid information. The main parameters of the grid should not change
        ! during the calculation, so we should only need to dump it once.
        if (present(mesh)) then
          call index_dump_lxyz(mesh%idx, mesh%np_part_global, restart%pwd, restart%mpi_grp, &
            restart%namespace, ierr)
          if (ierr /= 0) then
            message(1) = "Unable to write index map to '"//trim(restart%pwd)//"'."
            call messages_fatal(1)
          end if

          call mesh_write_fingerprint(mesh, restart%pwd, "grid", restart%mpi_grp, namespace, ierr)
          if (ierr /= 0) then
            message(1) = "Unable to write mesh fingerprint to '"//trim(restart%pwd)//"/grid'."
            call messages_fatal(1)
          end if
        end if
        
      end if

    case (RESTART_TYPE_LOAD)
      if (.not. dir_exists) then
        ierr = 1
        restart%skip = .true.

        message(1) = "Could not find '"//trim(restart%pwd)//"' directory for restart."
        message(2) = "No restart information will be read."
        call messages_warning(2)

      else
        message(1) = "Info: "//trim(tag)//" restart information will be read from '"//trim(restart%pwd)//"'."
        call messages_info(1)

        if (present(mesh)) then
          call mesh_check_dump_compatibility(mesh, restart%pwd, "grid", restart%namespace, &
            restart%mpi_grp, grid_changed, grid_reordered, restart%map, ierr)

          ! Check whether an error occurred. In this case we cannot read.
          if (ierr /= 0) then
            if (ierr == -1) then
              message(1) = "Unable to check mesh compatibility: unable to read mesh fingerprint"
              message(2) = "in '"//trim(restart%pwd)//"'."
            else if (ierr > 0) then
              message(1) = "Mesh from current calculation is not compatible with mesh found in"
              message(2) = "'"//trim(restart%pwd)//"'."
            end if
            message(3) = "No restart information will be read."
            call messages_warning(3)
            ierr = 1
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
              ierr = 1
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
    if(restart%mpi_grp%size > 1) &
      call MPI_Barrier(restart%mpi_grp%comm, mpi_err)
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
        call io_rm(trim(restart%pwd)//"/loading", restart%namespace)
      case (RESTART_TYPE_DUMP)
        call io_rm(trim(restart%pwd)//"/dumping", restart%namespace)
        message(1) = "Info: Finished writing information to '"//trim(restart%dir)//"'."
      end select
      call messages_info(1)
    end if

    restart%type = 0
    restart%data_type = 0
    restart%skip = .true.
    SAFE_DEALLOCATE_A(restart%map)
    restart%has_mesh = .false.
    nullify(restart%mc)
    
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
    character(len=MAX_PATH_LEN) :: restart_dir

    PUSH_SUB(restart_dir)

    restart_dir = io_workpath(restart%pwd, restart%namespace)

    POP_SUB(restart_dir)
  end function restart_dir


  ! ---------------------------------------------------------
  !> Change the restart directory to dirname, where "dirname" 
  !! is a subdirectory of the base restart directory.
  subroutine restart_open_dir(restart, dirname, ierr)
    type(restart_t),  intent(inout) :: restart
    character(len=*), intent(in)    :: dirname
    integer,          intent(out)   :: ierr

    PUSH_SUB(restart_open_dir)

    ASSERT(.not. restart%skip)

    ierr = 0

    select case (restart%type)
    case (RESTART_TYPE_DUMP)
      call restart_mkdir(restart, dirname)
    case (RESTART_TYPE_LOAD)
      if (.not. loct_dir_exists(trim(restart%dir)//"/"//trim(dirname))) then
        ierr = 1
      end if
    end select

    if (ierr == 0) then
      if (index(dirname, '/', .true.) == len_trim(dirname)) then
        restart%pwd = trim(restart%dir)//"/"//dirname(1:len_trim(dirname)-1)
      else
        restart%pwd = trim(restart%dir)//"/"//trim(dirname)
      end if
    end if

    POP_SUB(restart_open_dir)
  end subroutine restart_open_dir


  ! ---------------------------------------------------------
  !> Change back to the base directory. To be called after restart_open_dir.
  subroutine restart_close_dir(restart)
    type(restart_t),  intent(inout) :: restart

    PUSH_SUB(restart_close_dir)

    ASSERT(.not. restart%skip)

    restart%pwd = restart%dir

    POP_SUB(restart_close_dir)
  end subroutine restart_close_dir


  ! ---------------------------------------------------------
  !> Make directory "dirname" inside the current restart directory.
  subroutine restart_mkdir(restart, dirname)
    type(restart_t),  intent(in) :: restart
    character(len=*), intent(in) :: dirname

    PUSH_SUB(restart_mkdir)

    ASSERT(.not. restart%skip)

    ASSERT (restart%type == RESTART_TYPE_DUMP)

    call io_mkdir(trim(restart%pwd)//"/"//trim(dirname), restart%namespace, parents=.true.)

    POP_SUB(restart_mkdir)
  end subroutine restart_mkdir


  ! ---------------------------------------------------------  
  !> Remove directory or file "name" that is located inside the current restart directory.
  subroutine restart_rm(restart, name)
    type(restart_t),  intent(in) :: restart
    character(len=*), intent(in) :: name

    ASSERT(.not. restart%skip)
    ASSERT(restart%type == RESTART_TYPE_DUMP)

    PUSH_SUB(restart_rm)

    call io_rm(trim(restart%pwd)//"/"//trim(name), restart%namespace)

    POP_SUB(restart_rm)
  end subroutine restart_rm


  ! ---------------------------------------------------------
  !> Open file "filename" found inside the current restart
  !! directory. Depending on the type of restart, the file
  !! will be used for reading  (RESTART_TYPE_LOAD), or for
  !! writing (RESTART_TYPE_DUMP).
  !! Note that this function only be used to open *formatted* files and that
  !! unformatted data should be accessed using the restart_*_binary subroutines.
  function restart_open(restart, filename, status, position, silent)
    type(restart_t),            intent(in) :: restart
    character(len=*),           intent(in) :: filename
    character(len=*), optional, intent(in) :: status
    character(len=*), optional, intent(in) :: position
    logical,          optional, intent(in) :: silent
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
      call messages_fatal(1)
    end select

    if (present(status)) status_ = status

    restart_open = io_open(trim(restart%pwd)//"/"//trim(filename), restart%namespace, &
      action=trim(action), status=trim(status_), &
      die=die, position=position, form="formatted", grp=restart%mpi_grp)

    if (restart_open < 0 .and. .not. optional_default(silent, .false.)) then    
      message(1) = "Unable to open file '"//trim(restart%pwd)//"/"//trim(filename)//"'."
      call messages_warning(1)
    end if

    POP_SUB(restart_open)
  end function restart_open


  ! ---------------------------------------------------------
  subroutine restart_write(restart, iunit, lines, nlines, ierr)
    type(restart_t),  intent(in)  :: restart
    integer,          intent(in)  :: iunit
    character(len=*), intent(in)  :: lines(:)
    integer,          intent(in)  :: nlines
    integer,          intent(out) :: ierr

    integer :: iline

    PUSH_SUB(restart_write)

    if (iunit > 0) then
      ierr = 0
      if (mpi_grp_is_root(restart%mpi_grp)) then
        do iline = 1, nlines
          write(iunit,"(a)") trim(lines(iline))
        end do
      end if
    else
      ierr = 1
    end if

    POP_SUB(restart_write)
  end subroutine restart_write


  ! ---------------------------------------------------------
  subroutine restart_read(restart, iunit, lines, nlines, ierr)
    type(restart_t),  intent(in)  :: restart
    integer,          intent(in)  :: iunit
    character(len=*), intent(out) :: lines(:)
    integer,          intent(in)  :: nlines
    integer,          intent(out) :: ierr

    PUSH_SUB(restart_read)

    call iopar_read(restart%mpi_grp, iunit, lines, nlines, ierr)

    POP_SUB(restart_read)
  end subroutine restart_read


  ! ---------------------------------------------------------
  !> Close a file previously opened with restart_open.
  subroutine restart_close(restart, iunit)
    type(restart_t), intent(in)    :: restart
    integer,         intent(inout) :: iunit

    PUSH_SUB(restart_close)

    if (iunit > 0) call io_close(iunit, restart%mpi_grp)

    POP_SUB(restart_close)
  end subroutine restart_close


  ! ---------------------------------------------------------
  !> Returns true if the restart information should neither be
  !! read nor written. This might happen because the user
  !! chose not to write any restart information, or because
  !! the restart information is not available for reading.
  logical pure function restart_skip(restart)
    type(restart_t), intent(in) :: restart

    restart_skip = restart%skip .or. restart_has_flag(restart, RESTART_FLAG_SKIP)

  end function restart_skip


  ! ---------------------------------------------------------
  !> Returns true if...
  logical pure function restart_has_flag(restart, flag)
    type(restart_t), intent(in) :: restart
    integer,         intent(in) :: flag

    restart_has_flag = bitand(info(restart%data_type)%flags, flag) /= 0

  end function restart_has_flag


  ! ---------------------------------------------------------
  !> Returns true if the restart was from a different order of mesh points
  logical pure function restart_has_map(restart)
    type(restart_t), intent(in) :: restart

    restart_has_map = allocated(restart%map)

  end function restart_has_map


  ! ---------------------------------------------------------
  !> Returns true if...
  logical pure function restart_are_basedirs_equal(type1, type2)
    integer, intent(in) :: type1
    integer, intent(in) :: type2

    restart_are_basedirs_equal = trim(info(type1)%basedir) == trim(info(type2)%basedir)

  end function restart_are_basedirs_equal


#include "undef.F90"
#include "real.F90"
#include "restart_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "restart_inc.F90"

end module restart_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
