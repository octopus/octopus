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

module restart_m
  use curvilinear_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use io_function_m
  use io_m
  use loct_m
  use loct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use mesh_init_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_dim_m
  use states_calc_m
  use string_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                  &
    restart_init,            &
    clean_stop,              &
    read_free_states,        &
    restart_write,           &
    restart_read,            &
    restart_read_ob_intf,    &
    restart_read_ob_central, &
    restart_format,          &
    restart_dir,             &
    restart_look_and_read,   &
    restart_read_user_def_orbitals,    &
    drestart_write_function, &
    zrestart_write_function, &
    drestart_read_function,  &
    zrestart_read_function,  & 
    drestart_write_lr_rho,   &
    zrestart_write_lr_rho,   &
    drestart_read_lr_rho,    &
    zrestart_read_lr_rho
  

  integer, parameter :: &
    RESTART_NETCDF = 2, &
    RESTART_BINARY = 3

  integer           :: restart_format
  character(len=32) :: restart_dir

contains

  ! returns true if a file named stop exists
  function clean_stop()
    logical clean_stop, file_exists

    clean_stop = .false.
    inquire(file='stop', exist=file_exists)
    if(file_exists) then
      message(1) = 'Clean STOP'
      call write_warning(1)
      clean_stop = .true.
#if defined(HAVE_MPI)
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
      if(mpi_grp_is_root(mpi_world)) call loct_rm('stop')
    end if

    return
  end function clean_stop


  ! ---------------------------------------------------------
  ! read restart format information
  subroutine restart_init

    integer :: default, parsed

    call push_sub('restart.restart_init')

    !%Variable RestartFileFormat
    !%Type integer
    !%Default restart_netcdf
    !%Section Execution::IO
    !%Description
    !% Determines in which format the restart files should be written. The default
    !% is binary files (restart_binary).
    !%Option restart_netcdf 2
    !% NetCDF (platform independent) format. This requires the NETCDF library.
    !%Option restart_binary 3
    !% Octopus Binary Format, the new and more flexible binary format
    !% of Octopus. It is faster and produces smaller files than NetCDF
    !% and it is platform-independent.
    !%End

    default = RESTART_BINARY

    call loct_parse_int(datasets_check('RestartFileFormat'), default, parsed)
#ifndef HAVE_NETCDF
    if (parsed == RESTART_NETCDF) then 
      write(message(1),'(a)') 'Error: Octopus was compiled without NetCDF support but'
      write(message(2),'(a)') 'NetCDF restart files were requested.'
      call write_fatal(2)
    end if
#endif
    if(.not.varinfo_valid_option('RestartFileFormat', parsed)) call input_error('RestartFileFormat')

    select case(parsed) 
    case(RESTART_NETCDF)
      restart_format = io_function_fill_how("NETCDF")
    case(RESTART_BINARY)
      restart_format = io_function_fill_how("Binary")
    end select

    !%Variable RestartDir
    !%Type string
    !%Default ''
    !%Section Execution::IO
    !%Description
    !% When Octopus reads restart files, e. g. when running a time-propagation
    !% after a ground-state calculation, these files will be read from
    !% <tt>&lt;RestartDir&gt/</tt>. Usually, <tt>RestartDir</tt> is
    !% <tt>TmpDir</tt> but in a transport calculation, the output of
    !% a periodic dataset is required to calculate the extended ground state.
    !%End
    call loct_parse_string(datasets_check('RestartDir'), tmpdir, restart_dir)
    ! Append "/" if necessary.
    if(scan(restart_dir, '/', .true.).lt.1) then
      restart_dir = trim(restart_dir)//'/'
    end if

    call pop_sub()
  end subroutine restart_init


  ! ---------------------------------------------------------
  subroutine restart_look_and_read(st, gr, geo, is_complex, specify_dir)
    type(states_t),             intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    type(geometry_t),           intent(in)    :: geo
    logical, optional,          intent(in)    :: is_complex
    character(len=*), optional, intent(in)    :: specify_dir

    integer :: kpoints, dim, nst, ierr
    character(len=80) dir

    call push_sub('restart.restart_look_and_read')

    if(present(specify_dir)) then
       dir = specify_dir
    else
       dir = restart_dir
    endif

    !check how many wfs we have
    call states_look(trim(dir)//'gs', gr%mesh%mpi_grp, kpoints, dim, nst, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(dir)//'gs".'
      call write_fatal(1)
    end if

    st%nst    = nst
    st%st_end = nst
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_DEALLOCATE_P(st%occ)

    if (present(is_complex)) then
      if ( is_complex ) then 
        call states_allocate_wfns(st, gr%mesh, M_CMPLX)
      else 
        call states_allocate_wfns(st, gr%mesh, M_REAL)
      end if
    else
      ! allow states_allocate_wfns to decide for itself whether complex or real needed
      call states_allocate_wfns(st, gr%mesh)
    endif

    ALLOCATE(st%eigenval(st%nst, st%d%nik), st%nst*st%d%nik)
    ALLOCATE(st%occ(st%nst, st%d%nik), st%nst*st%d%nik)

    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%spin(3, st%nst, st%d%nik), st%nst*st%d%nik*3)
      st%spin = M_ZERO
    end if
    st%eigenval = huge(REAL_PRECISION)
    st%occ      = M_ZERO

    ! load wave-functions
    call restart_read(trim(dir)//'gs', st, gr, geo, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(dir)//"gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    call pop_sub()
  end subroutine restart_look_and_read


  ! ---------------------------------------------------------
  subroutine restart_write(dir, st, gr, ierr, iter, lr)
    character(len=*),  intent(in)  :: dir
    type(states_t),    intent(in)  :: st
    type(grid_t),      intent(in)  :: gr
    integer,           intent(out) :: ierr
    integer, optional, intent(in)  :: iter
    !if this next argument is present, the lr wfs are stored instead of the gs wfs
    type(lr_t), optional, intent(in)  :: lr 

    integer :: iunit, iunit2, iunit_mesh, iunit_states, err, ik, ist, idim, i
    character(len=80) :: filename, mformat
    logical :: wfns_are_associated, lr_wfns_are_associated

    call push_sub('restart.restart_write')

    wfns_are_associated = (associated(st%dpsi) .and. st%wfs_type == M_REAL) .or. &
         (associated(st%zpsi) .and. st%wfs_type == M_CMPLX)
    ASSERT(wfns_are_associated)

    if(present(lr)) then 
      lr_wfns_are_associated = (associated(lr%ddl_psi) .and. st%wfs_type == M_REAL) .or. &
                (associated(lr%zdl_psi) .and. st%wfs_type == M_CMPLX)
      ASSERT(lr_wfns_are_associated)
    endif

    mformat = '(f20.14,a,f20.14,a,4(f20.14,a),i10.10,3(a,i8))'
    ierr = 0

    call block_signals()
    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(dir, is_tmp=.true.)

      iunit = io_open(trim(dir)//'/wfns', action='write', is_tmp=.true.)
      write(iunit,'(a)') '#     #kpoint            #st            #dim    filename'
      write(iunit,'(a)') '%Wavefunctions'

      iunit2 = io_open(trim(dir)//'/occs', action='write', is_tmp=.true.)
      write(iunit2,'(a)') '# occupations | eigenvalue[a.u.] | K-Points | K-Weights | filename | ik | ist | idim'
      write(iunit2,'(a)') '%Occupations_Eigenvalues_K-Points'

      iunit_mesh = io_open(trim(dir)//'/mesh', action='write', is_tmp=.true.)
      write(iunit_mesh,'(a)') '# This file contains the necessary information to generate the'
      write(iunit_mesh,'(a)') '# mesh with which the functions in this directory were calculated,'
      write(iunit_mesh,'(a)') '# except for the geometry of the system.'
      call curvilinear_dump(gr%cv, iunit_mesh)
      call simul_box_dump(gr%sb, iunit_mesh)
      call mesh_dump(gr%mesh, iunit_mesh)
      call io_close(iunit_mesh)

      iunit_mesh = io_open(trim(dir)//'/Lxyz', action='write', is_tmp=.true.)
      call mesh_Lxyz_dump(gr%mesh, iunit_mesh)
      call io_close(iunit_mesh)      

      iunit_states = io_open(trim(dir)//'/states', action='write', is_tmp=.true.)
      call states_dump(st, iunit_states)
      call io_close(iunit_states)      
    end if

    i = 1
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          write(filename,'(i10.10)') i

          if(mpi_grp_is_root(mpi_world)) then
            write(unit=iunit,  fmt=*) ik, ' | ', ist, ' | ', idim, ' | "', trim(filename), '"'
            write(unit=iunit2, fmt=mformat) st%occ(ist,ik), ' | ', st%eigenval(ist, ik), ' | ',   &
                 st%d%kpoints(1,ik), ' | ', st%d%kpoints(2,ik), ' | ', st%d%kpoints(3,ik) , ' | ', &
                 st%d%kweights(ik), ' | ', i, ' | ', ik, ' | ', ist, ' | ', idim
          end if

          if(st%st_start <= ist .and. ist <= st%st_end) then
            if( .not. present(lr) ) then 
              if(st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
                if (st%wfs_type == M_REAL) then
                  call drestart_write_function(dir, filename, gr, st%dpsi(:, idim, ist, ik), err, size(st%dpsi,1))
                else
                  call zrestart_write_function(dir, filename, gr, st%zpsi(:, idim, ist, ik), err, size(st%zpsi,1))
                end if
                if(err == 0) ierr = ierr + 1
              end if
            else
              if (st%wfs_type == M_REAL) then
                call drestart_write_function(dir, filename, gr, lr%ddl_psi(:, idim, ist, ik), err, size(st%dpsi,1))
              else
                call zrestart_write_function(dir, filename, gr, lr%zdl_psi(:, idim, ist, ik), err, size(st%zpsi,1))
              end if
              if(err == 0) ierr = ierr + 1
            end if
          end if
#if defined(HAVE_MPI)
          call MPI_Barrier(MPI_COMM_WORLD, mpi_err) ! now we all wait
#endif
          i = i + 1
        end do
      end do
    end do

    ! do NOT use st%lnst here as it is not (st%st_end - st%st_start + 1)
    if(ierr == st%d%kpt%nlocal*(st%st_end - st%st_start + 1)*st%d%dim) ierr = 0 ! All OK

    if(mpi_grp_is_root(mpi_world)) then
      write(iunit,'(a)') '%'
      if(present(iter)) write(iunit,'(a,i7)') 'Iter = ', iter
      write(iunit2, '(a)') '%'

      call io_close(iunit)
      call io_close(iunit2)
    end if

#if defined(HAVE_MPI)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err) ! Since some processors did more than others...
#endif
    call unblock_signals()

    call pop_sub()
  end subroutine restart_write


  ! ---------------------------------------------------------
  ! Reads the interface regions of the ground state wavefunctions of
  ! an open boundaries calculation.
  ! Returns (in ierr)
  ! <0 => Fatal error,
  ! =0 => read all wave-functions,
  ! >0 => could only read x wavefunctions.
  subroutine restart_read_ob_intf(dir, st, gr, ierr)
    character(len=*), intent(in)    :: dir
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in)    :: gr
    integer,          intent(out)   :: ierr

    integer              :: io_wfns, io_occs, io_mesh, np, lead_np, i, err
    integer              :: ik, ist, idim, tnp
    character            :: char
    character(len=256)   :: line, filename
    CMPLX, allocatable   :: tmp(:)
    type(mesh_t)         :: gs_mesh
    type(mpi_grp_t)      :: mpi_grp

    call push_sub('restart.restart_read_ob_intf')

    write(message(1), '(a,i5)') 'Info: Reading ground state interface wave functions'
    call write_info(1)

    mpi_grp = mpi_world
    ! Sanity check.
    ASSERT(associated(st%ob_intf_psi))

    ierr = 0

    ! Read the mesh of the ground state calculation.
    io_mesh = io_open(trim(dir)//'/mesh', action='read', status='old', die=.false., is_tmp=.true.)
    if(io_mesh.lt.0) then
      ierr = -1
      call io_close(io_mesh, grp=mpi_grp)
      call pop_sub()
      return
    end if
      
    gs_mesh%sb => gr%mesh%sb
    call mesh_init_from_file(gs_mesh, io_mesh)

    io_wfns  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(io_wfns.lt.0) then
      ierr = -1
    end if

    io_occs = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(io_occs.lt.0) then
      if(io_wfns.gt.0) call io_close(io_wfns, grp=mpi_grp)
      ierr = -1
    end if

    if(ierr.ne.0) then
      write(message(1),'(a)') 'Could not load ground state interface wave functions.'
      call write_info(1)
      call messages_print_stress(stdout)
      call pop_sub()
      return
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, io_wfns, line, err); call iopar_read(mpi_grp, io_wfns, line, err)
    call iopar_read(mpi_grp, io_occs, line, err); call iopar_read(mpi_grp, io_occs, line, err)

    lead_np = gr%mesh%lead_unit_cell(LEFT)%np
    np = gr%intf(LEFT)%np
    ALLOCATE(tmp(gr%mesh%np+2*lead_np), gr%mesh%np+2*lead_np)
    
    do
      call iopar_read(mpi_grp, io_wfns, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char.eq.'%') exit

      call iopar_backspace(mpi_grp, io_wfns)

      call iopar_read(mpi_grp, io_wfns, line, err)
      read(line, *) ik, char, ist, char, idim, char, filename

      call iopar_read(mpi_grp, io_occs, line, err)
      read(line, *) st%occ(ist, ik), char, st%eigenval(ist, ik)

      if(ist.ge.st%st_start .and. ist.le.st%st_end .and. &
         ik.ge.st%d%kpt%start .and. ik.le.st%d%kpt%end) then
        ! Open boundaries imply complex wavefunctions.
        call zrestart_read_function(dir, filename, gs_mesh, tmp, err)
        ! Here we use the fact that transport is in x-direction (x is the index
        ! running slowest).
        ! Outer block.
        tnp = gr%mesh%np+2*lead_np
        st%ob_intf_psi(1:np, OUTER, idim, ist, ik, LEFT)  = tmp(lead_np-np+1:lead_np)
        st%ob_intf_psi(1:np, OUTER, idim, ist, ik, RIGHT) = tmp(tnp-lead_np+1:tnp-lead_np+np)
        ! Inner block.
        st%ob_intf_psi(1:np, INNER, idim, ist, ik, LEFT)  = tmp(lead_np+1:lead_np+np)
        st%ob_intf_psi(1:np, INNER, idim, ist, ik, RIGHT) = tmp(tnp-lead_np-np+1:tnp-lead_np)
        if(err.le.0) then
          ierr = ierr + 1
        end if
      end if
    end do
    
#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
      ierr = err
    end if
    if(st%d%kpt%parallel) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif

    if(ierr.eq.0) then
      ierr = -1
      call messages_print_stress(stdout, 'Loading restart information')
      message(1) = 'No interface wave functions could be read.'
      call write_info(1)
      call messages_print_stress(stdout)
    else
      if(ierr.eq.st%nst*st%d%nik*st%d%dim) then
        ierr = 0
      else
        call messages_print_stress(stdout, 'Loading restart information')
        write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' interface wave functions out of ', &
          st%nst*st%d%nik*st%d%dim, ' could be read.'
        call write_info(1)
        call messages_print_stress(stdout)
      end if
    end if

    call io_close(io_wfns, mpi_grp)
    call io_close(io_occs, mpi_grp)

    SAFE_DEALLOCATE_A(tmp)

    call pop_sub()
  end subroutine restart_read_ob_intf


  ! ---------------------------------------------------------
  ! Reads the central regions of the ground state wavefunctions of
  ! an open boundaries calculation.
  ! Returns (in ierr)
  ! <0 => Fatal error,
  ! =0 => read all wave-functions,
  ! >0 => could only read x wavefunctions.
  subroutine restart_read_ob_central(dir, st, gr, ierr)
    character(len=*), intent(in)    :: dir
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in)    :: gr
    integer,          intent(out)   :: ierr

    logical              :: psi_allocated
    integer              :: io_wfns, io_occs, io_mesh, lead_np, i, err
    integer              :: ik, ist, idim
    character            :: char
    character(len=256)   :: line, filename
    CMPLX, allocatable   :: tmp(:)
    type(mesh_t)         :: gs_mesh
    type(mpi_grp_t)      :: mpi_grp
    FLOAT                :: k_x, k_y, k_z

    call push_sub('restart.restart_read_ob_central')

    write(message(1), '(a,i5)') 'Info: Loading restart information'
    call write_info(1)

    mpi_grp = mpi_world

    ! Sanity check.
    psi_allocated = (associated(st%dpsi) .and. st%wfs_type.eq.M_REAL) .or. &
      (associated(st%zpsi) .and. st%wfs_type.eq.M_CMPLX)
    ASSERT(psi_allocated)

    ierr = 0

    ! Read the mesh of the ground state calculation.
    io_mesh = io_open(trim(dir)//'/mesh', action='read', status='old', die=.false., is_tmp=.true.)
    if(io_mesh.lt.0) then
      ierr = -1
      call io_close(io_mesh, grp=mpi_grp)
      call pop_sub()
      return
    end if
    gs_mesh%sb => gr%mesh%sb
    call mesh_init_from_file(gs_mesh, io_mesh)

    io_wfns  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp=.true., grp=mpi_grp)
    if(io_wfns.lt.0) then
      ierr = -1
    end if

    io_occs = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp = .true., grp = mpi_grp)
    if(io_occs.lt.0) then
      if(io_wfns.gt.0) call io_close(io_wfns, grp=mpi_grp)
      ierr = -1
    end if

    if(ierr.ne.0) then
      write(message(1),'(a)') 'Could not load any previous restart information.'
      call write_info(1)
      call messages_print_stress(stdout)
      call pop_sub()
      return
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, io_wfns, line, err); call iopar_read(mpi_grp, io_wfns, line, err)
    call iopar_read(mpi_grp, io_occs, line, err); call iopar_read(mpi_grp, io_occs, line, err)

    lead_np = gr%mesh%lead_unit_cell(LEFT)%np
    ALLOCATE(tmp(gr%mesh%np+2*lead_np), gr%mesh%np+2*lead_np)
    
    do
      call iopar_read(mpi_grp, io_wfns, line, i)
      if(i.ne.0) exit
      read(line, '(a)') char
      if(char.eq.'%') exit

      call iopar_backspace(mpi_grp, io_wfns)

      call iopar_read(mpi_grp, io_wfns, line, err)
      read(line, *) ik, char, ist, char, idim, char, filename

      call iopar_read(mpi_grp, io_occs, line, err)
      read(line, *) st%occ(ist, ik), char, st%eigenval(ist, ik), char, &
                    k_x, char, k_y, char, k_z, char, st%d%kweights(ik)
      st%d%kpoints(:, ik) = (/k_x, k_y, k_z/)

      if(ist.ge.st%st_start .and. ist.le.st%st_end .and. &
         ik.ge.st%d%kpt%start .and. ik.le.st%d%kpt%end) then
        ! Open boundaries imply complex wavefunctions.
        call zrestart_read_function(dir, filename, gs_mesh, tmp, err)
        ! Here we use the fact that transport is in x-direction (x is the index
        ! running slowest).
        st%zpsi(1:gr%mesh%np, idim, ist, ik) = tmp(lead_np+1:gr%mesh%np+lead_np)
        if(err.le.0) then
          ierr = ierr + 1
        end if
      end if
    end do
    
#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
      ierr = err
    end if
    if(st%d%kpt%parallel) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif
    if(ierr.eq.0) then
      ierr = -1
      call messages_print_stress(stdout, 'Loading restart information')
      message(1) = 'No files could be read. No restart information can be used.'
      call write_info(1)
      call messages_print_stress(stdout)
    else
      if(ierr.eq.st%nst*st%d%nik*st%d%dim) then
        ierr = 0
      else
        call messages_print_stress(stdout, 'Loading restart information')
        write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' files out of ', &
          st%nst*st%d%nik*st%d%dim, ' could be read.'
        call write_info(1)
        call messages_print_stress(stdout)
      end if
    end if

    call io_close(io_wfns, mpi_grp)
    call io_close(io_occs, mpi_grp)

    SAFE_DEALLOCATE_A(tmp)

    call pop_sub()
  end subroutine restart_read_ob_central


  ! ---------------------------------------------------------
  ! returns
  ! <0 => Fatal error
  ! =0 => read all wave-functions
  ! >0 => could only read x wavefunctions
  subroutine restart_read(dir, st, gr, geo, ierr, iter, lr)
    character(len=*),  intent(in)  :: dir
    type(states_t), intent(inout)  :: st
    type(grid_t),      intent(in)  :: gr
    type(geometry_t),  intent(in)  :: geo
    integer,           intent(out) :: ierr
    integer, optional, intent(inout) :: iter
    !if this next argument is present, the lr wfs are read instead of the gs wfs
    type(lr_t), optional, intent(inout)  :: lr 

    integer              :: iunit, iunit2, iunit_mesh, err, ik, ist, idim, i
    character(len=12)    :: filename
    character(len=1)     :: char
    logical, allocatable :: filled(:, :, :)
    character(len=256)   :: line
    character(len=50)    :: str

    FLOAT, allocatable   :: dphi(:)
    CMPLX, allocatable   :: zphi(:)
    type(mesh_t)         :: old_mesh
    type(curvilinear_t)   :: old_cv
    type(simul_box_t)    :: old_sb
    logical              :: mesh_change, full_interpolation, gs_allocated, lr_allocated

    call push_sub('restart.restart_read')

    if(.not. present(lr)) then 
      write(message(1), '(a,i5)') 'Info: Loading restart information'
    else
      write(message(1), '(a,i5)') 'Info: Loading restart information for linear response'
    end if
    call write_info(1)

    ! sanity check
    gs_allocated = (associated(st%dpsi) .and. st%wfs_type == M_REAL) .or. &
      (associated(st%zpsi) .and. st%wfs_type == M_CMPLX)
    ASSERT(gs_allocated)

    if(present(lr)) then 
      lr_allocated = (associated(lr%ddl_psi) .and. st%wfs_type == M_REAL) .or. &
        (associated(lr%zdl_psi) .and. st%wfs_type == M_CMPLX)
      ASSERT(lr_allocated)
    endif

    ierr = 0

    ! open files to read
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp = .true., grp = gr%mesh%mpi_grp)
    if(iunit < 0) then
      ierr = -1
    end if

    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp = .true., grp = gr%mesh%mpi_grp)
    if(iunit2 < 0) then
      if(iunit > 0) call io_close(iunit, grp = gr%mesh%mpi_grp)
      ierr = -1
    end if

    if(ierr.ne.0) then
      write(message(1),'(a)') 'Could not load any previous restart information.'
      call write_info(1)
      call messages_print_stress(stdout)
      call pop_sub()
      return
    end if

    ! Reads out the previous mesh info.
    call read_previous_mesh()

    ! now we really start
    ALLOCATE(filled(st%d%dim, st%st_start:st%st_end, st%d%nik), st%d%dim*st%lnst*st%d%nik)
    filled = .false.

    ! Skip two lines.
    call iopar_read(gr%mesh%mpi_grp, iunit,  line, err); call iopar_read(gr%mesh%mpi_grp, iunit,  line, err)
    call iopar_read(gr%mesh%mpi_grp, iunit2, line, err); call iopar_read(gr%mesh%mpi_grp, iunit2, line, err)

    do
      call iopar_read(gr%mesh%mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit

      call iopar_backspace(gr%mesh%mpi_grp, iunit)

      call iopar_read(gr%mesh%mpi_grp, iunit, line, err)
      read(line, *) ik, char, ist, char, idim, char, filename
      if(index_is_wrong()) then
        call iopar_read(gr%mesh%mpi_grp, iunit2, line, err)
        cycle
      end if

      call iopar_read(gr%mesh%mpi_grp, iunit2, line, err)
      read(line, *) st%occ(ist, ik), char, st%eigenval(ist, ik)

      if(ist >= st%st_start .and. ist <= st%st_end .and. &
         st%d%kpt%start <= ik .and. st%d%kpt%end >= ik) then
        if(.not.mesh_change) then
          if( .not. present(lr) ) then 
            if (st%wfs_type == M_REAL) then
              call drestart_read_function(dir, filename, gr%mesh, st%dpsi(:, idim, ist, ik), err)
            else
              call zrestart_read_function(dir, filename, gr%mesh, st%zpsi(:, idim, ist, ik), err)
            end if
          else
            if (st%wfs_type == M_REAL) then
              call drestart_read_function(dir, filename, gr%mesh, lr%ddl_psi(:, idim, ist, ik), err)
            else
              call zrestart_read_function(dir, filename, gr%mesh, lr%zdl_psi(:, idim, ist, ik), err)
            end if
          end if
          if(err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          end if
        else
          if (st%wfs_type == M_REAL) then
            call drestart_read_function(dir, filename, old_mesh, dphi, err)
            if( .not. present(lr) ) then 
              call dmf_interpolate(old_mesh, gr%mesh, full_interpolation, dphi, st%dpsi(:, idim, ist, ik))
            else
              call dmf_interpolate(old_mesh, gr%mesh, full_interpolation, dphi, lr%ddl_psi(:, idim, ist, ik))
            end if
          else
            call zrestart_read_function(dir, filename, old_mesh, zphi, err)
            if( .not. present(lr) ) then 
              call zmf_interpolate(old_mesh, gr%mesh, full_interpolation, zphi, st%zpsi(:, idim, ist, ik))
            else
              call zmf_interpolate(old_mesh, gr%mesh, full_interpolation, zphi, lr%zdl_psi(:, idim, ist, ik))
            end if
          end if
          if(err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          end if
        endif
      end if
    end do

    if(present(iter)) then
      call iopar_read(gr%mesh%mpi_grp, iunit, line, err)
      read(line, *) filename, filename, iter
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
      ierr = err
    end if
    if(st%d%kpt%parallel) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif

    call fill()
    if(ierr == 0) then
      ierr = -1 ! no files read
      if(.not. present(lr)) then 
        write(str, '(a,i5)') 'Loading restart information'
      else
        write(str, '(a,i5)') 'Loading restart information for linear response'
      end if
      call messages_print_stress(stdout, trim(str))
      write(message(1),'(a)') 'No files could be read. No restart information can be used.'
      call write_info(1)
      call messages_print_stress(stdout)
    else
      ! Everything o.k.
      if(ierr == st%nst*st%d%nik*st%d%dim) then
        ierr = 0
      else
        if(.not. present(lr)) then 
          write(str, '(a,i5)') 'Loading restart information'
        else
          write(str, '(a,i5)') 'Loading restart information for linear response'
        end if
        call messages_print_stress(stdout, trim(str))
        write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' files out of ', &
          st%nst*st%d%nik*st%d%dim, ' could be read.'
        call write_info(1)
        call messages_print_stress(stdout)
      endif
    end if

    SAFE_DEALLOCATE_A(filled)
    call io_close(iunit, grp = gr%mesh%mpi_grp)
    call io_close(iunit2, grp = gr%mesh%mpi_grp)

    if(mesh_change) call interpolation_end()

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine interpolation_init
      call mesh_init_stage_1(old_mesh, old_sb, old_cv, gr%der%n_ghost)
      call mesh_init_stage_2(old_mesh, old_sb, geo, old_cv, gr%stencil)
      call mesh_init_stage_3(old_mesh)

      if (states_are_real(st)) then
        ALLOCATE(dphi(old_mesh%np_global), old_mesh%np_global)
      else
        ALLOCATE(zphi(old_mesh%np_global), old_mesh%np_global)
      end if
    end subroutine interpolation_init


    ! ---------------------------------------------------------
    subroutine interpolation_end
      call mesh_end(old_mesh)
      if (states_are_real(st)) then
        SAFE_DEALLOCATE_A(dphi)
      else
        SAFE_DEALLOCATE_A(zphi)
      end if
    end subroutine interpolation_end


    ! ---------------------------------------------------------
    subroutine read_previous_mesh

      mesh_change = .false.
      full_interpolation = .true.

#if defined(HAVE_MPI)
      return ! For the moment, avoid the complications of parallel stuff.
#endif
      if(present(iter)) return ! No interpolation, in case we are in the td part.

      iunit_mesh  = io_open(trim(dir)//'/mesh', action='read', status='old', die=.false., grp = gr%mesh%mpi_grp)
      if(iunit_mesh < 0) return

      read(iunit_mesh, *); read(iunit_mesh, *); read(iunit_mesh, *)
      call curvilinear_init_from_file(old_cv, iunit_mesh)
      call simul_box_init_from_file(old_sb, iunit_mesh)
      call io_close(iunit_mesh, grp = gr%mesh%mpi_grp)

      if( .not. curvilinear_is_eq(old_cv, gr%cv) ) then
        mesh_change = .true.
        return
      end if
      if( .not. simul_box_is_eq(old_sb, gr%sb) ) then
        mesh_change = .true.
        ! First, check whether the spacings are the same.
        if(old_sb%h .app. gr%sb%h) full_interpolation = .false.
      end if

      if(mesh_change) then
        if(.not.full_interpolation) then
          write(message(1),'(a)') 'The functions stored in "restart/gs" were obtained with' 
          write(message(2),'(a)') 'a different simulation box. The possible missing regions will be'
          write(message(3),'(a)') 'padded with zeros.'
          call write_info(3)
        else
          write(message(1),'(a)') 'The functions stored in "restart/gs" were obtained with' 
          write(message(2),'(a)') 'a different mesh. The values of the functions for the current'
          write(message(3),'(a)') 'calculations will be interpolated/extrapolated.'
          call write_info(3)
        end if
        call interpolation_init()
      endif
    
    end subroutine read_previous_mesh


    ! ---------------------------------------------------------
    subroutine fill() ! Put random function in orbitals that could not be read.
      integer :: first

      do ik = st%d%kpt%start, st%d%kpt%end
        first = st%st_end + 1 ! we have to orthogonalize all states after first

        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            if(filled(idim, ist, ik)) then
              cycle
            else
              if(first > ist) first = ist
            end if

            call states_generate_random(st, gr%mesh, ist, ist)
            st%occ(ist, ik) = M_ZERO
          end do
        end do
      end do

      if(first.ne.st%st_end+1) call states_orthogonalize(st, gr%mesh, first)

    end subroutine fill


    ! ---------------------------------------------------------
    logical function index_is_wrong() ! .true. if the index (idim, ist, ik) is not present in st structure...
      if(idim > st%d%dim .or. idim < 1 .or.   &
           ist   > st%nst   .or. ist  < 1 .or.   &
           ik    > st%d%nik .or. ik   < 1) then
        index_is_wrong = .true.
      else
        index_is_wrong = .false.
      end if
    end function index_is_wrong

  end subroutine restart_read


  ! ---------------------------------------------------------
  ! When doing an open boundary calculation Octopus needs the
  ! unscattered states in order to solve the Lippmann-Schwinger
  ! equation to obtain extended eigenstates.
  ! Since we only need the occupied states we just read them.
  subroutine read_free_states(st, gr)
    type(states_t),       intent(inout) :: st
    type(grid_t), target, intent(in)    :: gr
    
    integer                    :: k, ik, ist, idim, jk, err, wfns, occs
    integer                    :: np, ip, lead_nr(2, MAX_DIM)
    character(len=256)         :: line, fname, filename, restart_dir, chars
    character                  :: char
    FLOAT                      :: occ, eval, k_x, k_y, k_z, w_k, w_sum
    type(simul_box_t), pointer :: sb
    type(mesh_t), pointer      :: m_lead, m_center
    CMPLX                      :: phase
    CMPLX, allocatable         :: tmp(:, :)
    type(mpi_grp_t)            :: mpi_grp

    call push_sub('restart.states_read_free_states')

    sb       => gr%sb
    m_lead   => gr%mesh%lead_unit_cell(LEFT)
    m_center => gr%mesh
    lead_nr(1, :) = m_lead%idx%nr(1, :)+m_lead%idx%enlarge
    lead_nr(1, :) = m_lead%idx%nr(2, :)-m_lead%idx%enlarge

    mpi_grp = mpi_world

    np = m_lead%np
    ALLOCATE(tmp(np, st%d%dim), np*st%d%dim)
    restart_dir = trim(sb%lead_restart_dir(LEFT))//'/gs'

    wfns = io_open(trim(restart_dir)//'/wfns', action='read', is_tmp=.true., grp=mpi_grp)
    if(wfns.lt.0) then
      message(1) = 'Could not read '//trim(restart_dir)//'/wfns.'
      call write_fatal(1)
    end if
    occs = io_open(trim(restart_dir)//'/occs', action='read', is_tmp=.true., grp=mpi_grp)
    if(occs.lt.0) then
      message(1) = 'Could not read '//trim(restart_dir)//'/occs.'
      call write_fatal(1)
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, wfns, line, err); call iopar_read(mpi_grp, wfns, line, err)
    call iopar_read(mpi_grp, occs, line, err); call iopar_read(mpi_grp, occs, line, err)

    jk  = 0 ! reset counter for k-points
    st%d%kpoints(:, :) = M_ZERO
    st%d%kweights(:) = M_ZERO

    st%ob_rho = M_ZERO
    call mpi_grp_copy(m_lead%mpi_grp, gr%mesh%mpi_grp)
    do
      ! Check for end of file. Check only one of the two files assuming
      ! they are written correctly, i. e. of same length.
      call iopar_read(mpi_grp, wfns, line, err)
      read(line, '(a)') char
      if(char.eq.'%') then
        exit
      end if

      call iopar_backspace(mpi_grp, wfns)

      call iopar_read(mpi_grp, wfns, line, err)
      read(line, *) ik, char, ist, char, idim, char, fname

      call iopar_read(mpi_grp, occs, line, err)
      read(line, *) occ, char, eval, char, k_x, char, k_y, char, k_z, char, w_k, char, &
        chars, char, ik, char, ist, char, idim
      ! FIXME for more than 1 state
      if(occ.gt.CNST(1e-5)) then
        ! count the occupied k-points (with idim==1)
        if(idim.eq.1) jk = jk + 1

        st%d%kpoints(:, jk) = (/k_x, k_y, k_z/)
        st%d%kweights(jk) = w_k
        st%occ(ist, jk) = occ
        ! if not in the corresponding node cycle
        if(jk < st%d%kpt%start .or. jk > st%d%kpt%end) cycle
      else ! we do not consider empty states
        cycle
      end if

      call zrestart_read_function(trim(sb%lead_restart_dir(LEFT))//'/gs', fname, m_lead, tmp(:, idim), err)

      call lead_dens_accum()


      ! This loop replicates the lead wavefunction into the central region.
      ! It only works in this compact form because the transport-direction (x) is
      ! the index running slowest.          
      do k = 1, gr%sb%n_ucells
        st%zphi((k-1)*np+1:k*np, idim, 1, jk) = tmp(1:np, idim)
      end do

      ! Apply phase.
      do ip = 1, gr%mesh%np
        phase = exp(-M_zI*sum(gr%mesh%x(ip, 1:MAX_DIM)*st%d%kpoints(1:MAX_DIM, jk)))
        st%zphi(ip, idim, 1, jk) = phase*st%zphi(ip, idim, 1, jk)
      end do

      ! For debugging: write phi in gnuplot format to files.
      ! Only the z=0 plane is written, so mainly useful for 1D and 2D
      ! debugging.
      if(in_debug_mode) then
        write(filename, '(a,i3.3,a,i4.4,a,i1.1)') 'phi-', jk, '-', ist, '-', idim
        select case(calc_dim)
        case(1)
          call zoutput_function(output_axis_x, 'debug/open_boundaries', filename, &
            m_center, sb, st%zphi(:, idim, ist, jk), M_ONE, err, is_tmp=.false.)
        case(2, 3)
          call zoutput_function(output_plane_z, 'debug/open_boundaries', filename, &
            m_center, sb, st%zphi(:, idim, ist, jk), M_ONE, err, is_tmp=.false.)
        end select
      end if

    end do ! Loop over all free states.

    ! renormalize weigths
    w_sum = sum(st%d%kweights(:))
    st%d%kweights(:) = st%d%kweights(:)/w_sum
    st%qtot = M_ZERO
    do ist = 1, st%nst
      st%qtot = st%qtot + sum(st%occ(ist, 1:st%d%nik) * st%d%kweights(1:st%d%nik))
    end do

    call io_close(wfns); call io_close(occs)
    SAFE_DEALLOCATE_A(tmp)

    message(1) = "Info: Sucessfully initialized free states from '"// &
      trim(sb%lead_restart_dir(LEFT))//"/gs'"
    write(message(2),'(a,i3,a)') 'Info:', jk, ' occupied states read by program.'
    call write_info(2)

    call pop_sub()

  contains

    subroutine lead_dens_accum()
      integer :: il

      call push_sub('restart.lead_dens_accum')
      !FIXME no spinors yet
      do il = 1, NLEADS
        do ip = 1, np
          st%ob_rho(ip, idim, il) = st%ob_rho(ip, idim, il) + w_k*occ* &
            (real(tmp(ip, idim), REAL_PRECISION)**2 + aimag(tmp(ip, idim))**2)
        end do
        select case(st%d%ispin)
        case(SPINORS)
          message(1) = "restart.lead_dens_accum() does not work with spinors yet!"
          call write_fatal(1)
!          if(k_idim.eq.2) then
!            do ip = 1, np
!              c = w_k*occ*tmp(ip, 1)*conjg(tmp(ip, 2))
!              st%ob_rho(ip, 3, il) = st%ob_rho(ip, 3, il) + real(c, REAL_PRECISION)
!              st%ob_rho(ip, 4, il) = st%ob_rho(ip, 4, il) + aimag(c)
!            end do
!          end if
        end select
      end do

      call pop_sub()
    end subroutine lead_dens_accum
  end subroutine read_free_states

  ! ---------------------------------------------------------
  ! the routine reads formulas for user-defined wavefunctions 
  ! from the input file and fills the respective orbitals
  subroutine restart_read_user_def_orbitals(mesh, st)
    type(mesh_t),      intent(in) :: mesh
    type(states_t), intent(inout) :: st    

    type(block_t) :: blk
    integer   :: ip, id, is, ik, nstates, state_from, ierr, ncols
    integer   :: ib, idim, inst, inik, normalize
    FLOAT     :: x(MAX_DIM), r, psi_re, psi_im
    character(len=150) :: filename

    integer, parameter ::      &
      state_from_formula  = 1, &
      state_from_file     = 0, &
      normalize_yes       = 1, &
      normalize_no        = 0

    call push_sub('restart.restart_read_user_def_orbitals')

    !%Variable UserDefinedStates
    !%Type block
    !%Section States
    !%Description
    !% Instead of using the ground state as initial state for
    !% time propagations it might be interesting in some cases 
    !% to specify alternative states. Similarly to user-defined
    !% potentials, this block allows you to specify formulas for
    !% the orbitals at t=0.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | formula | "exp(-r^2)*exp(-i*0.2*x)" | normalize_yes
    !% <br>%</tt>
    !%
    !% The first column specifies the component of the spinor, 
    !% the second column the number of the state and the third 
    !% contains kpoint and spin quantum numbers. Column four
    !% indicates that column five should be interpreted as a formula
    !% for the corresponding orbital.
    !% 
    !% Alternatively, if column four states file the state will
    !% be read from the file given in column five.
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | file | "/path/to/file" | normalize_no
    !% <br>%</tt>
    !% 
    !% Octopus reads first the ground-state orbitals from
    !% the restart/gs directory. Only the states that are
    !% specified in the above block will be overwritten with
    !% the given analytical expression for the orbital.
    !%
    !% The sixth (optional) column indicates if Octopus should renormalize the orbital
    !% or not. The default, i. e. no sixth column given, is to renormalize.
    !%
    !%Option file 0
    !% Read initial orbital from file
    !%Option formula 1
    !% Calculate initial orbital by given analytic expression
    !%Option normalize_yes 1
    !% Normalize orbitals (default)
    !%Option normalize_no 0
    !% Do not normalize orbitals
    !%End
    if(loct_parse_block(datasets_check('UserDefinedStates'), blk) == 0) then

      call messages_print_stress(stdout, trim('Substitution of orbitals'))

      ! find out how many lines (i.e. states) the block has
      nstates = loct_parse_block_n(blk)

      ! read all lines
      do ib = 1, nstates
        ! Check that number of columns is five or six.
        ncols = loct_parse_block_cols(blk, ib-1)
        if(ncols.lt.5.or.ncols.gt.6) then
          message(1) = 'Each line in the UserDefinedStates block must have'
          message(2) = 'five or six columns.'
          call write_fatal(2)
        end if

        call loct_parse_block_int(blk, ib-1, 0, idim)
        call loct_parse_block_int(blk, ib-1, 1, inst)
        call loct_parse_block_int(blk, ib-1, 2, inik)

        ! Calculate from expression or read from file?
        call loct_parse_block_int(blk, ib-1, 3, state_from)

        ! loop over all states
        do id = 1, st%d%dim
          do is = 1, st%nst
            do ik = 1, st%d%nik

              ! does the block entry match and is this node responsible?
              if(.not.(id.eq.idim .and. is.eq.inst .and. ik.eq.inik    &
                .and. st%st_start.le.is .and. st%st_end.ge.is          &
                .and. st%d%kpt%start.le.ik .and. st%d%kpt%end.ge.ik) ) cycle

              select case(state_from)

              case(state_from_formula)
                ! parse formula string
                call loct_parse_block_string(                            &
                  blk, ib-1, 4, st%user_def_states(id, is, ik))

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with the expression:'
                write(message(3), '(2a)') '  ',trim(st%user_def_states(id, is, ik))
                call write_info(3)

                ! convert to C string
                call conv_to_C_string(st%user_def_states(id, is, ik))

                ! fill states with user defined formulas
                do ip = 1, mesh%np
                  x = mesh%x(ip, :)
                  r = sqrt(sum(x(:)**2))

                  ! parse user defined expressions
                  call loct_parse_expression(psi_re, psi_im, mesh%sb%dim, x, r, M_ZERO, st%user_def_states(id, is, ik))
                  ! fill state
                  st%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
                end do

              case(state_from_file)
                ! The input format can be coded in column four now. As it is
                ! not used now, we just say "file".
                ! Read the filename.
                call loct_parse_block_string(blk, ib-1, 4, filename)

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with data from file:'
                write(message(3), '(2a)') '  ',trim(filename)
                call write_info(3)

                ! finally read the state
                call zinput_function(filename, mesh, st%zpsi(:, id, is, ik), ierr, .true.)

              case default
                message(1) = 'Wrong entry in UserDefinedStates, column 4.'
                message(2) = 'You may state "formula" or "file" here.'
                call write_fatal(2)
              end select

              ! normalize orbital
              if(loct_parse_block_cols(blk, ib-1).eq.6) then
                call loct_parse_block_int(blk, ib-1, 5, normalize)
              else
                normalize = 1
              end if
              select case(normalize)
              case(normalize_no)
              case(normalize_yes)
                call zstates_normalize_orbital(mesh, st%d%dim, st%zpsi(:,:, is, ik))
              case default
                message(1) = 'The sixth column in UserDefinedStates may either be'
                message(2) = '"normalize_yes" or "normalize_no"'
                call write_fatal(2)
              end select
            end do
          end do
        end do

      end do

      call loct_parse_block_end(blk)
      call messages_print_stress(stdout)

    else
      message(1) = '"UserDefinesStates" has to be specified as block.'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine restart_read_user_def_orbitals

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
