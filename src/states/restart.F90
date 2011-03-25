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
  use batch_m
  use curvilinear_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use io_binary_m
  use io_function_m
  use kpoints_m
  use lalg_basic_m
  use linear_response_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_batch_m
  use mesh_function_m
  use mesh_init_m
  use messages_m
  use mpi_m
  use ob_green_m
  use parser_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use string_m
  use unit_m
  use unit_system_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::                  &
    restart_init,            &
    clean_stop,              &
    read_free_states,        &
    restart_write,           &
    restart_read,            &
    restart_get_ob_intf,     &
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
  
  logical           :: restart_write_files
  integer           :: restart_format
  character(len=32) :: restart_dir

  type(profile_t), save :: prof_read, prof_write

contains

  ! returns true if a file named stop exists
  function clean_stop()
    logical clean_stop, file_exists

    PUSH_SUB(clean_stop)

    clean_stop = .false.
    inquire(file='stop', exist=file_exists)
    if(file_exists) then
      message(1) = 'Clean STOP'
      call messages_warning(1)
      clean_stop = .true.

      if(mpi_grp_is_root(mpi_world)) call loct_rm('stop')
    end if

    POP_SUB(clean_stop)
  end function clean_stop


  ! ---------------------------------------------------------
  ! read restart format information
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
      message(1) = 'Restart information will not be written'
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
    if(scan(restart_dir, '/', .true.).lt.1) then
      restart_dir = trim(restart_dir)//'/'
    end if

    POP_SUB(restart_init)
  end subroutine restart_init


  ! ---------------------------------------------------------
  subroutine restart_look_and_read(st, gr, geo, is_complex, specify_dir, exact)
    type(states_t),             intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    type(geometry_t),           intent(in)    :: geo
    logical, optional,          intent(in)    :: is_complex
    character(len=*), optional, intent(in)    :: specify_dir
    logical,          optional, intent(in)    :: exact ! if .true. we need all the wavefunctions and on the exact grid

    integer :: kpoints, dim, nst, ierr
    character(len=80) dir

    PUSH_SUB(restart_look_and_read)

    if(present(specify_dir)) then
      dir = specify_dir
    else
      dir = restart_dir
    endif

    !check how many wfs we have
    call states_look(trim(dir)//GS_DIR, gr%mesh%mpi_grp, kpoints, dim, nst, ierr)
    if(ierr .ne. 0) then
      message(1) = 'Could not properly read wavefunctions from "'//trim(dir)//GS_DIR//'".'
      call messages_fatal(1)
    end if
    
    ! FIXME: This wrong, one cannot just change the number of states
    ! without updating the internal structures.
    st%nst    = nst
    st%st_end = nst
    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_DEALLOCATE_P(st%occ)

    if (present(is_complex)) then
      if ( is_complex ) then 
        call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)
      else 
        call states_allocate_wfns(st, gr%mesh, TYPE_FLOAT)
      end if
    else
      ! allow states_allocate_wfns to decide for itself whether complex or real needed
      call states_allocate_wfns(st, gr%mesh)
    endif

    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(     st%occ(1:st%nst, 1:st%d%nik))

    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
      st%spin = M_ZERO
    end if
    st%eigenval = M_HUGE
    st%occ      = M_ZERO

    ! load wavefunctions
    call restart_read(trim(dir)//GS_DIR, st, gr, geo, ierr, exact)

    POP_SUB(restart_look_and_read)
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

    integer :: iunit, iunit2, iunit_mesh, iunit_states, err, ik, idir, ist, idim, itot
    character(len=80) :: filename
    logical :: wfns_are_associated, lr_wfns_are_associated
    FLOAT   :: kpoint(1:MAX_DIM)

    PUSH_SUB(restart_write)

    if(.not. restart_write_files) then
      ierr = 0
    POP_SUB(restart_write)
      return
    end if

    call profiling_in(prof_write, "RESTART_WRITE")

    wfns_are_associated = (associated(st%dpsi) .and. states_are_real(st)) .or. &
      (associated(st%zpsi) .and. states_are_complex(st))
    ASSERT(wfns_are_associated)

    if(present(lr)) then 
      lr_wfns_are_associated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_wfns_are_associated)
    endif

    ierr = 0

    call block_signals()

    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(dir, is_tmp=.true.)
    end if

#ifdef HAVE_MPI
    !we need a barrier to wait for the directory to be created
    call MPI_Barrier(st%dom_st_kpt_mpi_grp%comm, mpi_err)
#endif

    if(mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
      iunit = io_open(trim(dir)//'/wfns', action='write', is_tmp=.true.)
      write(iunit,'(a)') '#     #k-point            #st            #dim    filename'
      write(iunit,'(a)') '%Wavefunctions'

      iunit2 = io_open(trim(dir)//'/occs', action='write', is_tmp=.true.)
      write(iunit2,'(a)') '# occupations | eigenvalue[a.u.] | k-points | k-weights | filename | ik | ist | idim'
      write(iunit2,'(a)') '%Occupations_Eigenvalues_K-Points'

      iunit_mesh = io_open(trim(dir)//'/mesh', action='write', is_tmp=.true.)
      write(iunit_mesh,'(a)') '# This file contains the necessary information to generate the'
      write(iunit_mesh,'(a)') '# mesh with which the functions in this directory were calculated,'
      write(iunit_mesh,'(a)') '# except for the geometry of the system.'

      call simul_box_dump(gr%sb, iunit_mesh)
      call mesh_dump(gr%mesh, iunit_mesh)
      call io_close(iunit_mesh)

      call mesh_write_fingerprint(gr%mesh, trim(dir)//'/grid')

      ! write the lxyz array
      if(gr%mesh%sb%box_shape /= HYPERCUBE) then
        ASSERT(associated(gr%mesh%idx%lxyz))
        call io_binary_write(trim(dir)//'/lxyz.obf', gr%mesh%np_part_global*gr%sb%dim, gr%mesh%idx%lxyz, ierr)
      end if

      iunit_states = io_open(trim(dir)//'/states', action='write', is_tmp=.true.)
      call states_dump(st, iunit_states)
      call io_close(iunit_states)
    end if

    itot = 1
    do ik = 1, st%d%nik
      kpoint = M_ZERO
      kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

      do ist = 1, st%nst
        do idim = 1, st%d%dim
          write(filename,'(i10.10)') itot

          if(mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
            write(iunit, '(i8,a,i8,a,i8,3a)') ik, ' | ', ist, ' | ', idim, ' | "', trim(filename), '"'
            write(iunit2, '(e20.14,a,e20.14,a)', advance='no') st%occ(ist,ik), ' | ', st%eigenval(ist, ik), ' | '
            do idir = 1, gr%sb%dim
              write(iunit2, '(e21.14,a)', advance='no') kpoint(idir), ' | '
            enddo
            write(iunit2, '(e21.14,a,i10.10,3(a,i8))') st%d%kweights(ik), ' | ', itot, ' | ', ik, ' | ', ist, ' | ', idim
          end if

          if(st%st_start <= ist .and. ist <= st%st_end) then
            if( .not. present(lr) ) then 
              if(st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
                if (states_are_real(st)) then
                  call drestart_write_function(dir, filename, gr, st%dpsi(:, idim, ist, ik), err, size(st%dpsi,1))
                else
                  call zrestart_write_function(dir, filename, gr, st%zpsi(:, idim, ist, ik), err, size(st%zpsi,1))
                end if
                if(err == 0) ierr = ierr + 1
              end if
            else
              if (states_are_real(st)) then
                call drestart_write_function(dir, filename, gr, lr%ddl_psi(:, idim, ist, ik), err, size(st%dpsi,1))
              else
                call zrestart_write_function(dir, filename, gr, lr%zdl_psi(:, idim, ist, ik), err, size(st%zpsi,1))
              end if
              if(err == 0) ierr = ierr + 1
            end if
          end if

          itot = itot + 1
        end do
      end do
    end do

    ! do NOT use st%lnst here as it is not (st%st_end - st%st_start + 1)
    if(ierr == st%d%kpt%nlocal * (st%st_end - st%st_start + 1) * st%d%dim) ierr = 0 ! All OK

    if(mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
      write(iunit,'(a)') '%'
      if(present(iter)) write(iunit,'(a,i7)') 'Iter = ', iter
      write(iunit2, '(a)') '%'

      call io_close(iunit)
      call io_close(iunit2)
    end if

    call unblock_signals()

    call profiling_out(prof_write)
    POP_SUB(restart_write)
  end subroutine restart_write


  ! ---------------------------------------------------------
  ! Reads the interface regions of the wavefunctions
  subroutine restart_get_ob_intf(st, gr)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in) :: gr

    integer              :: ik, ist, idim, il, ip

    PUSH_SUB(restart_get_ob_intf)

    write(message(1), '(a,i5)') 'Info: Reading ground-state interface wavefunctions.'
    call messages_info(1)

    ! Sanity check.
    do il = 1, NLEADS
      ASSERT(associated(st%ob_lead(il)%intf_psi))
      ASSERT(il.le.2) ! FIXME: wrong if non-transport calculation
    end do

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          do il = 1, NLEADS
            forall(ip = 1:gr%intf(il)%np_uc)
              ! FIXME: this will probably fail for MeshBlockSize > 1
              st%ob_lead(il)%intf_psi(ip, idim, ist, ik) = st%zpsi(gr%intf(il)%index(ip), idim, ist, ik)
            end forall
          end do
        end do
      end do
    end do

    POP_SUB(restart_get_ob_intf)
  end subroutine restart_get_ob_intf


  ! ---------------------------------------------------------
  ! returns
  ! <0 => Fatal error
  ! =0 => read all wavefunctions
  ! >0 => could only read x wavefunctions
  subroutine restart_read(dir, st, gr, geo, ierr, read_occ, iter, lr, exact)
    character(len=*),     intent(in)    :: dir
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in)    :: gr
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(out)   :: ierr
    logical,    optional, intent(in)    :: read_occ !< should I read the occupations
    integer,    optional, intent(inout) :: iter
    type(lr_t), optional, intent(inout) :: lr       !< if present, the lr wfs are read instead of the gs wfs
    logical,    optional, intent(in)    :: exact    !< if .true. we need all the wavefunctions and on the same grid

    integer              :: wfns_file, occ_file, err, ik, ist, idir, idim, int 
    integer              :: read_np, read_np_part, read_ierr, ip, xx(1:MAX_DIM), iread, nread
    character(len=12)    :: filename
    character(len=1)     :: char
    logical, allocatable :: filled(:, :, :)
    character(len=256)   :: line
    character(len=50)    :: str
    integer, allocatable :: read_lxyz(:, :), map(:)

    FLOAT                :: my_occ, flt
    logical              :: read_occ_, gs_allocated, lr_allocated, grid_changed, grid_reordered
    logical              :: exact_
    FLOAT, pointer       :: dpsi(:, :, :, :)
    CMPLX, pointer       :: zpsi(:, :, :, :)
    character(len=256), allocatable :: restart_file(:, :, :)
    logical,            allocatable :: restart_file_present(:, :, :)

    PUSH_SUB(restart_read)

    call profiling_in(prof_read, "RESTART_READ")

    exact_ = optional_default(exact, .false.)

    if(.not. present(lr)) then 
      write(message(1), '(a,i5)') 'Info: Loading restart information.'
    else
      write(message(1), '(a,i5)') 'Info: Loading restart information for linear response.'
    end if
    call messages_info(1)

    ! If one restarts a GS calculation changing the %Occupations block, one
    ! cannot read the occupations, otherwise these overwrite the ones from
    ! the input file. The same is true when reading a linear-response restart.
    read_occ_ = optional_default(read_occ, .true.)

    ! sanity check
    gs_allocated = (associated(st%dpsi) .and. states_are_real(st)) .or. &
      (associated(st%zpsi) .and. states_are_complex(st))
    ASSERT(gs_allocated)

    if(present(lr)) then 
      lr_allocated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_allocated)
    endif

    ierr = 0

    ! open files to read
    wfns_file  = io_open(trim(dir)//'/wfns', action='read', &
      status='old', die=.false., is_tmp = .true., grp = st%dom_st_kpt_mpi_grp)
    if(wfns_file < 0) ierr = -1

    occ_file = io_open(trim(dir)//'/occs', action='read', &
      status='old', die=.false., is_tmp = .true., grp = st%dom_st_kpt_mpi_grp)
    if(occ_file < 0) ierr = -1

    ! now read the mesh information
    call mesh_read_fingerprint(gr%mesh, trim(dir)//'/grid', read_np_part, read_np)

    ! For the moment we continue reading if we receive -1 so we can
    ! read old restart files that do not have a fingerprint file.
    ! if (read_np < 0) ierr = -1

    if (read_np > 0 .and. gr%mesh%sb%box_shape /= HYPERCUBE) then

      grid_changed = .true.

      ! perhaps only the order of the points changed, this can only
      ! happen if the number of points is the same and no points maps
      ! to zero (this is checked below)
      grid_reordered = (read_np == gr%mesh%np_global)

      ! the grid is different, so we read the coordinates.
      SAFE_ALLOCATE(read_lxyz(1:read_np_part, 1:gr%mesh%sb%dim))
      ASSERT(associated(gr%mesh%idx%lxyz))
      call io_binary_read(trim(dir)//'/lxyz.obf', read_np_part*gr%mesh%sb%dim, read_lxyz, read_ierr)

      ! and generate the map
      SAFE_ALLOCATE(map(1:read_np))

      do ip = 1, read_np
        xx = 0
        xx(1:gr%mesh%sb%dim) = read_lxyz(ip, 1:gr%mesh%sb%dim)
        if(any(xx(1:gr%mesh%sb%dim) < gr%mesh%idx%nr(1, 1:gr%mesh%sb%dim)) .or. &
          any(xx(1:gr%mesh%sb%dim) > gr%mesh%idx%nr(2, 1:gr%mesh%sb%dim))) then
          map(ip) = 0
          grid_reordered = .false.
        else
          map(ip) = gr%mesh%idx%lxyz_inv(xx(1), xx(2), xx(3))
          if(map(ip) > gr%mesh%np_global) map(ip) = 0
        end if
      end do

      SAFE_DEALLOCATE_A(read_lxyz)

      if(grid_reordered) then
        message(1) = 'Octopus is attempting to restart from mesh with a different order of points.'
        call messages_warning(1)
      else
        message(1) = 'Octopus is attempting to restart from a different mesh.'
        call messages_warning(1)
        if(exact_) ierr = -1
      end if

    else
      grid_changed = .false.
      grid_reordered = .false.
    end if

    if(ierr .ne. 0) then
      if(wfns_file > 0) call io_close(wfns_file, grp = st%dom_st_kpt_mpi_grp)
      if(occ_file > 0) call io_close(occ_file, grp = st%dom_st_kpt_mpi_grp)
      if(exact_) call restart_fail()
      write(message(1),'(a)') 'Could not load restart information.'
      call messages_info(1)
      call messages_print_stress(stdout)
      call profiling_out(prof_read)
      POP_SUB(restart_read)
      return
    end if

    if (states_are_real(st)) then
      nullify(zpsi)
      if(.not. present(lr)) then
        dpsi => st%dpsi
      else
        dpsi => lr%ddl_psi
      end if
    else
      nullify(dpsi)
      if(.not. present(lr)) then
        zpsi => st%zpsi
      else
        zpsi => lr%zdl_psi
      end if
    end if

    SAFE_ALLOCATE(restart_file(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    SAFE_ALLOCATE(restart_file_present(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    restart_file_present = .true.

    ! Skip two lines.
    call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file,  line, err)
    call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file,  line, err)

    call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, line, err)
    call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, line, err)

    do
      call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, line, int)
      read(line, '(a)') char
      if(int .ne. 0 .or. char == '%') exit

      call iopar_backspace(st%dom_st_kpt_mpi_grp, wfns_file)

      call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, line, err)
      read(line, *) ik, char, ist, char, idim, char, filename
      if(index_is_wrong()) then
        call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, line, err)
        cycle
      end if

      call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, line, err)
      if(.not. present(lr)) then ! do not read eigenvalues when reading linear response
        ! # occupations | eigenvalue[a.u.] | k-points | k-weights | filename | ik | ist | idim
        read(line, *) my_occ, char, st%eigenval(ist, ik), char, (flt, char, idir = 1, gr%sb%dim), st%d%kweights(ik)
        if(read_occ_) st%occ(ist, ik) = my_occ
      end if

      if(ist >= st%st_start .and. ist <= st%st_end .and. &
        st%d%kpt%start <= ik .and. st%d%kpt%end >= ik) then

        restart_file(idim, ist, ik) = trim(filename)
        restart_file_present(idim, ist, ik) = .true.
      end if
    end do

    !now we really read, to avoid serialisation of reads in the
    !parallel case (io_par works as a barrier)

    SAFE_ALLOCATE(filled(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    filled = .false.

    if(mpi_grp_is_root(mpi_world)) then
      iread = 1
      nread = st%lnst*st%d%kpt%nlocal*st%d%dim
      call loct_progress_bar(-1, nread)
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim

          if(.not. restart_file_present(idim, ist, ik)) cycle

          if (states_are_real(st)) then
            if (.not. grid_changed) then
              call drestart_read_function(dir, restart_file(idim, ist, ik), gr%mesh, dpsi(:, idim, ist, ik), err)
            else
              call drestart_read_function(dir, restart_file(idim, ist, ik), gr%mesh, dpsi(:, idim, ist, ik), err, map)
            end if
          else
            if (.not. grid_changed) then
              call zrestart_read_function(dir, restart_file(idim, ist, ik), gr%mesh, zpsi(:, idim, ist, ik), err)
            else
              call zrestart_read_function(dir, restart_file(idim, ist, ik), gr%mesh, zpsi(:, idim, ist, ik), err, map)
            end if
          end if
          
          if(err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          end if
          
          if(mpi_grp_is_root(mpi_world)) then
            call loct_progress_bar(iread, nread)
            INCR(iread, 1)
          end if

        end do
      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then
      write(stdout, '(1x)')
    end if

    if(present(iter)) then
      call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, line, err)
      read(line, *) filename, filename, iter
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%st_kpt_mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif

    call fill_random()

    if(ierr == 0) then
      ierr = -1 ! no files read
      if(.not. present(lr)) then 
        write(str, '(a,i5)') 'Loading restart information'
      else
        write(str, '(a,i5)') 'Loading restart information for linear response'
      end if
      call messages_print_stress(stdout, trim(str))
      write(message(1),'(a)') 'No files could be read. No restart information can be used.'
      call messages_info(1)
      call messages_print_stress(stdout)
    else if(ierr == st%nst * st%d%nik * st%d%dim) then
      ierr = 0
      write(message(1), '(a)') 'Info: Restart loading done.'
      call messages_info(1)
    else
      if(.not. present(lr)) then 
        write(str, '(a,i5)') 'Loading restart information.'
      else
        write(str, '(a,i5)') 'Loading restart information for linear response.'
      end if
      call messages_print_stress(stdout, trim(str))
      write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' files out of ', &
        st%nst * st%d%nik * st%d%dim, ' could be read.'
      call messages_info(1)
      call messages_print_stress(stdout)

      if(ierr .ne. 0 .and. exact_) call restart_fail()
    end if

    SAFE_DEALLOCATE_A(filled)
    SAFE_DEALLOCATE_A(map)
    call io_close(wfns_file, grp = st%dom_st_kpt_mpi_grp)
    call io_close(occ_file, grp = st%dom_st_kpt_mpi_grp)

    call profiling_out(prof_read)
    POP_SUB(restart_read)

  contains

    ! ---------------------------------------------------------
    subroutine fill_random() ! Put random function in orbitals that could not be read.
      PUSH_SUB(restart_read.fill_random)

      do ik = st%d%kpt%start, st%d%kpt%end

        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            if(filled(idim, ist, ik)) cycle

            call states_generate_random(st, gr%mesh, ist, ist)
            if(read_occ_) st%occ(ist, ik) = M_ZERO
          end do
        end do
      end do

      POP_SUB(restart_read.fill_random)
    end subroutine fill_random

    ! ---------------------------------------------------------

    logical function index_is_wrong() ! .true. if the index (idim, ist, ik) is not present in st structure...
      PUSH_SUB(restart_read.index_is_wrong)

      if(idim > st%d%dim .or. idim < 1 .or.   &
        ist   > st%nst   .or. ist  < 1 .or.   &
        ik    > st%d%nik .or. ik   < 1) then
        index_is_wrong = .true.
      else
        index_is_wrong = .false.
      end if

      POP_SUB(restart_read.index_is_wrong)
    end function index_is_wrong

    ! ---------------------------------------------------------

    subroutine restart_fail()
      message(1) = "Could not read KS orbitals from '"//trim(dir)
      message(2) = "Please run a ground-state calculation first!"
      call messages_fatal(2)
    end subroutine restart_fail

  end subroutine restart_read


  ! ---------------------------------------------------------
  ! When doing an open-boundary calculation Octopus needs the
  ! unscattered states in order to solve the Lippmann-Schwinger
  ! equation to obtain extended eigenstates.
  subroutine read_free_states(st, gr)
    type(states_t),       intent(inout) :: st
    type(grid_t), target, intent(in)    :: gr
    
    integer                    :: ik, ist, idim, counter, err, wfns, occs, il
    integer                    :: np, ip, ip_lead, idir, ix, iy, iz, ixlead(1:3)
    character(len=256)         :: line, fname, filename, restart_dir, chars
    character                  :: char
    FLOAT                      :: occ, eval, kpoint(1:MAX_DIM), w_k
    type(simul_box_t), pointer :: sb
    type(mesh_t), pointer      :: m_lead, m_center
    CMPLX                      :: phase
    CMPLX, allocatable         :: tmp(:, :)
    type(mpi_grp_t)            :: mpi_grp
    integer                    :: start(1:3), end(1:3), start_lead(1:3), end_lead(1:3)

    PUSH_SUB(read_free_states)

    sb       => gr%sb
    m_lead   => gr%ob_grid%lead(LEFT)%mesh
    m_center => gr%mesh

    mpi_grp = st%dom_st_kpt_mpi_grp

    np = m_lead%np
    SAFE_ALLOCATE(tmp(1:np, 1:st%d%dim))
    restart_dir = trim(gr%ob_grid%lead(LEFT)%info%restart_dir)//'/'// GS_DIR

    wfns = io_open(trim(restart_dir)//'/wfns', action='read', is_tmp=.true., grp=mpi_grp)
    if(wfns .lt. 0) then
      message(1) = 'Could not read '//trim(restart_dir)//'/wfns.'
      call messages_fatal(1)
    end if
    occs = io_open(trim(restart_dir)//'/occs', action='read', is_tmp=.true., grp=mpi_grp)
    if(occs .lt. 0) then
      message(1) = 'Could not read '//trim(restart_dir)//'/occs.'
      call messages_fatal(1)
    end if

    ! Skip two lines.
    call iopar_read(mpi_grp, wfns, line, err)
    call iopar_read(mpi_grp, wfns, line, err)

    call iopar_read(mpi_grp, occs, line, err)
    call iopar_read(mpi_grp, occs, line, err)

    counter  = 0 ! reset counter
    st%d%kweights(:) = M_ZERO

    forall(il = 1:NLEADS) st%ob_lead(il)%rho = M_ZERO
    call mpi_grp_copy(m_lead%mpi_grp, gr%mesh%mpi_grp)
    do
      ! Check for end of file. Check only one of the two files assuming
      ! they are written correctly, i.e. of same length.
      call iopar_read(mpi_grp, wfns, line, err)
      read(line, '(a)') char
      if(char .eq. '%') then
        exit
      end if

      call iopar_backspace(mpi_grp, wfns)

      call iopar_read(mpi_grp, wfns, line, err)
      read(line, *) ik, char, ist, char, idim, char, fname

      call iopar_read(mpi_grp, occs, line, err)
      !# occupations | eigenvalue[a.u.] | k-points | k-weights | filename | ik | ist | idim
      read(line, *) occ, char, eval, char, (kpoint(idir), char, idir = 1, gr%sb%dim), &
        w_k, char, chars, char, ik, char, ist, char, idim

      st%d%kweights(ik) = w_k
      st%occ(ist, ik) = occ
      counter = counter + 1
      ! if not in the corresponding node cycle
      if(ik < st%d%kpt%start .or. ik > st%d%kpt%end .or. ist < st%st_start .or. ist > st%st_end) cycle

      call zrestart_read_function(trim(restart_dir), fname, m_lead, tmp(:, idim), err)

      call lead_dens_accum()

      ! Copy the wave-function from the lead to the central
      ! region, using periodicity.
      start(1:3) = gr%mesh%idx%nr(1, 1:3) + gr%mesh%idx%enlarge(1:3)
      end(1:3) = gr%mesh%idx%nr(2, 1:3) - gr%mesh%idx%enlarge(1:3)

      start_lead(1:3) = gr%ob_grid%lead(LEFT)%mesh%idx%nr(1, 1:3) + gr%ob_grid%lead(LEFT)%mesh%idx%enlarge(1:3)
      end_lead(1:3) = gr%ob_grid%lead(LEFT)%mesh%idx%nr(2, 1:3) - gr%ob_grid%lead(LEFT)%mesh%idx%enlarge(1:3)

      ixlead(1) = start_lead(1)
      do ix = start(1), end(1)

        ixlead(2) = start_lead(2)
        do iy = start(2), end(2)
          
          ixlead(3) = start_lead(3)
          do iz = start(3), end(3)

            ip = gr%mesh%idx%lxyz_inv(ix, iy, iz)
            ip_lead = gr%ob_grid%lead(LEFT)%mesh%idx%lxyz_inv(ixlead(1), ixlead(2), ixlead(3))

            if(ip > 0 .and. ip <= gr%mesh%np .and. ip_lead > 0 .and. ip_lead <= gr%ob_grid%lead(LEFT)%mesh%np) then
              st%zphi(ip, idim, ist, ik) = tmp(ip_lead, idim)
            end if

            INCR(ixlead(3), 1)
            if(TRANS_DIR == 3 .and. ixlead(3) > end_lead(3)) ixlead(3) = start_lead(3)
          end do

          INCR(ixlead(2), 1)
          if(TRANS_DIR == 2 .and. ixlead(2) > end_lead(2)) ixlead(2) = start_lead(2)
        end do

        INCR(ixlead(1), 1)
        if(TRANS_DIR == 1 .and. ixlead(1) > end_lead(1)) ixlead(1) = start_lead(1)
      end do

      ! Apply phase. Here the kpoint read from the file is different
      ! from the one calculated by octopus, since the box size changed.
      ! The read kpoints are use for describing the free states of the
      ! open system, which are mapped 1 to 1 (periodic incoming to scattered)
      ! with the Lippmann-Schwinger equation.

      do ip = 1, gr%mesh%np
        phase = exp(-M_zI * sum(gr%mesh%x(ip, 1:gr%sb%dim)*kpoint(1:gr%sb%dim)))
        st%zphi(ip, idim, ist, ik) = phase * st%zphi(ip, idim, ist, ik)
      end do

      ! For debugging: write phi in gnuplot format to files.
      ! Only the z=0 plane is written, so mainly useful for 1D and 2D
      ! debugging.
      if(in_debug_mode) then
        write(filename, '(a,i3.3,a,i4.4,a,i1.1)') 'phi-', ik, '-', ist, '-', idim
        select case(gr%sb%dim)
        case(1)
          call zoutput_function(output_axis_x, 'debug/open_boundaries', filename, &
            m_center, st%zphi(:, idim, ist, ik), sqrt(units_out%length**(-sb%dim)), err, is_tmp=.false.)
        case(2, 3)
          call zoutput_function(output_plane_z, 'debug/open_boundaries', filename, &
            m_center, st%zphi(:, idim, ist, ik), sqrt(units_out%length**(-sb%dim)), err, is_tmp=.false.)
        end select
      end if

    end do ! Loop over all free states.

    call io_close(wfns)
    call io_close(occs)
    SAFE_DEALLOCATE_A(tmp)

    message(1) = "Info: Successfully initialized free states from '"//trim(restart_dir)//"'"
    write(message(2),'(a,i3,a)') 'Info:', counter, ' wave functions read by program.'
    call messages_info(2)

    POP_SUB(read_free_states)

  contains

    subroutine lead_dens_accum()
      integer :: il

      PUSH_SUB(read_free_states.lead_dens_accum)
      !FIXME no spinors yet
      do il = 1, NLEADS
        do ip = 1, np
          st%ob_lead(il)%rho(ip, idim) = st%ob_lead(il)%rho(ip, idim) + w_k * occ * &
            (real(tmp(ip, idim), REAL_PRECISION)**2 + aimag(tmp(ip, idim))**2)
       end do
        select case(st%d%ispin)
        case(SPINORS)
          message(1) = "restart.lead_dens_accum() does not work with spinors yet!"
          call messages_fatal(1)
!          if(k_idim.eq.2) then
!            do ip = 1, np
!              c = w_k*occ*tmp(ip, 1)*conjg(tmp(ip, 2))
!              st%ob_lead(il)%rho(ip, 3) = st%ob_lead(il)%rho(ip, 3) + real(c, REAL_PRECISION)
!              st%ob_lead(il)%rho(ip, 4) = st%ob_lead(il)%rho(ip, 4l) + aimag(c)
!            end do
!          end if
        end select
      end do

      POP_SUB(read_free_states.lead_dens_accum)
    end subroutine lead_dens_accum
  end subroutine read_free_states

  ! ---------------------------------------------------------
  ! the routine reads formulas for user-defined wavefunctions 
  ! from the input file and fills the respective orbitals
  subroutine restart_read_user_def_orbitals(mesh, st)
    type(mesh_t),      intent(in) :: mesh
    type(states_t), intent(inout) :: st    

    type(block_t) :: blk
    integer :: ip, id, is, ik, nstates, state_from, ierr, ncols
    integer :: ib, idim, inst, inik, normalize
    FLOAT :: xx(MAX_DIM), rr, psi_re, psi_im
    character(len=150) :: filename

    integer, parameter ::      &
      state_from_formula  = 1, &
      state_from_file     = 0, &
      normalize_yes       = 1, &
      normalize_no        = 0

    PUSH_SUB(restart_read_user_def_orbitals)

    !%Variable UserDefinedStates
    !%Type block
    !%Section States
    !%Description
    !% Instead of using the ground state as initial state for
    !% time-propagations it might be interesting in some cases 
    !% to specify alternate states. Like with user-defined
    !% potentials, this block allows you to specify formulas for
    !% the orbitals at <i>t</i>=0.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | formula | "exp(-r^2)*exp(-i*0.2*x)" | normalize_yes
    !% <br>%</tt>
    !%
    !% The first column specifies the component of the spinor, 
    !% the second column the number of the state and the third 
    !% contains <i>k</i>-point and spin quantum numbers. Column four
    !% indicates that column five should be interpreted as a formula
    !% for the corresponding orbital.
    !% 
    !% Alternatively, if column four states <tt>file</tt> the state will
    !% be read from the file given in column five.
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 | file | "/path/to/file" | normalize_no
    !% <br>%</tt>
    !% 
    !% Octopus reads first the ground-state orbitals from
    !% the <tt>restart/gs</tt> directory. Only the states that are
    !% specified in the above block will be overwritten with
    !% the given analytic expression for the orbital.
    !%
    !% The sixth (optional) column indicates whether <tt>Octopus</tt> should renormalize
    !% the orbital. The default (no sixth column given) is to renormalize.
    !%
    !%Option file 0
    !% Read initial orbital from file.
    !%Option formula 1
    !% Calculate initial orbital by given analytic expression.
    !%Option normalize_yes 1
    !% Normalize orbitals (default).
    !%Option normalize_no 0
    !% Do not normalize orbitals.
    !%End
    if(parse_block(datasets_check('UserDefinedStates'), blk) == 0) then

      call messages_print_stress(stdout, trim('Substitution of orbitals'))

      ! find out how many lines (i.e. states) the block has
      nstates = parse_block_n(blk)

      ! read all lines
      do ib = 1, nstates
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, ib - 1)
        if(ncols .lt. 5 .or. ncols .gt. 6) then
          message(1) = 'Each line in the UserDefinedStates block must have'
          message(2) = 'five or six columns.'
          call messages_fatal(2)
        end if

        call parse_block_integer(blk, ib - 1, 0, idim)
        call parse_block_integer(blk, ib - 1, 1, inst)
        call parse_block_integer(blk, ib - 1, 2, inik)

        ! Calculate from expression or read from file?
        call parse_block_integer(blk, ib - 1, 3, state_from)

        ! loop over all states
        do id = 1, st%d%dim
          do is = 1, st%nst
            do ik = 1, st%d%nik

              ! does the block entry match and is this node responsible?
              if(.not.(id .eq. idim .and. is .eq. inst .and. ik .eq. inik    &
                .and. st%st_start .le. is .and. st%st_end .ge. is          &
                .and. st%d%kpt%start .le. ik .and. st%d%kpt%end .ge. ik) ) cycle

              select case(state_from)

              case(state_from_formula)
                ! parse formula string
                call parse_block_string(                            &
                  blk, ib - 1, 4, st%user_def_states(id, is, ik))

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with the expression:'
                write(message(3), '(2a)') '  ',trim(st%user_def_states(id, is, ik))
                call messages_info(3)

                ! convert to C string
                call conv_to_C_string(st%user_def_states(id, is, ik))

                ! fill states with user-defined formulas
                do ip = 1, mesh%np
                  xx = mesh%x(ip, :)
                  rr = sqrt(sum(xx(:)**2))

                  ! parse user-defined expressions
                  call parse_expression(psi_re, psi_im, mesh%sb%dim, xx, rr, M_ZERO, st%user_def_states(id, is, ik))
                  ! fill state
                  st%zpsi(ip, id, is, ik) = psi_re + M_zI * psi_im
                end do

              case(state_from_file)
                ! The input format can be coded in column four now. As it is
                ! not used now, we just say "file".
                ! Read the filename.
                call parse_block_string(blk, ib - 1, 4, filename)

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with data from file:'
                write(message(3), '(2a)') '  ',trim(filename)
                call messages_info(3)

                ! finally read the state
                call zinput_function(filename, mesh, st%zpsi(:, id, is, ik), ierr, .true.)
                if (ierr > 0) then
                  message(1) = 'Could not read the file!'
                  write(message(2),'(a,i1)') 'Error code: ', ierr
                  call messages_fatal(2)
                end if

              case default
                message(1) = 'Wrong entry in UserDefinedStates, column 4.'
                message(2) = 'You may state "formula" or "file" here.'
                call messages_fatal(2)
              end select

              ! normalize orbital
              if(parse_block_cols(blk, ib - 1) .eq. 6) then
                call parse_block_integer(blk, ib - 1, 5, normalize)
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
                call messages_fatal(2)
              end select
            end do
          end do
        end do

      end do

      call parse_block_end(blk)
      call messages_print_stress(stdout)

    else
      message(1) = '"UserDefinedStates" has to be specified as block.'
      call messages_fatal(1)
    end if

    POP_SUB(restart_read_user_def_orbitals)
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
