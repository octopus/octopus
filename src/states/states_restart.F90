!! Copyright (C) 2002-2016 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

module states_restart_oct_m
  use global_oct_m
  use grid_oct_m
  use io_function_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use linear_response_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_oct_m
  use states_dim_oct_m
  use string_oct_m
  use types_oct_m

  implicit none

  type(profile_t), save :: prof_read, prof_write

  private
  public ::                         &
    states_look_and_load,           &
    states_dump,                    &
    states_load,                    &
    states_dump_rho,                &
    states_load_rho,                &
    states_read_user_def_orbitals,  &
    states_dump_spin,               &
    states_load_spin

contains

  ! ---------------------------------------------------------
  subroutine states_look_and_load(restart, parser, st, gr, is_complex)
    type(restart_t),            intent(in)    :: restart
    type(parser_t),             intent(in)    :: parser
    type(states_t),     target, intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    logical,          optional, intent(in)    :: is_complex

    integer :: kpoints, dim, nst, ierr
    FLOAT, pointer :: new_occ(:,:)

    PUSH_SUB(states_look_and_load)

    !check how many wfs we have
    call states_look(restart, kpoints, dim, nst, ierr)
    if(ierr /= 0) then
      message(1) = "Unable to read states information."
      call messages_fatal(1)
    end if

    if(st%parallel_in_states) then
      message(1) = "Internal error: cannot use states_look_and_load when parallel in states."
      call messages_fatal(1)
    end if

    ! Resize st%occ, retaining information there
    SAFE_ALLOCATE(new_occ(1:nst, 1:st%d%nik))
    new_occ(:,:) = M_ZERO
    new_occ(1:min(nst, st%nst),:) = st%occ(1:min(nst, st%nst),:)
    SAFE_DEALLOCATE_P(st%occ)
    st%occ => new_occ

    ! FIXME: This wrong, one cannot just change the number of states
    ! without updating the internal structures, in the case of parallel in states.
    ! I think it is right now without state parallelism.
    st%nst      = nst
    st%st_start = 1
    st%st_end   = nst
    st%lnst     = nst

    SAFE_DEALLOCATE_P(st%node)
    SAFE_ALLOCATE(st%node(1:st%nst))
    st%node(:)  = 0

    SAFE_DEALLOCATE_P(st%eigenval)
    SAFE_ALLOCATE(st%eigenval(1:st%nst, 1:st%d%nik))
    st%eigenval = huge(st%eigenval)

    if (present(is_complex)) then
      if ( is_complex ) then
        call states_allocate_wfns(st, gr%mesh, TYPE_CMPLX)
      else
        call states_allocate_wfns(st, gr%mesh, TYPE_FLOAT)
      end if
    else
      ! allow states_allocate_wfns to decide for itself whether complex or real needed
      call states_allocate_wfns(st, gr%mesh)
    end if

    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
      st%spin = M_ZERO
    end if

    ! load wavefunctions
    call states_load(restart, parser, st, gr, ierr)
    if(ierr /= 0) then
      message(1) = "Unable to read wavefunctions."
      call messages_fatal(1)
    end if

    POP_SUB(states_look_and_load)
  end subroutine states_look_and_load


  ! ---------------------------------------------------------
  subroutine states_dump(restart, st, gr, ierr, iter, lr, st_start_writing, verbose)
    type(restart_t),      intent(in)  :: restart
    type(states_t),       intent(in)  :: st
    type(grid_t),         intent(in)  :: gr
    integer,              intent(out) :: ierr
    integer,    optional, intent(in)  :: iter
    !> if this next argument is present, the lr wfs are stored instead of the gs wfs
    type(lr_t), optional, intent(in)  :: lr
    integer,    optional, intent(in)  :: st_start_writing
    logical,    optional, intent(in)  :: verbose

    integer :: iunit_wfns, iunit_occs, iunit_states
    integer :: err, err2(2), ik, idir, ist, idim, itot
    integer :: root(1:P_STRATEGY_MAX)
    character(len=MAX_PATH_LEN) :: filename
    character(len=300) :: lines(3)
    logical :: lr_wfns_are_associated, should_write, verbose_
    FLOAT   :: kpoint(1:MAX_DIM)
    FLOAT,  allocatable :: dpsi(:), rff_global(:)
    CMPLX,  allocatable :: zpsi(:), zff_global(:)

    PUSH_SUB(states_dump)

    verbose_ = optional_default(verbose, .true.)

    ierr = 0

    if(restart_skip(restart)) then
      POP_SUB(states_dump)
      return
    end if

    if(verbose_) then
      message(1) = "Info: Writing states."
      call print_date(trim(message(1))//' ')
    end if

    call profiling_in(prof_write, "RESTART_WRITE")

    if (present(lr)) then
      lr_wfns_are_associated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_wfns_are_associated)
    end if

    call restart_block_signals()


    iunit_states = restart_open(restart, 'states')
    write(lines(1), '(a20,1i10)')  'nst=                ', st%nst
    write(lines(2), '(a20,1i10)')  'dim=                ', st%d%dim
    write(lines(3), '(a20,1i10)')  'nik=                ', st%d%nik
    call restart_write(restart, iunit_states, lines, 3, err)
    if (err /= 0) ierr = ierr + 1
    call restart_close(restart, iunit_states)


    iunit_wfns = restart_open(restart, 'wfns')
    lines(1) = '#     #k-point            #st            #dim    filename'
    if (states_are_real(st)) then
      lines(2) = '%Real_Wavefunctions'
    else
      lines(2) = '%Complex_Wavefunctions'
    end if
    call restart_write(restart, iunit_wfns, lines, 2, err)
    if (err /= 0) ierr = ierr + 2


    iunit_occs = restart_open(restart, 'occs')
    lines(1) = '# occupations | eigenvalue[a.u.] | Im(eigenvalue) [a.u.] | k-points | k-weights | filename | ik | ist | idim'
    lines(2) = '%Occupations_Eigenvalues_K-Points'
    call restart_write(restart, iunit_occs, lines, 2, err)
    if (err /= 0) ierr = ierr + 4


    if(states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np))
      SAFE_ALLOCATE(rff_global(1:gr%mesh%np_global))
    else
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
      SAFE_ALLOCATE(zff_global(1:gr%mesh%np_global))
    end if

    itot = 1
    root = -1
    err2 = 0
    do ik = 1, st%d%nik
      kpoint = M_ZERO
      kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik), absolute_coordinates = .true.)

      do ist = 1, st%nst
        do idim = 1, st%d%dim
          root(P_STRATEGY_DOMAINS) = mod(itot - 1, gr%mesh%mpi_grp%size)
          write(filename,'(i10.10)') itot
          
          write(lines(1), '(i8,a,i8,a,i8,3a)') ik, ' | ', ist, ' | ', idim, ' | "', trim(filename), '"'
          call restart_write(restart, iunit_wfns, lines, 1, err)
          if (err /= 0) err2(1) = err2(1) + 1

          write(lines(1), '(e21.14,a,e21.14)') st%occ(ist,ik), ' | ', st%eigenval(ist, ik)
          write(lines(1), '(a,a,e21.14)') trim(lines(1)), ' | ', CNST(0.0)
          do idir = 1, gr%sb%dim
            write(lines(1), '(a,a,e21.14)') trim(lines(1)), ' | ', kpoint(idir)
          end do
          write(lines(1), '(a,a,e21.14,a,i10.10,3(a,i8))') trim(lines(1)), &
            ' | ', st%d%kweights(ik), ' | ', itot, ' | ', ik, ' | ', ist, ' | ', idim
          call restart_write(restart, iunit_occs, lines, 1, err)
          if (err /= 0) err2(1) = err2(1) + 1

          should_write = st%st_start <= ist .and. ist <= st%st_end
          if (should_write .and. present(st_start_writing)) then
            if (ist < st_start_writing) should_write = .false.
          end if

          if (should_write) then
            if (.not. present(lr)) then
              if(st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
                if (states_are_real(st)) then
                  call states_get_state(st, gr%mesh, idim, ist, ik, dpsi)
                  call drestart_write_mesh_function(restart, filename, gr%mesh, dpsi, err, root = root)
                else
                  call states_get_state(st, gr%mesh, idim, ist, ik, zpsi)
                  call zrestart_write_mesh_function(restart, filename, gr%mesh, zpsi, err, root = root)
                end if
              else
                err = 0
              end if
            else
              if(st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
                if (states_are_real(st)) then
                  call drestart_write_mesh_function(restart, filename, gr%mesh, &
                    lr%ddl_psi(:, idim, ist, ik), err, root = root)
                else
                  call zrestart_write_mesh_function(restart, filename, gr%mesh, &
                    lr%zdl_psi(:, idim, ist, ik), err, root = root)
                end if
              else
                err = 0
              end if
            end if
            if (err /= 0) err2(2) = err2(2) + 1
          end if

          itot = itot + 1
        end do ! st%d%dim
      end do ! st%nst
    end do ! st%d%nik
    if (err2(1) /= 0) ierr = ierr + 8
    if (err2(2) /= 0) ierr = ierr + 16

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(rff_global)
    SAFE_DEALLOCATE_A(zff_global)

    lines(1) = '%'
    call restart_write(restart, iunit_occs, lines, 1, err)
    if (err /= 0) ierr = ierr + 32
    call restart_write(restart, iunit_wfns, lines, 1, err) 
    if (err /= 0) ierr = ierr + 64
    if (present(iter)) then
      write(lines(1),'(a,i7)') 'Iter = ', iter
      call restart_write(restart, iunit_wfns, lines, 1, err)
      if (err /= 0) ierr = ierr + 128
    end if

    call restart_close(restart, iunit_wfns)
    call restart_close(restart, iunit_occs)

    if(verbose_) then
      message(1) = "Info: Finished writing states."
      call print_date(trim(message(1))//' ')
    end if

    call restart_unblock_signals()

    call profiling_out(prof_write)
    POP_SUB(states_dump)
    return
  end subroutine states_dump


  ! ---------------------------------------------------------
  !> returns in ierr:
  !! <0 => Fatal error, or nothing read
  !! =0 => read all wavefunctions
  !! >0 => could only read ierr wavefunctions
  subroutine states_load(restart, parser, st, gr, ierr, iter, lr, lowest_missing, label, verbose)
    type(restart_t),            intent(in)    :: restart
    type(parser_t),             intent(in)    :: parser
    type(states_t),             intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    integer,                    intent(out)   :: ierr
    integer,          optional, intent(out)   :: iter
    type(lr_t),       optional, intent(inout) :: lr       !< if present, the lr wfs are read instead of the gs wfs
    integer,          optional, intent(out)   :: lowest_missing(:, :) !< all states below this one were read successfully
    character(len=*), optional, intent(in)    :: label
    logical,          optional, intent(in)    :: verbose

    integer              :: states_file, wfns_file, occ_file, err, ik, ist, idir, idim
    integer              :: idone, iread, ntodo
    character(len=12)    :: filename
    character(len=1)     :: char
    logical, allocatable :: filled(:, :, :)
    character(len=256)   :: lines(3), label_
    character(len=50)    :: str

    FLOAT                :: my_occ, imev, my_kweight
    logical              :: read_occ, lr_allocated, verbose_
    logical              :: integral_occs
    FLOAT, allocatable   :: dpsi(:)
    CMPLX, allocatable   :: zpsi(:), zpsiL(:)
    character(len=256), allocatable :: restart_file(:, :, :)
    logical,            allocatable :: restart_file_present(:, :, :)
    FLOAT                :: kpoint(MAX_DIM), read_kpoint(MAX_DIM)

#if defined(HAVE_MPI)
    integer              :: iread_tmp
    integer, allocatable :: lowest_missing_tmp(:, :)
#endif
    
    PUSH_SUB(states_load)

    ierr = 0

    ! make sure these intent(out)`s are initialized no matter what
    if (present(lowest_missing)) lowest_missing = 1
    if (present(iter)) iter = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load)
      return
    end if

    call profiling_in(prof_read, "RESTART_READ")

    verbose_ = optional_default(verbose, .true.)

    if (present(label)) then
      label_ = trim(label)
    else
      if (present(lr)) then
        label_ = " for linear response"
      else
        label_ = ""
      end if
    end if

    message(1) = 'Info: Reading states'
    if (len(trim(label_)) > 0) then
      message(1) = trim(message(1)) // trim(label_)
    end if
    message(1) = trim(message(1)) // "."
    if(verbose_) call print_date(trim(message(1))//' ')

    if(.not. present(lr)) then
      st%fromScratch = .false. ! obviously, we are using restart info
    end if

    ! If one restarts a GS calculation changing the %Occupations block, one
    ! cannot read the occupations, otherwise these overwrite the ones from
    ! the input file. restart_fixed_occ makes that we do use the ones in the file.
    integral_occs = .true. ! only used if restart_fixed_occ
    if (st%restart_fixed_occ) then
      read_occ = .true.
      st%fixed_occ = .true.
    else
      read_occ = .not. st%fixed_occ
    end if

    if(.not. present(lr)) then
      st%eigenval(:, :) = M_ZERO
      ! to be filled in from reading afterward
    end if

    if (.not. present(lr) .and. read_occ) then
      st%occ(:, :) = M_ZERO
      ! to be filled in from reading afterward
    end if

    ! sanity check
    if (present(lr)) then
      lr_allocated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_allocated)
    end if

    states_file  = restart_open(restart, 'states')
    ! sanity check on spin/k-points. Example file 'states':
    ! nst=                         2
    ! dim=                         1
    ! nik=                         2
    call restart_read(restart, states_file, lines, 3, err)
    if (err /= 0) then
      ierr = ierr - 2
    else
      read(lines(2), *) str, idim
      read(lines(3), *) str, ik
      if(idim == 2 .and. st%d%dim == 1) then
        write(message(1),'(a)') 'Incompatible restart information: saved calculation is spinors, this one is not.'
        call messages_warning(1)
        ierr = ierr - 2**2
      end if
      if(idim == 1 .and. st%d%dim == 2) then
        write(message(1),'(a)') 'Incompatible restart information: this calculation is spinors, saved one is not.'
        call messages_warning(1)
        ierr = ierr - 2**3
      end if
      if(ik < st%d%nik) then
        write(message(1),'(a)') 'Incompatible restart information: not enough k-points.'
        write(message(2),'(2(a,i6))') 'Expected ', st%d%nik, ' > Read ', ik
        call messages_warning(2)
      end if
      ! We will check that they are the right k-points later, so we do not need to put a specific error here.
    end if
    call restart_close(restart, states_file)


    ! open files to read
    wfns_file  = restart_open(restart, 'wfns')
    occ_file = restart_open(restart, 'occs')
    call restart_read(restart, wfns_file, lines, 2, err)
    if (err /= 0) then
      ierr = ierr - 2**5
    else if (states_are_real(st)) then
      read(lines(2), '(a)') str
      if (str(2:8) == 'Complex') then
        message(1) = "Cannot read real states from complex wavefunctions."
        call messages_warning(1)
        ierr = ierr - 2**6
      else if (str(2:5) /= 'Real') then 
        message(1) = "Restart file 'wfns' does not specify real/complex; cannot check compatibility."
        call messages_warning(1)
      end if
    end if
    ! complex can be restarted from real, so there is no problem.

    ! Skip two lines.
    call restart_read(restart, occ_file, lines, 2, err)
    if (err /= 0) ierr = ierr - 2**7

    ! If any error occured up to this point then it is not worth continuing,
    ! as there something fundamentally wrong with the restart files
    if (ierr /= 0) then
      call restart_close(restart, wfns_file)
      call restart_close(restart, occ_file)
      call profiling_out(prof_read)
      POP_SUB(states_load)
      return
    end if

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np))
    else
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
    end if

    SAFE_ALLOCATE(restart_file(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    SAFE_ALLOCATE(restart_file_present(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    restart_file_present = .false.

    ! Next we read the list of states from the files. 
    ! Errors in reading the information of a specific state from the files are ignored
    ! at this point, because later we will skip reading the wavefunction of that state.
    do
      call restart_read(restart, wfns_file, lines, 1, err)
      if (err == 0) then
        read(lines(1), '(a)') char
        if (char == '%') then
          !We reached the end of the file
          exit
        else
          read(lines(1), *) ik, char, ist, char, idim, char, filename
        end if
      end if

      if (err /= 0 .or. index_is_wrong()) then
        call restart_read(restart, occ_file, lines, 1, err)
        cycle
      end if

      if (ist >= st%st_start .and. ist <= st%st_end .and. &
           st%d%kpt%start <= ik .and. st%d%kpt%end >= ik) then
          
        restart_file(idim, ist, ik) = trim(filename)
        restart_file_present(idim, ist, ik) = .true.
      end if

      call restart_read(restart, occ_file, lines, 1, err)
      if (.not. present(lr)) then ! do not read eigenvalues or occupations when reading linear response
        ! # occupations | eigenvalue[a.u.] | Im(eigenvalue) [a.u.] | k-points | k-weights | filename | ik | ist | idim

        if (err == 0) then
          read(lines(1), *) my_occ, char, st%eigenval(ist, ik), char, imev, char, &
               (read_kpoint(idir), char, idir = 1, gr%sb%dim), my_kweight
          ! we do not want to read the k-weights, we have already set them appropriately
        else
          ! There is a problem with this states information, so we skip it.
          restart_file_present(idim, ist, ik) = .false.
          cycle
        end if

        kpoint(1:gr%sb%dim) = &
          kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik), absolute_coordinates = .true.)
        ! FIXME: maybe should ignore ik and just try to match actual vector k-points?
        if (any(abs(kpoint(1:gr%sb%dim) - read_kpoint(1:gr%sb%dim)) > CNST(1e-12))) then
          ! write only once for each k-point so as not to be too verbose
          if (ist == 1) then
            write(message(1),'(a,i6)') 'Incompatible restart information: k-point mismatch for ik ', ik
            write(message(2),'(a,99f18.12)') '  Expected : ', kpoint(1:gr%sb%dim)
            write(message(3),'(a,99f18.12)') '  Read     : ', read_kpoint(1:gr%sb%dim)
            call messages_warning(3)
          end if
          restart_file_present(idim, ist, ik) = .false.
        end if

        if (read_occ) then
          st%occ(ist, ik) = my_occ
          integral_occs = integral_occs .and. &
               abs((st%occ(ist, ik) - st%smear%el_per_state) * st%occ(ist, ik))  <=  M_EPSILON
        end if
      end if
    end do

    if (present(iter)) then
      call restart_read(restart, wfns_file, lines, 1, err)
      if (err /= 0) then
        ierr = ierr - 2**8
      else
        read(lines(1), *) filename, filename, iter
      end if
    end if

    call restart_close(restart, wfns_file)
    call restart_close(restart, occ_file)

    if (st%restart_fixed_occ) then
      ! reset to overwrite whatever smearing may have been set earlier
      call smear_init(st%smear, parser, st%d%ispin, fixed_occ = .true., integral_occs = integral_occs, kpoints = gr%sb%kpoints)
    end if


    ! Now we read the wavefunctions. At this point we need to have all the information from the
    ! states, occs, and wfns files in order to avoid serialisation of reads, as restart_read
    ! works as a barrier.

    SAFE_ALLOCATE(filled(1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
    filled = .false.

    if (present(lowest_missing)) lowest_missing = st%nst + 1

    iread = 0
    if (mpi_grp_is_root(mpi_world) .and. verbose_) then
      idone = 0
      ntodo = st%lnst*st%d%kpt%nlocal*st%d%dim
      call loct_progress_bar(-1, ntodo)
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim

          if (.not. restart_file_present(idim, ist, ik)) then
            if (present(lowest_missing)) &
              lowest_missing(idim, ik) = min(lowest_missing(idim, ik), ist)
            cycle
          end if

          if (states_are_real(st)) then
            call drestart_read_mesh_function(restart, restart_file(idim, ist, ik), gr%mesh, dpsi, err)
          else
            call zrestart_read_mesh_function(restart, restart_file(idim, ist, ik), gr%mesh, zpsi, err)
          end if

          if(states_are_real(st)) then
            if(.not. present(lr)) then
              call states_set_state(st, gr%mesh, idim, ist, ik, dpsi)
            else
              call lalg_copy(gr%mesh%np, dpsi, lr%ddl_psi(:, idim, ist, ik))
            end if
          else
            if(.not. present(lr)) then
              call states_set_state(st, gr%mesh, idim, ist, ik, zpsi)
            else
              call lalg_copy(gr%mesh%np, zpsi, lr%zdl_psi(:, idim, ist, ik))
            end if
          end if


          if (err == 0) then
            filled(idim, ist, ik) = .true.
            iread = iread + 1
          else if (present(lowest_missing)) then
            lowest_missing(idim, ik) = min(lowest_missing(idim, ik), ist)
          end if

          if (mpi_grp_is_root(mpi_world) .and. verbose_) then
            idone = idone + 1
            call loct_progress_bar(idone, ntodo)
          end if

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zpsiL)
    SAFE_DEALLOCATE_A(restart_file)
    SAFE_DEALLOCATE_A(restart_file_present)

    if(mpi_grp_is_root(mpi_world) .and. verbose_) then
      call messages_new_line()
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      iread_tmp = iread
      call MPI_Allreduce(iread_tmp, iread, 1, MPI_INTEGER, MPI_SUM, st%st_kpt_mpi_grp%comm, mpi_err)
    end if

    if(st%d%kpt%parallel) then
      ! make sure all tasks know lowest_missing over all k-points
      if(present(lowest_missing)) then
        SAFE_ALLOCATE(lowest_missing_tmp(1:st%d%dim, 1:st%d%nik))
        lowest_missing_tmp = lowest_missing
        call MPI_Allreduce(lowest_missing_tmp, lowest_missing, st%d%dim*st%d%nik, &
          MPI_INTEGER, MPI_MIN, st%d%kpt%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(lowest_missing_tmp)
      end if
    end if
#endif

    if (.not. present(lr)) call fill_random()
    ! it is better to initialize lr wfns to zero

    SAFE_DEALLOCATE_A(filled)

    if (ierr == 0 .and. iread /= st%nst * st%d%nik * st%d%dim) then
      if(iread > 0) then
        ierr = iread
      else
        ierr = -1
      end if
      ! otherwise ierr = 0 would mean either all was read correctly, or nothing at all was read!

      if(.not. present(lr)) then
        write(str, '(a,i5)') 'Reading states.'
      else
        write(str, '(a,i5)') 'Reading states information for linear response.'
      end if
      call messages_print_stress(stdout, trim(str))
      write(message(1),'(a,i6,a,i6,a)') 'Only ', iread,' files out of ', &
           st%nst * st%d%nik * st%d%dim, ' could be read.'
      call messages_info(1)
      call messages_print_stress(stdout)
    end if

    message(1) = 'Info: States reading done.'
    if(verbose_) call print_date(trim(message(1))//' ')

    call profiling_out(prof_read)
    POP_SUB(states_load)

  contains

    ! ---------------------------------------------------------
    subroutine fill_random() !< Put random function in orbitals that could not be read.
      PUSH_SUB(states_load.fill_random)

      do ik = st%d%kpt%start, st%d%kpt%end

        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            if(filled(idim, ist, ik)) cycle

            call states_generate_random(st, gr%mesh, gr%sb, ist, ist, ik, ik)
          end do
        end do
      end do

      POP_SUB(states_load.fill_random)
    end subroutine fill_random

    ! ---------------------------------------------------------

    logical function index_is_wrong() !< .true. if the index (idim, ist, ik) is not present in st structure...
      PUSH_SUB(states_load.index_is_wrong)

      if(idim > st%d%dim .or. idim < 1 .or.   &
        ist   > st%nst   .or. ist  < 1 .or.   &
        ik    > st%d%nik .or. ik   < 1) then
        index_is_wrong = .true.
      else
        index_is_wrong = .false.
      end if

      POP_SUB(states_load.index_is_wrong)
    end function index_is_wrong

  end subroutine states_load


  subroutine states_dump_rho(restart, st, gr, ierr, iter)
    type(restart_t),      intent(in)    :: restart
    type(states_t),       intent(in)    :: st
    type(grid_t),         intent(in)    :: gr
    integer,              intent(out)   :: ierr
    integer,    optional, intent(in)    :: iter

    integer :: iunit, isp, err, err2(2)
    character(len=80) :: filename
    character(len=300) :: lines(2)
    FLOAT, pointer :: rho(:), rho_fine(:)

    PUSH_SUB(states_dump_rho)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(states_dump_rho)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing density restart."
      call messages_info(1)
    end if

    call restart_block_signals()

    !write the densities
    iunit = restart_open(restart, 'density')
    lines(1) = '#     #spin    #nspin    filename'
    lines(2) = '%densities'
    call restart_write(restart, iunit, lines, 2, err)
    if (err /= 0) ierr = ierr + 1

    if(gr%have_fine_mesh) then
      SAFE_ALLOCATE(rho(1:gr%mesh%np))
      SAFE_ALLOCATE(rho_fine(1:gr%fine%mesh%np))
    end if

    err2 = 0
    do isp = 1, st%d%nspin
      if(st%d%nspin==1) then
        write(filename, fmt='(a)') 'density'
      else
        write(filename, fmt='(a,i1)') 'density-sp', isp
      end if
      write(lines(1), '(i8,a,i8,a)') isp, ' | ', st%d%nspin, ' | "'//trim(adjustl(filename))//'"'
      call restart_write(restart, iunit, lines, 1, err)
      if (err /= 0) err2(1) = err2(1) + 1

      if(gr%have_fine_mesh)then
        rho_fine(1:gr%fine%mesh%np) = st%rho(1:gr%fine%mesh%np,isp)
        call dmultigrid_fine2coarse(gr%fine%tt, gr%fine%der, gr%mesh, rho_fine, rho, INJECTION)
        call drestart_write_mesh_function(restart, filename, gr%mesh, rho, err)
      else
        call drestart_write_mesh_function(restart, filename, gr%mesh, st%rho(:,isp), err)
      end if
      if (err /= 0) err2(2) = err2(2) + 1

    end do
    if (err2(1) /= 0) ierr = ierr + 2
    if (err2(2) /= 0) ierr = ierr + 4

    if(gr%have_fine_mesh)then
      SAFE_DEALLOCATE_P(rho)
      SAFE_DEALLOCATE_P(rho_fine)
    end if

    lines(1) = '%'
    call restart_write(restart, iunit, lines, 1, err)
    if (err /= 0) ierr = ierr + 8
    if (present(iter)) then
      write(lines(1),'(a,i7)') 'Iter = ', iter
      call restart_write(restart, iunit, lines, 1, err)
      if (err /= 0) ierr = ierr + 16
    end if
    call restart_close(restart, iunit)

    call restart_unblock_signals()

    if (debug%info) then
      message(1) = "Debug: Writing density restart done."
      call messages_info(1)
    end if

    POP_SUB(states_dump_rho)
  end subroutine states_dump_rho


  ! ---------------------------------------------------------
  subroutine states_load_rho(restart, st, gr, ierr)
    type(restart_t), intent(in)    :: restart
    type(states_t),  intent(inout) :: st
    type(grid_t),    intent(in)    :: gr
    integer,         intent(out)   :: ierr

    integer              :: err, err2, isp
    character(len=12)    :: filename
    FLOAT, allocatable   :: rho_coarse(:)

    PUSH_SUB(states_load_rho)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load_rho)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading density restart."
      call messages_info(1)
    end if

    ! skip for now, since we know what the files are going to be called
    !read the densities
!    iunit_rho = io_open(trim(dir)//'/density', action='write')
!    call iopar_read(st%dom_st_kpt_mpi_grp, iunit_rho, line, err)
!    call iopar_read(st%dom_st_kpt_mpi_grp, iunit_rho, line, err)
!   we could read the iteration 'iter' too, not sure if that is useful.

    if(gr%have_fine_mesh) then
      SAFE_ALLOCATE(rho_coarse(1:gr%mesh%np_part))
    end if

    err2 = 0
    do isp = 1, st%d%nspin
      if(st%d%nspin==1) then
        write(filename, fmt='(a)') 'density'
      else
        write(filename, fmt='(a,i1)') 'density-sp', isp
      end if
!      if(mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
!        read(iunit_rho, '(i8,a,i8,a)') isp, ' | ', st%d%nspin, ' | "'//trim(adjustl(filename))//'"'
!      end if
      if(gr%have_fine_mesh)then
        call drestart_read_mesh_function(restart, filename, gr%mesh, rho_coarse, err)
        call dmultigrid_coarse2fine(gr%fine%tt, gr%der, gr%fine%mesh, rho_coarse, st%rho(:,isp), order = 2)
      else
        call drestart_read_mesh_function(restart, filename, gr%mesh, st%rho(:,isp), err)
      end if
      if (err /= 0) err2 = err2 + 1

    end do
    if (err2 /= 0) ierr = ierr + 1

    if(gr%have_fine_mesh)then
      SAFE_DEALLOCATE_A(rho_coarse)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading density restart done."
      call messages_info(1)
    end if

    POP_SUB(states_load_rho)
  end subroutine states_load_rho

  ! ---------------------------------------------------------
  !> the routine reads formulas for user-defined wavefunctions
  !! from the input file and fills the respective orbitals
  subroutine states_read_user_def_orbitals(mesh, st)
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(inout) :: st

    type(block_t) :: blk
    integer :: ip, id, is, ik, nstates, state_from, ierr, ncols
    integer :: ib, idim, inst, inik, normalize
    FLOAT :: xx(MAX_DIM), rr, psi_re, psi_im
    character(len=150) :: filename
    CMPLX, allocatable :: zpsi(:, :)

    integer, parameter ::           &
      STATE_FROM_FORMULA  = 1,      &
      STATE_FROM_FILE     = -10010, &
      NORMALIZE_YES       = 1,      &
      NORMALIZE_NO        = 0

    PUSH_SUB(states_read_user_def_orbitals)

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
    !%Option file -10010
    !% Read initial orbital from file. 
    !% Accepted file formats, detected by extension: obf, ncdf and csv (real only).
    !%Option formula 1
    !% Calculate initial orbital by given analytic expression.
    !%Option normalize_yes 1
    !% Normalize orbitals (default).
    !%Option normalize_no 0
    !% Do not normalize orbitals.
    !%End
    if(parse_block('UserDefinedStates', blk) == 0) then

      call messages_print_stress(stdout, trim('Substitution of orbitals'))

      ! find out how many lines (i.e. states) the block has
      nstates = parse_block_n(blk)

      SAFE_ALLOCATE(zpsi(1:mesh%np, 1:st%d%dim))

      ! read all lines
      do ib = 1, nstates
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, ib - 1)
        if(ncols  <  5 .or. ncols > 6) then
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
              if(.not.(id  ==  idim .and. is  ==  inst .and. ik  ==  inik    &
                .and. st%st_start  <=  is .and. st%st_end >= is          &
                .and. st%d%kpt%start  <=  ik .and. st%d%kpt%end >= ik) ) cycle

              select case(state_from)

              case(STATE_FROM_FORMULA)
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
                  zpsi(ip, 1) = psi_re + M_zI * psi_im
                end do

              case(STATE_FROM_FILE)
                ! The input format can be coded in column four now. As it is
                ! not used now, we just say "file".
                ! Read the filename.
                call parse_block_string(blk, ib - 1, 4, filename)

                write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
                write(message(2), '(2a)') '  with data from file:'
                write(message(3), '(2a)') '  ',trim(filename)
                call messages_info(3)

                ! finally read the state
                call zio_function_input(filename, mesh, zpsi(:, 1), ierr)
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

              call states_set_state(st, mesh, id, is, ik, zpsi(:, 1))

              ! normalize orbital
              if(parse_block_cols(blk, ib - 1)  ==  6) then
                call parse_block_integer(blk, ib - 1, 5, normalize)
              else
                normalize = 1
              end if
              select case(normalize)
              case(NORMALIZE_NO)
              case(NORMALIZE_YES)
                call states_get_state(st, mesh, is, ik, zpsi)
                call zmf_normalize(mesh, st%d%dim, zpsi)
                call states_set_state(st, mesh, is, ik, zpsi)
              case default
                message(1) = 'The sixth column in UserDefinedStates may either be'
                message(2) = '"normalize_yes" or "normalize_no"'
                call messages_fatal(2)
              end select

            end do
          end do
        end do

      end do

      SAFE_DEALLOCATE_A(zpsi)

      call parse_block_end(blk)
      call messages_print_stress(stdout)

    else
      message(1) = "'UserDefinedStates' has to be specified as block."
      call messages_fatal(1)
    end if

    POP_SUB(states_read_user_def_orbitals)
  end subroutine states_read_user_def_orbitals

  ! ---------------------------------------------------------
  subroutine states_dump_spin(restart, st, gr, ierr)
    type(restart_t),      intent(in)  :: restart
    type(states_t),       intent(in)  :: st
    type(grid_t),         intent(in)  :: gr
    integer,              intent(out) :: ierr

    integer :: iunit_spin
    integer :: err, err2(2), ik, idir, ist
    character(len=300) :: lines(3)

    PUSH_SUB(states_dump_spin)

    ierr = 0

    if(restart_skip(restart)) then
      POP_SUB(states_dump_spin)
      return
    end if

    call profiling_in(prof_write, "RESTART_WRITE_SPIN")

    call restart_block_signals()

    iunit_spin = restart_open(restart, 'spin')
    lines(1) = '#     #k-point            #st       #spin(x) spin(y) spin(z)'
    call restart_write(restart, iunit_spin, lines, 1, err)
    if (err /= 0) ierr = ierr + 1

    err2 = 0
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        write(lines(1), '(i8,a,i8,3(a,f18.12))') ik, ' | ', ist, ' | ', &
                            st%spin(1,ist,ik), ' | ', st%spin(2,ist,ik),' | ', st%spin(3,ist,ik)
        call restart_write(restart, iunit_spin, lines, 1, err)
        if (err /= 0) err2(1) = err2(1) + 1
      end do ! st%nst
    end do ! st%d%nik
    lines(1) = '%'
    call restart_write(restart, iunit_spin, lines, 1, err)
    if (err2(1) /= 0) ierr = ierr + 8
    if (err2(2) /= 0) ierr = ierr + 16

    call restart_close(restart, iunit_spin)

    call restart_unblock_signals()

    call profiling_out(prof_write)
    POP_SUB(states_dump_spin)
  end subroutine states_dump_spin




  ! ---------------------------------------------------------
  !> returns in ierr:
  !! <0 => Fatal error, or nothing read
  !! =0 => read all wavefunctions
  !! >0 => could only read ierr wavefunctions
  subroutine states_load_spin(restart, st, gr, ierr)
    type(restart_t),            intent(in)    :: restart
    type(states_t),             intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    integer,                    intent(out)   :: ierr

    integer              :: spin_file, err, ik, ist
    character(len=256)   :: lines(3)
    FLOAT                :: spin(3)
    character(len=1)     :: char

    
    PUSH_SUB(states_load_spin)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load_spin)
      return
    end if

    call profiling_in(prof_read, "RESTART_READ_SPIN")

    ! open files to read
    spin_file = restart_open(restart, 'spin')
    ! Skip two lines.
    call restart_read(restart, spin_file, lines, 1, err)
    if (err /= 0) ierr = ierr - 2**7

    ! If any error occured up to this point then it is not worth continuing,
    ! as there something fundamentally wrong with the restart files
    if (ierr /= 0) then
      call restart_close(restart, spin_file)
      call profiling_out(prof_read)
      POP_SUB(states_load_spin)
      return
    end if

    ! Next we read the list of states from the files. 
    ! Errors in reading the information of a specific state from the files are ignored
    ! at this point, because later we will skip reading the wavefunction of that state.
    do
      call restart_read(restart, spin_file, lines, 1, err)
      read(lines(1), '(a)') char
      if (char == '%') then
          !We reached the end of the file
          exit
        else
        read(lines(1), *) ik, char, ist, char,  spin(1), char, spin(2), char, spin(3)
      end if

      if (err /= 0) cycle

      st%spin(1:3, ist, ik) = spin(1:3)
    end do

    call restart_close(restart, spin_file)

    call profiling_out(prof_read)
    POP_SUB(states_load_spin)
 end subroutine states_load_spin

end module states_restart_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
