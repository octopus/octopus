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

module states_restart_m
  use datasets_m
  use global_m
  use grid_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_basic_m
  use linear_response_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multigrid_m
  use ob_interface_m
  use parser_m
  use profiling_m
  use restart_m
  use simul_box_m
  use smear_m
  use states_m
  use states_calc_m
  use states_dim_m
  use states_io_m
  use string_m
  use unit_m
  use unit_system_m
  use types_m

  implicit none

  type(profile_t), save :: prof_read, prof_write

  private
  public ::                         &
    states_look_and_load,           &
    states_dump,                    &
    states_load,                    &
    states_dump_rho,                &
    states_load_rho,                &
    states_load_free_states,        &
    states_read_user_def_orbitals

contains

  ! ---------------------------------------------------------
  subroutine states_look_and_load(restart, st, gr, is_complex)
    type(restart_t),            intent(inout) :: restart
    type(states_t),     target, intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    logical,          optional, intent(in)    :: is_complex

    integer :: kpoints, dim, nst, ierr
    FLOAT, pointer :: new_occ(:,:)

    PUSH_SUB(states_look_and_load)

    !check how many wfs we have
    call states_look(restart, gr%mesh%mpi_grp, kpoints, dim, nst, ierr)
    if(ierr /= 0) then
      message(1) = 'Could not read wavefunctions.'
      call messages_fatal(1)
    end if

    ! Resize st%occ, retaining information there
    SAFE_ALLOCATE(new_occ(1:nst, 1:st%d%nik))
    new_occ(:,:) = M_ZERO
    new_occ(1:min(nst, st%nst),:) = st%occ(1:min(nst, st%nst),:)
    SAFE_DEALLOCATE_P(st%occ)
    st%occ => new_occ

    ! FIXME: This wrong, one cannot just change the number of states
    ! without updating the internal structures.
    st%nst    = nst
    st%st_end = nst
    SAFE_DEALLOCATE_P(st%zeigenval%Re)
    nullify(st%eigenval)
    
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

    SAFE_ALLOCATE(st%zeigenval%Re(1:st%nst, 1:st%d%nik))
    st%eigenval => st%zeigenval%Re

    if(st%d%ispin == SPINORS) then
      SAFE_ALLOCATE(st%spin(1:3, 1:st%nst, 1:st%d%nik))
      st%spin = M_ZERO
    end if
    st%eigenval = M_HUGE

    ! load wavefunctions
    call states_load(restart, st, gr, ierr)

    POP_SUB(states_look_and_load)
  end subroutine states_look_and_load


  ! ---------------------------------------------------------
  subroutine states_dump(restart, st, gr, ierr, iter, lr, st_start_writing)
    type(restart_t),      intent(in)  :: restart
    type(states_t),       intent(in)  :: st
    type(grid_t),         intent(in)  :: gr
    integer,              intent(out) :: ierr
    integer,    optional, intent(in)  :: iter
    !> if this next argument is present, the lr wfs are stored instead of the gs wfs
    type(lr_t), optional, intent(in)  :: lr
    integer,    optional, intent(in)  :: st_start_writing

    integer :: iunit_wfns, iunit_occs, iunit_states
    integer :: err, ik, idir, ist, idim, itot
    character(len=80) :: filename, filename1 
    logical :: lr_wfns_are_associated, should_write, cmplxscl
    FLOAT   :: kpoint(1:MAX_DIM)
    FLOAT,  allocatable :: dpsi(:)
    CMPLX,  allocatable :: zpsi(:)

    PUSH_SUB(states_dump)

    ierr = 0

    if(restart_skip(restart)) then
      POP_SUB(states_dump)
      return
    end if

    message(1) = "Info: Writing states."
    call print_date(trim(message(1))//' ')

    call profiling_in(prof_write, "RESTART_WRITE")

    cmplxscl = .false.
    if (associated(st%zeigenval%Im)) cmplxscl = .true.

    if (present(lr)) then
      lr_wfns_are_associated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_wfns_are_associated)
    end if

    call restart_block_signals()

    iunit_wfns = restart_open(restart, 'wfns')
    if (iunit_wfns < 0) ierr = 1
    iunit_occs = restart_open(restart, 'occs')
    if (iunit_occs < 0) ierr = 1
    iunit_states = restart_open(restart, 'states')
    if (iunit_states < 0) ierr = 1

    if (mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
      if (iunit_wfns > 0) then
        write(iunit_wfns,'(a)') '#     #k-point            #st            #dim    filename'
        if (states_are_real(st)) then
          write(iunit_wfns,'(a)') '%Real_Wavefunctions'
        else
          write(iunit_wfns,'(a)') '%Complex_Wavefunctions'
        end if
      end if

      if (iunit_occs > 0) then
        write(iunit_occs,'(2a)') '# occupations | eigenvalue[a.u.] | Im(eigenvalue) [a.u.] ', &
             '| k-points | k-weights | filename | ik | ist | idim'
        write(iunit_occs,'(a)') '%Occupations_Eigenvalues_K-Points'
      end if

      if (iunit_states > 0) then
        write(iunit_states, '(a20,1i10)')  'nst=                ', st%nst
        write(iunit_states, '(a20,1i10)')  'dim=                ', st%d%dim
        write(iunit_states, '(a20,1i10)')  'nik=                ', st%d%nik
      end if
    end if
    if (iunit_states > 0) call restart_close(restart, iunit_states)

    if(states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np))
    else
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
    end if


    itot = 1
    do ik = 1, st%d%nik
      kpoint = M_ZERO
      kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik), absolute_coordinates = .true.)

      do ist = 1, st%nst
        do idim = 1, st%d%dim
          write(filename,'(i10.10)') itot
          if (st%have_left_states) filename1 = 'L'//trim(filename) !cmplxscl

          if (mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
            if (iunit_wfns > 0) then
              write(iunit_wfns, '(i8,a,i8,a,i8,3a)') ik, ' | ', ist, ' | ', idim, ' | "', trim(filename), '"'
            end if

            if (iunit_occs > 0) then
              write(iunit_occs, '(e21.14,a,e21.14,a)', advance='no') st%occ(ist,ik), ' | ', st%eigenval(ist, ik)
              if (cmplxscl) then
                write(iunit_occs, '(a,e21.14)', advance='no') ' | ', st%zeigenval%Im(ist, ik)
              else
                write(iunit_occs, '(a,e21.14)', advance='no') ' | ', CNST(0.0)
              end if
              write(iunit_occs, '(a)', advance='no')  ' | '
              do idir = 1, gr%sb%dim
                write(iunit_occs, '(e21.14,a)', advance='no') kpoint(idir), ' | '
              end do
              write(iunit_occs, '(e21.14,a,i10.10,3(a,i8))') st%d%kweights(ik), ' | ', itot, ' | ', ik, ' | ', ist, ' | ', idim
            end if
          end if

          should_write = st%st_start <= ist .and. ist <= st%st_end
          if (should_write .and. present(st_start_writing)) then
            if (ist < st_start_writing) then
              should_write = .false.
              ierr = ierr + 1
            end if
          end if

          if (should_write) then
            if (.not. present(lr)) then
              if(st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
                if (states_are_real(st)) then
                  call states_get_state(st, gr%mesh, idim, ist, ik, dpsi)
                  call drestart_write_function(restart, filename, gr%mesh, dpsi, err)
                else
                  call states_get_state(st, gr%mesh, idim, ist, ik, zpsi)
                  call zrestart_write_function(restart, filename, gr%mesh, zpsi, err)
                  if(st%have_left_states) then!cmplxscl
                    call states_get_state(st, gr%mesh, idim, ist, ik, zpsi, left = .true.)
                    call zrestart_write_function(restart, filename1, gr%mesh, zpsi, err)
                  end if
                end if
                if (err == 0) then
                  ierr = ierr + 1
                else
                  message(1) = "Unable to write state wavefunction to '"//trim(filename)//"'."
                  call messages_warning(1)
                end if
              end if
            else
              if (states_are_real(st)) then
                call drestart_write_function(restart, filename, gr%mesh, lr%ddl_psi(:, idim, ist, ik), err)
              else
                call zrestart_write_function(restart, filename, gr%mesh, lr%zdl_psi(:, idim, ist, ik), err)
              end if
              if (err == 0) then
                ierr = ierr + 1
              else
                message(1) = "Unable to write state wavefunction to '"//trim(filename)//"'."
                call messages_warning(1)
              end if
            end if
          end if

          itot = itot + 1
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)

    ! do NOT use st%lnst here as it is not (st%st_end - st%st_start + 1)
    if (ierr == st%d%kpt%nlocal * (st%st_end - st%st_start + 1) * st%d%dim) ierr = 0 ! All OK

    if (mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
      if (iunit_wfns > 0) then
        write(iunit_wfns,'(a)') '%'
        if(present(iter)) write(iunit_wfns,'(a,i7)') 'Iter = ', iter
      end if
      if (iunit_occs > 0) write(iunit_occs, '(a)') '%'
    end if
    if (iunit_wfns > 0) call restart_close(restart, iunit_wfns)
    if (iunit_occs > 0) call restart_close(restart, iunit_occs)

    message(1) = "Info: Finished writing states."
    call print_date(trim(message(1))//' ')

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
  subroutine states_load(restart, st, gr, ierr, iter, lr, lowest_missing, read_left, label, verbose, kpts_dont_matter)
    type(restart_t),            intent(inout) :: restart
    type(states_t),             intent(inout) :: st
    type(grid_t),               intent(in)    :: gr
    integer,                    intent(out)   :: ierr
    integer,          optional, intent(out)   :: iter
    type(lr_t),       optional, intent(inout) :: lr       !< if present, the lr wfs are read instead of the gs wfs
    integer,          optional, intent(out)   :: lowest_missing(:, :) !< all states below this one were read successfully
    logical,          optional, intent(in)    :: read_left !< if .true. read left states (default is .false.)
    character(len=*), optional, intent(in)    :: label
    logical,          optional, intent(in)    :: verbose
    logical,          optional, intent(in)    :: kpts_dont_matter !< don`t check k-points match with current calculation
      !! for td transport; they have been set to zero in simul_box_init, and won`t match regardless due to box size

    integer              :: states_file, wfns_file, occ_file, err, ik, ist, idir, idim
    integer              :: iread, nread
    character(len=12)    :: filename
    character(len=1)     :: char
    logical, allocatable :: filled(:, :, :)
    character(len=256)   :: lines(3), label_, mod_time, occ_filename
    character(len=50)    :: str

    FLOAT                :: my_occ, imev
    logical              :: read_occ, lr_allocated, verbose_
    logical              :: integral_occs, cmplxscl, read_left_
    FLOAT, allocatable   :: dpsi(:)
    CMPLX, allocatable   :: zpsi(:), zpsiL(:)
    character(len=256), allocatable :: restart_file(:, :, :)
    logical,            allocatable :: restart_file_present(:, :, :)
    FLOAT                :: kpoint(MAX_DIM), read_kpoint(MAX_DIM)

    PUSH_SUB(states_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load)
      return
    end if

    call profiling_in(prof_read, "RESTART_READ")

    verbose_ = optional_default(verbose, .true.)

    cmplxscl = .false.
    if (associated(st%zeigenval%Im)) cmplxscl = .true.
    
    read_left_ = optional_default(read_left, .false.)
    if (read_left_) then
       ASSERT(st%have_left_states)
    end if

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
      if (cmplxscl) st%zeigenval%Im(:, :) = M_ZERO
      ! to be filled in from reading afterward
    endif

    if (.not. present(lr) .and. read_occ) then
      st%occ(:, :) = M_ZERO
      ! to be filled in from reading afterward
    endif

    ! sanity check
    if (present(lr)) then
      lr_allocated = (associated(lr%ddl_psi) .and. states_are_real(st)) .or. &
        (associated(lr%zdl_psi) .and. states_are_complex(st))
      ASSERT(lr_allocated)
    endif

    ! make sure these intent(out)`s are initialized no matter what
    if (present(lowest_missing)) lowest_missing = 1
    if (present(iter)) iter = 0

    ! open files to read
    states_file  = restart_open(restart, 'states')
    if (states_file < 0) ierr = -1

    wfns_file  = restart_open(restart, 'wfns')
    if (wfns_file < 0) ierr = -1

    occ_file = restart_open(restart, 'occs')
    if (occ_file < 0) ierr = -1

    ! sanity check on spin/k-points. Example file 'states':
    ! nst=                         2
    ! dim=                         1
    ! nik=                         2
    if (ierr == 0) call iopar_read(st%dom_st_kpt_mpi_grp, states_file, lines, 3, ierr)
    if (ierr == 0) then
      read(lines(2), *) str, idim
      read(lines(3), *) str, ik
      if(idim == 2 .and. st%d%dim == 1) then
        write(message(1),'(a)') 'Incompatible restart information: saved calculation is spinors, this one is not.'
        call messages_warning(1)
        ierr = -1
      endif
      if(idim == 1 .and. st%d%dim == 2) then
        write(message(1),'(a)') 'Incompatible restart information: this calculation is spinors, saved one is not.'
        call messages_warning(1)
        ierr = -1
      endif
      if(ik /= st%d%nik) then
        write(message(1),'(a)') 'Incompatible restart information: wrong number of k-points.'
        write(message(2),'(2(a,i6))') 'Expected ', st%d%nik, '; Read ', ik
        call messages_warning(2)
        ierr = -1
      endif
      ! FIXME: we could restart anyway if we can check that existing k-points are used correctly
    end if
    if (states_file > 0) call restart_close(restart, states_file)

    if (ierr == 0) call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, lines, 2, ierr)
    if (ierr == 0 .and. states_are_real(st)) then
      read(lines(2), '(a)') str
      if (str(2:8) == 'Complex') then
        message(1) = "Cannot read real states from complex wavefunctions."
        call messages_warning(1)
        ierr = -2
      else if (str(2:5) /= 'Real') then 
        message(1) = "Restart file 'wfns' does not specify real/complex; cannot check compatibility."
        call messages_warning(1)
      end if
    end if
    ! complex can be restarted from real, so there is no problem.

    ! Skip two lines.
    if (ierr == 0) call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, lines, 2, ierr)

    ! If any error occured up to this point then it is not worth continuing,
    ! as there something fundamentally wrong the restart files
    if (ierr /= 0) then
      if (wfns_file > 0) call restart_close(restart, wfns_file)
      if (occ_file  > 0) call restart_close(restart, occ_file)
      call profiling_out(prof_read)
      POP_SUB(states_load)
      return
    end if

    if (states_are_real(st)) then
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np))
    else
      SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
      if(read_left_) then
        SAFE_ALLOCATE(zpsiL(1:gr%mesh%np))
      end if
    end if

    SAFE_ALLOCATE(restart_file(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    SAFE_ALLOCATE(restart_file_present(1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    restart_file_present = .false.

    ! Next we read the list of states from the files. 
    ! Errors in reading the information of a specific state from the files are ignored
    ! at this point, because later we will skip reading the wavefunction of that state.
    do
      call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, lines, 1, err)
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
        call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, lines, 1, err)
        cycle
      end if

      if (ist >= st%st_start .and. ist <= st%st_end .and. &
           st%d%kpt%start <= ik .and. st%d%kpt%end >= ik) then
          
        restart_file(idim, ist, ik) = trim(filename)
        restart_file_present(idim, ist, ik) = .true.
      end if

      call iopar_read(st%dom_st_kpt_mpi_grp, occ_file, lines, 1, err)
      if (.not. present(lr)) then ! do not read eigenvalues or occupations when reading linear response
        ! # occupations | eigenvalue[a.u.] | Im(eigenvalue) [a.u.] | k-points | k-weights | filename | ik | ist | idim

        if (err == 0) then
          read(lines(1), *) my_occ, char, st%eigenval(ist, ik), char, imev, char, &
               (read_kpoint(idir), char, idir = 1, gr%sb%dim), st%d%kweights(ik)
        else
          ! There is a problem with this states information, so we skip it.
          restart_file_present(idim, ist, ik) = .false.
          cycle
        end if

        kpoint(1:gr%sb%dim) = &
          kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik), absolute_coordinates = .true.)
        ! FIXME: maybe should ignore ik and just try to match actual vector k-points?
        if ((.not. optional_default(kpts_dont_matter, .false.)) .and. &
          any(abs(kpoint(1:gr%sb%dim) - read_kpoint(1:gr%sb%dim)) > CNST(1e-12))) then
          ! write only once for each k-point so as not to be too verbose
          if (ist == 1) then
            write(message(1),'(a,i6)') 'Incompatible restart information: k-point mismatch for ik ', ik
            write(message(2),'(a,99f18.12)') '  Expected : ', kpoint(1:gr%sb%dim)
            write(message(3),'(a,99f18.12)') '  Read     : ', read_kpoint(1:gr%sb%dim)
            call messages_warning(3)
          end if
          restart_file_present(idim, ist, ik) = .false.
        end if

        if (cmplxscl) st%zeigenval%Im(ist, ik) = imev

        if (read_occ) then
          st%occ(ist, ik) = my_occ
          integral_occs = integral_occs .and. &
               abs((st%occ(ist, ik) - st%smear%el_per_state) * st%occ(ist, ik))  <=  M_EPSILON
        end if
      end if
    end do

    if (present(iter)) then
      call iopar_read(st%dom_st_kpt_mpi_grp, wfns_file, lines, 1, err)
      if (err == 0) read(lines(1), *) filename, filename, iter
    end if

    call restart_close(restart, wfns_file)
    call restart_close(restart, occ_file)

    if (st%restart_fixed_occ) then
      ! reset to overwrite whatever smearing may have been set earlier
      call smear_init(st%smear, st%d%ispin, fixed_occ = .true., integral_occs = integral_occs, kpoints = gr%sb%kpoints)
    endif


    ! Now we read the wavefunction. At this point we must have all the information from the
    ! states, occs, and wfns files in order to avoid serialisation of reads, as iopar_read
    ! works as a barrier.

    SAFE_ALLOCATE(filled(1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
    filled = .false.

    if (present(lowest_missing)) lowest_missing = st%nst + 1

    if (mpi_grp_is_root(mpi_world) .and. verbose_) then
      iread = 1
      nread = st%lnst*st%d%kpt%nlocal*st%d%dim
      call loct_progress_bar(-1, nread)
    end if

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim

          if (.not. restart_file_present(idim, ist, ik)) then
            if (present(lowest_missing)) &
              lowest_missing(idim, ik) = min(lowest_missing(idim, ik), ist)
            cycle
          endif

          if (states_are_real(st)) then
            call drestart_read_function(restart, restart_file(idim, ist, ik), gr%mesh, dpsi, err)
          else
            call zrestart_read_function(restart, restart_file(idim, ist, ik), gr%mesh, zpsi, err)
            if (read_left_) call zrestart_read_function(restart, 'L'//restart_file(idim, ist, ik), gr%mesh, zpsiL, err)  
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
              if(st%have_left_states) then
                if(read_left_) then
                  call states_set_state(st, gr%mesh, idim, ist, ik, zpsiL, left = .true.)
                else
                  call states_set_state(st, gr%mesh, idim, ist, ik, zpsi, left = .true.)
                end if
              end if  
            else
              call lalg_copy(gr%mesh%np, zpsi, lr%zdl_psi(:, idim, ist, ik))
            end if
          end if


          if (err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          else if (present(lowest_missing)) then
            lowest_missing(idim, ik) = min(lowest_missing(idim, ik), ist)
          end if

          if (mpi_grp_is_root(mpi_world) .and. verbose_) then
            call loct_progress_bar(iread, nread)
            INCR(iread, 1)
          end if

        end do
      end do
    end do

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zpsiL)

    if(mpi_grp_is_root(mpi_world) .and. verbose_) then
      call messages_new_line()
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%st_kpt_mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif

    if (.not. present(lr)) call fill_random()
    ! it is better to initialize lr wfns to zero

    SAFE_DEALLOCATE_A(filled)

    if (ierr == 0) then
      ierr = -1 ! no files read
      if(.not. present(lr)) then
        write(str, '(a,i5)') 'Reading states.'
      else
        write(str, '(a,i5)') 'Reading states for linear response.'
      end if
      call messages_print_stress(stdout, trim(str))
      message(1) = 'No files could be read. No states restart information can be used.'
      call messages_info(1)
      call messages_print_stress(stdout)
    else if (ierr == st%nst * st%d%nik * st%d%dim) then
      ierr = 0
      message(1) = 'Info: States reading done.'
      if(verbose_) call print_date(trim(message(1))//' ')
    else
      if(.not. present(lr)) then
        write(str, '(a,i5)') 'Reading states.'
      else
        write(str, '(a,i5)') 'Reading states information for linear response.'
      end if
      call messages_print_stress(stdout, trim(str))
      write(message(1),'(a,i6,a,i6,a)') 'Only ', ierr,' files out of ', &
        st%nst * st%d%nik * st%d%dim, ' could be read.'
      call messages_info(1)
      call messages_print_stress(stdout)

    end if

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

            call states_generate_random(st, gr%mesh, ist, ist)
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

    ! ---------------------------------------------------------

    subroutine load_fail()
      PUSH_SUB(states_load.load_fail)

      message(1) = "Could not read orbitals from restart directory"
      message(2) = "Please run a ground-state calculation first!"
      call messages_fatal(2)

      POP_SUB(states_load.load_fail)
    end subroutine load_fail

  end subroutine states_load


  subroutine states_dump_rho(restart, st, gr, ierr, iter)
    type(restart_t),      intent(in)    :: restart
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in)    :: gr
    integer,              intent(out)   :: ierr
    integer,    optional, intent(in)    :: iter

    integer :: iunit, isp, err
    character(len=80) :: filename
    FLOAT, pointer :: rho(:)
    CMPLX, pointer :: zrho(:), zrho_fine(:)

    PUSH_SUB(states_dump_rho)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(states_dump_rho)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Writing density restart."
      call messages_info(1)
    end if

    !write the densities
    iunit = restart_open(restart, 'density')
    if (ierr == 0) then
      if (mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
        write(iunit,'(a)') '#     #spin    #nspin    filename'
        write(iunit,'(a)') '%densities'
      end if
    else
      ierr = 1
    end if

    if(gr%have_fine_mesh) then
      if(st%cmplxscl%space) then
        SAFE_ALLOCATE(zrho(1:gr%mesh%np))
        SAFE_ALLOCATE(zrho_fine(1:gr%fine%mesh%np))
      else
        SAFE_ALLOCATE(rho(1:gr%mesh%np))
      end if
    end if

    if (ierr == 0) then
      do isp = 1, st%d%nspin
        if(st%d%nspin==1) then
          write(filename, fmt='(a)') 'density'
        else
          write(filename, fmt='(a,i1)') 'density-sp', isp
        endif
        if (iunit > 0 .and. mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
          write(iunit, '(i8,a,i8,a)') isp, ' | ', st%d%nspin, ' | "'//trim(adjustl(filename))//'"'
        end if

        if(gr%have_fine_mesh)then
          if(st%cmplxscl%space) then
            zrho_fine(:) = st%zrho%Re(:,isp) + M_zI*st%zrho%Im(:,isp)
            call zmultigrid_fine2coarse(gr%fine%tt, gr%fine%der, gr%mesh, zrho_fine, zrho, INJECTION)
            call zrestart_write_function(restart, filename, gr%mesh, zrho, err)
          else
            call dmultigrid_fine2coarse(gr%fine%tt, gr%fine%der, gr%mesh, st%rho(:,isp), rho, INJECTION)
            call drestart_write_function(restart, filename, gr%mesh, rho, err)
          end if
        else
          if(st%cmplxscl%space) then
            call zrestart_write_function(restart, filename, gr%mesh, st%zrho%Re(:,isp)+M_zI*st%zrho%Im(:,isp), err)
          else
            call drestart_write_function(restart, filename, gr%mesh, st%rho(:,isp), err)
          end if
        end if
        if (err /= 0) then
          message(1) = "Unsuccessful write of '"//trim(filename)//"'."
          call messages_warning(1)
          ierr = ierr + 1
        end if

      end do
    end if
    if(gr%have_fine_mesh)then
      SAFE_DEALLOCATE_P(rho)
      if(st%cmplxscl%space) then
        SAFE_DEALLOCATE_P(zrho)
        SAFE_DEALLOCATE_P(zrho_fine)
      endif
    end if

    if (iunit > 0 .and. mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
      write(iunit,'(a)') '%'
      if(present(iter)) write(iunit,'(a,i7)') 'Iter = ', iter
    end if
    if (iunit > 0) call restart_close(restart, iunit)

    if (in_debug_mode) then
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

    integer              :: err, isp
    character(len=12)    :: filename
    FLOAT, allocatable   :: rho_coarse(:)
    CMPLX, allocatable   :: zrho(:), zrho_coarse(:)

    PUSH_SUB(states_load_rho)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load_rho)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading density restart."
      call messages_info(1)
    end if

    ! skip for now, since we know what the files are going to be called
    !read the densities
!    iunit_rho = io_open(trim(dir)//'/density', action='write', is_tmp=.true.)
!    call iopar_read(st%dom_st_kpt_mpi_grp, iunit_rho, line, err)
!    call iopar_read(st%dom_st_kpt_mpi_grp, iunit_rho, line, err)
!   we could read the iteration 'iter' too, not sure if that is useful.

    if(gr%have_fine_mesh) then
      if(st%cmplxscl%space) then
        SAFE_ALLOCATE(zrho(1:gr%fine%mesh%np))
        SAFE_ALLOCATE(zrho_coarse(1:gr%mesh%np_part))
      else
        SAFE_ALLOCATE(rho_coarse(1:gr%mesh%np_part))
      end if
    end if
    do isp = 1, st%d%nspin
      if(st%d%nspin==1) then
        write(filename, fmt='(a)') 'density'
      else
        write(filename, fmt='(a,i1)') 'density-sp', isp
      endif
!      if(mpi_grp_is_root(st%dom_st_kpt_mpi_grp)) then
!        read(iunit_rho, '(i8,a,i8,a)') isp, ' | ', st%d%nspin, ' | "'//trim(adjustl(filename))//'"'
!      end if
      if(gr%have_fine_mesh)then
        if(st%cmplxscl%space) then
          call zrestart_read_function(restart, filename, gr%mesh, zrho_coarse, err)
          call zmultigrid_coarse2fine(gr%fine%tt, gr%der, gr%fine%mesh, zrho_coarse, zrho, order = 2)
          st%zrho%Re(:,isp) =  real(zrho, REAL_PRECISION)
          st%zrho%Im(:,isp) = aimag(zrho)
        else
          call drestart_read_function(restart, filename, gr%mesh, rho_coarse, err)
          call dmultigrid_coarse2fine(gr%fine%tt, gr%der, gr%fine%mesh, rho_coarse, st%rho(:,isp), order = 2)
        end if
      else
        if(st%cmplxscl%space) then
          call zrestart_read_function(restart, filename, gr%mesh, zrho, err)
          st%zrho%Re(:,isp) =  real(zrho, REAL_PRECISION)
          st%zrho%Im(:,isp) = aimag(zrho)
        else
          call drestart_read_function(restart, filename, gr%mesh, st%rho(:,isp), err)
        end if
      end if
      if (err /= 0) then
        message(1) = "Unsuccessful read of '"//trim(filename)//"'."
        call messages_warning(1)
        ierr = ierr + 1
      endif
    end do

    if(gr%have_fine_mesh)then
      if(st%cmplxscl%space) then
         SAFE_DEALLOCATE_A(zrho_coarse)
      else
        SAFE_DEALLOCATE_A(rho_coarse)
      endif
    end if
    if(st%cmplxscl%space) then
      SAFE_DEALLOCATE_A(zrho)
    endif

    if (in_debug_mode) then
      message(1) = "Debug: Reading density restart done."
      call messages_info(1)
    end if

    POP_SUB(states_load_rho)
  end subroutine states_load_rho


  ! ---------------------------------------------------------
  !> When doing an open-boundary calculation Octopus needs the
  !! unscattered states in order to solve the Lippmann-Schwinger
  !! equation to obtain extended eigenstates.
  subroutine states_load_free_states(restart, st, gr, ierr)
    type(restart_t),      intent(inout) :: restart
    type(states_t),       intent(inout) :: st
    type(grid_t), target, intent(inout) :: gr
    integer,              intent(out)   :: ierr

    integer                    :: ik, ist, idim, counter, wfns, occs, il, err
    integer                    :: np, ip, idir
    character(len=256)         :: lines(2), fname, filename, chars
    character                  :: char
    FLOAT                      :: occ, eval, imeval, kpoint(1:MAX_DIM), w_k
    type(simul_box_t), pointer :: sb
    type(mesh_t), pointer      :: m_lead, m_center
    CMPLX                      :: phase
    CMPLX, allocatable         :: tmp(:, :)
    type(mpi_grp_t)            :: mpi_grp
    integer                    :: start(1:3), end(1:3), start_lead(1:3), end_lead(1:3)

    PUSH_SUB(states_load_free_states)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(states_load_free_states)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading free states restart."
      call messages_info(1)
    end if

    sb       => gr%sb
    m_lead   => gr%ob_grid%lead(LEFT)%mesh
    m_center => gr%mesh

    mpi_grp = st%dom_st_kpt_mpi_grp

    np = m_lead%np
    SAFE_ALLOCATE(tmp(1:np, 1:st%d%dim))

    wfns = restart_open(restart, 'wfns')
    if (wfns <  0) ierr = -1

    occs = restart_open(restart, 'occs')
    if (occs <  0) ierr = -1

    ! Skip two lines.
    if (ierr == 0) call iopar_read(mpi_grp, wfns, lines, 2, ierr)
    if (ierr == 0) call iopar_read(mpi_grp, occs, lines, 2, ierr)

    counter  = 0 ! reset counter
    st%d%kweights(:) = M_ZERO

    forall(il = 1:NLEADS) st%ob_lead(il)%rho = M_ZERO
    call mpi_grp_copy(m_lead%mpi_grp, gr%mesh%mpi_grp)

    do
      ! Check for end of file. Check only one of the two files assuming
      ! they are written correctly, i.e. of same length.
      if (ierr == 0) call iopar_read(mpi_grp, wfns, lines, 1, ierr)
      if (ierr /= 0) exit
      read(lines(1), '(a)') char
      if(char  ==  '%') then
        exit
      end if
      read(lines(1), *) ik, char, ist, char, idim, char, fname

      call iopar_read(mpi_grp, occs, lines, 1, ierr)
      if (ierr /= 0) exit
      !# occupations | eigenvalue[a.u.] | k-points | k-weights | filename | ik | ist | idim
      read(lines(1), *) occ, char, eval, char, imeval, char, (kpoint(idir), char, idir = 1, gr%sb%dim), &
        w_k, char, chars, char, ik, char, ist, char, idim

      ! we need the kpoints from the periodic run for the scattering states
      ! so overwrite the "false" kpoints of the finite center
      call kpoints_set_point(gr%sb%kpoints, ik, kpoint(1:gr%sb%dim))
      st%d%kweights(ik) = w_k
      st%occ(ist, ik) = occ
      counter = counter + 1
      ! if not in the corresponding node cycle
      if(ik < st%d%kpt%start .or. ik > st%d%kpt%end .or. ist < st%st_start .or. ist > st%st_end) cycle

      call zrestart_read_function(restart, fname, m_lead, tmp(:, idim), ierr)
      if (ierr /= 0) exit

      call lead_dens_accum()


      ! Copy the wave-function from the lead to the central region, using periodicity.
      start(1:3) = m_center%idx%nr(1, 1:3) + m_center%idx%enlarge(1:3)
      end(1:3) = m_center%idx%nr(2, 1:3) - m_center%idx%enlarge(1:3)

      start_lead(1:3) = m_lead%idx%nr(1, 1:3) + m_lead%idx%enlarge(1:3)
      end_lead(1:3) = m_lead%idx%nr(2, 1:3) - m_lead%idx%enlarge(1:3)

      st%zphi(:, idim, ist, ik)  = M_z0
      call zmf_add(m_lead, start_lead, end_lead, tmp(:, idim), m_center, start, end, &
                    st%zphi(:, idim, ist, ik), TRANS_DIR)

      ! Apply phase. Here the kpoint read from the file is different
      ! from the one calculated by octopus, since the box size changed.
      ! The read kpoints are used for describing the free states of the
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
          call zio_function_output(C_OUTPUT_HOW_AXIS_X, 'debug/open_boundaries', filename, &
            m_center, st%zphi(:, idim, ist, ik), sqrt(units_out%length**(-sb%dim)), err, is_tmp=.false.)
        case(2, 3)
          call zio_function_output(C_OUTPUT_HOW_PLANE_Z, 'debug/open_boundaries', filename, &
            m_center, st%zphi(:, idim, ist, ik), sqrt(units_out%length**(-sb%dim)), err, is_tmp=.false.)
        end select
        if (err /= 0) then
          !Failure to write debug information should not affect rest of the routine, so we just print a warning.
          message(1) = "Unable to write debug information to '"//"debug/open_boundaries/"//trim(filename)//"'."
          call messages_warning(1)
        end if
      end if

    end do ! Loop over all free states.

    if (wfns > 0) call restart_close(restart, wfns)
    if (occs > 0) call restart_close(restart, occs)
    SAFE_DEALLOCATE_A(tmp)

    if (ierr == 0) then
      message(1) = "Info: Successfully initialized free states"
      write(message(2),'(a,i3,a)') 'Info:', counter, ' wave functions read by program.'
      call messages_info(2)
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading free states restart done."
      call messages_info(1)
    end if

    POP_SUB(states_load_free_states)
  contains

    subroutine lead_dens_accum()
      integer :: il

      PUSH_SUB(states_load_free_states.lead_dens_accum)
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
!          if(k_idim == 2) then
!            do ip = 1, np
!              c = w_k*occ*tmp(ip, 1)*conjg(tmp(ip, 2))
!              st%ob_lead(il)%rho(ip, 3) = st%ob_lead(il)%rho(ip, 3) + real(c, REAL_PRECISION)
!              st%ob_lead(il)%rho(ip, 4) = st%ob_lead(il)%rho(ip, 4l) + aimag(c)
!            end do
!          end if
        end select
      end do

      POP_SUB(states_load_free_states.lead_dens_accum)
    end subroutine lead_dens_accum
  end subroutine states_load_free_states


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

    integer, parameter ::      &
      state_from_formula  = 1, &
      state_from_file     = 0, &
      normalize_yes       = 1, &
      normalize_no        = 0

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
    !%Option file 0
    !% Read initial orbital from file. 
    !% Accepted file formats: obf, ncdf and csv.
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
                  zpsi(ip, 1) = psi_re + M_zI * psi_im
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
                call zio_function_input(filename, mesh, zpsi(:, 1), ierr, .true.)
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
              case(normalize_no)
              case(normalize_yes)
                call states_get_state(st, mesh, is, ik, zpsi)
                call zstates_normalize_orbital(mesh, st%d%dim, zpsi)
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

end module states_restart_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
