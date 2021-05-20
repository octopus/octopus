!! Copyright (C) 2019 R. Jestaedt, F. Bonafe, H. Appel, A. Rubio
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

module states_mxll_restart_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use global_oct_m
  use grid_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use states_elec_restart_oct_m
  use states_mxll_oct_m
  use states_elec_dim_oct_m
  use string_oct_m
  use types_oct_m
  use unit_system_oct_m
  use unit_oct_m
    
  implicit none

  type(profile_t), save :: prof_read, prof_write

  private
  public ::                      &
    states_mxll_read_user_def,   &
    states_mxll_dump,            &
    states_mxll_load

contains

  !----------------------------------------------------------
  subroutine states_mxll_read_user_def(namespace, space, mesh, st, user_def_rs_state)
    type(namespace_t),   intent(in)    :: namespace
    type(space_t),       intent(in)    :: space
    type(mesh_t),        intent(inout) :: mesh
    type(states_mxll_t), intent(inout) :: st
    CMPLX,               intent(inout) :: user_def_rs_state(:,:)

    type(block_t)      :: blk
    integer            :: il, nlines, idim, ncols, ip, state_from, ierr, maxwell_field
    FLOAT              :: xx(space%dim), rr, e_value, dummy, b_value
    FLOAT, allocatable :: e_field(:), b_field(:)
    FLOAT, allocatable :: total_efield(:,:), total_bfield(:,:)
    CMPLX, allocatable :: rs_state_add(:), rs_state(:,:)
    character(len=150), pointer :: filename_e_field, filename_b_field
    character(1) :: cdim
    type(profile_t), save :: prof
    
    integer, parameter ::           &
      STATE_FROM_FORMULA  = 1,      &
      STATE_FROM_FILE     = -10010

    PUSH_SUB(states_mxll_read_user_def)

    call profiling_in(prof, 'STATES_MXLL_READ_USER_DEF')

    !%Variable UserDefinedInitialMaxwellStates
    !%Type block
    !%Section MaxwellStates
    !%Description
    !% The initial electromagnetic fields can be set by the user 
    !% with the <tt>UserDefinedMaxwellStates</tt> block variable.
    !% The electromagnetic fields have to fulfill the 
    !% Maxwells equations in vacuum.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellStates
    !% <br>&nbsp;&nbsp; 2 | formula | "magnetic_field" | "-1/P_c * sin(x)"
    !% <br>&nbsp;&nbsp; 3 | formula | "electric_field" | "   sin(x)      "
    !% <br>%</tt>
    !%
    !% The first column specifies the component of the dimension of
    !% the electric field and magnetic field. Column four
    !% indicates that column five and six should be interpreted 
    !% as a formula for the corresponding orbital. P_c is the 
    !% speed of light constant.
    !%
    !% Alternatively, if column four states <tt>file</tt> the 
    !% electric field and magnetic field will be read from 
    !% the files given in column five and six.
    !%
    !% <tt>%MUserDefinedMaxwellStates
    !% <br>&nbsp;&nbsp; 3 | file | electric_field | "/path/to/file_electric_field_of_dimension_3"
    !% <br>&nbsp;&nbsp; 2 | file | magnetic_field | "/path/to/file_magnetic_field_of_dimension_2"
    !% <br>%</tt>
    !%
    !%Option file -10010
    !% Read initial orbital from file. 
    !% Accepted file formats: obf, ncdf and csv.
    !%Option formula 1
    !% Calculate initial orbital by given analytic expression.
    !%Option electric_field 1
    !% This row defines the electric field component of the corresponding dimension
    !%Option magnetic_field 2
    !% This row defines the magnetic field component of the corresponding dimension
    !%End

    if(parse_block(namespace, 'UserDefinedInitialMaxwellStates', blk) == 0) then

      SAFE_ALLOCATE(rs_state_add(1:mesh%np_part))
      SAFE_ALLOCATE(rs_state(1:mesh%np, 1:3))

      ! Set electromagnetic field equal to zero in the whole simulation box.
      user_def_rs_state(:,:) = M_ZERO

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      write(message(1), '(a,i5)') 'Maxwell electromagnetic fields are added.'
      write(message(2), '(a,i5)') ''
      call messages_info(2)

    
      SAFE_ALLOCATE(total_efield(1:mesh%np, 1:3))
      SAFE_ALLOCATE(total_bfield(1:mesh%np, 1:3))
      total_efield = M_ZERO
      total_bfield = M_ZERO

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is five or six.
        ncols = parse_block_cols(blk, il - 1)
        if(ncols  <  4 .or. ncols > 4) then
          message(1) = 'Each line in the UserDefinedMaxwellStates block must have'
          message(2) = 'four columns.'
          call messages_fatal(2, namespace=namespace)
        end if

        call parse_block_integer(blk, il - 1, 0, idim)
        write(cdim,'(I1)')idim

        ! Calculate from expression or read from file?
        call parse_block_integer(blk, il - 1, 1, state_from)
    
        rs_state_add(:) = M_ZERO

        select case(state_from)

        case(STATE_FROM_FORMULA)

          ! parse formula string
          call parse_block_integer(blk, il - 1, 2, maxwell_field )
          if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__ELECTRIC_FIELD) then
            call parse_block_string( blk, il - 1, 3, st%user_def_e_field(idim))
            call messages_write("  E-field in dimension "//trim(cdim)//" : "//trim(st%user_def_e_field(idim)), fmt='(a,i1,2a)')
            call conv_to_C_string(st%user_def_e_field(idim))
          else if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__MAGNETIC_FIELD) then
            call parse_block_string( blk, il - 1, 3, st%user_def_b_field(idim))
            call messages_write("  B-field in dimension "//trim(cdim)//" : "//trim(st%user_def_b_field(idim)), fmt='(a,i1,2a)')
            call conv_to_C_string(st%user_def_b_field(idim))
          end if

          if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__ELECTRIC_FIELD) then

            ! fill Maxwell states with user-defined formulas
            do ip = 1, mesh%np
              xx = mesh%x(ip, :)
              rr = sqrt(sum(xx**2))
              ! parse user-defined expressions
              call parse_expression(e_value, dummy, st%dim, xx, rr, M_ZERO, &
                                    st%user_def_e_field(idim))
              total_efield(ip, idim) = total_efield(ip, idim) + e_value
            end do

          else if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__MAGNETIC_FIELD) then
            ! fill Maxwell states with user-defined formulas
            do ip = 1, mesh%np
              xx = mesh%x(ip, :)
              rr = sqrt(sum(xx**2))

              call parse_expression(b_value, dummy, st%dim, xx, rr, M_ZERO, &
                                      st%user_def_b_field(idim))
              total_bfield(ip, idim) = total_bfield(ip, idim) + b_value
            end do

          end if

        case(STATE_FROM_FILE)
          ! The input format can be coded in column four now. As it is
          ! not used now, we just say "file".
          ! Read the filename.
          call parse_block_integer(blk, il - 1, 2, maxwell_field )
          SAFE_ALLOCATE(e_field(1:mesh%np))
          SAFE_ALLOCATE(b_field(1:mesh%np))
          e_field = M_ZERO
          b_field = M_ZERO
          if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__ELECTRIC_FIELD) then
            call parse_block_string(blk, il - 1, 3, filename_e_field)
            call messages_write("  E-field in dimension "//trim(cdim)//" : "//trim(filename_e_field), fmt='(a,i1,2a)')
            call dio_function_input(filename_e_field, namespace, space, mesh, e_field(:), ierr)
            if (ierr > 0) then
              message(1) = 'Could not read the file!'
              write(message(2),'(a,i1)') 'Error code: ', ierr
              call messages_fatal(2, namespace=namespace)
            end if
            e_field = units_to_atomic(units_inp%energy/units_inp%length, e_field)
          else if (maxwell_field == OPTION__USERDEFINEDINITIALMAXWELLSTATES__MAGNETIC_FIELD) then
            call parse_block_string(blk, il - 1, 3, filename_b_field)
            call messages_write("  B-field in dimension "//trim(cdim)//" : "//trim(filename_b_field), fmt='(a,i1,2a)')
            call dio_function_input(filename_b_field, namespace, space, mesh, b_field(:), ierr)
            if (ierr > 0) then
              message(1) = 'Could not read the file!'
              write(message(2),'(a,i1)') 'Error code: ', ierr
              call messages_fatal(2, namespace=namespace)
            end if
            b_field = units_to_atomic(unit_one/units_inp%length**2, b_field)
          end if
          ! fill state
          call build_rs_vector(e_field(:), b_field(:), st%rs_sign, rs_state_add(:), &
                                       st%ep(ip), st%mu(ip))

          SAFE_DEALLOCATE_A(e_field)
          SAFE_DEALLOCATE_A(b_field)
        
          call lalg_axpy(mesh%np, M_ONE, rs_state_add, user_def_rs_state(:,idim))

        case default
          message(1) = 'Wrong entry in UserDefinedMaxwellStates, column 2.'
          message(2) = 'You may state "formula" or "file" here.'
          call messages_fatal(2, namespace=namespace)
        end select

      end do

      if(state_from == STATE_FROM_FORMULA) then
        do idim = 1, 3
          total_efield(:, idim) = units_to_atomic(units_inp%energy/units_inp%length, total_efield(:, idim))
          total_bfield(:, idim) = units_to_atomic(unit_one/(units_inp%length**2), total_bfield(:, idim)) 
        end do

       ! fill state
       call build_rs_state(total_efield, total_bfield, st%rs_sign, rs_state, mesh, st%ep, st%mu)
       call lalg_axpy(mesh%np, 3, M_ONE, rs_state, user_def_rs_state)
      end if

      
      SAFE_DEALLOCATE_A(total_efield)
      SAFE_DEALLOCATE_A(total_bfield)

      SAFE_DEALLOCATE_A(rs_state)
      SAFE_DEALLOCATE_A(rs_state_add)
      call parse_block_end(blk)
      !call messages_print_stress(stdout, namespace=namespace)

    else
      call messages_variable_is_block(namespace, 'UserDefineInitialdStates')
    end if

    call profiling_out(prof)

    POP_SUB(states_mxll_read_user_def)
  end subroutine states_mxll_read_user_def


  !----------------------------------------------------------
  subroutine states_mxll_dump(restart, st, space, mesh, zff, zff_dim, ierr, iter, st_start_writing, verbose)
    type(restart_t),      intent(in)  :: restart
    type(states_mxll_t),  intent(in)  :: st
    type(space_t),        intent(in)  :: space
    type(mesh_t),         intent(in)  :: mesh
    CMPLX,                intent(in)  :: zff(:,:)
    integer,              intent(in)  :: zff_dim
    integer,              intent(out) :: ierr
    integer,    optional, intent(in)  :: iter
    integer,    optional, intent(in)  :: st_start_writing
    logical,    optional, intent(in)  :: verbose

    integer :: iunit_wfns, iunit_states
    integer :: err, err2(2), ist, idim, itot
    integer :: root(1:P_STRATEGY_MAX)
    character(len=MAX_PATH_LEN) :: filename
    character(len=300) :: lines(3)
    logical :: should_write, verbose_

    PUSH_SUB(states_mxll_dump)

    verbose_ = optional_default(verbose, .true.)

    ierr = 0

    if(restart_skip(restart)) then
      POP_SUB(states_mxll_dump)
      return
    end if

    if(verbose_) then
      message(1) = "Info: Writing Maxwell states."
      call print_date(trim(message(1))//' ')
    end if

    call profiling_in(prof_write, "MAXWELL_RESTART_WRITE")

    call restart_block_signals()

    iunit_states = restart_open(restart, 'maxwell_states')
    write(lines(2), '(a20,1i10)')  'dim=                ', zff_dim
    call restart_write(restart, iunit_states, lines, 3, err)
    if (err /= 0) ierr = ierr + 1
    call restart_close(restart, iunit_states)

    iunit_wfns = restart_open(restart, 'wfns')
    lines(1) = '#     #dim    filename'
    lines(2) = '%RS States'
    call restart_write(restart, iunit_wfns, lines, 2, err)
    if (err /= 0) ierr = ierr + 2


    itot = 1
    root = 0
    err2 = 0
    ist = 1
    
    do idim = 1, zff_dim
       itot = itot + 1

       root(P_STRATEGY_DOMAINS) = mod(itot - 1, mesh%mpi_grp%size)
       write(filename,'(i10.10)') itot

       write(lines(1), '(i8,3a)') idim, ' | "', trim(filename), '"'
       call restart_write(restart, iunit_wfns, lines, 1, err)
       if (err /= 0) err2(1) = err2(1) + 1

       should_write = st%st_start <= ist .and. ist <= st%st_end
       if (should_write .and. present(st_start_writing)) then
          if (ist < st_start_writing) should_write = .false.
       end if

       if (should_write) then
          call zrestart_write_mesh_function(restart, space, filename, mesh, zff(:,idim), err, root = root)
          if (err /= 0) err2(2) = err2(2) + 1
       end if
    end do ! zff_dim

    if (err2(1) /= 0) ierr = ierr + 8
    if (err2(2) /= 0) ierr = ierr + 16

    lines(1) = '%'
    call restart_write(restart, iunit_wfns, lines, 1, err) 
    if (err /= 0) ierr = ierr + 64
    if (present(iter)) then
      write(lines(1),'(a,i7)') 'Iter = ', iter
      call restart_write(restart, iunit_wfns, lines, 1, err)
      if (err /= 0) ierr = ierr + 128
    end if

    call restart_close(restart, iunit_wfns)

    if(verbose_) then
      message(1) = "Info: Finished writing Maxwell states."
      call print_date(trim(message(1))//' ')
    end if

    call restart_unblock_signals()

    call profiling_out(prof_write)
    POP_SUB(states_mxll_dump)
    return

  end subroutine states_mxll_dump

  !----------------------------------------------------------
  subroutine states_mxll_load(restart, st, mesh, namespace, space, zff, zff_dim, ierr, iter, lowest_missing, label, verbose)
    type(restart_t),            intent(in)    :: restart
    type(states_mxll_t),        intent(inout) :: st
    type(mesh_t),               intent(in)    :: mesh
    type(namespace_t),          intent(in)    :: namespace
    type(space_t),              intent(in)    :: space
    CMPLX,                      intent(inout) :: zff(:,:)
    integer,                    intent(in)    :: zff_dim
    integer,                    intent(out)   :: ierr
    integer,          optional, intent(out)   :: iter
    integer,          optional, intent(out)   :: lowest_missing(:) !< all states below this one were read successfully
    character(len=*), optional, intent(in)    :: label
    logical,          optional, intent(in)    :: verbose

    integer              :: states_file, wfns_file, err, ist, idim, dim, mx_st_start, mx_st_end
    integer              :: idone, iread, ntodo
    character(len=12)    :: filename
    character(len=1)     :: char
    logical, allocatable :: filled(:, :)
    character(len=256)   :: lines(3), label_
    character(len=50)    :: str

    logical              :: verbose_
    character(len=256), allocatable :: restart_file(:, :)
    logical,            allocatable :: restart_file_present(:, :)

    PUSH_SUB(states_mxll_load)

    ierr = 0
    dim  = zff_dim

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
    end if

    message(1) = 'Info: Reading Maxwell states'
    if (len(trim(label_)) > 0) then
      message(1) = trim(message(1)) // trim(label_)
    end if
    message(1) = trim(message(1)) // "."
    if(verbose_) call print_date(trim(message(1))//' ')

    states_file  = restart_open(restart, 'maxwell_states')
    call restart_read(restart, states_file, lines, 3, err)
    if (err /= 0) then
      ierr = ierr - 2
    else
      read(lines(2), *) idim
    end if
    call restart_close(restart, states_file)

    ! open files to read
    wfns_file  = restart_open(restart, 'wfns')
    call restart_read(restart, wfns_file, lines, 2, err)
    if (err /= 0) then
      ierr = ierr - 2**5
    end if

    ! If any error occured up to this point then it is not worth continuing,
    ! as there something fundamentally wrong with the restart files
    if (ierr /= 0) then
      call restart_close(restart, wfns_file)
      call profiling_out(prof_read)
      POP_SUB(states_mxll_load)
      return
    end if

    SAFE_ALLOCATE(restart_file(1:zff_dim, st%st_start:st%st_end))
    SAFE_ALLOCATE(restart_file_present(1:zff_dim,st%st_start:st%st_end))
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
          read(lines(1), *) idim, char, filename
        end if
      end if

      if (ist >= st%st_start .and. ist <= st%st_end) then
        restart_file(idim, ist) = trim(filename)
        restart_file_present(idim, ist) = .true.
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

    ! Now we read the wavefunctions. At this point we need to have all the information from the
    ! states, and wfns files in order to avoid serialisation of reads, as restart_read
    ! works as a barrier.

    mx_st_start=st%st_start
    mx_st_end=st%st_end
    SAFE_ALLOCATE(filled(1:zff_dim,mx_st_start:mx_st_end))
    filled = .false.

    if (present(lowest_missing)) lowest_missing = st%nst + 1

    iread = 0
    if (mpi_grp_is_root(mpi_world) .and. verbose_) then
      idone = 0
      ntodo = st%lnst*zff_dim
      call loct_progress_bar(-1, ntodo)
    end if

    ist = 1
    do idim = 1, zff_dim
      if (.not. restart_file_present(idim, ist)) then
        if (present(lowest_missing)) &
          lowest_missing(idim) = min(lowest_missing(idim), ist)
        cycle
      endif

      call zrestart_read_mesh_function(restart, space, restart_file(idim, ist), mesh, &
        zff(:,idim), err)

      if (err == 0) then
        filled(idim, ist) = .true.
        iread = iread + 1
      else if (present(lowest_missing)) then
        lowest_missing(idim) = min(lowest_missing(idim), ist)
      end if

      if (mpi_grp_is_root(mpi_world) .and. verbose_) then
        idone = idone + 1
        call loct_progress_bar(idone, ntodo)
      end if

    end do
    
    SAFE_DEALLOCATE_A(restart_file)
    SAFE_DEALLOCATE_A(restart_file_present)
    SAFE_DEALLOCATE_A(filled)

    if(mpi_grp_is_root(mpi_world) .and. verbose_) then
      call messages_new_line()
    end if

    if (ierr == 0 .and. iread /= st%nst * zff_dim) then
      if(iread > 0) then
        ierr = iread
      else
        ierr = -1
      endif
      ! otherwise ierr = 0 would mean either all was read correctly, or nothing at all was read!

      write(str, '(a,i5)') 'Reading Maxwell states.'
      call messages_print_stress(stdout, trim(str), namespace=namespace)
      write(message(1),'(a,i6,a,i6,a)') 'Only ', iread,' files out of ', &
           st%nst * zff_dim, ' could be read.'
      call messages_info(1)
      call messages_print_stress(stdout, namespace=namespace)
    end if

    message(1) = 'Info: Maxwell states reading done.'
    if(verbose_) call print_date(trim(message(1))//' ')

    call profiling_out(prof_read)
    POP_SUB(states_mxll_load)

  end subroutine states_mxll_load

end module states_mxll_restart_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
