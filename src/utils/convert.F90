!! Copyright (C) 2013 J. Alberdi-Rodriguez, J. Jornet-Somoza
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

program oct_convert
  use calc_mode_m
  use command_line_m
  use datasets_m
  use fft_m
  use fftw_m
  use geometry_m
  use global_m
  use io_m
  use io_function_m
  use io_binary_m
  use loct_m
  use messages_m
  use mesh_m
  use mpi_m
  use output_m
  use parser_m
  use poisson_m
  use profiling_m
  use string_m
  use system_m
  use restart_m
  use unit_m
  use unit_system_m
  use utils_m

  implicit none

  character*256 :: config_str
  integer :: ierr
  
  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call global_init()
  call calc_mode_init()
  call messages_init()

  call datasets_init(1)
  call io_init()
  call profiling_init()
  call messages_experimental("oct-convert utility")

  call print_header()

  call messages_print_stress(stdout, "Convert mode")
  call messages_print_stress(stdout)

  call fft_all_init()
  call unit_system_init()

  call convert()

  call fft_all_end()
  call profiling_output()
  call profiling_end()
  call io_end()
  call print_date("Calculation ended on ")
  call datasets_end()
  call messages_end()
  call global_end()

contains

  ! -------------
  !> Reads an binary file and writes the equivalent files, 
  !! defined with OutputHow.
  !! This is a high-level interface that reads the input file and
  !! calls the proper function.
  subroutine convert()
    type(system_t) :: sys

    character(64)  :: basename, folder, ref_name, ref_folder, folder_default
    integer        :: c_start, c_end, c_step, c_start_default, length
    logical        :: iterate_folder, subtract_file, fourier_trans

    PUSH_SUB(convert)

    call calc_mode_set_parallelization(P_STRATEGY_STATES, default = .false.)
    call system_init(sys)

    message(1) = 'Info: Converting files'
    message(2) = ''
    call messages_info(2)

    !%Variable ConvertFilename
    !%Type string
    !%Default density
    !%Section Utilities::oct-convert
    !%Description
    !% Input filename. The original filename which is going to be converted in the format
    !% specified in <tt>OutputHow</tt>. It is going to convert various files, it should 
    !% only contain the beginning of the name. For instance, in the case of the restart 
    !% files it should be one space ' '.
    !%End
    call parse_string(datasets_check('ConvertFilename'), 'density', basename)
    if ( basename == " " ) basename = ""
    ! Delete the extension if present
    length = len_trim(basename)
    if ( length > 4) then
      if ( basename(length-3:length) == '.obf' ) then
        basename = trim(basename(1:length-4))
      end if
    end if


    !%Variable ConvertIterateFolder
    !%Type logical
    !%Default true
    !%Section Utilities::oct-convert
    !%Description
    !% This variable decides if a folder is going to be iterated or the 
    !% filename is going to be iterated.
    !%End
    call parse_logical(datasets_check('ConvertIterateFolder'), .true., iterate_folder)

    if (iterate_folder) then
      folder_default  = 'td.'
      c_start_default = 0
    else
      folder_default  = 'restart'
      c_start_default = 1
    end if
    
    !%Variable ConvertFolder
    !%Type string
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name where the input files are.
    !%End
    call parse_string(datasets_check('ConvertFolder'), folder_default, folder)
    call add_last_slash(folder)

    !%Variable ConvertStart
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-convert
    !%Description
    !% The starting number of the filename or folder.
    !%End
    call parse_integer(datasets_check('ConvertStart'), c_start_default, c_start)

    !%Variable ConvertEnd
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The last number of the filename or folder.
    !%End
    call parse_integer(datasets_check('ConvertEnd'), 1, c_end)

    !%Variable ConvertStep
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The padding between the filenames or folder.
    !%End
    call parse_integer(datasets_check('ConvertStep'), 1, c_step)

    !%Variable ConvertSubtractFilename
    !%Type string
    !%Default density
    !%Section Utilities::oct-convert
    !%Description
    !% Input filename. The file which is going to subtracted to rest of the files.
    !% specified in <tt>OutputHow</tt>.
    !%End
    call parse_string(datasets_check('ConvertSubtractFilename'), 'density', ref_name)
    if ( ref_name == " " ) ref_name = ""
    ! Delete the extension if present
    length = len_trim(ref_name)
    if ( length > 4) then
      if ( ref_name(length-3:length) == '.obf' ) then
        ref_name = trim(ref_name(1:length-4))
      end if
    end if
    
    !%Variable ConvertSubtract
    !%Type logical
    !%Default false
    !%Section Utilities::oct-convert
    !%Description
    !% Decides if a reference file is going to be subtracted.
    !%End
    call parse_logical(datasets_check('ConvertSubtract'), .false., subtract_file)

    !%Variable ConvertSubtractFolder
    !%Type string
    !%Default [blank]
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name which is going to be subtracted.
    !%End
    call parse_string(datasets_check('ConvertSubtractFolder'), '.', ref_folder)
    call add_last_slash(folder)
    
    !%Variable ConvertTransform
    !%Type logical
    !%Default false
    !%Section Utilities::oct-convert
    !%Description
    !% Decides if the input files are going to be transformed via Fourier Transform.
    !%End
    call parse_logical(datasets_check('ConvertTransform'), .false., fourier_trans)

    ! Compute Fourier transform 
    if (fourier_trans) then
      call convert_transform(sys%gr%mesh, sys%geo, basename, folder, &
         c_start, c_end, c_step, sys%outp, subtract_file, &
         ref_name, ref_folder)
    else
      call convert_low(sys%gr%mesh, sys%geo, basename, folder, &
         c_start, c_end, c_step, sys%outp, iterate_folder, &
         subtract_file, ref_name, ref_folder)
    end if

    call system_end(sys)

    POP_SUB(convert)
  end subroutine convert

  ! ---------------------------------------------------------
  !> Giving a range of input files, it writes the corresponding 
  !! output files
  subroutine convert_low(mesh, geo, basename, in_folder, c_start, c_end, c_step, outp, iterate_folder, & 
                                 subtract_file, ref_name, ref_folder)
    type(mesh_t)    , intent(in)    :: mesh
    type(geometry_t), intent(in)    :: geo
    character(len=*), intent(inout) :: basename       !< File name
    character(len=*), intent(in)    :: in_folder      !< Folder name
    integer,          intent(in)    :: c_start        !< The first file number
    integer,          intent(in)    :: c_end          !< The last file number
    integer,          intent(in)    :: c_step         !< The step between files
    type(output_t),   intent(in)    :: outp           !< Output objetct; Decides the kind, what and where to output
    logical,          intent(in)    :: iterate_folder !< If true, it iterates over the folders, keeping the filename fixed.
                                                      !! If false, it iterates over the filenames
    logical,          intent(in)    :: subtract_file  !< If true, it subtracts the density from the reference 
    character(len=*), intent(inout) :: ref_name       !< Reference file name 
    character(len=*), intent(inout) :: ref_folder     !< Reference folder name

    type(restart_t)    :: restart
    integer            :: ierr, ii
    character(64)      :: filename, out_name, folder, frmt
    FLOAT, allocatable :: read_ff(:), read_rff(:), pot(:)

    PUSH_SUB(convert_low)

    SAFE_ALLOCATE(read_ff(1:mesh%np))
    SAFE_ALLOCATE(read_rff(1:mesh%np))
    SAFE_ALLOCATE(pot(1:mesh%np))
    read_rff(:) = M_ZERO
   
    write(message(1),'(5a,i5,a,i5,a,i5)') "Converting '", trim(in_folder), "/", trim(basename), &
         "' from ", c_start, " to ", c_end, " every ", c_step
    call messages_info(1)
 
    if (subtract_file) then
      write(message(1),'(a,a,a,a)') "Reading ref-file from ", trim(ref_folder), trim(ref_name),".obf"
      call restart_init(restart, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mesh%mpi_grp, &
                      ierr, dir=trim(ref_folder))
      ! FIXME: why only real functions? Please generalize.
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, trim(ref_name), mesh, read_rff, ierr)
        call restart_end(restart)
      else
        write(message(1),'(2a)') "Failed to read from ref-file ", trim(ref_name)
        write(message(2), '(2a)') "from folder ", trim(ref_folder)
        call messages_fatal(2)
      endif
    end if

    call loct_progress_bar(-1, c_end-c_start)
    do ii = c_start, c_end, c_step
      if (iterate_folder) then
        ! Delete the last / and add the corresponding folder number
        write(folder,'(a,i0.7,a)') in_folder(1:len_trim(in_folder)-1),ii,"/"
        write(filename, '(a,a,a)') trim(folder), trim(basename), ".obf"
        out_name = trim(basename)
      else
        folder = in_folder
        if ( c_start /= c_end ) then
          ! Here, we are only considering 10 character long filenames.
          ! Subtract the initial part given at 'ConvertFilename' from the format and pad
          ! with zeros.
          write(frmt,'(a,i0,a)')"(a,i0.",10-len_trim(basename),")"
          write(filename, fmt=trim(frmt)) trim(basename), ii
          write(out_name, '(a)') trim(filename)
          write(filename, '(a,a,a,a)') trim(folder),"/", trim(out_name),".obf"
        else 
          ! Assuming filename is given complete in the 'ConvertFilename'
          write(filename, '(a,a,a,a)') trim(folder),"/", trim(basename),".obf"
          write(out_name, '(a)') trim(basename)
        end if
      end if

      ! Read the obf file
      call restart_init(restart, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mesh%mpi_grp, &
                      ierr, dir=trim(folder))
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, trim(out_name), mesh, read_ff, ierr)
        call restart_end(restart)
      endif

      if (ierr /= 0) then
        write(message(1), '(a,a)') "Error reading the file ", filename
        write(message(2), '(a,i4)') "Error code: ",ierr
        write(message(3), '(a)') "Skipping...."
        call messages_warning(3)
        cycle
      end if
      if (subtract_file) then
        read_ff(:) = read_ff(:) - read_rff(:) 
        write(out_name, '(a,a)') trim(out_name),"-ref"
      end if
      ! Write the corresponding output
      call dio_function_output(outp%how, &
        trim(folder), trim(out_name), mesh, read_ff, units_out%length**(-mesh%sb%dim), ierr, geo = geo)
      
      if (iand(outp%what, C_OUTPUT_POTENTIAL) /= 0) then
        write(out_name, '(a)') "potential"
        call dpoisson_solve(psolver, pot, read_ff)
        call dio_function_output(outp%how, &
             trim(folder), trim(out_name), mesh, pot, units_out%energy, ierr, geo = geo)
      end if
      call loct_progress_bar(ii-c_start, c_end-c_start) 
    end do
    call restart_end(restart)
    
    SAFE_DEALLOCATE_A(read_ff)
    SAFE_DEALLOCATE_A(read_rff)
    SAFE_DEALLOCATE_A(pot)
    POP_SUB(convert_low)
  end subroutine convert_low

  ! ---------------------------------------------------------
  !> Giving a range of input files, it computes the Fourier transform
  !! of the file.
  subroutine convert_transform(mesh, geo, basename, in_folder, c_start, c_end, c_step, outp, & 
       subtract_file, ref_name, ref_folder)
    type(mesh_t)    , intent(in)    :: mesh
    type(geometry_t), intent(in)    :: geo
    character(len=*), intent(inout) :: basename       !< File name
    character(len=*), intent(in)    :: in_folder      !< Folder name
    integer,          intent(in)    :: c_start        !< The first file number
    integer,          intent(in)    :: c_end          !< The last file number
    integer,          intent(in)    :: c_step         !< The step between files
    type(output_t),   intent(in)    :: outp           !< Output object; Decides the kind, what and where to output
    logical,          intent(in)    :: subtract_file  !< If true, it subtracts the density from the reference 
    character(len=*), intent(inout) :: ref_name       !< Reference file name 
    character(len=*), intent(inout) :: ref_folder     !< Reference folder name

    integer             :: ierr, ii, i_space, i_time, nn(1:3), optimize_parity(1:3)
    integer             :: i_energy, e_end, e_start, e_step, e_point, no_e
    logical             :: optimize(1:3)
    character(64)       :: filename, ref_filename, folder
    FLOAT               :: fdefault, w_max
    FLOAT, allocatable  :: read_ft(:), read_rff(:), read_point(:), write_point(:,:)
    CMPLX, allocatable  :: out_fft(:)
    type(fft_t)         :: fft

    FLOAT   :: start_time          !< start time for the transform
    integer :: time_steps          !< number of time steps
    FLOAT   :: dt                  !< step in time mesh
    FLOAT   :: dw                  !< step in energy mesh
    FLOAT   :: max_energy          !< maximum of energy mesh
    FLOAT   :: min_energy          !< minimum of energy mesh

    PUSH_SUB(convert_transform)

    ! set default time_step as dt from TD section
    fdefault = M_ZERO
    call parse_float(datasets_check('TDTimeStep'), fdefault, dt, unit = units_inp%time)
    if (dt <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
      write(message(2),'(a)') 'Input: TDTimeStep reset to 0. Check input file'
      call messages_info(2)
    end if

    e_point = (c_end - c_start) / c_step + 1
    nn(1) = e_point
    nn(2) = 1
    nn(3) = 1
    start_time = M_ZERO
    SAFE_ALLOCATE(read_ft(1:e_point+1))
    SAFE_ALLOCATE(out_fft(1:e_point+1))
    SAFE_ALLOCATE(read_point(1:1))
    SAFE_ALLOCATE(read_rff(1:mesh%np))

    !%Variable ConvertEnergyMin
    !%Type float
    !%Default 0.
    !%Section Utilities::oct-convert
    !%Description
    !% The starting number of the filename.
    !%End
    call parse_float(datasets_check('ConvertEnergyMin'), M_ZERO, min_energy, units_inp%energy)

    ! Calculate the limits in frequency space.
    start_time = c_start * dt
    dt = dt * c_step
    time_steps = (c_end - c_start)/c_step 
    w_max = M_TWO * M_PI / dt 

    !%Variable ConvertEnergyMax
    !%Type float
    !%Default w_max
    !%Section Utilities::oct-convert
    !%Description
    !% The last number of the filename.
    !%End
    fdefault = units_from_atomic(units_inp%energy, w_max)
    call parse_float(datasets_check('ConvertEnergyMax'),fdefault, max_energy, units_inp%energy)
    if( (max_energy > w_max)) then
      write(message(1),'(a,f12.7)')'Impossible to set ConvertEnergyMax to ', &
           units_from_atomic(units_inp%energy, max_energy)
      write(message(2),'(a)')'ConvertEnergyMax is too large.'
      write(message(3),'(a,f12.7,a)')'ConvertEnergyMax reset to ', &
           units_from_atomic(units_inp%energy, w_max),'[' // trim(units_abbrev(units_out%energy)) // ']'
      call messages_info(3)
      max_energy = w_max
    end if

    optimize = .false.
    optimize_parity = -1
    call fft_init(fft, nn, 1, FFT_REAL, FFTLIB_FFTW, optimize, optimize_parity)
    dw = M_TWO*M_PI / (dt * time_steps)
    e_start = int(min_energy / dw)
    e_end   = int(max_energy / dw)
    no_e    = e_end - e_start + 1
    e_step  = 1
    write(message(1),'(a,1x,i0.7,a,f12.7,a,i0.7,a,f12.7,a)')'Frequency index:',e_start,'(',&
         units_from_atomic(units_out%energy, e_start * dw),')-',e_end,'(',units_from_atomic(units_out%energy, e_end * dw),')' 
    write(message(2),'(a,f12.7,a)')'Frequency Step, dw:  ', units_from_atomic(units_out%energy, dw), &
         '[' // trim(units_abbrev(units_out%energy)) // ']'
    call messages_info(2)

    if (subtract_file) then
      write(ref_filename, '(a,a,a,a)') trim(ref_folder),"/", trim(ref_name),".obf"
      write(message(1),'(a,a)') "Reading ref-file from ", trim(ref_filename)
      call io_binary_read(trim(ref_filename), mesh%np, read_rff, ierr)
    end if
    
    !For each mesh point, open density file and read corresponding point.  
    if (mpi_world%rank == 0) call loct_progress_bar(-1, mesh%np)
    SAFE_ALLOCATE(write_point(1:mesh%np,e_point+1))

    ! Space
    do i_space = 1, mesh%np
      ! Time
      e_point = 0
      do i_time = c_start, c_end, c_step
         e_point = e_point + 1
        ! Here, we always iterate folders
        ! Delete the last / and add the corresponding folder number
        write(folder,'(a,i0.7,a)') in_folder(1:len_trim(in_folder)-1),i_time,"/"
        write(filename, '(a,a,a,a)') trim(outp%iter_dir), trim(folder), trim(basename), ".obf"
        if (mpi_world%size > 1) then
          ii = mesh%vp%local(mesh%vp%xlocal + i_space - 1)
        else
          ii = i_space
        end if
        ! Read the obf files, only one point per file
        call io_binary_read(trim(filename), 1, read_point(1:1), ierr, offset = ii-1)
        if (subtract_file) then
          read_ft(e_point) =  read_point(1) - read_rff(ii)
        else
          read_ft(e_point) = read_point(1)
        end if
        if (ierr /= 0) then
          write(message(1), '(a,a,2i10)') "Error reading the file ", trim(filename), ii, i_time
          write(message(2), '(a)') "Skipping...."
          call messages_warning(2)
          cycle
        end if
      end do ! Time

      call fftw_execute_dft(fft%planf, read_ft(1), out_fft(1))

      ! save densities
      do i_energy = e_start+1, e_end+1, e_step
        write_point(i_space, i_energy) = DBLE(out_fft(i_energy))
      end do ! Energy

      if (mod(i_space, 100) == 0 .and. mpi_world%rank == 0) then
        call loct_progress_bar(i_space-1, mesh%np) 
      end if
    end do ! Space

    if (mpi_world%rank == 0) call loct_progress_bar(mesh%np, mesh%np) 
#ifdef HAVE_MPI
    call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif
    
    ! write the output files
    do i_energy = e_start+1, e_end+1, e_step
      call io_mkdir('wd.general')
      write(filename,'(a14,i0.7,a1)')'wd.general/wd.',i_energy-1,'/'
      write(message(1),'(a,a,f12.7,a,1x,i7,a)')trim(filename),' w =', &
           units_from_atomic(units_out%energy,(i_energy-1) * dw), & 
           '[' // trim(units_abbrev(units_out%energy)) // ']'
      call messages_info(1)
      call io_mkdir(trim(filename))
      call dio_function_output(outp%how, trim(filename), & 
           trim('density'), mesh, write_point(:, i_energy), units_out%length**(-mesh%sb%dim), ierr, geo = geo)
    end do
    
    SAFE_DEALLOCATE_A(write_point)
    SAFE_DEALLOCATE_A(read_ft)
    SAFE_DEALLOCATE_A(out_fft)
    SAFE_DEALLOCATE_A(read_point)
    SAFE_DEALLOCATE_A(read_rff)

    POP_SUB(convert_transform)
  end subroutine convert_transform
end program

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
