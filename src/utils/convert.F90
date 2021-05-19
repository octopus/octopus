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

!> This utility runs in parallel and can be used for post-processing of the results of Output.

#include "global.h"

program oct_convert
  use batch_oct_m
  use calc_mode_par_oct_m
  use command_line_oct_m
  use electrons_oct_m
  use fft_oct_m
  use fftw_params_oct_m
  use global_oct_m
  use io_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use ions_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use loct_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use spectrum_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m

  implicit none

  character(len=256) :: config_str
  integer :: ierr
  
  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call global_init()
  call calc_mode_par_init()

  call parser_init()
  
  call messages_init()

  call io_init()
  call profiling_init(global_namespace)
  call messages_experimental("oct-convert utility")

  call print_header()

  call messages_print_stress(stdout, "Convert mode")
  call messages_print_stress(stdout)

  call restart_module_init(global_namespace)
  call fft_all_init(global_namespace)
  call unit_system_init(global_namespace)

  call convert()

  call fft_all_end()
  call profiling_end(global_namespace)
  call io_end()
  call print_date("Calculation ended on ")
  call messages_end()

  call parser_end()

  call global_end()

contains

  ! -------------
  !> Reads an binary file and writes the equivalent files, 
  !! defined with OutputFormat.
  !! This is a high-level interface that reads the input file and
  !! calls the proper function.
  subroutine convert()
    type(electrons_t), pointer :: sys

    character(MAX_PATH_LEN)  :: basename, folder, ref_name, ref_folder, folder_default
    integer                  :: c_start, c_end, c_step, c_start_default, length, c_how
    logical                  :: iterate_folder, subtract_file
    integer, parameter       :: CONVERT_FORMAT = 1, FOURIER_TRANSFORM = 2, OPERATION = 3 

    PUSH_SUB(convert)

    call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
    sys => electrons_t(global_namespace)
    call sys%init_parallelization(mpi_world)

    message(1) = 'Info: Converting files'
    message(2) = ''
    call messages_info(2)

    !%Variable ConvertFilename
    !%Type string
    !%Default "density"
    !%Section Utilities::oct-convert
    !%Description
    !% Input filename. The original filename which is going to be converted in the format
    !% specified in <tt>OutputFormat</tt>. It is going to convert various files, it should 
    !% only contain the beginning of the name. For instance, in the case of the restart 
    !% files, it should be one space ' '.
    !%End
    call parse_variable(global_namespace, 'ConvertFilename', 'density', basename)
    if (basename == " ") basename = ""
    ! Delete the extension if present
    length = len_trim(basename)
    if (length > 4) then
      if (basename(length-3:length) == '.obf') then
        basename = trim(basename(1:length-4))
      end if
    end if

    !%Variable ConvertHow
    !%Type integer
    !%Default convert_format
    !%Section Utilities::oct-convert
    !%Description
    !% Select how the mesh function will be converted.
    !%Option format 1
    !% The format of the mesh function will be convert from the binary file.obf.
    !% The format of the output function is set by OutputHow variable.
    !%Option fourier_transform 2
    !% A fourier transform of the mesh function will be computed.
    !% It requieres that ConvertStart and ConvertEnd have to be set.
    !%Option operation 3
    !% Convert utility will generate a new mesh function constructed by linear 
    !% combination of scalar function of different mesh functions,
    !%End
    call parse_variable(global_namespace, 'ConvertHow', CONVERT_FORMAT, c_how)

    !%Variable ConvertIterateFolder
    !%Type logical
    !%Default true
    !%Section Utilities::oct-convert
    !%Description
    !% This variable decides if a folder is going to be iterated or the 
    !% filename is going to be iterated.
    !%End
    call parse_variable(global_namespace, 'ConvertIterateFolder', .true., iterate_folder)

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
    !% The folder name where the input files are. The default is
    !% <tt>td.</tt> if <tt>ConvertIterateFolder = true</tt>, otherwise <tt>restart</tt>.
    !%End
    call parse_variable(global_namespace, 'ConvertFolder', folder_default, folder)
    call add_last_slash(folder)

    !%Variable ConvertStart
    !%Type integer
    !%Section Utilities::oct-convert
    !%Description
    !% The starting number of the filename or folder.
    !% Default is 0 if <tt>ConvertIterateFolder = true</tt>, otherwise 1.
    !%End
    call parse_variable(global_namespace, 'ConvertStart', c_start_default, c_start)

    !%Variable ConvertEnd
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The last number of the filename or folder.
    !%End
    call parse_variable(global_namespace, 'ConvertEnd', 1, c_end)

    !%Variable ConvertStep
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The padding between the filenames or folder.
    !%End
    call parse_variable(global_namespace, 'ConvertStep', 1, c_step)

    !%Variable ConvertSubtractFilename
    !%Type string
    !%Default density
    !%Section Utilities::oct-convert
    !%Description
    !% Input filename. The file which is going to subtracted to rest of the files.
    !%End
    call parse_variable(global_namespace, 'ConvertSubtractFilename', 'density', ref_name)
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
    call parse_variable(global_namespace, 'ConvertSubtract', .false., subtract_file)

    !%Variable ConvertSubtractFolder
    !%Type string
    !%Default .
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name which is going to be subtracted.
    !%End
    call parse_variable(global_namespace, 'ConvertSubtractFolder', '.', ref_folder)
    call add_last_slash(folder)
    
    select case (c_how)
    CASE(OPERATION)
      call convert_operate(sys%gr%mesh, global_namespace, sys%space, sys%ions, sys%mc, sys%outp)

    CASE(FOURIER_TRANSFORM)
      ! Compute Fourier transform 
      call convert_transform(sys%gr%mesh, global_namespace, sys%space, sys%ions, sys%mc, sys%kpoints, basename, folder, &
         c_start, c_end, c_step, sys%outp, subtract_file, &
         ref_name, ref_folder)

    CASE(CONVERT_FORMAT)
      call convert_low(sys%gr%mesh, global_namespace, sys%space, sys%ions, sys%hm%psolver, sys%mc, basename, folder, &
         c_start, c_end, c_step, sys%outp, iterate_folder, &
         subtract_file, ref_name, ref_folder)
    end select

    SAFE_DEALLOCATE_P(sys)

    POP_SUB(convert)
  end subroutine convert

  ! ---------------------------------------------------------
  !> Giving a range of input files, it writes the corresponding 
  !! output files
  subroutine convert_low(mesh, namespace, space, ions, psolver, mc, basename, in_folder, c_start, c_end, c_step, outp, &
    iterate_folder, subtract_file, ref_name, ref_folder)
    type(mesh_t),      intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(ions_t),      intent(in)    :: ions
    type(poisson_t),   intent(in)    :: psolver
    type(multicomm_t), intent(in)    :: mc
    character(len=*),  intent(inout) :: basename       !< File name
    character(len=*),  intent(in)    :: in_folder      !< Folder name
    integer,           intent(in)    :: c_start        !< The first file number
    integer,           intent(in)    :: c_end          !< The last file number
    integer,           intent(in)    :: c_step         !< The step between files
    type(output_t),    intent(in)    :: outp           !< Output object; Decides the kind, what and where to output
    logical,           intent(in)    :: iterate_folder !< If true, it iterates over the folders, keeping the filename fixed.
                                                       !! If false, it iterates over the filenames
    logical,           intent(in)    :: subtract_file  !< If true, it subtracts the density from the reference 
    character(len=*),  intent(inout) :: ref_name       !< Reference file name 
    character(len=*),  intent(inout) :: ref_folder     !< Reference folder name

    type(restart_t)          :: restart
    integer                  :: ierr, ii, folder_index, output_i
    character(MAX_PATH_LEN)  :: filename, out_name, folder, frmt, restart_folder
    FLOAT, allocatable       :: read_ff(:), read_rff(:), pot(:)

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
      call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mc, ierr, &
                        dir=trim(ref_folder), mesh = mesh)
      ! FIXME: why only real functions? Please generalize.
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, space, trim(ref_name), mesh, read_rff, ierr)
        call restart_end(restart)
      else
        write(message(1),'(2a)') "Failed to read from ref-file ", trim(ref_name)
        write(message(2), '(2a)') "from folder ", trim(ref_folder)
        call messages_fatal(2)
      end if
    end if

    ! Initialize the restart directory from <tt>ConvertFolder</tt> value.
    ! This directory has to have the files 'grid' and 'lxyz.obf'
    ! and the files that are going to be converged, must be inside this folder
    if (iterate_folder) then
      ! Delete the last / and find the previous /, if any
      folder = in_folder(1:len_trim(in_folder)-1)
      folder_index = index(folder, '/', .true.)
      restart_folder = folder(1:folder_index)
    else 
      restart_folder = in_folder
    end if
    call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mc, ierr, &
                      dir=trim(restart_folder), mesh = mesh)
    call loct_progress_bar(-1, c_end-c_start)
    do ii = c_start, c_end, c_step
      if (iterate_folder) then
        ! Delete the last / and add the corresponding folder number
        write(folder,'(a,i0.7,a)') in_folder(folder_index+1:len_trim(in_folder)-1),ii,"/"
        write(filename, '(a,a,a)') trim(folder), trim(basename)
        out_name = trim(basename)
      else
        folder = ""
        if (c_start /= c_end) then
          ! Here, we are only considering 10 character long filenames.
          ! Subtract the initial part given at 'ConvertFilename' from the format and pad
          ! with zeros.
          write(frmt,'(a,i0,a)')"(a,i0.",10-len_trim(basename),")"
          write(filename, fmt=trim(frmt)) trim(basename), ii
          write(out_name, '(a)') trim(filename)
        else 
          ! Assuming filename is given complete in the 'ConvertFilename'
          write(filename, '(a,a,a,a)') trim(folder),"/", trim(basename)
          filename = basename
          write(out_name, '(a)') trim(basename)
        end if
      end if

      ! Read the obf file
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, space, trim(filename), mesh, read_ff, ierr)
      end if

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
      do output_i = lbound(outp%how, 1), ubound(outp%how, 1)
        if (outp%how(output_i) /= 0) then
          call dio_function_output(outp%how(output_i), trim(restart_folder)//trim(folder), & 
           trim(out_name), namespace, space, mesh, read_ff, units_out%length**(-space%dim), ierr, ions = ions)
        end if
      end do
      if (outp%what(OPTION__OUTPUT__POTENTIAL)) then
        write(out_name, '(a)') "potential"
        call dpoisson_solve(psolver, pot, read_ff)
        call dio_function_output(outp%how(OPTION__OUTPUT__POTENTIAL), trim(restart_folder)//trim(folder), &
             trim(out_name), namespace, space, mesh, pot, units_out%energy, ierr, ions = ions)
      end if
      call loct_progress_bar(ii-c_start, c_end-c_start) 
      ! It does not matter if the current write has failed for the next iteration
      ierr = 0
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
  subroutine convert_transform(mesh, namespace, space, ions, mc, kpoints, basename, in_folder, c_start, c_end, c_step, outp, & 
       subtract_file, ref_name, ref_folder)
    type(mesh_t)    ,  intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(ions_t),      intent(in)    :: ions
    type(multicomm_t), intent(in)    :: mc
    type(kpoints_t),   intent(in)    :: kpoints
    character(len=*),  intent(inout) :: basename       !< File name
    character(len=*),  intent(in)    :: in_folder      !< Folder name
    integer,           intent(in)    :: c_start        !< The first file number
    integer,           intent(in)    :: c_end          !< The last file number
    integer,           intent(in)    :: c_step         !< The step between files
    type(output_t),    intent(in)    :: outp           !< Output object; Decides the kind, what and where to output
    logical,           intent(in)    :: subtract_file  !< If true, it subtracts the density from the reference 
    character(len=*),  intent(inout) :: ref_name       !< Reference file name 
    character(len=*),  intent(inout) :: ref_folder     !< Reference folder name

    integer                 :: ierr, i_space, i_time, nn(1:3), optimize_parity(1:3), wd_info, output_i
    integer                 :: i_energy, e_end, e_start, e_point, chunk_size, read_count, t_point
    logical                 :: optimize(1:3)
    integer                 :: folder_index
    character(MAX_PATH_LEN) :: filename, folder, restart_folder
    FLOAT                   :: fdefault, w_max
    FLOAT, allocatable      :: read_ft(:), read_rff(:), point_tmp(:,:)

    integer, parameter      :: FAST_FOURIER = 1, STANDARD_FOURIER = 2
    type(kick_t)            :: kick
    integer                 :: ft_method
    type(spectrum_t)        :: spectrum
    type(batch_t)           :: tdrho_b, wdrho_b 
    FLOAT, allocatable      :: tdrho_a(:,:,:), wdrho_a(:,:,:)
    type(fft_t)             :: fft
    type(profile_t), save   :: prof_fftw, prof_io

    type(restart_t)         :: restart
    CMPLX, allocatable      :: out_fft(:)

    FLOAT   :: start_time          !< start time for the transform
    integer :: time_steps          !< number of time steps
    FLOAT   :: dt                  !< step in time mesh
    FLOAT   :: dw                  !< step in energy mesh
    FLOAT   :: max_energy          !< maximum of energy mesh
    FLOAT   :: min_energy          !< minimum of energy mesh

    PUSH_SUB(convert_transform)

    ! set default time_step as dt from TD section
    fdefault = M_ZERO
    call parse_variable(namespace, 'TDTimeStep', fdefault, dt, unit = units_inp%time)
    if (dt <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
      write(message(2),'(a)') 'Input: TDTimeStep reset to 0. Check input file'
      call messages_info(2)
    end if

    call io_mkdir('wd.general', namespace)
    wd_info = io_open('wd.general/wd.info', global_namespace, action='write')
    call messages_print_stress(wd_info, "Fourier Transform Options")

    !%Variable ConvertEnergyMin
    !%Type float
    !%Default 0.0
    !%Section Utilities::oct-convert
    !%Description
    !% Minimum energy to output from Fourier transform.
    !%End
    call parse_variable(namespace, 'ConvertEnergyMin', M_ZERO, min_energy, units_inp%energy)
    call messages_print_var_value(wd_info, 'ConvertEnergyMin', min_energy, unit = units_out%energy)

    !%Variable ConvertReadSize
    !%Type integer
    !%Default mesh%np
    !%Section Utilities::oct-convert
    !%Description
    !% How many points are read at once. For the parallel run this has not been
    !% yet tested, so it should be one. For the serial run, a number
    !% of 100-1000 will speed-up the execution time by this factor.
    !%End
    call parse_variable(namespace, 'ConvertReadSize', mesh%np, chunk_size)
    call messages_print_var_value(wd_info, 'ConvertReadSize', chunk_size)
    !Check that default value is set when ConvertReadSize = 0
    if ( chunk_size == 0) chunk_size = mesh%np
    ! Parallel version just work in domains and chunk_size equal to mesh%np 
    if(mesh%mpi_grp%size > 1 .and. chunk_size /= mesh%np) then
      write(message(1),*)'Incompatible value for ConvertReadSize and Parallelizaion in Domains'
      write(message(2),*)'Use the default value for ConvertReadSize (or set it to 0)'
      call messages_fatal(2)
    end if
    
    ! Calculate the limits in frequency space.
    start_time = c_start * dt
    dt         = dt * c_step
    time_steps = (c_end - c_start) / c_step 
    dw         = M_TWO * M_PI / (dt * time_steps)
    w_max      = M_TWO * M_PI / dt 

    !%Variable ConvertEnergyMax
    !%Type float
    !%Default w_max
    !%Section Utilities::oct-convert
    !%Description
    !% Maximum energy to output from Fourier transform.
    !%End
    fdefault = units_from_atomic(units_inp%energy, w_max)
    call parse_variable(namespace, 'ConvertEnergyMax',fdefault, max_energy, units_inp%energy)
    if (max_energy > w_max) then
      write(message(1),'(a,f12.7)')'Impossible to set ConvertEnergyMax to ', &
           units_from_atomic(units_inp%energy, max_energy)
      write(message(2),'(a)')'ConvertEnergyMax is too large.'
      write(message(3),'(a,f12.7,a)')'ConvertEnergyMax reset to ', &
           units_from_atomic(units_inp%energy, w_max),'[' // trim(units_abbrev(units_out%energy)) // ']'
      call messages_info(3)
      max_energy = w_max
    end if
    call messages_print_var_value(wd_info, 'ConvertEnergyMax', max_energy, unit = units_out%energy)

    !%Variable ConvertFTMethod
    !%Type integer
    !%Default FAST_FOURIER
    !%Section Utilities::oct-convert
    !%Description
    !% Describes the method used to perform the Fourier Transform
    !%Option fast_fourier 1
    !% Uses Fast Fourier Transform as implemented in the external library.
    !%Option standard_fourier 2
    !% Uses polinomial approach to the computation of discrete Fourier Transform.
    !% It uses the same variable described in how to obtain spectrum from
    !% a time-propagation calculation. 
    !%End
    call parse_variable(namespace, 'ConvertFTMethod', 1, ft_method)
    call messages_print_var_option(wd_info, 'ConvertFTMethod', ft_method)

    !TODO: check if e_point can be used instead of e_point+1
    SAFE_ALLOCATE(read_ft(0:time_steps))
    read_ft = M_ZERO
    SAFE_ALLOCATE(read_rff(1:mesh%np))

    select case(ft_method)
      case (FAST_FOURIER)
        nn(1) = time_steps + 1
        nn(2) = 1
        nn(3) = 1
        SAFE_ALLOCATE(out_fft(0:time_steps))
        optimize = .false.
        optimize_parity = -1
        call fft_init(fft, nn, 1, FFT_REAL, FFTLIB_FFTW, optimize, optimize_parity)
      case (STANDARD_FOURIER)
        !%Variable ConvertEnergyStep
        !%Type float
        !%Default <math>2 \pi / T</math>, where <math>T</math> is the total propagation time
        !%Section Utilities::oct-convert
        !%Description
        !% Energy step to output from Fourier transform.
        !% Sampling rate for the Fourier transform. If you supply a number equal or smaller than zero, then
        !% the sampling rate will be <math>2 \pi / T</math>, where <math>T</math> is the total propagation time.
        !%End
        fdefault = M_TWO * M_PI / (dt * time_steps)
        call parse_variable(namespace, 'ConvertEnergyStep',fdefault, dw, units_inp%energy)
        if (dw <= M_ZERO) dw = M_TWO * M_PI / (dt * time_steps)
        
        call spectrum_init(spectrum, namespace, dw, w_max)
        ! Manually setting already defined variables on spectrum.
        spectrum%start_time = c_start * dt
        spectrum%end_time = c_end * dt 
        spectrum%energy_step =  dw
        spectrum%max_energy = max_energy
        SAFE_ALLOCATE(tdrho_a(0:time_steps, 1, 1))
        SAFE_ALLOCATE(wdrho_a(0:time_steps, 1, 1))
    end select
    call messages_print_var_value(wd_info, 'ConvertEnergyStep', dw, unit = units_out%energy)

    !TODO: set system variable common for all the program in 
    !      order to use call kick_init(kick, sy%st%d%nspin, sys%space%dim, sys%ions%periodic_dim)
    call kick_init(kick, namespace, space, kpoints, 1)

    e_start = nint(min_energy / dw)
    e_end   = nint(max_energy / dw)
    write(message(1),'(a,1x,i0.7,a,f12.7,a,i0.7,a,f12.7,a)')'Frequency index:',e_start,'(',&
         units_from_atomic(units_out%energy, e_start * dw),')-',e_end,'(',units_from_atomic(units_out%energy, e_end * dw),')' 
    write(message(2),'(a,f12.7,a)')'Frequency Step, dw:  ', units_from_atomic(units_out%energy, dw), &
         '[' // trim(units_abbrev(units_out%energy)) // ']'
    call messages_info(2)

    if (subtract_file) then
      write(message(1),'(a,a,a,a)') "Reading ref-file from ", trim(ref_folder), trim(ref_name),".obf"
      call messages_info(1)
      call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mc, ierr, &
                        dir=trim(ref_folder), mesh = mesh)
      ! FIXME: why only real functions? Please generalize.
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, space, trim(ref_name), mesh, read_rff, ierr)
        call restart_end(restart)
      else
        write(message(1),'(2a)') "Failed to read from ref-file ", trim(ref_name)
        write(message(2), '(2a)') "from folder ", trim(ref_folder)
        call messages_fatal(2)
      end if
    end if
    
    call messages_print_stress(wd_info, "File Information")
    do i_energy = e_start, e_end
      write(filename,'(a14,i0.7,a1)')'wd.general/wd.',i_energy,'/'
      write(message(1),'(a,a,f12.7,a,1x,i7,a)')trim(filename),' w =', &
           units_from_atomic(units_out%energy,(i_energy) * dw), & 
           '[' // trim(units_abbrev(units_out%energy)) // ']'
      if (mpi_world%rank == 0) write(wd_info,'(a)') message(1)
      call messages_info(1)
      call io_mkdir(trim(filename), namespace)
    end do
    call io_close(wd_info)

    if(mesh%parallel_in_domains) then
      ! Delete the last / and find the previous /, if any
      folder = in_folder(1:len_trim(in_folder)-1)
      folder_index = index(folder, '/', .true.)
      restart_folder = folder(1:folder_index)
      call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mc, ierr, &
                        dir=trim(restart_folder), mesh = mesh)
    end if
   
    !For each mesh point, open density file and read corresponding point.  
    if (mpi_world%rank == 0) call loct_progress_bar(-1, mesh%np)

    SAFE_ALLOCATE(point_tmp(1:chunk_size, 0:time_steps))
    read_count = 0
    ! Space
    do i_space = 1, mesh%np
      ! Time
      t_point = 0
      do i_time = c_start, c_end, c_step
        
        if (mesh%parallel_in_domains .and. i_space == 1) then
          write(folder,'(a,i0.7,a)') in_folder(folder_index+1:len_trim(in_folder)-1),i_time,"/"
          write(filename, '(a,a,a)') trim(folder), trim(basename)
          call drestart_read_mesh_function(restart, space, trim(filename), mesh, point_tmp(:, t_point), ierr)
        else
          ! Here, we always iterate folders
          ! Delete the last / and add the corresponding folder number
          write(folder,'(a,i0.7,a)') in_folder(1:len_trim(in_folder)-1),i_time,"/"
          write(filename, '(a,a,a,a)') trim(folder), trim(basename), ".obf"
          if (mod(i_space-1, chunk_size) == 0) then
            call profiling_in(prof_io,"READING")
            !TODO: check for any error on the whole file before reading by parts.
            call io_binary_read(trim(filename), chunk_size, point_tmp(1:chunk_size, t_point), ierr, offset = i_space-1)
            call profiling_out(prof_io)
            if (i_time == c_start) read_count = 0
          end if
        end if

        call profiling_out(prof_io)

        if (ierr /= 0 .and. i_space == 1) then
          write(message(1), '(a,a,2i10)') "Error reading the file ", trim(filename), i_space, i_time
          write(message(2), '(a)') "Skipping...."
          write(message(3), '(a,i0)') "Error :", ierr
          call messages_warning(3)
          cycle
        end if

        if (i_time == c_start) read_count = read_count + 1
        if (subtract_file) then
          read_ft(t_point) =  point_tmp(read_count, t_point) - read_rff(i_space)
        else
          read_ft(t_point) = point_tmp(read_count, t_point)
        end if

        t_point = t_point + 1
      end do ! Time

      select case (ft_method)
      case (FAST_FOURIER)
        call profiling_in(prof_fftw, "CONVERT_FFTW")
        call dfft_forward(fft, read_ft, out_fft)
        call profiling_out(prof_fftw)
        ! Should the value be multiplied by dt ??? as in standard discrete Fourier Transform ?
        point_tmp(read_count, 0:time_steps) = AIMAG(out_fft(0:time_steps)) * dt
      case (STANDARD_FOURIER)
        tdrho_a(0:time_steps, 1, 1) = read_ft(0:time_steps)
        call batch_init(tdrho_b, 1, 1, 1, tdrho_a)
        call batch_init(wdrho_b, 1, 1, 1, wdrho_a)
        call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, c_start + 1, c_start + time_steps + 1, & 
                                  kick%time, dt, tdrho_b)
        call spectrum_fourier_transform(spectrum%method, spectrum%transform, spectrum%noise, &
              c_start + 1, c_start + time_steps + 1, kick%time, dt, tdrho_b, min_energy, max_energy, &
              spectrum%energy_step, wdrho_b)
        call tdrho_b%end()
        call wdrho_b%end()
        do e_point = e_start, e_end
          point_tmp(read_count, e_point) = - wdrho_a(e_point, 1, 1)
        end do
      end select

      if (mod(i_space-1, 1000) == 0 .and. mpi_world%rank == 0) then
        call loct_progress_bar(i_space-1, mesh%np) 
      end if
      
      !print out wd densities from (ii-chunksize,ii] if running in serial
      if (mesh%mpi_grp%size == 1) then
        if (mod(i_space, chunk_size) == 0) then
          write(message(1),'(a)') ""
          write(message(2),'(a,i0)') "Writing binary output: step ", i_space/chunk_size
          call messages_info(2)
          do i_energy = e_start, e_end
            write(filename,'(a14,i0.7,a12)')'wd.general/wd.',i_energy,'/density.obf'
            ! If it is the first time entering here, write the header. But, only once
            if (i_space == chunk_size) &
               !call write_header(trim(filename), mesh%np_global, ierr)
                call dwrite_header(trim(filename),mesh%np_global, ierr)
            call io_binary_write(trim(filename), chunk_size, point_tmp(1:chunk_size, i_energy), ierr, &
                nohead = .true.)
          end do
        end if
      end if
    end do ! Space

    if (mpi_world%rank == 0) call loct_progress_bar(mesh%np, mesh%np) 
#ifdef HAVE_MPI
    call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif

    if(mesh%parallel_in_domains) then
      do i_energy = e_start, e_end
        write(filename,'(a14,i0.7,a1)')'wd.general/wd.',i_energy,'/'
        do output_i = lbound(outp%how, 1), ubound(outp%how, 1)
          if (outp%how(output_i) /= 0) then
            call dio_function_output(0_8, trim(filename), & 
            trim('density'), namespace, space, mesh, point_tmp(:, i_energy), &
            units_out%length**(-space%dim), ierr, ions = ions)
          end if
        end do
      end do
      call restart_end(restart)
    else
      ! write the output files
      if (any(outp%how /= OPTION__OUTPUTFORMAT__BINARY)) then
        do i_energy = e_start, e_end
          write(filename,'(a14,i0.7,a1)')'wd.general/wd.',i_energy,'/'
          call io_binary_read(trim(filename)//'density.obf', mesh%np, read_rff, ierr)
          do output_i = lbound(outp%how, 1), ubound(outp%how, 1)
            if ((outp%how(output_i) /= 0) .and. (outp%how(output_i) /= OPTION__OUTPUTFORMAT__BINARY)) then
              call dio_function_output(outp%how(output_i), trim(filename), & 
                trim('density'), namespace, space, mesh, read_rff, &
                units_out%length**(-space%dim), ierr, ions = ions)
            end if
          end do
        end do
      end if
    end if
    
    SAFE_DEALLOCATE_A(point_tmp)
    SAFE_DEALLOCATE_A(read_ft)
    SAFE_DEALLOCATE_A(read_rff)

    select case (ft_method)
    case (FAST_FOURIER)
      SAFE_DEALLOCATE_A(out_fft)
    case (STANDARD_FOURIER)
      call kick_end(kick)
      SAFE_DEALLOCATE_A(tdrho_a)
      SAFE_DEALLOCATE_A(wdrho_a)
    end select

    POP_SUB(convert_transform)
  end subroutine convert_transform
  ! ---------------------------------------------------------
  !> Given a set of mesh function operations it computes a  
  !! a resulting mesh function from linear combination of them.
  subroutine convert_operate(mesh, namespace, space, ions, mc, outp)
    type(mesh_t),      intent(in)   :: mesh
    type(namespace_t), intent(in)   :: namespace
    type(space_t),     intent(in)   :: space
    type(ions_t),      intent(in)   :: ions
    type(multicomm_t), intent(in)   :: mc
    type(output_t)  ,  intent(in)   :: outp           !< Output object; Decides the kind, what and where to output

    integer             :: ierr, ip, i_op, length, n_operations, output_i
    type(block_t)       :: blk
    type(restart_t)     :: restart 
    type(unit_t)        :: units
    FLOAT               :: f_re, f_im
    FLOAT, allocatable  :: tmp_ff(:), scalar_ff(:)

    character(len=200) :: var, scalar_expression
    character(len=MAX_PATH_LEN) :: folder, filename, out_folder, out_filename

    PUSH_SUB(convert_operate)

    !%Variable ConvertScalarOperation
    !%Type block
    !%Section Utilities::oct-convert
    !%Description
    !% This variable is used to generate a new mesh function as a linear combination
    !% different mesh function having the same mesh. Each row defines an operation for
    !% for a single mesh function. 
    !% The format of the block is the following: <br>
    !% 'variable name' | 'folder' | 'file' | 'operation'
    !%End
    ! First, find out if there is a ConvertScalarOperation block.
    n_operations = 0
    if(parse_block(namespace, 'ConvertScalarOperation', blk) == 0) then
      n_operations = parse_block_n(blk)
    end if

    if (n_operations == 0) then
      write(message(1),'(a)')'No operations found. Check the input file'
      call messages_fatal(1)
    end if

    !%Variable ConvertOutputFolder
    !%Type string
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name where the output files will be write. The default is
    !% <tt>convert</tt>. 
    !%End
    call parse_variable(namespace, 'ConvertOutputFolder', "convert", out_folder)
    call add_last_slash(out_folder)
    call io_mkdir(out_folder, namespace)

    !%Variable ConvertOutputFilename
    !%Type string
    !%Default "density"
    !%Section Utilities::oct-convert
    !%Description
    !% Output filename. The name of the file in which the converted mesh function will be 
    !% written in the format specified in <tt>OutputFormat</tt>. 
    !%End
    call parse_variable(namespace, 'ConvertOutputFilename', 'density', out_filename)

    SAFE_ALLOCATE(tmp_ff(1:mesh%np))
    SAFE_ALLOCATE(scalar_ff(1:mesh%np))
    scalar_ff = M_ZERO

    do i_op = 1, n_operations
      !read variable name
      call parse_block_string(blk, i_op-1, 0, var)
      !read folder path
      call parse_block_string(blk, i_op-1, 1, folder)
      call add_last_slash(folder)
      !read file
      call parse_block_string(blk, i_op-1, 2, filename)
      ! Delete the extension if present
      length = len_trim(filename)
      if (length > 4) then
        if (filename(length-3:length) == '.obf') then
          filename = trim(filename(1:length-4))
        end if
      end if
      ! FIXME: why only real functions? Please generalize.
      ! TODO: check if mesh function are real or complex.
      call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, mc, ierr, &
                        dir=trim(folder), mesh = mesh, exact=.true.)
      if(ierr == 0) then
        call drestart_read_mesh_function(restart, space, trim(filename), mesh, tmp_ff, ierr)
      else
        write(message(1),'(2a)') "Failed to read from file ", trim(filename)
        write(message(2), '(2a)') "from folder ", trim(folder)
        call messages_fatal(2)
      end if
      !read scalar expression
      call parse_block_string(blk, i_op-1, 3, scalar_expression)

      do ip = 1, mesh%np
        call parse_expression(f_re, f_im, trim(var), TOFLOAT(tmp_ff(ip)), trim(scalar_expression))
        !TODO: implement use of complex functions. 
        scalar_ff(ip) = scalar_ff(ip) + f_re
      end do
      
      call restart_end(restart)

    end do

    call parse_block_end(blk)

#ifdef HAVE_MPI
    call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif
    ! Write the corresponding output
    !TODO: add variable ConvertFunctionType to select the type(density, wfs, potential, ...) 
    !      and units of the conversion.
    units = units_out%length**(-space%dim)
    do output_i = lbound(outp%how, 1), ubound(outp%how, 1)
      if (outp%how(output_i) /= 0) then
        call dio_function_output(outp%how(output_i), trim(out_folder), trim(out_filename), namespace, space, mesh, & 
          scalar_ff, units, ierr, ions = ions)
      end if
    end do

    SAFE_DEALLOCATE_A(tmp_ff)
    SAFE_DEALLOCATE_A(scalar_ff)

    POP_SUB(convert_operate)
  end subroutine convert_operate
end program

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
