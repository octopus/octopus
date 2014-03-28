!! Copyright (C) 2013 J. Alberdi-Rodriguez
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
  use command_line_m
  use datasets_m
  use derivatives_m
  use fft_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use io_binary_m
  use loct_m
  use messages_m
  use mesh_m
  use mpi_m
  use multicomm_m
  use output_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_calc_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  character*256 :: config_str
  integer :: test_type
  integer :: conv_mode
  integer :: ierr
  
  integer, parameter ::              &
       CONV_FROM_BINARY          =   1
  
  
  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call global_init()
  call messages_init()

  !%Variable ConvertMode
  !%Type integer
  !%Default from binary to dx
  !%Section Utilities::oct-convert
  !%Description
  !% Decides what kind of conversion should be performed.
  !%Option conv_from_binary 1
  !% Reads an obf file and writes an equivalent (in the same
  !% folder) in the format given by <tt>OutputHow</tt>.
  !%End
  call parse_integer('ConvertMode', CONV_FROM_BINARY, conv_mode)
  call datasets_init(conv_mode)  

  call io_init()
  call profiling_init()

  call print_header()

  call messages_print_stress(stdout, "Convert mode")
  call messages_print_var_option(stdout, "ConvertMode", conv_mode)
  call messages_print_stress(stdout)

  call fft_all_init()
  call unit_system_init()

  select case(conv_mode)
  case(CONV_FROM_BINARY)
    call convert()
  end select

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

    character(64)  :: basename, folder, refname, ref_folder, folder_default
    integer        :: c_start, c_end, c_step, c_start_default
    logical        :: iterate_folder, subtract_file

    PUSH_SUB(convert)

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
    !% specified in <tt>OutputHow</tt>.
    !%End
    call parse_string(datasets_check('ConvertFilename'), 'density', basename)
    if ( basename == " " ) basename = ""

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

    !%Variable ConvertStart
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-convert
    !%Description
    !% The starting number of the filename.
    !%End
    call parse_integer(datasets_check('ConvertStart'), c_start_default, c_start)

    !%Variable ConvertEnd
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The last number of the filename.
    !%End
    call parse_integer(datasets_check('ConvertEnd'), 1, c_end)

    !%Variable ConvertStep
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-convert
    !%Description
    !% The padding between the filenames.
    !%End
    call parse_integer(datasets_check('ConvertStep'), 1, c_step)

    !%Variable ConvertRefFileName
    !%Type string
    !%Default density
    !%Section Utilities::oct-convert
    !%Description
    !% Input filename. The original filename which is going to convert in the formats 
    !% specified in <tt>OutputHow</tt>.
    !%End
    call parse_string(datasets_check('ConvertRefFileName'), 'density', refname)
    if ( refname == " " ) refname = ""

    !%Variable ConvertSubtractFile
    !%Type logical
    !%Default false
    !%Section Utilities::oct-convert
    !%Description
    !% The reference file that is going to be used to subtract from.
    !%End
    call parse_logical(datasets_check('ConvertSubtractFile'), .false., subtract_file)

    !%Variable ConvertSubtractFolder
    !%Type string
    !%Default [blank]
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name which is going to be subtracted.
    !%End
    call parse_string(datasets_check('ConvertSubtractFolder'), ' ', ref_folder)
    if ( ref_folder == " " ) ref_folder = ""

    call convert_low(sys%gr%mesh, sys%geo, basename, folder, &
         c_start, c_end, c_step, sys%outp%how, sys%outp%what, iterate_folder, &
         subtract_file, refname, ref_folder )

    call system_end(sys)

    POP_SUB(convert)
  end subroutine convert

  ! ---------------------------------------------------------
  !> Giving a range of input files, it writes the corresponding 
  !! output files
  subroutine convert_low(mesh, geo, basename, folder, c_start, c_end, c_step, how, what, iterate_folder, & 
                                 subtract_file, ref_name, ref_folder)
    type(mesh_t)    , intent(in)    :: mesh
    type(geometry_t), intent(in)    :: geo
    character(len=*), intent(inout) :: basename       !< File name
    character(len=*), intent(inout) :: folder         !< Folder name
    integer,          intent(in)    :: c_start        !< The first file number
    integer,          intent(in)    :: c_end          !< The last file number
    integer,          intent(in)    :: c_step         !< The step between files
    integer,          intent(in)    :: how            !< Decides the kind of the output
    integer,          intent(in)    :: what           !< Decides what is going to be written
    logical,          intent(in)    :: iterate_folder !< If true, it iterates over the folders, keeping the filename fixed.
                                                      !! If false, it iterates over the filenames 
    character(len=*), intent(inout) :: ref_name       !< Reference file name 
    character(len=*), intent(inout) :: ref_folder     !< Reference folder name
    logical,          intent(in)    :: subtract_file  !< If true, it subtracts the density from the reference

    integer            :: ierr, ii
    character(64)      :: filename, out_name, ref_filename
    FLOAT, allocatable :: read_ff(:), read_rff(:), pot(:)

    PUSH_SUB(convert_low)

    SAFE_ALLOCATE(read_ff(1:mesh%np))
    SAFE_ALLOCATE(read_rff(1:mesh%np))
    SAFE_ALLOCATE(pot(1:mesh%np))
    read_rff(:) = M_ZERO
   
    write(message(1),'(5a,i5,a,i5,a,i5)') "Converting '", trim(folder), "//", trim(basename), &
         "' from ", c_start, " to ", c_end, " every ", c_step
    call messages_info(1)
 
    if (subtract_file) then
      write(ref_filename, '(a,a,a)') trim(ref_folder), trim(ref_name),".obf"
      call io_binary_read(trim(ref_filename), mesh%np, read_rff, ierr)
    endif

    call loct_progress_bar(-1, c_end-c_start) 
    do ii = c_start, c_end, c_step
      if (iterate_folder) then
        write(folder,'(a,i0.7,a)') "td.",ii,"/"
        write(filename, '(a,a,a)') trim(folder), trim(basename), ".obf"
        out_name = trim(basename)
      else
        write(filename, '(a,a,a,a)') trim(folder),"/", trim(basename),".obf"
        write(out_name, '(a)') trim(basename)
      end if

      ! Read the obf file
      call io_binary_read(trim(filename), mesh%np, read_ff, ierr)

      if (ierr /= 0) then
        write(message(1), '(a,a)') "Error reading the file ", filename
        write(message(2), '(a)') "Skipping...."
        call messages_warning(2)
      end if
      if (subtract_file) read_ff(:) = read_ff(:) - read_rff(:) 
      if (subtract_file) write(out_name, '(a,a)') trim(out_name),"-ref"
      ! Write the corresponding output
      call dio_function_output(how, &
        trim(folder), trim(out_name), mesh, read_ff, units_out%length, ierr, geo = geo)
      
      if (iand(what, C_OUTPUT_POTENTIAL) /= 0) then
        write(out_name, '(a)') "potential"
        call dpoisson_solve(psolver, pot, read_ff)
        call dio_function_output(how, &
             trim(folder), trim(out_name), mesh, pot, units_out%length, ierr, geo = geo)
      end if
      call loct_progress_bar(ii-c_start, c_end-c_start) 
    end do
    
    SAFE_DEALLOCATE_A(read_ff)
    SAFE_DEALLOCATE_A(read_rff)
    SAFE_DEALLOCATE_A(pot)
    POP_SUB(convert_low)
  end subroutine convert_low
  
end program

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
