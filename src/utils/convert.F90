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
  use calc_mode_m
  use command_line_m
  use datasets_m
  use derivatives_m
  use fft_m
  use global_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use loct_m
  use messages_m
  use mpi_m
  use multicomm_m
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

  if(no_datasets > 1) then
    message(1) = 'Info: Multi-Dataset Mode'
    message(2) = 'Info: Running dataset "'//trim(current_label)//'"'
    call messages_info(2, stress = .true.)
  end if

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

    call io_function_convert(sys%gr%mesh, sys%geo, basename, folder, &
         c_start, c_end, c_step, sys%outp%how, iterate_folder, &
         subtract_file, refname, ref_folder )

    call system_end(sys)

    POP_SUB(convert)
  end subroutine convert

end program

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
