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

#include "global.h"

program xyzanim
  use command_line_oct_m
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  integer :: ierr
  
  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_xyz_anim()
  call getopt_end()

  call parser_init()
  
  call messages_init()
  call io_init()
  call unit_system_init(global_namespace)

  call generate_xyz_anim()

  call io_end()
  call messages_end()

  call parser_end()
  call global_end()

contains

  subroutine generate_xyz_anim()

    character(len=256) :: coords_file, comment
    integer :: ierr, sampling, i, coords_unit, iter, j, record_length
    logical :: multifiles
    FLOAT :: time
    type(ions_t),     pointer :: ions

    ! Sets the filenames
    coords_file = 'td.general/coordinates'

    !%Variable AnimationSampling
    !%Type integer
    !%Default 100
    !%Section Utilities::oct-xyz-anim
    !%Description
    !% Sampling rate of the animation. The animation will be constructed using
    !% the iteration numbers that are multiples of <tt>AnimationSampling<tt>.
    !%End
    call parse_variable(global_namespace, 'AnimationSampling', 100, sampling)
    if(sampling < 1) then
      message(1) = 'Sampling rate (AnimationSampling) should be bigger than 0'
      call messages_fatal(1)
    end if

    !%Variable AnimationMultiFiles
    !%Type logical
    !%Default false
    !%Section Utilities::oct-xyz-anim
    !%Description
    !% If true, each iteration written will be in a separate file.
    !%End
    call parse_variable(global_namespace, 'AnimationMultiFiles', .false., multifiles)

    ions => ions_t(global_namespace)

    record_length = 100 + ions%space%dim*ions%natoms*3*20

    ! Opens the coordinates file
    coords_unit = io_open(coords_file, ions%namespace, action='read', recl = record_length)

    call io_skip_header(coords_unit)
    ierr = 0
    do while(ierr == 0)
      read(unit = coords_unit, iostat = ierr, fmt = *) iter, time, &
        ((ions%pos(j, i), j = 1, ions%space%dim), i = 1, ions%natoms)
      ions%pos = units_to_atomic(units_out%length, ions%pos)
      if(mod(iter, sampling) == 0) then
        write(comment, '(i10,f20.6)') iter, time
        if(.not.multifiles)then
          call io_mkdir('td.general', global_namespace)
          call ions%write_xyz('td.general/movie', append = .true., comment = trim(comment))
        else
          call io_mkdir('td.general/movie/', ions%namespace)
          write(coords_file,'(i7.7)')iter
          call ions%write_xyz('td.general/movie/geo-' + trim(coords_file), append = .false.)
        end if
      end if
    end do

    SAFE_DEALLOCATE_P(ions)

    call io_close(coords_unit)

  end subroutine generate_xyz_anim

end program xyzanim

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
