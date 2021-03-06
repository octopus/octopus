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
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use space_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  character(len=256) :: coords_file, comment
  integer :: ierr, sampling, i, coords_unit, iter, j, record_length
  logical :: multifiles
  FLOAT :: time
  type(geometry_t)  :: geo
  type(space_t)     :: space
  
  ! Initialize stuff
  call global_init(is_serial = .true.)

  call getopt_init(ierr)
  if(ierr == 0) call getopt_xyz_anim()
  call getopt_end()

  call parser_init()
  
  call messages_init()
  call io_init()
  call unit_system_init(global_namespace)

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

  call space_init(space, global_namespace)
  call geometry_init(geo, global_namespace, space)

  record_length = 100 + geo%space%dim*geo%natoms*3*20

  ! Opens the coordinates file
  coords_unit = io_open(coords_file, global_namespace, action='read', recl = record_length)

  call io_skip_header(coords_unit)
  ierr = 0
  do while(ierr == 0)
    read(unit = coords_unit, iostat = ierr, fmt = *) iter, time, &
      ((geo%atom(i)%x(j), j = 1, geo%space%dim), i = 1, geo%natoms)
      do i = 1, geo%natoms
        do j = 1, geo%space%dim
          geo%atom(i)%x(j)=units_to_atomic(units_out%length, geo%atom(i)%x(j))
        end do
      end do
    if(mod(iter, sampling) == 0) then
      write(comment, '(i10,f20.6)') iter, time
      if(.not.multifiles)then
        call io_mkdir('td.general', global_namespace)
        call geometry_write_xyz(geo, 'td.general/movie', global_namespace, &
          append = .true., comment = trim(comment))
      else
        call io_mkdir('td.general/movie/', global_namespace)
        write(coords_file,'(i7.7)')iter
        call geometry_write_xyz(geo,'td.general/movie/geo-' + trim(coords_file), global_namespace, &
          append = .false.)
      end if
    end if
  end do

  call geometry_end(geo)

  call io_close(coords_unit)

  call io_end()
  call messages_end()

  call parser_end()
  call global_end()

end program xyzanim

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
