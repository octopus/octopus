!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module pulpo
  use io
  use global
  use lib_oct

contains
  subroutine pulpo_print()
    character(len=256) :: filename

    ! some white space
    message(1) = ''; message(2) = ''
    call write_info(2)

    call loct_printRecipe(trim(conf%share), filename)
    call io_dump_file(stdout, filename)
    call write_info(2)    
    call io_dump_file(stdout, trim(conf%share)//"/recipes/disclaimer.txt")
    call write_info(2)

  end subroutine pulpo_print
end module pulpo
