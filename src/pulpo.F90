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
  use global
  use liboct

contains
  subroutine pulpo_print()
    character(len=2) :: lang
    
    ! some white space
    message(1) = ''; message(2) = ''
    call write_info(2)

    call oct_parse_str('RecipeLang', 'en', lang)
    call lowcase(lang)
    select case(lang)
    case('en')
      call oct_printRecipe(1)
    case('es')
      call oct_printRecipe(2)
    case default
      write(message(1), '(a,a,a)') "Language '", trim(lang), "' is not recognized"
      message(2) = "  RecipeLang = en | es"
      call write_fatal(2)
    end select

    ! DISCLAIMER
    call oct_printRecipe(0)

  end subroutine pulpo_print
end module pulpo
