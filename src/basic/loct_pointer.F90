!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!!
!! $Id$

#include "global.h"

module loct_pointer_m

  implicit none

  private
  public :: loct_pointer_copy

  interface loct_pointer_copy
    module procedure sloct_pointer_copy_1
    module procedure sloct_pointer_copy_2
    module procedure sloct_pointer_copy_3
    module procedure sloct_pointer_copy_4
    module procedure dloct_pointer_copy_1
    module procedure dloct_pointer_copy_2
    module procedure dloct_pointer_copy_3
    module procedure dloct_pointer_copy_4
    module procedure cloct_pointer_copy_1
    module procedure cloct_pointer_copy_2
    module procedure cloct_pointer_copy_3
    module procedure cloct_pointer_copy_4
    module procedure zloct_pointer_copy_1
    module procedure zloct_pointer_copy_2
    module procedure zloct_pointer_copy_3
    module procedure zloct_pointer_copy_4
    module procedure iloct_pointer_copy_1
    module procedure iloct_pointer_copy_2
    module procedure iloct_pointer_copy_3
    module procedure iloct_pointer_copy_4
    module procedure aloct_pointer_copy_1
    module procedure aloct_pointer_copy_2
    module procedure aloct_pointer_copy_3
    module procedure aloct_pointer_copy_4
    module procedure lloct_pointer_copy_1
    module procedure lloct_pointer_copy_2
    module procedure lloct_pointer_copy_3
    module procedure lloct_pointer_copy_4
  end interface loct_pointer_copy

contains

#  define TYPE real(4)
#  define SUBNAME(x) s ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE real(8)
#  define SUBNAME(x) d ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE complex(4)
#  define SUBNAME(x) c ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE complex(8)
#  define SUBNAME(x) z ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  integer
#  define SUBNAME(x) i ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  character(len=*)
#  define SUBNAME(x) a ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

#  define TYPE  logical
#  define SUBNAME(x) l ## x
#  include "loct_pointer_inc.F90"
#  undef SUBNAME
#  undef TYPE

end module loct_pointer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
