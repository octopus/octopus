!! Copyright (C) 2019 N. Tancogne-Dejean
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
!! along with st program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module hamiltonian_abst_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    hamiltonian_abst_t,                &
    hamiltonian_are_complex,               &
    hamiltonian_are_real,                  &
    hamiltonian_set_complex

  type, abstract :: hamiltonian_abst_t
    private

  contains

    procedure(unpack),     deferred :: unpack  
  end type hamiltonian_abst_t

  abstract interface

    subroutine unpack(st, copy)
      import hamiltonian_abst_t
      class(hamiltonian_abst_t), intent(inout) :: st
      logical, optional,    intent(in)    :: copy
    end subroutine unpack
  end interface

contains


end module hamiltonian_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
