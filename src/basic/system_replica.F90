!! Copyright (C) 2020 Heiko Appel
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

module system_replica_oct_m
  implicit none

  private
  public ::                   &
     system_replica_t


  integer, public, parameter ::        &
    UNIFORM_REPLICA             =  1,  &
    GAUSS_REPLICA               =  2,  &
    INPUT_REPLICA               =  3

  type :: system_replica_t
    integer :: n_replicas = 0
    logical :: is_replica = .false.
    integer :: replica_distribution = 1
    FLOAT   :: width = CNST(0.01)
  end type system_replica_t

contains

end module system_replica_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
