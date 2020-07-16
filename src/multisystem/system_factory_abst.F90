!! Copyright (C) 2020 M. Oliveira
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

module system_factory_abst_oct_m
  use namespace_oct_m
  use system_oct_m
  use system_replica_oct_m
  implicit none

  private
  public ::                        &
    system_factory_abst_t

  type, abstract :: system_factory_abst_t
  contains
    procedure(system_factory_abst_create), deferred :: create
  end type system_factory_abst_t

  abstract interface
    function system_factory_abst_create(this, namespace, name, type, system_replica) result(system)
      import :: system_factory_abst_t
      import system_t
      import namespace_t
      import system_replica_t
      class(system_factory_abst_t), intent(in) :: this
      type(namespace_t),            intent(in) :: namespace
      character(len=*),             intent(in) :: name
      integer,                      intent(in) :: type
      type(system_replica_t),       intent(inout) :: system_replica
      class(system_t),              pointer    :: system
    end function system_factory_abst_create
  end interface

end module system_factory_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
