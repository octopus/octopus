!! Copyright (C) 2005-2006 Heiko Appel
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

module td_transport_m
  use messages_m

  implicit none

  private
  public :: &
    td_transport_run


contains

  ! ---------------------------------------------------------
  subroutine td_transport_run()
    call push_sub('td_transport.td_transport_run')

    message(1) = 'Time dependent quantum transport not yet implemented.'
    call write_fatal(1)

    call pop_sub()
  end subroutine td_transport_run
end module td_transport_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
