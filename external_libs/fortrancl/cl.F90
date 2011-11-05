!! Copyright (C) 2010 X. Andrade
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
!! $Id: cl.F90 3587 2007-11-22 16:43:00Z xavier $

module cl

  ! This is the module that should be used by the users.

  use cl_types_m
  use cl_buffer_m
  use cl_command_queue_m
  use cl_constants_m
  use cl_context_m
  use cl_device_m
  use cl_kernel_m
  use cl_platform_m
  use cl_program_m

  implicit none 

end module cl

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
