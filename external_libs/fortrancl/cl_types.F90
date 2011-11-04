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

#include "config_F90.h"

module cl_types_m
  implicit none 
  
  type :: cl_platform_id
    private 
    integer, pointer :: p 
  end type cl_platform_id
  
  type :: cl_device_id
    private 
    integer, pointer :: p 
  end type cl_device_id

  type :: cl_context
    private 
    integer, pointer :: p 
  end type cl_context

  type :: cl_command_queue
    private 
    integer, pointer :: p 
  end type cl_command_queue

  type :: cl_mem
    private 
    integer, pointer :: p 
  end type cl_mem

  type :: cl_program
    private 
    integer, pointer :: p 
  end type cl_program

  type :: cl_kernel
    private 
    integer, pointer :: p 
  end type cl_kernel

  type :: cl_event
    private 
    integer, pointer :: p 
  end type cl_event

  type :: cl_sampler
    private 
    integer, pointer :: p 
  end type cl_sampler

end module cl_types_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
