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

#include "config_F90.h"

module output
  use liboct
  use units
  use mesh
  
implicit none

type output_type
  logical :: what(6)
  integer :: how    ! how to output

  integer :: iter   ! output every iter
  
  integer :: wfs(32) ! which wfs to output
end type output_type

integer, parameter :: &
     output_potential = 1, &
     output_density   = 2, &
     output_wfs       = 3, &
     output_ELF       = 4, &
     output_geometry  = 5, &
     output_something = 6   ! this one should be the last

integer, parameter ::     &
     output_plane_x = 1,  &
     output_plane_y = 2,  &
     output_plane_z = 4,  &
     output_xyz     = 8,  &
     output_dx      = 16

contains

subroutine output_init(outp)
  type(output_type), intent(out) :: outp

  integer :: i
  logical :: l
  character(len=80) :: owf

  call oct_parse_logical("OutputKSPotential", .false., outp%what(output_potential))
  call oct_parse_logical("OutputDensity",     .false., outp%what(output_density))
  call oct_parse_logical("OutputWfs",         .false., outp%what(output_wfs))
  call oct_parse_logical("OutputELF",         .false., outp%what(output_elf))
  call oct_parse_logical("OutputGeometry",    .false., outp%what(output_geometry))
  
  outp%what(output_something) = .false.
  do i = 1, output_something - 1
    outp%what(output_something) = outp%what(output_something).or.outp%what(i)
  end do

  if(outp%what(output_wfs)) then
    call oct_parse_str("OutputWfsNumber", "1-1024", owf)
    call oct_wfs_list(C_string(owf), outp%wfs)
  end if

  if(outp%what(output_something)) then

#if defined(THREE_D)
    call oct_parse_logical("OutputPlaneX", .false., l)
    if(l) outp%how = ior(outp%how, output_plane_x)
    call oct_parse_logical("OutputPlaneY", .false., l)
    if(l) outp%how = ior(outp%how, output_plane_y)
    call oct_parse_logical("OutputPlaneZ", .false., l)
    if(l) outp%how = ior(outp%how, output_plane_z)
    call oct_parse_logical("OutputXYZ",    .false., l)
    if(l) outp%how = ior(outp%how, output_xyz)
    call oct_parse_logical("OutputDX", .false., l)
    if(l) outp%how = ior(outp%how, output_dx)
#endif

    call oct_parse_int(C_string("OutputEvery"), 1000, outp%iter)
  end if
  
end subroutine output_init

#include "undef.F90"
#include "real.F90"
#include "out_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "out_inc.F90"
#include "undef.F90"

end module output
