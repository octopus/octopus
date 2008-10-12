!! Copyright (C) 2007 X. Andrade
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
!! $Id: submesh.F90 2781 2007-03-23 10:58:32Z lorenzen $
 
#include "global.h"
  
module zoltan_m
  private
  
  public :: zoltan_partition
  
  integer, public ::    &
       GEOMETRIC  = 2,  &
       GRAPH      = 3,  &
       HYPERGRAPH = 4

  interface
    subroutine zoltan_partition(method, sbdim, np_global, np_part_global, x_global, xedges, edges, ipart, part, fcomm)
      integer, intent(in)    :: method
      integer, intent(in)    :: sbdim
      integer, intent(in)    :: np_global
      integer, intent(in)    :: np_part_global
      FLOAT,   intent(in)    :: x_global
      integer, intent(in)    :: xedges
      integer, intent(in)    :: edges
      integer, intent(in)    :: ipart
      integer, intent(out)   :: part
      integer, intent(in)    :: fcomm
    end subroutine zoltan_partition
  end interface

end module zoltan_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
