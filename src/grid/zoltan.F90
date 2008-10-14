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
  
  integer, public, parameter ::    &
       RCB        = 2,             &
       GRAPH      = 3,             &
       HYPERGRAPH = 4

  interface
    subroutine zoltan_partition(method, sbdim, np_global, np_part_global, x_global, estart, xedges, edges, ipart, part, fcomm)
      integer, intent(in)    :: method  
      integer, intent(in)    :: sbdim             ! the dimension of the space
      integer, intent(in)    :: np_global         ! the number of points to partition
      integer, intent(in)    :: np_part_global    ! the number of points including boundaries (required by x_global)
      FLOAT,   intent(in)    :: x_global          ! the coordinates of the points
      integer, intent(in)    :: estart            ! (local) where do the points for this partition start
      integer, intent(in)    :: xedges            ! (local) stores the position of each point in the array of the edges
      integer, intent(in)    :: edges             ! (local) the array of edges
      integer, intent(in)    :: ipart             ! (local) the index of the current partition
      integer, intent(inout) :: part              ! marks to which partition belongs each point
      integer, intent(in)    :: fcomm             ! the communicator
    end subroutine zoltan_partition
  end interface

end module zoltan_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
