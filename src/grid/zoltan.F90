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
  use global_m
  use messages_m
  
  private

#ifdef HAVE_MPI  
  public ::                         &
       zoltan_partition
#endif

  public ::                         &
       zoltan_method_info,          &
       zoltan_method_is_geometric 
  
  ! this values have to match with the ones defined in zoltan_low.c
  integer, public, parameter ::    &
       RCB        = 2,             &
       RIB        = 3,             &
       HSFC       = 4,             &
       REFTREE    = 5,             &
       GRAPH      = 6,             &
       HYPERGRAPH = 7

  !%Variable MeshPartition
  !%Type integer
  !%Section Execution::Parallelization
  !%Description
  !% Decides which algorithm is used to partition the mesh. By default,
  !% <tt>graph</tt> partitioning is used for 8 or more partitions, and <tt>rcb</tt> for fewer.
  !%Option rcb 2
  !% Recursive coordinate bisection partitioning.
  !%Option rib 3
  !% Recursive inertial bisection partitioning.
  !%Option hsfc 4
  !% Hilbert space-filling curve partioning.
  !%Option reftree 5
  !% Refinement-tree-based partitioning.
  !%Option graph 6
  !% Graph partitioning.
  !%Option hypergraph 7
  !% Hypergraph partitioning.
  !%End

#ifdef HAVE_MPI
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
#endif

contains

  subroutine zoltan_method_info(method)
    integer, intent(in) :: method

    PUSH_SUB(zoltan_method_info)

    message(1) = 'Info: Using Zoltan to partition the mesh.'
  
    select case(method)
    case(RCB)
      message(2) = '      Recursive coordinate bisection partitioning.'
    case(RIB)
      message(2) = '      Recursive inertial bisection partitioning.'
    case(HSFC)
      message(2) = '      Hilbert space-filling curve partitioning.'
    case(REFTREE)
      message(2) = '      Refinement-tree-based partitioning.'
    case(GRAPH)
      message(2) = '      Graph partition algorithm.'
    case(HYPERGRAPH)
      message(2) = '      Hypergraph partition algorithm.'
    end select
    message(3) = ''
    call messages_info(3)

    POP_SUB(zoltan_method_info)
  end subroutine zoltan_method_info

  logical function zoltan_method_is_geometric(method) result(geometric)
    integer, intent(in) :: method

    PUSH_SUB(zoltan_method_is_geometric)

    geometric = method == RCB .or. method == RIB .or. method == HSFC .or. method == REFTREE

    POP_SUB(zoltan_method_is_geometric)
  end function zoltan_method_is_geometric

end module zoltan_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
