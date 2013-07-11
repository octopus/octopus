!! Copyright (C) 2013 M. Oliveira
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
!! $Id$

#include "global.h"

! -----------------------------------------------------------------------
!> This module contains interfaces for METIS and PARMETIS routines
! -----------------------------------------------------------------------

module metis_m
  implicit none

  public

  ! Options codes (copied from metis.h)
  integer, parameter ::             &
       METIS_OPTION_PTYPE     = 1,  &
       METIS_OPTION_OBJTYPE   = 2,  &
       METIS_OPTION_CTYPE     = 3,  &
       METIS_OPTION_IPTYPE    = 4,  &
       METIS_OPTION_RTYPE     = 5,  &
       METIS_OPTION_DBGLVL    = 6,  &
       METIS_OPTION_NITER     = 7,  &
       METIS_OPTION_NCUTS     = 8,  &
       METIS_OPTION_SEED      = 9,  &
       METIS_OPTION_NO2HOP    = 10, &
       METIS_OPTION_MINCONN   = 11, &
       METIS_OPTION_CONTIG    = 12, &
       METIS_OPTION_COMPRESS  = 13, &
       METIS_OPTION_CCORDER   = 14, &
       METIS_OPTION_PFACTOR   = 15, &
       METIS_OPTION_NSEPS     = 16, &
       METIS_OPTION_UFACTOR   = 17, &
       METIS_OPTION_NUMBERING = 18, &
       METIS_OPTION_HELP      = 19, &
       METIS_OPTION_TPWGTS    = 20, &
       METIS_OPTION_NCOMMON   = 21, &
       METIS_OPTION_NOOUTPUT  = 22, &
       METIS_OPTION_BALANCE   = 23, &
       METIS_OPTION_GTYPE     = 24, &
       METIS_OPTION_UBVEC     = 25

  interface 

#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
    subroutine oct_metis_setdefaultoptions(options)
      integer, intent(inout) :: options
    end subroutine oct_metis_setdefaultoptions

    subroutine oct_metis_partgraphrecursive(nvtxs, ncon, xadj, adjncy, nparts, tpwgts, ubvec, options, objval, part)
      integer,     intent(in)  :: nvtxs      !< The number of vertices in the graph.
      integer,     intent(in)  :: ncon       !< The number of balancing constraints. It should be at least 1.
      integer,     intent(in)  :: xadj       !< The adjacency structure of the graph.
      integer,     intent(in)  :: adjncy     !< The adjacency structure of the graph.
      integer,     intent(in)  :: nparts     !< The number of parts to partition the graph.
      REAL_SINGLE, intent(in)  :: tpwgts     !< This is an array of size nparts x ncon that specifies the desired 
                                             !! weight for each partition and constraint.
      REAL_SINGLE, intent(in)  :: ubvec      !< This is an array of size ncon that specifies the allowed load imbalance 
                                             !! tolerance for each constraint.
      integer,     intent(in)  :: options    !< This is the array of options.
      integer,     intent(out) :: objval     !< Upon successful completion, this variable stores the edge-cut or the total 
                                             !! communication volume of the partitioning solution.
      integer,     intent(out) :: part       !< This is a vector of size nvtxs that upon successful completion stores the 
                                             !! partition vector of the graph.
    end subroutine oct_metis_partgraphrecursive

    subroutine oct_metis_partgraphkway(nvtxs, ncon, xadj, adjncy, nparts, tpwgts, ubvec, options, objval, part)
      integer,     intent(in)  :: nvtxs      !< The number of vertices in the graph.
      integer,     intent(in)  :: ncon       !< The number of balancing constraints. It should be at least 1.
      integer,     intent(in)  :: xadj       !< The adjacency structure of the graph.
      integer,     intent(in)  :: adjncy     !< The adjacency structure of the graph.
      integer,     intent(in)  :: nparts     !< The number of parts to partition the graph.
      REAL_SINGLE, intent(in)  :: tpwgts     !< This is an array of size nparts x ncon that specifies the desired weight for 
                                             !! each partition and constraint.
      REAL_SINGLE, intent(in)  :: ubvec      !< This is an array of size ncon that specifies the allowed load imbalance 
                                             !! tolerance for each constraint.
      integer,     intent(in)  :: options    !< This is the array of options.
      integer,     intent(out) :: objval     !< Upon successful completion, this variable stores the edge-cut or the total 
                                             !! communication volume of the partitioning solution.
      integer,     intent(out) :: part       !< This is a vector of size nvtxs that upon successful completion stores the 
                                             !! partition vector of the graph.
    end subroutine oct_metis_partgraphkway

#endif
#if defined(HAVE_PARMETIS)

    subroutine oct_parmetis_v3_partkway(vtxdist, xadj, adjncy, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm)
      integer,     intent(in)  :: vtxdist    !< This array describes how the vertices of the graph are distributed among 
                                             !! the processors. Its contents are identical for every processor.
      integer,     intent(in)  :: xadj       !< These store the (local) adjacency structure of the graph at each processor.
      integer,     intent(in)  :: adjncy     !< These store the (local) adjacency structure of the graph at each processor.
      integer,     intent(in)  :: ncon       !< This is used to specify the number of weights that each vertex has. It is 
                                             !! also the number of balance constraints that must be satisfied.
      integer,     intent(in)  :: nparts     !< This is used to specify the number of sub-domains that are desired. Note 
                                             !! that the number of sub-domains is independent of the number of processors 
                                             !! that call this routine.
      REAL_SINGLE, intent(in)  :: tpwgts     !< An array of size ncon × nparts that is used to specify the fraction of vertex 
                                             !! weight that should be distributed to each sub-domain for each balance 
                                             !! constraint. If all of the sub-domains are to be of the same size for every 
                                             !! vertex weight, then each of the ncon × nparts elements should be set to
                                             !! a value of 1/nparts. If ncon is greater than 1, the target sub-domain weights
                                             !! for each sub-domain are stored contiguously (similar to the vwgt array). 
                                             !! Note that the sum of all of the tpwgts for a give vertex weight should be one.
      REAL_SINGLE, intent(in)  :: ubvec      !< An array of size ncon that is used to specify the imbalance tolerance for each 
                                             !! vertex weight, with 1 being perfect balance and nparts being perfect imbalance.
                                             !! A value of 1.05 for each of the ncon weights is recommended.
      integer,     intent(in)  :: options    !<  This is an array of integers that is used to pass additional parameters for 
                                             !! the routine. The first element (i.e., options[0]) can take either the value of 
                                             !! 0 or 1. If it is 0, then the default values are used.
      integer,     intent(out) :: edgecut    !< Upon successful completion, the number of edges that are cut by the 
                                             !! partitioning is written to this parameter.
      integer,     intent(out) :: part       !<  This is an array of size equal to the number of locally-stored vertices. Upon 
                                             !! successful completion the partition vector of the locally-stored vertices is 
                                             !! written to this array.
      integer,     intent(in)  :: comm       !< This is a pointer to the MPI communicator of the processes that call PARMETIS. 
    end subroutine oct_parmetis_v3_partkway

#endif
  end interface

end module metis_m
