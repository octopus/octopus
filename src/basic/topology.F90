!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!> Given a number n_index of indices and their ranges index_range(1:n_index),
!! we divide the n_nodes in groups and create communicators for each group.
!! Each group is associated with one index. The min_range indicates the minimum
!! number of elements in each processor. For example, given min_range(1) = 25000,
!! the algorithm will try to put at least 25000 points in each processor
!!
!! Example:
!! Given 3 indices with ranges
!!   index_range(1) = 100000, (number of points in the mesh)
!!   index_range(2) = 15,     (number of states)
!!   index_range(3) = 2,      (number of k-points)
!! and 12 processors, we could get
!!   mc%group_sizes = (2, 3, 2)
!! which means that
!! * space is divided in 2 domains per state,
!! * the states are divided in 3 groups, i.e. 5 states per processor, and
!! * the whole setting is duplicated because of the 2 k-points.
!!
!! To perform collective operations (like a reduce), you can use the communicators
!! provided in mc%group_comm(:). For example, to sum over states, the communicator
!! to use is mc%group_comm(P_STRATEGY_STATES)
!!
!! You can use the routine multicomm_strategy_is_parallel to know if a certain
!! index is parallelized.

  module topology_m
    use datasets_m
    use global_m
    use io_m
    use loct_m
    use messages_m
    use mpi_m
    use parser_m
    use profiling_m
    use utils_m

    implicit none

    private

    public ::                          &
      topology_t,                      &
      topology_init,                   &
      topology_groups_are_equal,       &
      topology_end

    type topology_t
      integer          :: ng
      integer          :: maxgsize
      integer, pointer :: distance(:, :)
      integer, pointer :: groups(:,:)
      integer, pointer :: gsize(:)
    end type topology_t

  contains
    !-----------------------------------------------------------------
    !> this routine tries to guess the distribution of the processors
    !! we got, currently only checks processes that are running in the
    !! same node

    subroutine topology_init(this, base_grp)
      type(topology_t), intent(out)   :: this
      type(mpi_grp_t),  intent(inout) :: base_grp

#ifdef HAVE_MPI
      character(len=256) :: my_name, its_name
      integer :: ir,  wsize, ig

      PUSH_SUB(topology_init)

      wsize = base_grp%size

      !get the system name
      call loct_sysname(my_name)

      SAFE_ALLOCATE(this%distance(1:wsize, 1:wsize))

      do ir = 1, wsize

        if(ir - 1 == base_grp%rank) then

          call MPI_Bcast(my_name, 256, MPI_CHARACTER, ir - 1, base_grp%comm, mpi_err)

          this%distance(ir, base_grp%rank + 1) = 0

        else

          call MPI_Bcast(its_name, 256, MPI_CHARACTER, ir - 1, base_grp%comm, mpi_err)

          if(my_name == its_name) then
            this%distance(ir, base_grp%rank + 1) = 1
          else
            this%distance(ir, base_grp%rank + 1) = 2
          end if

        end if

      end do

      do ir = 1, wsize
        call MPI_Bcast(this%distance(1, ir), wsize, MPI_INTEGER, ir - 1, base_grp%comm, mpi_err)
      end do

      !classify processors in groups

      SAFE_ALLOCATE(this%groups(1:wsize, 1:wsize))
      SAFE_ALLOCATE(this%gsize(1:wsize))

      ! put the first node in the first group
      this%ng = 1
      this%gsize = 0
      this%groups = 0
      this%groups(1, 1) = 1
      this%gsize(1) = 1

      ! check the other processors
      do ir = 2, wsize
        do ig = 1, wsize
          ! if the group has elements
          if(this%gsize(ig) > 0) then
            ! check the distance
            if(this%distance(ir, this%groups(1, ig)) == 1) then
              ! if it is close enough
              ! add it to the group
              this%gsize(ig) = this%gsize(ig) + 1
              this%groups(this%gsize(ig), ig) = ir
              exit
            end if
          else
            ! if the group is empty 
            ! set this processor as the head of a group
            this%gsize(ig) = this%gsize(ig) + 1
            this%groups(1, ig) = ir
            this%ng = this%ng + 1
            exit
          end if
        end do
      end do

      ASSERT(sum(this%gsize(1:this%ng)) == base_grp%size)

      this%maxgsize = maxval(this%gsize(1:this%ng))

      !convert to mpi ranks
      this%groups = this%groups - 1

      POP_SUB(topology_init)
#endif
    end subroutine topology_init


    ! ---------------------------------------------------------
    logical function topology_groups_are_equal(this) result(are_equal)
      type(topology_t), intent(in) :: this

      PUSH_SUB(topology_groups_are_equal)
      are_equal = all(this%gsize(2:this%ng) == this%maxgsize)

      POP_SUB(topology_groups_are_equal)
    end function topology_groups_are_equal


    ! ---------------------------------------------------------
    subroutine topology_end(this)
      type(topology_t), intent(inout) :: this

      PUSH_SUB(topology_end)

      SAFE_DEALLOCATE_P(this%groups)
      SAFE_DEALLOCATE_P(this%distance)
      SAFE_DEALLOCATE_P(this%gsize)
      
      POP_SUB(topology_end)
    end subroutine topology_end

end module topology_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
