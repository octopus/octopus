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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module states_elec_parallel_oct_m
  use batch_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use states_abst_oct_m
  use states_elec_oct_m

  implicit none

  private

  public ::                                    &
    states_elec_parallel_blacs_blocksize,           &
    states_elec_parallel_remote_access_start,       &
    states_elec_parallel_remote_access_stop,        &
    states_elec_parallel_get_block,                 &
    states_elec_parallel_release_block,             &
    states_elec_parallel_gather

  interface states_elec_parallel_gather
    module procedure dstates_elec_parallel_gather_1, zstates_elec_parallel_gather_1
    module procedure dstates_elec_parallel_gather_3, zstates_elec_parallel_gather_3
  end interface states_elec_parallel_gather

  type(profile_t), save:: prof_gather
  
contains

  subroutine states_elec_parallel_blacs_blocksize(st, namespace, mesh, blocksize, total_np)
    type(states_elec_t),  intent(in)    :: st
    type(namespace_t),    intent(in)    :: namespace
    type(mesh_t),         intent(in)    :: mesh
    integer,              intent(out)   :: blocksize(2)
    integer,              intent(out)   :: total_np

    PUSH_SUB(states_elec_parallel_blacs_blocksize)

#ifdef HAVE_SCALAPACK
    ! We need to select the block size of the decomposition. This is
    ! tricky, since not all processors have the same number of
    ! points.
    !
    ! What we do for now is to use the maximum of the number of
    ! points and we set to zero the remaining points.

    if(.not. st%scalapack_compatible) then
      message(1) = "Attempt to use ScaLAPACK when processes have not been distributed in compatible layout."
      message(2) = "You need to set ScaLAPACKCompatible = yes in the input file and re-run."
      call messages_fatal(2, only_root_writes = .true., namespace=namespace)
    end if
    
    if (mesh%parallel_in_domains) then
      blocksize(1) = maxval(mesh%vp%np_local_vec) + (st%d%dim - 1) * &
       maxval(mesh%vp%np_local_vec + mesh%vp%np_bndry + mesh%vp%np_ghost)
    else
      blocksize(1) = mesh%np + (st%d%dim - 1)*mesh%np_part
    end if

    if (st%parallel_in_states) then
      blocksize(2) = maxval(st%dist%num)
    else
      blocksize(2) = st%nst
    end if

    total_np = blocksize(1)*st%dom_st_proc_grid%nprow


    ASSERT(st%d%dim*mesh%np_part >= blocksize(1))
#else
    blocksize(1) = 0
    blocksize(2) = 0
    total_np = 0
#endif

    POP_SUB(states_elec_parallel_blacs_blocksize)
  end subroutine states_elec_parallel_blacs_blocksize

  ! ------------------------------------------------------------

  subroutine states_elec_parallel_remote_access_start(this)
    type(states_elec_t),       intent(inout) :: this
    
    integer :: ib, iqn
    
    PUSH_SUB(states_elec_parallel_remote_access_start)

    ASSERT(associated(this%group%psib))

    SAFE_ALLOCATE(this%group%rma_win(1:this%group%nblocks, 1:this%d%nik))
    
    do iqn = this%d%kpt%start, this%d%kpt%end
      do ib = 1, this%group%nblocks
        if(this%group%block_is_local(ib, iqn)) then
          call this%group%psib(ib, iqn)%remote_access_start(this%mpi_grp, this%group%rma_win(ib, iqn))
        else
#ifdef HAVE_MPI2
          ! create an empty window
          call MPI_Win_create(0, int(0, MPI_ADDRESS_KIND), 1, &
            MPI_INFO_NULL, this%mpi_grp%comm, this%group%rma_win(ib, iqn), mpi_err)
#endif
        end if
      end do
    end do
    
    POP_SUB(states_elec_parallel_remote_access_start)
  end subroutine states_elec_parallel_remote_access_start

    ! ---------------------------------------------------------

  subroutine states_elec_parallel_remote_access_stop(this)
    type(states_elec_t),       intent(inout) :: this
    
    integer :: ib, iqn
    
    PUSH_SUB(states_elec_parallel_remote_access_stop)

    ASSERT(associated(this%group%psib))
    
    do iqn = this%d%kpt%start, this%d%kpt%end
      do ib = 1, this%group%nblocks
        if(this%group%block_is_local(ib, iqn)) then
          call this%group%psib(ib, iqn)%remote_access_stop(this%group%rma_win(ib, iqn))
        else
#ifdef HAVE_MPI2
          call MPI_Win_free(this%group%rma_win(ib, iqn), mpi_err)
#endif
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(this%group%rma_win)
    
    POP_SUB(states_elec_parallel_remote_access_stop)
  end subroutine states_elec_parallel_remote_access_stop

  ! --------------------------------------

  subroutine states_elec_parallel_get_block(this, mesh, ib, iqn, psib)
    type(states_elec_t), target, intent(in) :: this
    type(mesh_t),                intent(in) :: mesh
    integer,                     intent(in) :: ib
    integer,                     intent(in) :: iqn
    type(batch_t),               pointer    :: psib

    type(profile_t), save :: prof
    
    PUSH_SUB(states_elec_parallel_get_block)

    call profiling_in(prof, "STATES_GET_BLOCK")
    
    if(this%group%block_is_local(ib, iqn)) then
      psib => this%group%psib(ib, iqn)
    else
      SAFE_ALLOCATE(psib)
      call batch_init(psib, this%d%dim, this%group%block_size(ib))

      if(states_are_real(this)) then
        call psib%dallocate(this%group%block_range(ib, 1), this%group%block_range(ib, 2), mesh%np_part)
      else
        call psib%zallocate(this%group%block_range(ib, 1), this%group%block_range(ib, 2), mesh%np_part)
      end if
      
      call psib%do_pack(copy = .false.)
      
#ifdef HAVE_MPI2
      call MPI_Win_lock(MPI_LOCK_SHARED, this%group%block_node(ib), 0, this%group%rma_win(ib, iqn),  mpi_err)

      if(states_are_real(this)) then
        call MPI_Get(psib%pack%dpsi(1, 1), product(psib%pack%size), MPI_FLOAT, &
          this%group%block_node(ib), int(0, MPI_ADDRESS_KIND), product(psib%pack%size), MPI_FLOAT, &
          this%group%rma_win(ib, iqn), mpi_err)
      else
        call MPI_Get(psib%pack%zpsi(1, 1), product(psib%pack%size), MPI_CMPLX, &
          this%group%block_node(ib), int(0, MPI_ADDRESS_KIND), product(psib%pack%size), MPI_CMPLX, &
          this%group%rma_win(ib, iqn), mpi_err)
      end if
        
      call MPI_Win_unlock(this%group%block_node(ib), this%group%rma_win(ib, iqn),  mpi_err)
#endif
    end if

    call profiling_out(prof)
    
    POP_SUB(states_elec_parallel_get_block)
  end subroutine states_elec_parallel_get_block

  ! --------------------------------------

  subroutine states_elec_parallel_release_block(this, ib, iqn, psib)
    type(states_elec_t), target, intent(in) :: this
    integer,                     intent(in) :: ib
    integer,                     intent(in) :: iqn
    type(batch_t),               pointer    :: psib

    PUSH_SUB(states_elec_parallel_release_block)

    if(this%group%block_is_local(ib, iqn)) then
      nullify(psib)
    else
      call psib%end
      SAFE_DEALLOCATE_P(psib)
    end if
    
    POP_SUB(states_elec_parallel_release_block)
  end subroutine states_elec_parallel_release_block

  ! --------------------------------------
  
#include "undef.F90"
#include "real.F90"
#include "states_elec_parallel_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_elec_parallel_inc.F90"
#include "undef.F90"
  
  
end module states_elec_parallel_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
