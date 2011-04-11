!! Copyright (C) 2005-2010 Florian Lorenzen, Heiko Appel, X. Andrade
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

module boundaries_m
  use batch_m
  use global_m
  use messages_m
  use mesh_m
  use mpi_m
  use mpi_debug_m
  use opencl_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use subarray_m

  implicit none
  
  private

  type boundaries_t
    type(mesh_t), pointer :: mesh
    integer          :: nper            !< the number of points that correpond to pbc
    integer, pointer :: per_points(:)   !< (1:nper) the list of points that correspond to pbc 
    integer, pointer :: per_map(:)      !< (1:nper) the inner point that corresponds to each pbc point
#ifdef HAVE_MPI
    integer, pointer :: nsend(:)
    integer, pointer :: nrecv(:)
    integer, pointer :: dsend_type(:)
    integer, pointer :: zsend_type(:)
    integer, pointer :: drecv_type(:)
    integer, pointer :: zrecv_type(:)
#endif
  end type boundaries_t

  public ::                        &
    boundaries_t,                  &
    boundaries_init,               &
    boundaries_end

#if defined(HAVE_MPI)
  public ::                        &
    pv_handle_batch_t,             &
    dvec_ghost_update,             &
    zvec_ghost_update,             &
    dghost_update_batch_start,     &
    zghost_update_batch_start,     &
    dghost_update_batch_finish,    &
    zghost_update_batch_finish

  integer :: SEND = 1, RECV = 2

  type pv_handle_batch_t
    private
    type(batch_t)        :: ghost_send
    integer,     pointer :: requests(:)
    integer              :: nnb
    ! these are needed for CL
    FLOAT, pointer       :: drecv_buffer(:)
    CMPLX, pointer       :: zrecv_buffer(:)
    FLOAT, pointer       :: dsend_buffer(:)
    CMPLX, pointer       :: zsend_buffer(:)
    type(batch_t),   pointer :: v_local
    type(pv_t),      pointer :: vp
  end type pv_handle_batch_t

  type(profile_t), save :: prof_start
  type(profile_t), save :: prof_wait
  type(profile_t), save :: prof_update
  
#endif
contains
  
  ! ---------------------------------------------------------
  subroutine boundaries_init(this, mesh)
    type(boundaries_t),   intent(out)   :: this
    type(mesh_t), target, intent(in)    :: mesh

    integer :: sp, ip, ip_inner, iper, ip_global
#ifdef HAVE_MPI
    integer :: ip_inner_global, ipart, nblocks
    integer, allocatable :: recv_rem_points(:, :)
    integer :: nper_recv
    integer :: maxmax
    integer, allocatable :: recv_points(:, :), blocklengths(:), offsets(:)
    integer, allocatable :: send_points(:, :)
    integer, allocatable :: send_buffer(:)
    integer :: bsize, status(MPI_STATUS_SIZE)
#endif

    PUSH_SUB(boundaries_init)

    this%mesh => mesh

    nullify(this%per_points)
    nullify(this%per_map)

    if (simul_box_is_periodic(mesh%sb)) then

      sp = mesh%np
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) sp = mesh%np + mesh%vp%np_ghost(mesh%vp%partno)
#endif

      !count the number of points that are periodic
      this%nper = 0
#ifdef HAVE_MPI
      nper_recv = 0
#endif
      do ip = sp + 1, mesh%np_part

        ip_global = ip

#ifdef HAVE_MPI
        !translate to a global point
        if(mesh%parallel_in_domains) ip_global = mesh%vp%bndry(ip - sp - 1 + mesh%vp%xbndry(mesh%vp%partno))
#endif

        ip_inner = mesh_periodic_point(mesh, ip_global)

#ifdef HAVE_MPI
        !translate back to a local point
        if(mesh%parallel_in_domains) ip_inner = vec_global2local(mesh%vp, ip_inner, mesh%vp%partno)
#endif
        
        ! If the point is the periodic of another point, is not zero
        ! (this might happen in the parallel case) and is inside the
        ! grid then we have to copy it from the grid points.  
        !
        ! If the point index is larger than mesh%np then it is the
        ! periodic copy of a point that is zero, so we don`t count it
        ! as it will be initialized to zero anyway. For different
        ! mixed boundary conditions the last check should be removed.
        !
        if(ip /= ip_inner .and. ip_inner /= 0 .and. ip_inner <= mesh%np) then 
          this%nper = this%nper + 1
#ifdef HAVE_MPI
        else if(mesh%parallel_in_domains .and. ip /= ip_inner) then
          nper_recv = nper_recv + 1
#endif
        end if
      end do

      SAFE_ALLOCATE(this%per_points(1:this%nper))
      SAFE_ALLOCATE(this%per_map(1:this%nper))

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then
        SAFE_ALLOCATE(recv_points(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(recv_rem_points(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(this%nrecv(1:mesh%vp%npart))
        this%nrecv = 0
      end if
#endif

      iper = 0
      do ip = sp + 1, mesh%np_part

        ip_global = ip

        !translate to a global point
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) ip_global = mesh%vp%bndry(ip - sp - 1 + mesh%vp%xbndry(mesh%vp%partno))
#endif

        ip_inner = mesh_periodic_point(mesh, ip_global)
        
        !translate to local (and keep a copy of the global)
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) then
          ip_inner_global = ip_inner
          ip_inner = vec_global2local(mesh%vp, ip_inner, mesh%vp%partno)
        end if
#endif

        if(ip /= ip_inner .and. ip_inner /= 0 .and. ip_inner <= mesh%np) then
          iper = iper + 1
          this%per_points(iper) = ip
          this%per_map(iper) = ip_inner

#ifdef HAVE_MPI
        else if(mesh%parallel_in_domains .and. ip /= ip_inner) then ! the point is in another node
          ! find in which paritition it is
          do ipart = 1, mesh%vp%npart
            if(ipart == mesh%vp%partno) cycle

            ip_inner = vec_global2local(mesh%vp, ip_inner_global, ipart)
            
            if(ip_inner /= 0) then
              if(ip_inner <= mesh%vp%np_local(ipart)) then
                ! count the points to receive from each node
                this%nrecv(ipart) = this%nrecv(ipart) + 1
                ! and store the number of the point
                recv_points(this%nrecv(ipart), ipart) = ip
                ! and where it is in the other partition
                recv_rem_points(this%nrecv(ipart), ipart) = ip_inner

                ASSERT(mesh%vp%rank /= ipart - 1) ! if we are here, the point must be in another node
              
                exit
              end if
            end if
            
          end do
#endif
        end if
      end do

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then

        ! first we allocate the buffer to be able to use MPI_Bsend
        bsize = mesh%vp%npart - 1 + nper_recv + MPI_BSEND_OVERHEAD*2*(mesh%vp%npart - 1)
        SAFE_ALLOCATE(send_buffer(1:bsize))
        call MPI_Buffer_attach(send_buffer(1), bsize*4, mpi_err)

        ! Now we communicate to each node the points they will have to
        ! send us. Probably this could be done without communication,
        ! but this way it seems simpler to implement.
        
        ! We send the number of points we expect to receive.
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Bsend(this%nrecv(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, mpi_err)
        end do

        ! And we receive it
        SAFE_ALLOCATE(this%nsend(1:mesh%vp%npart))
        this%nsend = 0
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Recv(this%nsend(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, status, mpi_err)
        end do

        ! Now we send the indices of the points
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. this%nrecv(ipart) == 0) cycle
          call MPI_Bsend(recv_rem_points(1, ipart), this%nrecv(ipart), MPI_INTEGER, ipart - 1, 1, mesh%vp%comm, mpi_err)
        end do

        SAFE_ALLOCATE(send_points(1:maxval(this%nsend), 1:mesh%vp%npart))

        ! And we receive them
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. this%nsend(ipart) == 0) cycle
          call MPI_Recv(send_points(1, ipart), this%nsend(ipart), MPI_INTEGER, &
               ipart - 1, 1, mesh%vp%comm, status, mpi_err)
        end do

        ! we no longer need this
        SAFE_DEALLOCATE_A(recv_rem_points)

        ! Now we have all the indices required locally, so we can
        ! build the mpi datatypes

        SAFE_ALLOCATE(this%dsend_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(this%zsend_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(this%drecv_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(this%zrecv_type(1:mesh%vp%npart))

        maxmax = max(maxval(this%nsend), maxval(this%nrecv))

        SAFE_ALLOCATE(blocklengths(1:maxmax))
        SAFE_ALLOCATE(offsets(1:maxmax))

        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle

          if(this%nsend(ipart) > 0) then

            ASSERT(all(send_points(1:this%nsend(ipart), ipart) <= mesh%np))

            ! MPI indices start from zero
            send_points(1:this%nsend(ipart), ipart) = send_points(1:this%nsend(ipart), ipart) - 1

            call get_blocks(this%nsend(ipart), send_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_FLOAT, this%dsend_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_CMPLX, this%zsend_type(ipart), mpi_err)
            call MPI_Type_commit(this%dsend_type(ipart), mpi_err)
            call MPI_Type_commit(this%zsend_type(ipart), mpi_err)

          end if
          
          if(this%nrecv(ipart) > 0) then
            ASSERT(all(recv_points(1:this%nrecv(ipart), ipart) <= mesh%np_part))
            ASSERT(all(recv_points(1:this%nrecv(ipart), ipart) > mesh%np))

            recv_points(1:this%nrecv(ipart), ipart) = recv_points(1:this%nrecv(ipart), ipart) - 1

            call get_blocks(this%nrecv(ipart), recv_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_FLOAT, this%drecv_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_CMPLX, this%zrecv_type(ipart), mpi_err)
            call MPI_Type_commit(this%drecv_type(ipart), mpi_err)
            call MPI_Type_commit(this%zrecv_type(ipart), mpi_err)

          end if

        end do

        call MPI_Buffer_detach(send_buffer(1), bsize, mpi_err)

      end if
#endif
      ASSERT(iper == this%nper)

    end if

    POP_SUB(boundaries_init)
  end subroutine boundaries_init

  ! ---------------------------------------------------------
  
  subroutine boundaries_end(this)
    type(boundaries_t),  intent(out)   :: this

#ifdef HAVE_MPI
    integer :: ipart
#endif

    PUSH_SUB(boundaries_end)

    if(simul_box_is_periodic(this%mesh%sb)) then
#ifdef HAVE_MPI
      if(this%mesh%parallel_in_domains) then    
        
        ASSERT(associated(this%nsend))
        ASSERT(associated(this%nrecv))
        
        do ipart = 1, this%mesh%vp%npart
          if(this%nsend(ipart) /= 0) then 
            call MPI_Type_free(this%dsend_type(ipart), mpi_err)
            call MPI_Type_free(this%zsend_type(ipart), mpi_err)
          end if
          if(this%nrecv(ipart) /= 0) then 
            call MPI_Type_free(this%drecv_type(ipart), mpi_err)
            call MPI_Type_free(this%zrecv_type(ipart), mpi_err)
          end if
        end do
        
        SAFE_DEALLOCATE_P(this%dsend_type)
        SAFE_DEALLOCATE_P(this%zsend_type)
        SAFE_DEALLOCATE_P(this%drecv_type)
        SAFE_DEALLOCATE_P(this%zrecv_type)
        SAFE_DEALLOCATE_P(this%nsend)
        SAFE_DEALLOCATE_P(this%nrecv)
      end if
#endif
      
      SAFE_DEALLOCATE_P(this%per_points)
      SAFE_DEALLOCATE_P(this%per_map)
    end if

    POP_SUB(boundaries_end)
  end subroutine boundaries_end

#if defined(HAVE_MPI)

#include "undef.F90"
#include "complex.F90"
#include "boundaries_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "boundaries_inc.F90"

#endif
end module boundaries_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
