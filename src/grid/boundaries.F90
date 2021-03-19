!! Copyright (C) 2005-2010 Florian Lorenzen, Heiko Appel, X. Andrade
!! Copyright (C) 2021 Sebastian Ohlmann
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

module boundaries_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use global_oct_m
  use index_oct_m
  use math_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use multiresolution_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use partition_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none
  
  private

  type boundaries_t
    private
    type(mesh_t), pointer :: mesh
    logical              :: periodic = .false.       !< some boundaries are to be treated periodic
    logical              :: fully_periodic = .false. !< all boundaries are to be treated periodic
    integer              :: nper             !< the number of points that correspond to pbc
    integer, allocatable :: per_points(:, :) !< (1:2, 1:nper) the list of points that correspond to pbc 
    integer, allocatable :: per_send(:, :)
    integer, allocatable :: per_recv(:, :)
    integer, allocatable :: nsend(:)
    integer, allocatable :: nrecv(:)
    type(accel_mem_t)    :: buff_per_points
    type(accel_mem_t)    :: buff_per_send
    type(accel_mem_t)    :: buff_per_recv
    type(accel_mem_t)    :: buff_nsend
    type(accel_mem_t)    :: buff_nrecv
    logical, public      :: spiralBC           !< set .true. when SpiralBoundaryCondition are set in the input file
    logical, public      :: spiral             !< set .true. after first time step IF spiralBC == .true. (see td_run in td.F90)
    FLOAT,   public      :: spiral_q(MAX_DIM)
  end type boundaries_t

  public ::                        &
    boundaries_t,                  &
    boundaries_nullify,            &
    boundaries_init,               &
    boundaries_end,                &
    boundaries_set

  public ::                        &
    pv_handle_batch_t,             &
    dvec_ghost_update,             &
    zvec_ghost_update,             &
    dghost_update_batch_start,     &
    zghost_update_batch_start,     &
    dghost_update_batch_finish,    &
    zghost_update_batch_finish

  integer, parameter, public ::    &
    POINT_BOUNDARY = 1,            &
    POINT_INNER    = 2

  type pv_handle_batch_t
    private
    type(batch_t)        :: ghost_send   !< batch for sending data; it is packed into this one
    integer, allocatable :: requests(:)
    integer              :: nnb
    ! these are needed for CL
    FLOAT, pointer       :: drecv_buffer(:)
    CMPLX, pointer       :: zrecv_buffer(:)
    FLOAT, pointer       :: dsend_buffer(:)
    CMPLX, pointer       :: zsend_buffer(:)
    type(batch_t),   pointer :: v_local
    type(pv_t),      pointer :: vp
  end type pv_handle_batch_t

  interface boundaries_set
    module procedure boundaries_set_batch
    module procedure dboundaries_set_single
    module procedure zboundaries_set_single
  end interface boundaries_set

contains
  
  ! ---------------------------------------------------------
  elemental subroutine boundaries_nullify(this)
    type(boundaries_t), intent(out) :: this

    nullify(this%mesh)
    this%nper = 0
    call accel_mem_nullify(this%buff_per_points)
    call accel_mem_nullify(this%buff_per_send)
    call accel_mem_nullify(this%buff_per_recv)
    call accel_mem_nullify(this%buff_nsend)
    call accel_mem_nullify(this%buff_nrecv)
    this%spiralBC = .false.
    this%spiral = .false.
    this%spiral_q(1:MAX_DIM) = M_ZERO

  end subroutine boundaries_nullify

  ! ---------------------------------------------------------
  subroutine boundaries_init(this, namespace, space, mesh)
    type(boundaries_t),   intent(inout) :: this
    type(namespace_t),    intent(in)    :: namespace
    type(space_t),        intent(in)    :: space
    type(mesh_t), target, intent(in)    :: mesh

    integer :: sp, ip, ip_inner, iper, idir
    integer :: ip_inner_global, ipart
    integer, allocatable :: recv_rem_points(:, :), points(:), part(:), points_local(:)
    integer :: nper_recv, iper_recv
#ifdef HAVE_MPI
    integer, allocatable :: send_buffer(:)
    integer :: bsize, status(MPI_STATUS_SIZE)
#endif
    type(block_t) :: blk

    PUSH_SUB(boundaries_init)

    this%mesh => mesh

    this%periodic = space%is_periodic()
    this%fully_periodic = space%periodic_dim == space%dim

    if (space%is_periodic()) then

      !%Variable SpiralBoundaryCondition
      !%Type logical
      !%Default no
      !%Section Mesh
      !%Description
      !% (Experimental) If set to yes, Octopus will apply spin-spiral boundary conditions.
      !% The momentum of the spin spiral is defined by the variable 
      !% <tt>TDMomentumTransfer</tt> 
      !%End
      call parse_variable(namespace, 'SpiralBoundaryCondition', .false., this%spiralBC)
      if(this%spiralBC) then
        if(parse_is_defined(namespace, 'TDMomentumTransfer')) then
          if(parse_block(namespace, 'TDMomentumTransfer', blk)==0) then
            do idir = 1, MAX_DIM
             call parse_block_float(blk, 0, idir - 1, this%spiral_q(idir))
             this%spiral_q(idir) = units_to_atomic(unit_one / units_inp%length, this%spiral_q(idir))
            end do
            call messages_experimental("SpiralBoundaryCondition")
          else 
            message(1) = "TDMomentumTransfer must be defined if SpiralBoundaryCondition=yes"
            call messages_fatal(1, namespace=namespace)
          end if
        else
         message(1) = "TDMomentumTransfer must be defined if SpiralBoundaryCondition=yes"
         call messages_fatal(1, namespace=namespace) 
        end if
      end if

      sp = mesh%np
      if(mesh%parallel_in_domains) sp = mesh%np + mesh%vp%np_ghost

      ! count the number of points that are periodic
      this%nper = 0
      nper_recv = 0
      do ip = sp + 1, mesh%np_part
        ip_inner_global = mesh_periodic_point(mesh, space, ip)
        ip_inner = mesh_global2local(mesh, ip_inner_global)

        ! it is the same point, can happen for mixed periodicity
        if (ip == ip_inner) cycle
        ! the point maps to the boundary, can happen for mixed periodicity
        ! in this case the point is already set to zero, so we can ignore it
        ! for different mixed boundary conditions, we would need to be careful here
        if (ip_inner_global > mesh%np_global) cycle
        ! now check if point is local or if it needs to be communicated
        if (ip_inner /= 0 .and. ip_inner <= mesh%np) then
          this%nper = this%nper + 1
        else
          nper_recv = nper_recv + 1
        end if
      end do
      if(.not.mesh%parallel_in_domains) then
        ASSERT(nper_recv == 0)
      end if

      SAFE_ALLOCATE(this%per_points(1:2, 1:this%nper))

      if(mesh%parallel_in_domains) then
        SAFE_ALLOCATE(this%per_recv(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(this%nrecv(1:mesh%vp%npart))
        SAFE_ALLOCATE(recv_rem_points(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(points(1:nper_recv))
        SAFE_ALLOCATE(points_local(1:nper_recv))
        SAFE_ALLOCATE(part(1:nper_recv))
        this%nrecv = 0
      end if

      iper = 0
      iper_recv = 0
      do ip = sp + 1, mesh%np_part
        ip_inner_global = mesh_periodic_point(mesh, space, ip)
        ip_inner = mesh_global2local(mesh, ip_inner_global)
        
        ! it is the same point, can happen for mixed periodicity
        if (ip == ip_inner) cycle
        ! the point maps to the boundary, can happen for mixed periodicity
        ! in this case the point is already set to zero, so we can ignore it
        if (ip_inner_global > mesh%np_global) cycle
        ! now check if point is local or if it needs to be communicated
        if (ip_inner /= 0 .and. ip_inner <= mesh%np) then
          iper = iper + 1
          this%per_points(POINT_BOUNDARY, iper) = ip
          this%per_points(POINT_INNER, iper) = ip_inner
        else
          ! this can only happen if parallel in domain
          ! the point is on another node
          iper_recv = iper_recv + 1
          points(iper_recv) = ip_inner_global
          points_local(iper_recv) = ip
        end if
      end do
      if(mesh%parallel_in_domains) then
        ! find the points in the other partitions
        call partition_get_partition_number(mesh%partition, nper_recv, &
          points, part)
        do iper_recv = 1, nper_recv
          ipart = part(iper_recv)
          ! count the points to receive from each node
          this%nrecv(ipart) = this%nrecv(ipart) + 1
          ! and store the number of the point
          this%per_recv(this%nrecv(ipart), ipart) = points_local(iper_recv)
          ! and its global index
          recv_rem_points(this%nrecv(ipart), ipart) = points(iper_recv)
        end do
      end if

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
          call MPI_Bsend(this%nrecv(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%mpi_grp%comm, mpi_err)
        end do

        ! And we receive it
        SAFE_ALLOCATE(this%nsend(1:mesh%vp%npart))
        this%nsend = 0
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Recv(this%nsend(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%mpi_grp%comm, status, mpi_err)
        end do

        ! Now we send the indices of the points
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. this%nrecv(ipart) == 0) cycle
          call MPI_Bsend(recv_rem_points(1, ipart), this%nrecv(ipart), MPI_INTEGER, ipart - 1, 1, mesh%mpi_grp%comm, mpi_err)
        end do

        SAFE_ALLOCATE(this%per_send(1:maxval(this%nsend), 1:mesh%vp%npart))

        ! And we receive them
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. this%nsend(ipart) == 0) cycle
          call MPI_Recv(this%per_send(1, ipart), this%nsend(ipart), MPI_INTEGER, &
               ipart - 1, 1, mesh%mpi_grp%comm, status, mpi_err)
          ! get local index of the points to send
          do ip = 1, this%nsend(ipart)
            this%per_send(ip, ipart) = mesh_global2local(mesh, this%per_send(ip, ipart))
            ! make sure we have local points here
            ASSERT(this%per_send(ip, ipart) > 0)
            ASSERT(this%per_send(ip, ipart) <= mesh%np)
          end do
        end do

        ! we no longer need this
        SAFE_DEALLOCATE_A(recv_rem_points)

        call MPI_Buffer_detach(send_buffer(1), bsize, mpi_err)
        SAFE_DEALLOCATE_A(send_buffer)
        
        SAFE_DEALLOCATE_A(recv_rem_points)
        SAFE_DEALLOCATE_A(points)
        SAFE_DEALLOCATE_A(points_local)
        SAFE_DEALLOCATE_A(part)
      end if
#endif

      if(accel_is_enabled()) then
        call accel_create_buffer(this%buff_per_points, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, 2*this%nper)
        call accel_write_buffer(this%buff_per_points, 2*this%nper, this%per_points)

        if(mesh%parallel_in_domains) then
          call accel_create_buffer(this%buff_per_send, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, product(ubound(this%per_send)))
          call accel_write_buffer(this%buff_per_send, product(ubound(this%per_send)), this%per_send)

          call accel_create_buffer(this%buff_per_recv, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, product(ubound(this%per_recv)))
          call accel_write_buffer(this%buff_per_recv, product(ubound(this%per_recv)), this%per_recv)

          call accel_create_buffer(this%buff_nsend, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, mesh%vp%npart)
          call accel_write_buffer(this%buff_nsend, mesh%vp%npart, this%nsend)

          call accel_create_buffer(this%buff_nrecv, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, mesh%vp%npart)
          call accel_write_buffer(this%buff_nrecv, mesh%vp%npart, this%nrecv)
        end if
      end if

    end if

    POP_SUB(boundaries_init)
  end subroutine boundaries_init

  ! ---------------------------------------------------------
  
  subroutine boundaries_end(this)
    type(boundaries_t),  intent(inout) :: this

    PUSH_SUB(boundaries_end)

    if (this%periodic) then
      if(this%mesh%parallel_in_domains) then    
        
        ASSERT(allocated(this%nsend))
        ASSERT(allocated(this%nrecv))
        
        SAFE_DEALLOCATE_A(this%per_send)
        SAFE_DEALLOCATE_A(this%per_recv)
        SAFE_DEALLOCATE_A(this%nsend)
        SAFE_DEALLOCATE_A(this%nrecv)

        if(accel_is_enabled()) then
          call accel_release_buffer(this%buff_per_send)
          call accel_release_buffer(this%buff_per_recv)
          call accel_release_buffer(this%buff_nsend)
          call accel_release_buffer(this%buff_nrecv)
        end if
      end if

      if(accel_is_enabled()) call accel_release_buffer(this%buff_per_points)

      SAFE_DEALLOCATE_A(this%per_points)
    end if

    POP_SUB(boundaries_end)
  end subroutine boundaries_end

  ! -------------------------------------------------------

  subroutine boundaries_set_batch(this, ffb, phase_correction)
    type(boundaries_t), intent(in)    :: this
    class(batch_t),     intent(inout) :: ffb
    CMPLX, optional,    intent(in)    :: phase_correction(:)

    PUSH_SUB(boundaries_set_batch)
    
    if(ffb%type() == TYPE_FLOAT) then 
      call dboundaries_set_batch(this, ffb, phase_correction)
    else if(ffb%type() == TYPE_CMPLX) then 
      call zboundaries_set_batch(this, ffb, phase_correction)
    else
      ASSERT(.false.)
     end if
     
     POP_SUB(boundaries_set_batch)
   end subroutine boundaries_set_batch

#include "undef.F90"
#include "complex.F90"
#include "boundaries_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "boundaries_inc.F90"

end module boundaries_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
