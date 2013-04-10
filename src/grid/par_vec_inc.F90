!! Copyright (C) 2005-2006 Florian Lorenzen, Heiko Appel
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

!> Generally:
!! Xvec_gather and Xvec_scatter only consider inner points.
!! Xvec_scatter_bndry takes care of boundary points (there is
!! no Xvec_gather_bndry as they are only written and not read).
!! Xvec_scatter_all is Xvec_scatter followd by Xvec_scatter_bndry.

!! ---------------------------------------------------------
!! Scatters a vector v to all nodes in vp with respect to
!! to point -> node mapping in vp.
!! v_local has at least to be of size vp%np_local(vp%partno).
subroutine X(vec_scatter)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(in)  :: v(:)
  R_TYPE,     intent(out) :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Send buffer.

  PUSH_SUB(X(vec_scatter))
  call profiling_in(prof_scatter, "VEC_SCATTER")

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart < 2) then
    v_local(1:vp%np) = v(1:vp%np)
    POP_SUB(X(vec_scatter))
    return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal - 1

  SAFE_ALLOCATE(v_tmp(1:1))
  if(vp%rank == root) then
  ! Fill send buffer.
    SAFE_DEALLOCATE_A(v_tmp)
    SAFE_ALLOCATE(v_tmp(1:vp%np))

    ! Rearrange copy of v. All points of node r are in
    ! v_tmp(xlocal(r):xlocal(r)+np_local(r)-1).
    do i = 1, vp%np
      v_tmp(i) = v(vp%local(i))
    end do
  end if

  ! Careful: MPI rank numbers range from 0 to mpiv%numprocs-1
  ! But partition numbers from 1 to vp%npart with usually
  ! vp%npart = mpiv%numprocs.
  call mpi_debug_in(vp%comm, C_MPI_SCATTERV)
  call MPI_Scatterv(v_tmp(1), vp%np_local, displs(1), R_MPITYPE, v_local(1), &
                    vp%np_local(vp%partno), R_MPITYPE,              &
                    root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_SCATTERV)

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call profiling_out(prof_scatter)
  POP_SUB(X(vec_scatter))

end subroutine X(vec_scatter)


! ---------------------------------------------------------
!> Reverse operation of Xvec_scatter.
!! All v_locals from the nodes are packed together
!! into v on node root in correct order.
subroutine X(vec_gather)(vp, root, v, v_local)
  type(pv_t), intent(in)  :: vp
  integer,    intent(in)  :: root
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Receive buffer.

  PUSH_SUB(X(vec_gather))

  ! Skip the MPI call if domain parallelization is not used.
  if(vp%npart < 2) then
    v(1:vp%np) = v_local(1:vp%np)
    POP_SUB(X(vec_gather))
    return
  end if

  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_GATHERV)
  call MPI_Gatherv(v_local(1), vp%np_local(vp%partno), R_MPITYPE, v_tmp(1), &
                   vp%np_local, displs(1), R_MPITYPE,                        &
                   root, vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_GATHERV)

  ! Copy values from v_tmp to their original position in v.
  if(vp%rank == root) then
    do i = 1, vp%np
      v(vp%local(i)) = v_tmp(i)
    end do

  end if

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  POP_SUB(X(vec_gather))

end subroutine X(vec_gather)

! ---------------------------------------------------------
!> Like Xvec_gather but the result is gathered
!! on all nodes, i. e. v has to be a properly
!! allocated array on all nodes.
subroutine X(vec_allgather)(vp, v, v_local)
  type(pv_t), intent(in)  :: vp
  R_TYPE,     intent(out) :: v(:)
  R_TYPE,     intent(in)  :: v_local(:)

  integer              :: i         !< Counter.
  integer, allocatable :: displs(:) !< Displacements for scatter.
  R_TYPE,  allocatable :: v_tmp(:)  !< Receive buffer.

  PUSH_SUB(X(vec_allgather))
  call profiling_in(prof_allgather, "VEC_ALLGATHER")
  
  ! Unfortunately, vp%xlocal ist not quite the required
  ! displacement vector.
  SAFE_ALLOCATE(displs(1:vp%npart))
  displs = vp%xlocal-1

  SAFE_ALLOCATE(v_tmp(1:vp%np))

  call mpi_debug_in(vp%comm, C_MPI_ALLGATHERV)
  call MPI_Allgatherv(v_local(1), vp%np_local(vp%partno), R_MPITYPE, v_tmp(1), &
                      vp%np_local, displs(1), R_MPITYPE,                       &
                      vp%comm, mpi_err)
  call mpi_debug_out(vp%comm, C_MPI_ALLGATHERV)

  ! Copy values from v_tmp to their original position in v.
  do i = 1, vp%np
    v(vp%local(i)) = v_tmp(i)
  end do

  SAFE_DEALLOCATE_A(v_tmp)
  SAFE_DEALLOCATE_A(displs)

  call profiling_out(prof_allgather)
  POP_SUB(X(vec_allgather))

end subroutine X(vec_allgather)

!--------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
