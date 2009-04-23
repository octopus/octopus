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
!! $Id: submesh_inc.F90 2781 2007-03-23 10:58:32Z lorenzen $

R_TYPE function X(sm_integrate)(m, sm, f) result(res)
  type(mesh_t),      intent(in) :: m
  type(submesh_t),   intent(in) :: sm
  R_TYPE,            intent(in) :: f(:)

#if defined(HAVE_MPI)
  R_TYPE :: tmp
#endif

  if (m%use_curvilinear) then
    res = sum(f(1:sm%ns)*m%vol_pp(sm%jxyz(1:sm%ns)) )
  else
    res = sum(f(1:sm%ns))*m%vol_pp(1)
  end if

#if defined(HAVE_MPI)
  if(m%parallel_in_domains) then
    call MPI_Allreduce(res, tmp, 1, R_MPITYPE, MPI_SUM, m%vp%comm, mpi_err)
    res = tmp
  end if
#endif

end function X(sm_integrate)

#ifdef HAVE_MPI

subroutine X(submesh_comm_reduce)(this, sm, mesh, count, val)
  type(submesh_comm_t), intent(out)   :: this
  type(submesh_t),      intent(in)    :: sm
  type(mesh_t),         intent(in)    :: mesh
  integer,              intent(in)    :: count
  R_TYPE,               intent(inout) :: val(:)

  integer :: ipart, tag

  ASSERT(mesh%parallel_in_domains)

  tag = tagcounter
  tagcounter = tagcounter + 1

  SAFE_ALLOCATE(this%X(allval)(1:count, 1:mesh%vp%npart))

  this%X(allval)(1:count, 1:mesh%vp%npart) = M_ZERO

  SAFE_ALLOCATE(this%requests(1:2*mesh%vp%npart))

  this%nreq = 0

  do ipart = 1, mesh%vp%npart
    if(ipart == mesh%vp%partno .or. sm%psize(ipart) == 0) cycle
    this%nreq = this%nreq + 1
    call MPI_Irecv(this%X(allval)(:, ipart), count, R_MPITYPE, ipart - 1, tag, &
         mesh%mpi_grp%comm, this%requests(this%nreq), mpi_err)
  end do

  if(sm%has_points) then
    do ipart = 1, mesh%vp%npart
      if(ipart == mesh%vp%partno) cycle
      this%nreq = this%nreq + 1
      call MPI_Isend(val, count, R_MPITYPE, ipart - 1, tag, mesh%mpi_grp%comm, this%requests(this%nreq), mpi_err)
    end do
  end if

  this%X(allval)(1:count, mesh%vp%partno) = val(1:count)

end subroutine X(submesh_comm_reduce)

subroutine X(submesh_comm_finish)(this, sm, mesh, count, val)
  type(submesh_comm_t), intent(inout) :: this
  type(submesh_t),      intent(in)    :: sm
  type(mesh_t),         intent(in)    :: mesh
  integer,              intent(in)    :: count
  R_TYPE,               intent(inout) :: val(:)

  integer, allocatable :: statuses(:, :)
  integer :: ipart

  SAFE_ALLOCATE(statuses(1:MPI_STATUS_SIZE, 1:this%nreq))

  call MPI_Waitall(this%nreq, this%requests, statuses, mpi_err)

  val(1:count) = M_ZERO
  do ipart = 1, mesh%vp%npart
    val(1:count) = val(1:count) + this%X(allval)(1:count, ipart)
  end do
  
  SAFE_DEALLOCATE_A(statuses)
  SAFE_DEALLOCATE_P(this%X(allval))
  SAFE_DEALLOCATE_P(this%requests)

end subroutine X(submesh_comm_finish)
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
