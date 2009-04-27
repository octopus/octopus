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

subroutine X(dsubmesh_add_to_mesh)(this, sphi, phi, factor)
  type(submesh_t),  intent(in)    :: this
  FLOAT,            intent(in)    :: sphi(:)
  R_TYPE,           intent(inout) :: phi(:)
  R_TYPE, optional, intent(in)    :: factor

  integer :: is

  if(present(factor)) then
    forall(is = 1:this%ns) phi(this%jxyz(is)) = phi(this%jxyz(is)) + factor*sphi(is)
  else
    forall(is = 1:this%ns) phi(this%jxyz(is)) = phi(this%jxyz(is)) + sphi(is)
  end if

end subroutine X(dsubmesh_add_to_mesh)

R_TYPE function X(dsubmesh_to_mesh_dotp)(this, dim, sphi, phi) result(dotp)
  type(submesh_t),  intent(in)  :: this
  integer,          intent(in)  :: dim
  FLOAT,            intent(in)  :: sphi(:)
  R_TYPE,           intent(in)  :: phi(:, :)

  integer :: is, idim

  dotp = R_TOTYPE(M_ZERO)

  if(this%mesh%use_curvilinear) then
    do idim = 1, dim
      do is = 1, this%ns
        dotp = dotp + this%mesh%vol_pp(this%jxyz(is))*phi(this%jxyz(is), idim)*sphi(is)
      end do
    end do
  else
    do idim = 1, dim
      do is = 1, this%ns
        dotp = dotp + phi(this%jxyz(is), idim)*sphi(is)
      end do
    end do
    dotp = dotp*this%mesh%vol_pp(1)
  end if

end function X(dsubmesh_to_mesh_dotp)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
