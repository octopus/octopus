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

R_TYPE function X(sm_integrate)(mesh, sm, ff) result(res)
  type(mesh_t),      intent(in) :: mesh
  type(submesh_t),   intent(in) :: sm
  R_TYPE, optional,  intent(in) :: ff(:)

  PUSH_SUB(X(sm_integrate))

  ASSERT(present(ff) .or. sm%np .eq. 0)

  if(sm%np > 0) then
    if (mesh%use_curvilinear) then
      res = sum(ff(1:sm%np)*mesh%vol_pp(sm%map(1:sm%np)) )
    else
      res = sum(ff(1:sm%np))*mesh%vol_pp(1)
    end if
  else
    res = M_ZERO
  endif

  if(mesh%parallel_in_domains) call comm_allreduce(mesh%vp%comm, res)

  POP_SUB(X(sm_integrate))
end function X(sm_integrate)

!------------------------------------------------------------

subroutine X(dsubmesh_add_to_mesh)(this, sphi, phi, factor)
  type(submesh_t),  intent(in)    :: this
  FLOAT,            intent(in)    :: sphi(:)
  R_TYPE,           intent(inout) :: phi(:)
  R_TYPE, optional, intent(in)    :: factor

  integer :: is

  PUSH_SUB(X(dsubmesh_add_to_mesh))

  if(present(factor)) then
    forall(is = 1:this%np) phi(this%map(is)) = phi(this%map(is)) + factor*sphi(is)
  else
    forall(is = 1:this%np) phi(this%map(is)) = phi(this%map(is)) + sphi(is)
  end if

  POP_SUB(X(dsubmesh_add_to_mesh))
end subroutine X(dsubmesh_add_to_mesh)

!------------------------------------------------------------

R_TYPE function X(dsubmesh_to_mesh_dotp)(this, dim, sphi, phi, reduce) result(dotp)
  type(submesh_t),   intent(in) :: this
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: sphi(:)
  R_TYPE,            intent(in) :: phi(:, :)
  logical, optional, intent(in) :: reduce

  integer :: is, idim

  PUSH_SUB(X(dsubmesh_to mesh_dotp))

  dotp = R_TOTYPE(M_ZERO)

  if(this%mesh%use_curvilinear) then
    do idim = 1, dim
      do is = 1, this%np
        dotp = dotp + this%mesh%vol_pp(this%map(is))*phi(this%map(is), idim)*sphi(is)
      end do
    end do
  else
    do idim = 1, dim
      do is = 1, this%np
        dotp = dotp + phi(this%map(is), idim)*sphi(is)
      end do
    end do
    dotp = dotp*this%mesh%vol_pp(1)
  end if

  if(optional_default(reduce, .true.) .and. this%mesh%parallel_in_domains) call comm_allreduce(this%mesh%vp%comm, dotp)

  POP_SUB(X(dsubmesh_to mesh_dotp))
end function X(dsubmesh_to_mesh_dotp)

!------------------------------------------------------------

!> The following functions takes a batch of functions defined in
!! submesh (ss) and adds all of them to each of the mesh functions in
!! other batch (mm).  Each one is multiplied by a factor given by the
!! array 'factor'.
subroutine X(submesh_batch_add_matrix)(this, factor, ss, mm)
  type(submesh_t),  intent(in)    :: this
  R_TYPE,           intent(in)    :: factor(:, :)
  type(batch_t),    intent(in)    :: ss
  type(batch_t),    intent(inout) :: mm

  integer :: ist, jst, idim, jdim, is
  type(profile_t), save :: prof
  
  PUSH_SUB(X(submesh_batch_add_matrix))
  call profiling_in(prof, 'SUBMESH_ADD_MATRIX')

  !$omp parallel do private(ist, idim, jdim, jst, is, aa)
  do ist =  1, mm%nst
    do idim = 1, mm%dim
      jdim = min(idim, ss%dim)
      do jst = 1, ss%nst
        if(associated(ss%states(jst)%dpsi)) then
          forall(is = 1:this%np)
            mm%states(ist)%X(psi)(this%map(is), idim) = &
              mm%states(ist)%X(psi)(this%map(is), idim) + factor(jst, ist)*ss%states(jst)%dpsi(is, jdim)
          end forall
        else
          forall(is = 1:this%np)
            mm%states(ist)%X(psi)(this%map(is), idim) = &
              mm%states(ist)%X(psi)(this%map(is), idim) + factor(jst, ist)*ss%states(jst)%zpsi(is, jdim)
          end forall
        end if
      end do
    end do
  end do
  !$omp end parallel do
  
  call profiling_out(prof)
  POP_SUB(X(submesh_batch_add_matrix))
end subroutine X(submesh_batch_add_matrix)


!------------------------------------------------------------ 

!> The following function takes a batch of functions defined in
!! submesh (ss) and adds one of them to each of the mesh functions in
!! other batch (mm).  Each one is multiplied by a factor given by the
!! array factor.
subroutine X(submesh_batch_add)(this, ss, mm)
  type(submesh_t),  intent(in)    :: this
  type(batch_t),    intent(in)    :: ss
  type(batch_t),    intent(inout) :: mm

  integer :: ist, idim, jdim, is

  PUSH_SUB(X(submesh_batch_add))

  ASSERT(mm%nst == ss%nst)

  !$omp parallel do private(ist, idim, jdim, is)
  do ist =  1, mm%nst
    do idim = 1, mm%dim
      jdim = min(idim, ss%dim)

      if(associated(ss%states(ist)%dpsi)) then

        forall(is = 1:this%np)
          mm%states(ist)%X(psi)(this%map(is), idim) = &
            mm%states(ist)%X(psi)(this%map(is), idim) + ss%states(ist)%dpsi(is, jdim)
        end forall
        
      else
        
        forall(is = 1:this%np)
          mm%states(ist)%X(psi)(this%map(is), idim) = &
            mm%states(ist)%X(psi)(this%map(is), idim) + ss%states(ist)%zpsi(is, jdim)
        end forall
        
      end if
    end do
  end do
  !$omp end parallel do
  
  POP_SUB(X(submesh_batch_add))
end subroutine X(submesh_batch_add)

!----------------------------------------------------------------------------------

subroutine X(submesh_batch_dotp_matrix)(this, mm, ss, dot, reduce)
  type(submesh_t),   intent(in)    :: this
  type(batch_t),     intent(in)    :: ss
  type(batch_t),     intent(in)    :: mm
  R_TYPE,            intent(inout) :: dot(:, :)
  logical, optional, intent(in)    :: reduce

  integer :: ist, jst, idim, jdim, is
  R_TYPE :: dotp

  PUSH_SUB(X(submesh_batch_dotp_matrix))

  if(this%mesh%use_curvilinear) then

    do ist = 1, ss%nst
      do jst = 1, mm%nst
        dotp = R_TOTYPE(M_ZERO)
        do idim = 1, ss%dim
          jdim = min(idim, ss%dim)

          if(associated(ss%states(ist)%dpsi)) then

            do is = 1, this%np
              dotp = dotp + this%mesh%vol_pp(this%map(is))*&
                mm%states(jst)%X(psi)(this%map(is), idim)*&
                ss%states(ist)%dpsi(is, jdim)
            end do

          else

            do is = 1, this%np
              dotp = dotp + this%mesh%vol_pp(this%map(is))*&
                mm%states(jst)%X(psi)(this%map(is), idim)*&
                ss%states(ist)%zpsi(is, jdim)
            end do

          end if
        end do

        dot(ist, jst) = dotp
      end do
    end do
    
  else

    !$omp parallel do private(ist, jst, dotp, idim, jdim, is)
    do ist = 1, ss%nst
      do jst = 1, mm%nst
        dotp = R_TOTYPE(M_ZERO)

        do idim = 1, mm%dim
          jdim = min(idim, ss%dim)

          if(associated(ss%states(ist)%dpsi)) then
            do is = 1, this%np
              dotp = dotp + &
                mm%states(jst)%X(psi)(this%map(is), idim)*&
                ss%states(ist)%dpsi(is, jdim)
            end do
          else
            do is = 1, this%np
              dotp = dotp + &
                mm%states(jst)%X(psi)(this%map(is), idim)*&
                ss%states(ist)%zpsi(is, jdim)
            end do
          end if

        end do

        dot(jst, ist) = dotp*this%mesh%vol_pp(1)
      end do
    end do
    
  end if

#if defined(HAVE_MPI)
  if(optional_default(reduce, .true.) .and. this%mesh%parallel_in_domains) then
    call comm_allreduce(this%mesh%mpi_grp%comm, dot, dim = (/mm%nst, ss%nst/))
  end if
#endif

  POP_SUB(X(submesh_batch_dotp_matrix))
end subroutine X(submesh_batch_dotp_matrix)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
