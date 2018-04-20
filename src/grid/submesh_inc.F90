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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: submesh_inc.F90 15314 2016-04-30 08:40:18Z xavier $

!Here ff is a function in the submesh 
R_TYPE function X(sm_integrate)(mesh, sm, ff) result(res)
  type(mesh_t),      intent(in) :: mesh
  type(submesh_t),   intent(in) :: sm
  R_TYPE, optional,  intent(in) :: ff(:)

  integer :: is

  PUSH_SUB(X(sm_integrate))

  ASSERT(present(ff) .or. sm%np  ==  0)

  if(sm%np > 0) then
    if (mesh%use_curvilinear) then
      res = sum(ff(1:sm%np)*mesh%vol_pp(sm%map(1:sm%np)) )
    else
      res = R_TOTYPE(M_ZERO)
      !$omp parallel do reduction(+:res)
      do is = 1, sm%np
        res = res+ff(is)
      end do
      !$omp end parallel do
      res=res*mesh%volume_element
    end if
  else
    res = M_ZERO
  end if

  if(mesh%parallel_in_domains) then
    call profiling_in(C_PROFILING_SM_REDUCE, "SM_REDUCE")
    call comm_allreduce(mesh%vp%comm, res)
    call profiling_out(C_PROFILING_SM_REDUCE)
  end if 
 
  POP_SUB(X(sm_integrate))
end function X(sm_integrate)

!Here ff is a function expressed in mesh
R_TYPE function X(sm_integrate_frommesh)(mesh, sm, ff) result(res)
  type(mesh_t),      intent(in) :: mesh
  type(submesh_t),   intent(in) :: sm
  R_TYPE, optional,  intent(in) :: ff(:)

  PUSH_SUB(X(sm_integrate_frommesh))

  ASSERT(present(ff) .or. sm%np  ==  0)

  if(sm%np > 0) then
    if (mesh%use_curvilinear) then
      res = sum(ff(sm%map(1:sm%np))*mesh%vol_pp(sm%map(1:sm%np)) )
    else
      res = sum(ff(sm%map(1:sm%np)))*mesh%volume_element
    end if
  else
    res = M_ZERO
  end if

  if(mesh%parallel_in_domains) call comm_allreduce(mesh%vp%comm, res)

  POP_SUB(X(sm_integrate_frommesh))
end function X(sm_integrate_frommesh)


!------------------------------------------------------------

subroutine X(dsubmesh_add_to_mesh)(this, sphi, phi, factor)
  type(submesh_t),  intent(in)    :: this
  FLOAT,            intent(in)    :: sphi(:)
  R_TYPE,           intent(inout) :: phi(:)
  R_TYPE, optional, intent(in)    :: factor

  integer :: ip, m

  PUSH_SUB(X(dsubmesh_add_to_mesh))

  if(present(factor)) then
    !Loop unrolling inspired by BLAS axpy routine
    m = mod(this%np,4)
    do ip = 1, m
      phi(this%map(ip)) = phi(this%map(ip)) + factor*sphi(ip)
    end do
    if( this%np.lt.4) return
    do ip = m+1, this%np, 4
      phi(this%map(ip))   = phi(this%map(ip))   + factor*sphi(ip)
      phi(this%map(ip+1)) = phi(this%map(ip+1)) + factor*sphi(ip+1)
      phi(this%map(ip+2)) = phi(this%map(ip+2)) + factor*sphi(ip+2)
      phi(this%map(ip+3)) = phi(this%map(ip+3)) + factor*sphi(ip+3)
    end do
  else
    m = mod(this%np,4)
    do ip = 1, m
      phi(this%map(ip)) = phi(this%map(ip)) + sphi(ip)
    end do
    if( this%np.lt.4) return
    do ip = m+1, this%np, 4
      phi(this%map(ip))   = phi(this%map(ip))   + sphi(ip)
      phi(this%map(ip+1)) = phi(this%map(ip+1)) + sphi(ip+1)
      phi(this%map(ip+2)) = phi(this%map(ip+2)) + sphi(ip+2)
      phi(this%map(ip+3)) = phi(this%map(ip+3)) + sphi(ip+3)
    end do
  end if

  POP_SUB(X(dsubmesh_add_to_mesh))
end subroutine X(dsubmesh_add_to_mesh)


! ---------------------------------------------------------
subroutine X(submesh_copy_from_mesh)(this, phi, sphi, conjugate)
  type(submesh_t),  intent(in)    :: this
  R_TYPE,           intent(in)    :: phi(:)
  R_TYPE,           intent(inout) :: sphi(:)
  logical, optional,   intent(in) :: conjugate


  integer :: ip

  PUSH_SUB(X(submesh_copy_from_mesh))

  if(.not. optional_default(conjugate, .false.) ) then
    !$omp parallel do 
    do ip = 1,this%np
      sphi(ip) = phi(this%map(ip))
    end do
  else
    !$omp parallel do 
    do ip = 1,this%np
      sphi(ip) = R_CONJ(phi(this%map(ip)))
    end do
  end if

  POP_SUB(X(submesh_copy_from_mesh))
end subroutine X(submesh_copy_from_mesh)

! ---------------------------------------------------------
subroutine X(submesh_copy_from_mesh_batch)(this, psib, spsi)
  type(submesh_t),  intent(in)    :: this
  type(batch_t),    intent(in)    :: psib
  R_TYPE,           intent(inout) :: spsi(:,:)

  integer :: ip, ist, ii, m, ip_map
  type(profile_t), save :: prof

  call profiling_in(prof, "SM_CP_MESH_BATCH")
  PUSH_SUB(X(submesh_copy_from_mesh_batch))

  ASSERT(batch_status(psib)/= BATCH_CL_PACKED)

  select case(batch_status(psib))
    case(BATCH_NOT_PACKED)
      do ist = 1, psib%nst_linear
        !$omp parallel do
        do ip = 1,this%np
          spsi(ip,ist) = psib%states_linear(ist)%X(psi)(this%map(ip))
        end do
      end do
    case(BATCH_PACKED)
      m = mod(psib%nst_linear,4)
      !$omp parallel do private(ii,ip_map)
      do ip = 1, this%np
        ip_map = this%map(ip)
        do ii = 1, m
          spsi(ii,ip) = psib%pack%X(psi)(ii,ip_map)
        end do
        do ii = m+1, psib%nst_linear, 4
          spsi(ii,ip) = psib%pack%X(psi)(ii,ip_map)
          spsi(ii+1,ip) = psib%pack%X(psi)(ii+1,ip_map)
          spsi(ii+2,ip) = psib%pack%X(psi)(ii+2,ip_map)
          spsi(ii+3,ip) = psib%pack%X(psi)(ii+3,ip_map)
        end do
      end do
      !$omp end parallel do
  end select

  POP_SUB(X(submesh_copy_from_mesh_batch))
   call profiling_out(prof)
end subroutine X(submesh_copy_from_mesh_batch)

 
! ---------------------------------------------------------
!> this function returns the the norm of a vector
FLOAT function X(sm_nrm2)(sm, ff, reduce) result(nrm2)
  type(submesh_t),   intent(in) :: sm
  R_TYPE,            intent(in) :: ff(:)
  logical, optional, intent(in) :: reduce

  R_TYPE, allocatable :: ll(:)

  call profiling_in(C_PROFILING_SM_NRM2, "SM_NRM2")
  PUSH_SUB(X(sm_nrm2))

  if(sm%mesh%use_curvilinear) then
    SAFE_ALLOCATE(ll(1:sm%np))
    ll(1:sm%np) = ff(1:sm%np)*sqrt(sm%mesh%vol_pp(sm%map(1:sm%np)))
    nrm2 = lalg_nrm2(sm%np, ll)
    SAFE_DEALLOCATE_A(ll)
  else
    nrm2 = lalg_nrm2(sm%np, ff)
  end if

  nrm2 = nrm2*sqrt(sm%mesh%volume_element)

  if(sm%mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(C_PROFILING_SM_REDUCE, "SM_REDUCE")
    nrm2 = nrm2**2
    call comm_allreduce(sm%mesh%vp%comm, nrm2)
    nrm2 = sqrt(nrm2)
    call profiling_out(C_PROFILING_SM_REDUCE)
  end if

  POP_SUB(X(sm_nrm2))
  call profiling_out(C_PROFILING_SM_NRM2)

end function X(sm_nrm2)


!------------------------------------------------------------

R_TYPE function X(dsubmesh_to_mesh_dotp)(this, sphi, phi, reduce) result(dotp)
  type(submesh_t),   intent(in) :: this
  FLOAT,             intent(in) :: sphi(:)
  R_TYPE,            intent(in) :: phi(:)
  logical, optional, intent(in) :: reduce

  integer :: is, m, ip

  PUSH_SUB(X(dsubmesh_to_mesh_dotp))

  dotp = R_TOTYPE(M_ZERO)

  if(this%mesh%use_curvilinear) then
    do is = 1, this%np
      dotp = dotp + this%mesh%vol_pp(this%map(is))*phi(this%map(is))*sphi(is)
    end do
  else
    m = mod(this%np,4)
    do ip = 1, m
      dotp = dotp + phi(this%map(ip))*sphi(ip)
    end do
    if( this%np.lt.4) return
    do ip = m+1, this%np, 4
      dotp = dotp + phi(this%map(ip))*sphi(ip) &
                  + phi(this%map(ip+1))*sphi(ip+1) &
                  + phi(this%map(ip+2))*sphi(ip+2) &
                  + phi(this%map(ip+3))*sphi(ip+3)
    end do
    dotp = dotp*this%mesh%vol_pp(1)
  end if

  if(optional_default(reduce, .true.) .and. this%mesh%parallel_in_domains) &
    call comm_allreduce(this%mesh%vp%comm, dotp)

  POP_SUB(X(dsubmesh_to_mesh_dotp))
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

  !$omp parallel do private(ist, idim, jdim, jst, is)
  do ist =  1, min(mm%nst, ubound(factor, 2))
    do idim = 1, mm%dim
      ! FIXME: this line should instead be assert(mm%dim == ss%dim)!!
      jdim = min(idim, ss%dim)
      do jst = 1, ss%nst
        if(associated(ss%states(jst)%dpsi)) then
          forall(is = 1:this%np)
            mm%states(ist)%X(psi)(this%map(is), idim) = &
              mm%states(ist)%X(psi)(this%map(is), idim) + factor(jst, ist)*ss%states(jst)%dpsi(is, jdim)
          end forall
        else

#ifdef R_TCOMPLEX
          forall(is = 1:this%np)
            mm%states(ist)%X(psi)(this%map(is), idim) = &
              mm%states(ist)%X(psi)(this%map(is), idim) + factor(jst, ist)*ss%states(jst)%zpsi(is, jdim)
          end forall
#else
          message(1) = "Internal error: cannot call dsubmesh_batch_add_matrix with complex batch ss"
          call messages_fatal(1)
#endif

        end if
      end do
    end do
  end do
  !$omp end parallel do

  call profiling_count_operations(mm%nst*mm%dim*ss%nst*this%np*(R_ADD + R_MUL))
  
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

#ifdef R_TCOMPLEX        
        forall(is = 1:this%np)
          mm%states(ist)%X(psi)(this%map(is), idim) = &
            mm%states(ist)%X(psi)(this%map(is), idim) + ss%states(ist)%zpsi(is, jdim)
        end forall
#else
        message(1) = "Internal error: cannot call dsubmesh_batch_add with complex batch ss"
        call messages_fatal(1)
#endif
        
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
                R_CONJ(mm%states(jst)%X(psi)(this%map(is), idim))*&
                ss%states(ist)%dpsi(is, jdim)
            end do

          else

#ifdef R_TCOMPLEX
            do is = 1, this%np
              dotp = dotp + this%mesh%vol_pp(this%map(is))*&
                R_CONJ(mm%states(jst)%X(psi)(this%map(is), idim))*&
                ss%states(ist)%zpsi(is, jdim)
            end do
#else
            message(1) = "Internal error: cannot call dsubmesh_batch_dotp_matrix with complex batch ss"
            call messages_fatal(1)
#endif

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
                R_CONJ(mm%states(jst)%X(psi)(this%map(is), idim))*&
                ss%states(ist)%dpsi(is, jdim)
            end do
          else

#ifdef R_TCOMPLEX
            do is = 1, this%np
              dotp = dotp + &
                R_CONJ(mm%states(jst)%X(psi)(this%map(is), idim))*&
                ss%states(ist)%zpsi(is, jdim)
            end do
#else
            message(1) = "Internal error: cannot call dsubmesh_batch_dotp_matrix with complex batch ss"
            call messages_fatal(1)
#endif

          end if

        end do

        dot(jst, ist) = dotp*this%mesh%volume_element
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
