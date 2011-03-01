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

!------------------------------------------------------------------------------
! X(project_psi) calculates the action of a projector on the psi wavefunction.
! The result is summed up to ppsi
subroutine X(project_psi)(mesh, pj, npj, dim, psi, ppsi, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj(:)
  integer,           intent(in)    :: npj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(in)    :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: ik

  type(batch_t) :: psib, ppsib

  PUSH_SUB(X(project_psi))

  call batch_init(psib, dim, 1)
  call batch_add_state(psib, 1, psi)
  call batch_init(ppsib, dim, 1)
  call batch_add_state(ppsib, 1, ppsi)

  call X(project_psi_batch)(mesh, pj, npj, dim, psib, ppsib, ik)

  call batch_end(psib)
  call batch_end(ppsib)

  POP_SUB(X(project_psi))
end subroutine X(project_psi)


!------------------------------------------------------------------------------
subroutine X(project_psi_batch)(mesh, pj, npj, dim, psib, ppsib, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj(:)
  integer,           intent(in)    :: npj
  integer,           intent(in)    :: dim
  type(batch_t),     intent(in)    :: psib
  type(batch_t),     intent(inout) :: ppsib
  integer,           intent(in)    :: ik

  integer :: ipj, nreduce, ii, ns, idim, ll, mm, is, ist
  R_TYPE, allocatable :: reduce_buffer(:), lpsi(:, :)
  integer, allocatable :: ireduce(:, :, :, :)
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof

  PUSH_SUB(X(project_psi_batch))
  call profiling_in(prof, "VNLPSI")

  ! To optimize the application of the non-local operator in parallel,
  ! the projectors are applied in steps. First the <p|psi> is
  ! calculated for all projectors and the result is stored on an array
  ! (reduce_buffer). Then the array is reduced (as it is contiguous
  ! only one reduction is required). Finally |ppsi> += |p><p|psi> is
  ! calculated.

  ! generate the reduce buffer and related structures
  SAFE_ALLOCATE(ireduce(1:npj, 0:MAX_L, -MAX_L:MAX_L, 1:psib%nst))
  nreduce = 0
  
  ! count the number of elements in the reduce buffer
  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == M_NONE) cycle
      do ll = 0, pj(ipj)%lmax
        if (ll == pj(ipj)%lloc) cycle
        do mm = -ll, ll
          ireduce(ipj, ll, mm, ist) = 1 + nreduce
          nreduce = nreduce + pj(ipj)%nprojections
        end do
      end do
    end do
  end do

  ! Check whether we have or not "real" projectors. If we do not, return.
  if(nreduce == 0) then
    call profiling_out(prof)
    POP_SUB(X(project_psi_batch))
    return
  end if

  SAFE_ALLOCATE(reduce_buffer(1:nreduce))

  reduce_buffer = R_TOTYPE(M_ZERO)
  
  !$omp parallel private(ist, ipj, ns, lpsi, ll, mm, ii, idim, is)
  SAFE_ALLOCATE(lpsi(1:maxval(pj(1:npj)%sphere%ns), 1:dim))
  !$omp do
  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == M_NONE) cycle
      ns = pj(ipj)%sphere%ns
      if(ns < 1) cycle

      ! copy psi to the small spherical grid
      do idim = 1, dim
        if(associated(pj(ipj)%phase)) then
          forall (is = 1:ns) 
            lpsi(is, idim) = psib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim)*pj(ipj)%phase(is, ik)
          end forall
        else
          forall (is = 1:ns) lpsi(is, idim) = psib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim)
        end if
      end do

      ! apply the projectors for each angular momentum component
      do ll = 0, pj(ipj)%lmax
        if (ll == pj(ipj)%lloc) cycle
        do mm = -ll, ll

          ii = ireduce(ipj, ll, mm, ist)
          select case(pj(ipj)%type)
          case(M_KB)
            call X(kb_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(ll, mm), dim, lpsi(1:ns, 1:dim), reduce_buffer(ii:))
          case(M_RKB)
#ifdef R_TCOMPLEX
            if(ll /= 0) then
              call rkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%rkb_p(ll, mm), lpsi(1:ns, 1:dim), reduce_buffer(ii:))
            else
              call zkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(1, 1), dim, lpsi(1:ns, 1:dim), reduce_buffer(ii:))
            end if
#endif
          case(M_HGH)
            call X(hgh_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%hgh_p(ll, mm), dim, pj(ipj)%reltype, &
              lpsi(1:ns, 1:dim), reduce_buffer(ii:))
          end select
        end do ! mm
      end do ! ll

    end do ! ipj
  end do ! ist
  !$omp end do nowait
  SAFE_DEALLOCATE_A(lpsi)
  !$omp end parallel

  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "VNLPSI_REDUCE_BATCH")
    call comm_allreduce(mesh%mpi_grp%comm, reduce_buffer, dim = nreduce)
    call profiling_out(reduce_prof)
  end if

  ! calculate |ppsi> += |p><p|psi>
  !$omp parallel private(ist, ipj, ns, lpsi, ll, mm, ii, idim, is)
  SAFE_ALLOCATE(lpsi(1:maxval(pj(1:npj)%sphere%ns), 1:dim))
  !$omp do
  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == M_NONE) cycle

      ns = pj(ipj)%sphere%ns
      if(ns < 1) cycle

      lpsi(1:ns, 1:dim) = M_ZERO

      do ll = 0, pj(ipj)%lmax
        if (ll == pj(ipj)%lloc) cycle
        do mm = -ll, ll
          ii = ireduce(ipj, ll, mm, ist)

          select case(pj(ipj)%type)
          case(M_KB)
            call X(kb_project_ket)(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(ll, mm), dim, &
              reduce_buffer(ii:), lpsi(1:ns, 1:dim))
          case(M_RKB)
#ifdef R_TCOMPLEX
            if(ll /= 0) then
              call rkb_project_ket(pj(ipj)%rkb_p(ll, mm), reduce_buffer(ii:), lpsi(1:ns, 1:dim))
            else
              call zkb_project_ket(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(1, 1), dim, &
                reduce_buffer(ii:), lpsi(1:ns, 1:dim))
            end if
#endif
          case(M_HGH)
            call X(hgh_project_ket)(mesh, pj(ipj)%sphere, pj(ipj)%hgh_p(ll, mm), dim, &
              pj(ipj)%reltype, reduce_buffer(ii:), lpsi(1:ns, 1:dim))
          end select
          
        end do ! mm
      end do ! ll
    
      !put the result back in the complete grid
      do idim = 1, dim
        if(associated(pj(ipj)%phase)) then
          forall (is = 1:ns) 
            ppsib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim) = &
              ppsib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim) + &
              lpsi(is, idim)*conjg(pj(ipj)%phase(is, ik))
          end forall
        else
          forall (is = 1:ns) 
            ppsib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim) = &
              ppsib%states(ist)%X(psi)(pj(ipj)%sphere%jxyz(is), idim) + lpsi(is, idim)
          end forall
        end if
      end do

    end do ! ipj
  end do ! ist
  !$omp end do nowait
  SAFE_DEALLOCATE_A(lpsi)
  !$omp end parallel

  SAFE_DEALLOCATE_A(reduce_buffer)
  SAFE_DEALLOCATE_A(ireduce)

  call profiling_out(prof)
  POP_SUB(X(project_psi_batch))

end subroutine X(project_psi_batch)


!------------------------------------------------------------------------------
! X(projector_matrix_element) calculates <psia|projector|psib>
R_TYPE function X(projector_matrix_element)(pj, dim, ik, psia, psib) result(apb)
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  integer,           intent(in)    :: ik
  R_TYPE,            intent(in)    :: psia(:, :)  ! psia(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: psib(:, :)  ! psib(1:mesh%np, dim)

  integer ::  ns, idim, ll, mm, nc, is
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)
  R_TYPE, allocatable :: uvpsi(:, :, :, :)
#if defined(HAVE_MPI)
  integer :: size
  R_TYPE, allocatable :: uvpsi_tmp(:, :, :, :)
#endif
  type(mesh_t), pointer :: mesh

  PUSH_SUB(X(projector_matrix_element))

  ns = pj%sphere%ns

  ASSERT(associated(pj%sphere%mesh))
  mesh => pj%sphere%mesh

  SAFE_ALLOCATE(lpsi(1:ns, 1:dim))
  SAFE_ALLOCATE(plpsi(1:ns, 1:dim))

  do idim = 1, dim
    if(simul_box_is_periodic(mesh%sb)) then
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%jxyz(is), idim)*pj%phase(is, ik)
      end do
    else
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%jxyz(is), idim)
      end do
    end if
  end do

  if(pj%type == M_KB) then

    do idim = 1, dim
      if(simul_box_is_periodic(mesh%sb)) then
        do is = 1, ns
          plpsi(is, idim) = psia(pj%sphere%jxyz(is), idim)*pj%phase(is, ik)
        end do
      else
        do is = 1, ns
          plpsi(is, idim) = psia(pj%sphere%jxyz(is), idim)
        end do
      end if
    end do

    SAFE_ALLOCATE(uvpsi(1:pj%nprojections, 1:2, 0:pj%lmax, -pj%lmax:pj%lmax))
    uvpsi = R_TOTYPE(M_ZERO)

    ASSERT(associated(pj%kb_p))

    do ll = 0, pj%lmax
      if (ll == pj%lloc) cycle
      do mm = -ll, ll
        call X(kb_project_bra)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, lpsi, uvpsi(:, 1, ll, mm))
        call X(kb_project_bra)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, plpsi, uvpsi(:, 2, ll, mm))
      end do
    end do

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(uvpsi_tmp(1:pj%nprojections, 1:2, 0:pj%lmax, -pj%lmax:pj%lmax))
      uvpsi_tmp = R_TOTYPE(M_ZERO)
      size = pj%nprojections*2*(pj%lmax + 1)*(2*pj%lmax + 1)
      call MPI_Allreduce(uvpsi(1, 1, 0, -pj%lmax), uvpsi_tmp(1, 1, 0, -pj%lmax), size, &
        R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
      uvpsi = uvpsi_tmp
      SAFE_DEALLOCATE_A(uvpsi_tmp)
    end if
#endif

    apb = M_ZERO

    do ll = 0, pj%lmax
      if (ll == pj%lloc) cycle
      do mm = -ll, ll
        nc = pj%kb_p(ll, mm)%n_c
        call X(kb_mul_energies)(pj%kb_p(ll, mm), dim, uvpsi(1:nc, 1, ll, mm))
        apb = apb + sum(R_CONJ(uvpsi(1:nc, 2, ll, mm))*uvpsi(1:nc, 1, ll, mm))
      end do
    end do

    SAFE_DEALLOCATE_A(uvpsi)

   else

    call X(project_sphere)(mesh, pj, dim, lpsi, plpsi)

    apb = M_ZERO
    do idim = 1, dim
      if(simul_box_is_periodic(mesh%sb)) then
        do is = 1, ns
          plpsi(is, idim) = R_CONJ(psia(pj%sphere%jxyz(is), idim))*plpsi(is, idim)*conjg(pj%phase(is, ik))
        end do
      else
        do is = 1, ns
          plpsi(is, idim) = R_CONJ(psia(pj%sphere%jxyz(is), idim))*plpsi(is, idim)
        end do
      end if

      if(ns > 0) then
        apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:ns, idim))
      else
        apb = apb + X(sm_integrate)(mesh, pj%sphere)
      endif  
    end do

  end if

  SAFE_DEALLOCATE_A(lpsi)
  SAFE_DEALLOCATE_A(plpsi)

  POP_SUB(X(projector_matrix_element))
end function X(projector_matrix_element)

!------------------------------------------------------------------------------
subroutine X(project_sphere)(mesh, pj, dim, psi, ppsi)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:ns, dim)
  R_TYPE,            intent(out)   :: ppsi(:, :)  ! ppsi(1:ns, dim)

  integer :: ll, mm

  PUSH_SUB(X(project_sphere))

  ppsi = M_ZERO

  do ll = 0, pj%lmax
    if (ll == pj%lloc) cycle
    do mm = -ll, ll
      
      select case (pj%type)
      case (M_HGH)
        call X(hgh_project)(mesh, pj%sphere, pj%hgh_p(ll, mm), dim, psi, ppsi, pj%reltype)
      case (M_KB)
        call X(kb_project)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, psi, ppsi)
      case (M_RKB)
#ifdef R_TCOMPLEX
        if(ll /= 0) then
          call rkb_project(mesh, pj%sphere, pj%rkb_p(ll, mm), psi, ppsi)
        else
          call zkb_project(mesh, pj%sphere, pj%kb_p(1, 1), dim, psi, ppsi)
        end if
#endif
      end select
  
    end do
  end do

  POP_SUB(X(project_sphere))
end subroutine X(project_sphere)


!------------------------------------------------------------------------------
!This function calculates |cpsi> += [x, V_nl] |psi>
subroutine X(projector_commute_r)(pj, gr, dim, idir, ik, psi, cpsi)
  type(projector_t),     intent(in)     :: pj
  type(grid_t),          intent(in)     :: gr
  integer,               intent(in)     :: dim
  integer,               intent(in)     :: idir
  integer,               intent(in)     :: ik
  R_TYPE,                intent(in)     :: psi(:, :)
  R_TYPE,                intent(inout)  :: cpsi(:,:)

  integer ::  ns, idim
  R_TYPE, allocatable :: lpsi(:, :), pxlpsi(:,:), xplpsi(:,:)
  integer, pointer :: jxyz(:)
  FLOAT,   pointer :: smx(:, :)
  type(profile_t), save :: prof

  PUSH_SUB(X(projector_commute_r))
  call profiling_in(prof, "PROJ_COMMUTE")

  if(pj%type .ne. M_NONE) then

    ns = pj%sphere%ns
    jxyz => pj%sphere%jxyz
    smx => pj%sphere%x

    SAFE_ALLOCATE(  lpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(xplpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(pxlpsi(1:ns, 1:dim))

    if(simul_box_is_periodic(gr%mesh%sb)) then
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(jxyz(1:ns), idim)*pj%phase(1:ns, ik)
      end do
    else
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(jxyz(1:ns), idim)
      end do
    end if

    ! x V_nl |psi>
    call X(project_sphere)(gr%mesh, pj, dim, lpsi, xplpsi)
    do idim = 1, dim
      xplpsi(1:ns, idim) = smx(1:ns, idir) * xplpsi(1:ns, idim)
    end do

    ! V_nl x |psi>
    do idim = 1, dim
      lpsi(1:ns, idim) = smx(1:ns, idir) * lpsi(1:ns, idim)
    end do
    call X(project_sphere)(gr%mesh, pj, dim, lpsi, pxlpsi)
    
    ! |cpsi> += x V_nl |psi> - V_nl x |psi> 
    if(simul_box_is_periodic(gr%mesh%sb)) then
      do idim = 1, dim
        cpsi(jxyz(1:ns), idim) = cpsi(jxyz(1:ns), idim) + &
          (xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)) * R_CONJ(pj%phase(1:ns, ik))
      end do
    else
      do idim = 1, dim
        cpsi(jxyz(1:ns), idim) = cpsi(jxyz(1:ns), idim) + xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)
      end do
    end if

    SAFE_DEALLOCATE_A(lpsi)
    SAFE_DEALLOCATE_A(xplpsi)
    SAFE_DEALLOCATE_A(pxlpsi)
  end if
  call profiling_out(prof)
  POP_SUB(X(projector_commute_r))

end subroutine X(projector_commute_r)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
