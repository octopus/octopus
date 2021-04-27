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

!------------------------------------------------------------------------------
!> X(project_psi) calculates the action of a projector on the psi wavefunction.
!! The result is summed up to ppsi
subroutine X(project_psi)(mesh, bnd, pj, npj, dim, psi, ppsi, ik)
  type(mesh_t),       intent(in)    :: mesh
  type(boundaries_t), intent(in)    :: bnd
  type(projector_t),  intent(in)    :: pj(:)
  integer,            intent(in)    :: npj
  integer,            intent(in)    :: dim
  R_TYPE, contiguous, intent(inout)    :: psi(:, :)   !< (1:mesh%np, dim)
  R_TYPE, contiguous, intent(inout) :: ppsi(:, :)  !< (1:mesh%np, dim)
  integer,            intent(in)    :: ik

  type(wfs_elec_t) :: psib, ppsib

  PUSH_SUB(X(project_psi))

  call wfs_elec_init(psib, dim, 1, 1, psi, ik)
  call wfs_elec_init(ppsib, dim, 1, 1, ppsi, ik)

  call X(project_psi_batch)(mesh, bnd, pj, npj, dim, psib, ppsib)

  call psib%end()
  call ppsib%end()

  POP_SUB(X(project_psi))
end subroutine X(project_psi)


!------------------------------------------------------------------------------
!> To optimize the application of the non-local operator in parallel,
!! the projectors are applied in steps. 
!!
!! First the <p|psi> is
!! calculated for all projectors and the result is stored on an array
!! (reduce_buffer). Then the array is reduced (as it is contiguous
!! only one reduction is required). Finally |ppsi> += |p><p|psi> is
!! calculated.
subroutine X(project_psi_batch)(mesh, bnd, pj, npj, dim, psib, ppsib)
  type(mesh_t),      intent(in)    :: mesh
  type(boundaries_t),intent(in)    :: bnd
  type(projector_t), intent(in)    :: pj(:)
  integer,           intent(in)    :: npj
  integer,           intent(in)    :: dim
  type(wfs_elec_t),  intent(in)    :: psib
  type(wfs_elec_t),  intent(inout) :: ppsib

  integer :: ipj, nreduce, ii, ns, idim, ll, mm, is, ist, bind
  R_TYPE, allocatable :: reduce_buffer(:,:), lpsi(:, :), uvpsi(:,:,:)
  integer, allocatable :: ireduce(:, :, :, :)
  type(profile_t), save :: prof
  type(profile_t), save :: reduce_prof

  PUSH_SUB(X(project_psi_batch))
  call profiling_in(prof, TOSTRING(X(VNLPSI)))

  ASSERT(.not.bnd%spiral)
  ASSERT(psib%status() /= BATCH_DEVICE_PACKED)

  ! generate the reduce buffer and related structures
  SAFE_ALLOCATE(ireduce(1:npj, 0:MAX_L, -MAX_L:MAX_L, 1:psib%nst))
  nreduce = 0
  
  ! count the number of elements in the reduce buffer
  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == PROJ_NONE) cycle
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

  SAFE_ALLOCATE(reduce_buffer(1:dim, 1:nreduce))

  reduce_buffer = R_TOTYPE(M_ZERO)
  
  ! FIXME: restore openmp
  SAFE_ALLOCATE(lpsi(1:maxval(pj(1:npj)%sphere%np), 1:dim))

  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == PROJ_NONE) cycle
      ns = pj(ipj)%sphere%np
      if(ns < 1) cycle

      ! copy psi to the small spherical grid
      select case(psib%status())
      case(BATCH_NOT_PACKED)
        do idim = 1, dim
          bind = psib%ist_idim_to_linear((/ist, idim/))
          if (allocated(pj(ipj)%phase)) then
#ifdef R_TCOMPLEX
            do is = 1, ns
              lpsi(is, idim) = psib%zff_linear(pj(ipj)%sphere%map(is), bind)*pj(ipj)%phase(is, 1, psib%ik)
            end do
#else
            ! Phase not allowed for real batches
            ASSERT(.false.)
#endif
          else
            do is = 1, ns
              lpsi(is, idim) = psib%X(ff_linear)(pj(ipj)%sphere%map(is), bind)
            end do
          end if
        end do

      case(BATCH_PACKED)
        do idim = 1, dim
          bind = psib%ist_idim_to_linear((/ist, idim/))
          if (allocated(pj(ipj)%phase)) then
#ifdef R_TCOMPLEX
            do is = 1, ns
              lpsi(is, idim) = psib%zff_pack(bind, pj(ipj)%sphere%map(is))*pj(ipj)%phase(is, 1, psib%ik)
            end do
#else
            ! Phase not allowed for real batches
            ASSERT(.false.)
#endif  
          else
            do is = 1, ns
              lpsi(is, idim) = psib%X(ff_pack)(bind, pj(ipj)%sphere%map(is))
            end do
          end if
        end do

      end select

      ! apply the projectors for each angular momentum component
      do ll = 0, pj(ipj)%lmax
        if (ll == pj(ipj)%lloc) cycle
        do mm = -ll, ll

          ii = ireduce(ipj, ll, mm, ist)
          select case(pj(ipj)%type)
          case(PROJ_KB)
            call X(kb_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(ll, mm), dim, lpsi(1:ns, 1:dim), reduce_buffer(1:dim, ii:))
          case(PROJ_RKB)
#ifdef R_TCOMPLEX
            if(ll /= 0) then
              call rkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%rkb_p(ll, mm), lpsi(1:ns, 1:dim), reduce_buffer(1:dim, ii:))
            else
              call zkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(1, 1), dim, lpsi(1:ns, 1:dim), reduce_buffer(1:dim, ii:))
            end if
#endif
          case(PROJ_HGH)
            call X(hgh_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%hgh_p(ll, mm), dim, pj(ipj)%reltype, &
              lpsi(1:ns, 1:dim), reduce_buffer(1:dim, ii:))
          end select
        end do ! mm
      end do ! ll

    end do ! ipj
  end do ! ist

  SAFE_DEALLOCATE_A(lpsi)

  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, TOSTRING(X(VNLPSI_REDUCE_BATCH)))
    call mesh%allreduce(reduce_buffer)
    call profiling_out(reduce_prof)
  end if

  ! calculate |ppsi> += |p><p|psi>
  !$omp parallel private(ist, ipj, ns, lpsi, ll, mm, ii, idim, is, uvpsi, bind)
  SAFE_ALLOCATE(lpsi(1:maxval(pj(1:npj)%sphere%np), 1:dim))
  !$omp do
  do ist = 1, psib%nst
    do ipj = 1, npj
      if(pj(ipj)%type == PROJ_NONE) cycle

      ns = pj(ipj)%sphere%np
      if(ns < 1) cycle

      lpsi(1:ns, 1:dim) = M_ZERO

      do ll = 0, pj(ipj)%lmax
        if (ll == pj(ipj)%lloc) cycle
        do mm = -ll, ll
          ii = ireduce(ipj, ll, mm, ist)

          select case(pj(ipj)%type)
          case(PROJ_KB)
            call X(kb_project_ket)(pj(ipj)%kb_p(ll, mm), dim, reduce_buffer(1:dim, ii:), lpsi(1:ns, 1:dim))
          case(PROJ_RKB)
#ifdef R_TCOMPLEX
            if(ll /= 0) then
              call rkb_project_ket(pj(ipj)%rkb_p(ll, mm), reduce_buffer(1:dim, ii:), lpsi(1:ns, 1:dim))
            else
              call zkb_project_ket(pj(ipj)%kb_p(1, 1), dim, reduce_buffer(1:dim, ii:), lpsi(1:ns, 1:dim))
            end if
#endif
          end select
        end do ! mm
 
        if(pj(ipj)%type == PROJ_HGH) then
          SAFE_ALLOCATE(uvpsi(1:dim, 1:3, -ll:ll))
          do mm = -ll,ll
            ii = ireduce(ipj, ll, mm, ist)
            uvpsi(1:dim, 1:3, mm) = reduce_buffer(1:dim, ii:ii+2)
          end do
          call X(hgh_project_ket)(pj(ipj)%hgh_p(ll, :), ll, pj(ipj)%lmax, dim, &
              pj(ipj)%reltype, uvpsi(1:dim, 1:3, -ll:ll), lpsi(1:ns, 1:dim))
          SAFE_DEALLOCATE_A(uvpsi)
        end if
      end do ! ll

    !  print *, ll, lpsi(1, 1:dim)  
  
      !put the result back in the complete grid
      select case(psib%status())
      case(BATCH_NOT_PACKED)
        do idim = 1, dim
          bind = psib%ist_idim_to_linear((/ist, idim/))
          if (allocated(pj(ipj)%phase)) then
#ifdef R_TCOMPLEX
            do is = 1, ns
              ppsib%zff_linear(pj(ipj)%sphere%map(is), bind) = &
                ppsib%zff_linear(pj(ipj)%sphere%map(is), bind) + lpsi(is, idim)*conjg(pj(ipj)%phase(is, 1, psib%ik))
            end do
#else
            ! Phase not allowed for real batches
            ASSERT(.false.)
#endif
          else
            do is = 1, ns
              ppsib%X(ff_linear)(pj(ipj)%sphere%map(is), bind) = &
                ppsib%X(ff_linear)(pj(ipj)%sphere%map(is), bind) + lpsi(is, idim)
            end do
          end if
        end do

      case(BATCH_PACKED)
        do idim = 1, dim
          bind = psib%ist_idim_to_linear((/ist, idim/))
          if (allocated(pj(ipj)%phase)) then
#ifdef R_TCOMPLEX
            do is = 1, ns
              ppsib%zff_pack(bind, pj(ipj)%sphere%map(is)) = &
                ppsib%zff_pack(bind, pj(ipj)%sphere%map(is)) + lpsi(is, idim)*conjg(pj(ipj)%phase(is, 1, psib%ik))
            end do
#else
            ! Phase not allowed for real batches
            ASSERT(.false.)
#endif
          else
            do is = 1, ns
              ppsib%X(ff_pack)(bind, pj(ipj)%sphere%map(is)) = &
                ppsib%X(ff_pack)(bind, pj(ipj)%sphere%map(is)) + lpsi(is, idim)
            end do
          end if
        end do

      end select

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
!> X(projector_matrix_element) calculates <psia|projector|psib>
R_TYPE function X(projector_matrix_element)(pj, bnd, dim, ik, psia, psib) result(apb)
  type(projector_t), target, intent(in)    :: pj
  type(boundaries_t),        intent(in)    :: bnd
  integer,                   intent(in)    :: dim
  integer,                   intent(in)    :: ik
  R_TYPE,                    intent(in)    :: psia(:, :)  !< psia(1:mesh%np, dim)
  R_TYPE,                    intent(in)    :: psib(:, :)  !< psib(1:mesh%np, dim)

  integer ::  ns, idim, is
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)
  type(mesh_t), pointer :: mesh
  type(profile_t), save :: prof

  PUSH_SUB(X(projector_matrix_element))

  call profiling_in(prof, TOSTRING(X(PROJ_MAT_ELEM)))
  ASSERT(.not. bnd%spiral)

  ns = pj%sphere%np

  ASSERT(associated(pj%sphere%mesh))
  mesh => pj%sphere%mesh

  SAFE_ALLOCATE(lpsi(1:ns, 1:dim))
  SAFE_ALLOCATE(plpsi(1:ns, 1:dim))

  do idim = 1, dim
    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%map(is), idim)*pj%phase(is, 1, ik)
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
    else
      do is = 1, ns
        lpsi(is, idim) = psib(pj%sphere%map(is), idim)
      end do
    end if
  end do

  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi)

  apb = M_ZERO
  do idim = 1, dim
    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do is = 1, ns
        plpsi(is, idim) = R_CONJ(psia(pj%sphere%map(is), idim))*plpsi(is, idim)*conjg(pj%phase(is, 1, ik))
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
    else
      do is = 1, ns
        plpsi(is, idim) = R_CONJ(psia(pj%sphere%map(is), idim))*plpsi(is, idim)
      end do
    end if

    if(ns > 0) then
      apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:ns, idim))
    else
      apb = apb + X(sm_integrate)(mesh, pj%sphere)
    end if
  end do

  SAFE_DEALLOCATE_A(lpsi)
  SAFE_DEALLOCATE_A(plpsi)

  call profiling_out(prof)
  
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
      case (PROJ_KB)
        call X(kb_project)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, psi, ppsi)
      case (PROJ_RKB)
#ifdef R_TCOMPLEX
        if(ll /= 0) then
          call rkb_project(mesh, pj%sphere, pj%rkb_p(ll, mm), psi, ppsi)
        else
          call zkb_project(mesh, pj%sphere, pj%kb_p(1, 1), dim, psi, ppsi)
        end if
#endif
      end select
  
    end do
    if(pj%type == PROJ_HGH) then
      call X(hgh_project)(mesh, pj%sphere, pj%hgh_p(ll, :), ll, pj%lmax, dim, psi, ppsi, pj%reltype)
    end if
  end do

  POP_SUB(X(project_sphere))
end subroutine X(project_sphere)


!------------------------------------------------------------------------------
!> This function calculates |cpsi> += [x, V_nl] |psi>
subroutine X(projector_commute_r)(pj, mesh, bnd, dim, idir, ik, psi, cpsi)
  type(projector_t), target, intent(in)     :: pj
  type(mesh_t),              intent(in)     :: mesh
  type(boundaries_t),        intent(in)     :: bnd
  integer,                   intent(in)     :: dim
  integer,                   intent(in)     :: idir
  integer,                   intent(in)     :: ik
  R_TYPE,                    intent(in)     :: psi(:, :)
  R_TYPE,                    intent(inout)  :: cpsi(:,:)

  integer ::  ns, idim
  R_TYPE, allocatable :: lpsi(:, :), pxlpsi(:,:), xplpsi(:,:)
  integer, pointer :: map(:)
  FLOAT,   pointer :: smx(:, :)
  type(profile_t), save :: prof

  PUSH_SUB(X(projector_commute_r))
  call profiling_in(prof, TOSTRING(X(PROJ_COMMUTE)))

  ASSERT(.not. bnd%spiral)

  if(pj%type /= PROJ_NONE) then

    ns = pj%sphere%np
    map => pj%sphere%map
    smx => pj%sphere%x

    SAFE_ALLOCATE(  lpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(xplpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(pxlpsi(1:ns, 1:dim))

    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(map(1:ns), idim)*pj%phase(1:ns, 1, ik)
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
    else
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(map(1:ns), idim)
      end do
    end if

    ! x V_nl |psi>
    call X(project_sphere)(mesh, pj, dim, lpsi, xplpsi)
    do idim = 1, dim
      xplpsi(1:ns, idim) = smx(1:ns, idir) * xplpsi(1:ns, idim)
    end do

    ! V_nl x |psi>
    do idim = 1, dim
      lpsi(1:ns, idim) = smx(1:ns, idir) * lpsi(1:ns, idim)
    end do
    call X(project_sphere)(mesh, pj, dim, lpsi, pxlpsi)
    
    ! |cpsi> += x V_nl |psi> - V_nl x |psi> 
    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do idim = 1, dim
        cpsi(map(1:ns), idim) = cpsi(map(1:ns), idim) + &
          (xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)) * R_CONJ(pj%phase(1:ns, 1, ik))
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
    else
      do idim = 1, dim
        cpsi(map(1:ns), idim) = cpsi(map(1:ns), idim) + xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)
      end do
    end if

    SAFE_DEALLOCATE_A(lpsi)
    SAFE_DEALLOCATE_A(xplpsi)
    SAFE_DEALLOCATE_A(pxlpsi)
  end if
  call profiling_out(prof)
  POP_SUB(X(projector_commute_r))

end subroutine X(projector_commute_r)

!------------------------------------------------------------------------------
!> This function calculates |cpsi> += [x, V_nl] |psi>
subroutine X(projector_commute_r_allatoms_alldir)(pj, ions, mesh, dim, bnd, ik, psi, cpsi)
  type(projector_t), target, intent(in)     :: pj(:)
  type(ions_t),              intent(in)     :: ions
  type(mesh_t),              intent(in)     :: mesh
  integer,                   intent(in)     :: dim
  type(boundaries_t),        intent(in)     :: bnd
  integer,                   intent(in)     :: ik
  R_TYPE,                    intent(in)     :: psi(:, :)
  R_TYPE,                    intent(inout)  :: cpsi(:,:,:)

  integer ::  ns, idim, idir, iatom
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:), xlpsi(:,:), pxlpsi(:,:), xplpsi(:,:)
  integer, pointer :: map(:)
  FLOAT,   pointer :: smx(:, :)
  logical  :: phase
  type(profile_t), save :: prof

  PUSH_SUB(X(projector_commute_r_allatoms_alldir))
  call profiling_in(prof, TOSTRING(X(PROJ_COMMUTE_ALL)))

  ASSERT(.not. bnd%spiral)

  do iatom = 1, ions%natoms
    if(species_is_ps(ions%atom(iatom)%species) .and. pj(iatom)%type /= PROJ_NONE) then

      ns = pj(iatom)%sphere%np
      map => pj(iatom)%sphere%map
      smx => pj(iatom)%sphere%x

      SAFE_ALLOCATE(  lpsi(1:ns, 1:dim))
      SAFE_ALLOCATE( plpsi(1:ns, 1:dim))
      SAFE_ALLOCATE( xlpsi(1:ns, 1:dim))
      SAFE_ALLOCATE(xplpsi(1:ns, 1:dim))
      SAFE_ALLOCATE(pxlpsi(1:ns, 1:dim))

      phase = allocated(pj(iatom)%phase)

      if(phase) then
#ifdef R_TCOMPLEX
        do idim = 1, dim 
          lpsi(1:ns, idim) = psi(map(1:ns), idim)*pj(iatom)%phase(1:ns, 1, ik) 
        end do
#else
        ! Phase not allowed for real functions
        ASSERT(.false.)
#endif
      else
        do idim = 1, dim 
          lpsi(1:ns, idim) = psi(map(1:ns), idim)
        end do
      end if

      !V_nl |psi>
      call X(project_sphere)(mesh, pj(iatom), dim, lpsi, plpsi)

      do idir = 1, mesh%sb%dim 
        ! x V_nl |psi>
        do idim = 1, dim
          ! x V_nl |psi>
          xplpsi(1:ns, idim) = smx(1:ns, idir) * plpsi(1:ns, idim)
          ! x |psi>
          xlpsi(1:ns, idim) = smx(1:ns, idir) * lpsi(1:ns, idim)
        end do
        ! V_nl x |psi>
        call X(project_sphere)(mesh, pj(iatom), dim, xlpsi, pxlpsi)
     
        ! |cpsi> += x V_nl |psi> - V_nl x |psi> 
        if(phase) then
#ifdef R_TCOMPLEX
          do idim = 1, dim
           cpsi(map(1:ns), idir, idim) = cpsi(map(1:ns), idir, idim) + &
             (xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)) * R_CONJ(pj(iatom)%phase(1:ns, 1, ik))
          end do
#else
          ! Phase not allowed for real functions
          ASSERT(.false.)
#endif
        else
          do idim = 1, dim
           cpsi(map(1:ns), idir, idim) = cpsi(map(1:ns), idir, idim) + xplpsi(1:ns, idim) - pxlpsi(1:ns, idim)
          end do
        end if

      end do !idir

      SAFE_DEALLOCATE_A(lpsi)
      SAFE_DEALLOCATE_A(plpsi)
      SAFE_DEALLOCATE_A(xlpsi)
      SAFE_DEALLOCATE_A(xplpsi)
      SAFE_DEALLOCATE_A(pxlpsi)
    end if
  end do
  call profiling_out(prof)
  POP_SUB(X(projector_commute_r_allatoms_alldir))

end subroutine X(projector_commute_r_allatoms_alldir)

!------------------------------------------------------------------------------
!> This function calculates |cpsi> += r * V_nl  |psi>
subroutine X(r_project_psi)(pj, mesh, dim, ik, psi, cpsi)
  type(projector_t), target, intent(in)     :: pj
  type(mesh_t),              intent(in)     :: mesh
  integer,                   intent(in)     :: dim
  integer,                   intent(in)     :: ik
  R_TYPE,                    intent(in)     :: psi(:, :)
  R_TYPE,                    intent(inout)  :: cpsi(:,:,:)

  integer ::  ns, idim, sb_dim, isb_dim
#ifdef R_TCOMPLEX
  integer :: ip
#endif
  R_TYPE, allocatable :: lpsi(:, :), xplpsi(:, :), xplpsi_t(:, :, :)
  integer, pointer :: map(:)
  FLOAT,   pointer :: smx(:, :)
  type(profile_t), save :: prof

  PUSH_SUB(X(r_project_psi))
  call profiling_in(prof, TOSTRING(X(P_PROJECT_PSI)))

  sb_dim = mesh%sb%dim
  
  if(pj%type /= PROJ_NONE) then

    ns = pj%sphere%np
    map => pj%sphere%map
    smx => pj%sphere%x

    SAFE_ALLOCATE(  lpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(xplpsi(1:ns, 1:dim))
    SAFE_ALLOCATE(xplpsi_t(1:ns, 1:sb_dim+1, 1:dim))

    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(map(1:ns), idim)*pj%phase(1:ns, 1, ik)
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
    else
      do idim = 1, dim
        lpsi(1:ns, idim) = psi(map(1:ns), idim)
      end do
    end if

    ! x V_nl |psi>
    call X(project_sphere)(mesh, pj, dim, lpsi, xplpsi)
    do idim = 1, dim

       do isb_dim = 1,sb_dim
          xplpsi_t(1:ns, isb_dim, idim) = smx(1:ns, isb_dim) * xplpsi(1:ns, idim)
       end do

       xplpsi_t(1:ns, sb_dim+1, idim) = xplpsi(1:ns, idim)       
    end do

    ! |cpsi> += x V_nl |psi> 
    if (allocated(pj%phase)) then
#ifdef R_TCOMPLEX
      do idim = 1, dim
        do ip = 1, ns
          cpsi(map(ip), 1:sb_dim+1, idim) = cpsi(map(ip), 1:sb_dim+1, idim) + &
            xplpsi_t(ip, 1:sb_dim+1, idim) * R_CONJ(pj%phase(ip, 1, ik))
        end do
      end do
#else
      ! Phase not allowed for real functions
      ASSERT(.false.)
#endif
   else
      do idim = 1, dim
         cpsi(map(1:ns), 1:sb_dim+1, idim) = cpsi(map(1:ns), 1:sb_dim+1, idim) &
              + xplpsi_t(1:ns, 1:sb_dim+1, idim)
      end do
    end if

    SAFE_DEALLOCATE_A(lpsi)
    SAFE_DEALLOCATE_A(xplpsi)
    SAFE_DEALLOCATE_A(xplpsi_t)
  end if
  call profiling_out(prof)
  POP_SUB(X(r_project_psi))

end subroutine X(r_project_psi)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
