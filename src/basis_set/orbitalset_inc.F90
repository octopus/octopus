!! Copyright (C) 2017 N. Tancogne-Dejean 
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

subroutine X(orbitalset_get_coefficients)(os, ndim, psi, ik, has_phase, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  R_TYPE,               intent(in) :: psi(:,:)
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:)

  integer :: im, ip, idim, idim_orb
  type(profile_t), save :: prof, prof_reduce
  R_TYPE, allocatable :: spsi(:)

  call profiling_in(prof, "ORBSET_GET_COEFFICIENTS")

  PUSH_SUB(X(orbitalset_get_coefficients))

  do im = 1, os%norbs
    !If we need to add the phase, we explicitly do the operation using the sphere
    if(has_phase) then
#ifdef R_TCOMPLEX
      if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. os%submeshforperiodic) then
        do idim = 1, ndim
          idim_orb = min(idim,os%ndim)
          dot(idim,im) = zmf_dotp(os%sphere%mesh, os%eorb_mesh(1:os%sphere%mesh%np,idim_orb,im,ik),&
                              psi(1:os%sphere%mesh%np,idim), reduce=.false.)
        end do
      else
        do idim = 1, ndim
          idim_orb = min(idim,os%ndim)
          dot(idim, im) = submesh_to_mesh_dotp(os%sphere, os%eorb_submesh(1:os%sphere%np,idim_orb,im,ik),&
                              psi(1:os%sphere%mesh%np, idim), reduce=.false.)
        end do
      endif 
#endif
    else
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)
        dot(idim,im) = submesh_to_mesh_dotp(os%sphere, os%X(orb)(1:os%sphere%np,idim_orb,im),&
                               psi(1:os%sphere%mesh%np, idim), reduce=.false.)
      end do
    end if
  end do

  if(os%sphere%mesh%parallel_in_domains) then
    call profiling_in(prof_reduce, "ORBSET_GET_COEFF_REDUCE")
    call comm_allreduce(os%sphere%mesh%mpi_grp%comm, dot) 
    call profiling_out(prof_reduce)
  end if

  POP_SUB(X(orbitalset_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbitalset_get_coefficients)


subroutine X(orbitalset_get_coeff_batch)(os, ndim, psib, ik, has_phase, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(batch_t),        intent(in) :: psib
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:,:)

  integer :: ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: spsi(:,:), psi(:,:)

  call profiling_in(prof, "ORBSET_GET_COEFF_BATCH")

  PUSH_SUB(X(orbitalset_get_coeff_batch))

  SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
  do ist = 1, psib%nst
    call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
    call X(orbitalset_get_coefficients)(os, ndim, psi, ik, has_phase, dot(1:ndim,1:os%norbs,ist))
  end do
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(orbitalset_get_coeff_batch))
  call profiling_out(prof)
end subroutine X(orbitalset_get_coeff_batch)

subroutine X(orbitalset_add_to_psi)(os, ndim, psi, ik, has_phase, weight)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  R_TYPE,            intent(inout) :: psi(:,:)
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,               intent(in) :: weight(:,:)

  integer :: im, ip, idim, idim_orb
  type(profile_t), save :: prof

  call profiling_in(prof, "ORBSET_ADD_TO_PSI")

  PUSH_SUB(X(orbitalset_add_to_psi))

  do im = 1, os%norbs
    do idim = 1, ndim
      idim_orb = min(idim,os%ndim)
      !In case of phase, we have to apply the conjugate of the phase here
      if(has_phase) then
#ifdef R_TCOMPLEX
        if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. os%submeshforperiodic) then
          call lalg_axpy(os%sphere%mesh%np, weight(idim,im), os%eorb_mesh(1:os%sphere%mesh%np,idim_orb,im,ik), &
                                  psi(1:os%sphere%mesh%np,idim))
        else
          call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,idim_orb,im,ik), &
                                  psi(1:os%sphere%mesh%np,idim), weight(idim,im))
        endif
#endif
      else
        call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np,idim_orb, im), &
                                  psi(1:os%sphere%mesh%np,idim), weight(idim,im))
      end if
    end do
  end do
 
  POP_SUB(X(orbitalset_add_to_psi))
  call profiling_out(prof)
end subroutine X(orbitalset_add_to_psi)


subroutine X(orbitalset_add_to_batch)(os, ndim, psib, ik, has_phase, weight)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(batch_t),     intent(inout) :: psib
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,               intent(in) :: weight(:,:)

  integer :: ip, iorb, ii, ist, idim, bind, idim_orb
  integer :: idim1, idim2, idim3, idim4
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:), sorb(:)
  R_TYPE :: tmp
  integer :: block_size, size, sp, ep

  call profiling_in(prof, "ORBSET_ADD_TO_BATCH")

  PUSH_SUB(X(orbitalset_add_to_batch))

  ! This routine uses blocking to optimize cache usage.   
  block_size = hardware%X(block_size)

  call batch_pack_was_modified(psib)

  if(os%sphere%mesh%use_curvilinear .or. batch_status(psib) == BATCH_CL_PACKED) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      do iorb = 1, os%norbs
        do idim = 1, ndim
         idim_orb = min(idim,os%ndim)
         bind = batch_ist_idim_to_linear(psib, (/ist, idim/))
          !In case of phase, we have to apply the conjugate of the phase here
          if(has_phase) then
#ifdef R_TCOMPLEX
            if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. os%submeshforperiodic) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb,bind),os%eorb_mesh(1:os%sphere%mesh%np,idim_orb,iorb,ik), &
                                   psi(1:os%sphere%mesh%np,idim))
            else
              call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,idim_orb,iorb,ik), &
                                  psi(1:os%sphere%mesh%np,idim), weight(iorb,bind))
            endif
#endif
          else
            call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, idim_orb, iorb), &
                                    psi(1:os%sphere%mesh%np,idim), weight(iorb,bind))
          end if
        end do !idim
      end do ! iorb
      call batch_set_state(psib, ist, os%sphere%mesh%np, psi)
    end do
    SAFE_DEALLOCATE_A(psi)
    !
  else
    !
    select case(batch_status(psib))
    case(BATCH_NOT_PACKED)
      !
      if(has_phase) then
#ifdef R_TCOMPLEX
        if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. os%submeshforperiodic) then
          do ist = 1, psib%nst_linear
            idim = min(batch_linear_to_idim(psib,ist),os%ndim)
            do sp = 1, os%sphere%mesh%np, block_size
              ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
              do iorb = 1, os%norbs
                tmp =  weight(iorb,ist)
                forall(ip = sp:ep)
                psib%states_linear(ist)%zpsi(ip) = &
                   psib%states_linear(ist)%zpsi(ip) + tmp*os%eorb_mesh(ip,idim,iorb,ik)
                end forall
              end do
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(batch_linear_to_idim(psib,ist),os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb,ist)
              forall (ip = 1:os%sphere%np)
                sorb(ip) = sorb(ip) + os%eorb_submesh(ip,idim,iorb,ik)*tmp
              end forall
            end do
            do ip = 1,os%sphere%np
              psib%states_linear(ist)%zpsi(os%sphere%map(ip)) = &
                    psib%states_linear(ist)%zpsi(os%sphere%map(ip)) + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
        SAFE_ALLOCATE(sorb(1:os%sphere%np))
        do ist = 1, psib%nst_linear
          idim = min(batch_linear_to_idim(psib,ist),os%ndim)
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%X(orb)(ip,idim,iorb)*tmp
            end forall
          end do
          do ip = 1,os%sphere%np
            psib%states_linear(ist)%X(psi)(os%sphere%map(ip)) = &
                  psib%states_linear(ist)%X(psi)(os%sphere%map(ip)) + sorb(ip)
          end do
        end do
        SAFE_DEALLOCATE_A(sorb)
      end if

    case(BATCH_PACKED)
      !
      if(has_phase) then
#ifdef R_TCOMPLEX
        if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. os%submeshforperiodic) then
          do iorb = 1, os%norbs
            do sp = 1, os%sphere%mesh%np, block_size
              ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
              do ist = 1, psib%nst_linear - 4 + 1, 4
                idim1 = min(batch_linear_to_idim(psib,ist), os%ndim)
                idim2 = min(batch_linear_to_idim(psib,ist+1), os%ndim)
                idim3 = min(batch_linear_to_idim(psib,ist+2), os%ndim)
                idim4 = min(batch_linear_to_idim(psib,ist+3), os%ndim)

                do ip = sp,ep !1:os%sphere%mesh%np)
                  psib%pack%zpsi(ist,ip) = &
                     psib%pack%zpsi(ist,ip) + weight(iorb,ist)*os%eorb_mesh(ip,idim1,iorb,ik)
                  psib%pack%zpsi(ist+1,ip) = &
                     psib%pack%zpsi(ist+1,ip) + weight(iorb,ist+1)*os%eorb_mesh(ip,idim2,iorb,ik)
                  psib%pack%zpsi(ist+2,ip) = &
                     psib%pack%zpsi(ist+2,ip) + weight(iorb,ist+2)*os%eorb_mesh(ip,idim3,iorb,ik)
                  psib%pack%zpsi(ist+3,ip) = &
                    psib%pack%zpsi(ist+3,ip) + weight(iorb,ist+3)*os%eorb_mesh(ip,idim4,iorb,ik)
                end do !forall 
              end do

              do ist = ist, psib%nst_linear
                idim = min(batch_linear_to_idim(psib,ist), os%ndim)
                ! forall(ip = 1:os%sphere%mesh%np)
                do ip=sp,ep
                  psib%pack%zpsi(ist,ip) = &
                     psib%pack%zpsi(ist,ip) + weight(iorb,ist)*os%eorb_mesh(ip,idim,iorb,ik)
                end do! forall
              end do
            end do
          end do

        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(batch_linear_to_idim(psib,ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb,ist)
              forall (ip = 1:os%sphere%np)
                sorb(ip) = sorb(ip) + os%eorb_submesh(ip,idim,iorb,ik)*tmp
              end forall
            end do
            do ip = 1,os%sphere%np
              psib%pack%zpsi(ist,os%sphere%map(ip)) = psib%pack%zpsi(ist,os%sphere%map(ip)) &
                                + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
       SAFE_ALLOCATE(sorb(1:os%sphere%np))
       do ist = 1, psib%nst_linear
          idim = min(batch_linear_to_idim(psib,ist), os%ndim)
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%X(orb)(ip,idim,iorb)*tmp
            end forall
          end do
          do ip = 1,os%sphere%np
            psib%pack%X(psi)(ist,os%sphere%map(ip)) = psib%pack%X(psi)(ist,os%sphere%map(ip)) &
                              + sorb(ip)
          end do
        end do
        SAFE_DEALLOCATE_A(sorb)
      end if
    end select
  end if

  POP_SUB(X(orbitalset_add_to_batch))
  call profiling_out(prof)
end subroutine X(orbitalset_add_to_batch)


