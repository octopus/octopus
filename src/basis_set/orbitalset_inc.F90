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
  R_TYPE, allocatable :: spsi(:,:)

  call profiling_in(prof, TOSTRING(X(ORBSET_GET_COEFFICIENTS)))

  PUSH_SUB(X(orbitalset_get_coefficients))

  if(os%submesh) then
    SAFE_ALLOCATE(spsi(1:os%sphere%np, 1:ndim))
    do idim = 1, ndim
      do ip = 1, os%sphere%np
        spsi(ip,idim) = psi(os%sphere%map(ip), idim)
      end do
    end do
  end if

  !If we need to add the phase, we explicitly do the operation using the sphere
  if(has_phase) then
#ifdef R_TCOMPLEX
    do im = 1, os%norbs
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)
        if(.not. os%submesh) then
          dot(idim,im) = zmf_dotp(os%sphere%mesh, os%eorb_mesh(1:os%sphere%mesh%np, im, idim_orb, ik),&
                              psi(1:os%sphere%mesh%np, idim), reduce=.false.)
        else
          dot(idim, im) = zmf_dotp(os%sphere%mesh, os%eorb_submesh(1:os%sphere%np, idim_orb, im, ik),&
                              spsi(1:os%sphere%np, idim), reduce=.false., np=os%sphere%np)
        end if
      end do
    end do
#endif
  else
    do im = 1, os%norbs
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)
        if(.not. os%submesh) then
          dot(idim,im) = X(mf_dotp)(os%sphere%mesh, os%X(orb)(1:os%sphere%mesh%np, idim_orb, im),&
                              psi(1:os%sphere%mesh%np, idim), reduce=.false.)
        else
          dot(idim,im) = X(mf_dotp)(os%sphere%mesh, os%X(orb)(1:os%sphere%np, idim_orb, im),&
                               spsi(1:os%sphere%np, idim), reduce=.false., np=os%sphere%np)
        end if
      end do
    end do
  end if

  if(os%sphere%mesh%parallel_in_domains) then
    call profiling_in(prof_reduce, TOSTRING(X(ORBSET_GET_COEFF_REDUCE)))
    call os%sphere%mesh%allreduce(dot)
    call profiling_out(prof_reduce)
  end if

  SAFE_DEALLOCATE_A(spsi)

  POP_SUB(X(orbitalset_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbitalset_get_coefficients)


subroutine X(orbitalset_get_coeff_batch)(os, ndim, psib, dot)
  type(orbitalset_t),   intent(in) :: os
  integer,              intent(in) :: ndim
  type(wfs_elec_t),     intent(in) :: psib
  R_TYPE,            intent(inout) :: dot(:,:,:)

  integer :: ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:)

  call profiling_in(prof, TOSTRING(X(ORBSET_GET_COEFF_BATCH)))

  PUSH_SUB(X(orbitalset_get_coeff_batch))

  SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
  do ist = 1, psib%nst
    call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
    call X(orbitalset_get_coefficients)(os, ndim, psi, psib%ik, psib%has_phase, dot(1:ndim,1:os%norbs,ist))
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
  R_TYPE,               intent(in) :: weight(:,:) !(min(ndim,os%ndim), os%norbs)

  integer :: im, idim, idim_orb
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(ORBSET_ADD_TO_PSI)))

  PUSH_SUB(X(orbitalset_add_to_psi))

  if(has_phase) then
#ifdef R_TCOMPLEX
    do im = 1, os%norbs
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)
        if(.not. os%submesh) then
          call lalg_axpy(os%sphere%mesh%np, weight(idim, im), os%eorb_mesh(1:os%sphere%mesh%np,im,idim_orb,ik), &
                                  psi(1:os%sphere%mesh%np,idim))
        else
          call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,idim_orb,im,ik), &
                                  psi(1:os%sphere%mesh%np, idim), weight(idim, im))
        end if
      end do
    end do
#endif
  else
    do im = 1, os%norbs
      do idim = 1, ndim
        idim_orb = min(idim,os%ndim)

        if(.not. os%submesh) then
          call lalg_axpy(os%sphere%mesh%np, weight(idim, im), os%X(orb)(1:os%sphere%mesh%np, idim_orb, im), &
                                  psi(1:os%sphere%mesh%np, idim))
        else
          call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, idim_orb, im), &
                                    psi(1:os%sphere%mesh%np, idim), weight(idim, im))
        end if
      end do
    end do
  end if

  POP_SUB(X(orbitalset_add_to_psi))
  call profiling_out(prof)
end subroutine X(orbitalset_add_to_psi)


subroutine X(orbitalset_add_to_batch)(os, ndim, psib, weight)
  type(orbitalset_t),   intent(in)    :: os
  integer,              intent(in)    :: ndim
  type(wfs_elec_t),     intent(inout) :: psib
  R_TYPE,               intent(in)    :: weight(:,:) !(os%norbs, psib%nst_linear)

  integer :: ip, iorb, ist, idim, bind, idim_orb
  integer :: idim1, idim2, idim3, idim4
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:), sorb(:)
  R_TYPE :: tmp
  integer :: block_size, size, sp, ep

  call profiling_in(prof, TOSTRING(X(ORBSET_ADD_TO_BATCH)))

  PUSH_SUB(X(orbitalset_add_to_batch))

  ! This routine uses blocking to optimize cache usage.   
  block_size = hardware%X(block_size)

  if(os%sphere%mesh%use_curvilinear .or. psib%status() == BATCH_DEVICE_PACKED) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:ndim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      !In case of phase, we have to apply the conjugate of the phase here
      if(psib%has_phase) then
#ifdef R_TCOMPLEX
        do idim = 1, ndim
          idim_orb = min(idim,os%ndim)
          bind = psib%ist_idim_to_linear((/ist, idim/))
          do iorb = 1, os%norbs
            if(.not. os%submesh) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb, bind), os%eorb_mesh(1:os%sphere%mesh%np, iorb, idim_orb, psib%ik), &
                psi(1:os%sphere%mesh%np, idim))
            else
              call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np, idim_orb, iorb, psib%ik), &
                psi(1:os%sphere%mesh%np,idim), weight(iorb,bind))
            end if
          end do
        end do
#endif
      else
        do iorb = 1, os%norbs
          do idim = 1, ndim
            idim_orb = min(idim,os%ndim)
            bind = psib%ist_idim_to_linear((/ist, idim/))
            if(.not. os%submesh) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb, bind), os%X(orb)(1:os%sphere%mesh%np,idim_orb,iorb), &
                                  psi(1:os%sphere%mesh%np,idim))
            else
              call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, idim_orb, iorb), &
                                  psi(1:os%sphere%mesh%np,idim), weight(iorb, bind))
            end if
          end do
        end do
      end if
      call batch_set_state(psib, ist, os%sphere%mesh%np, psi)
    end do
    SAFE_DEALLOCATE_A(psi)
    !
  else
    !
    select case(psib%status())
    case(BATCH_NOT_PACKED)
      !
      if(psib%has_phase) then
#ifdef R_TCOMPLEX
        if(.not. os%submesh) then
          do sp = 1, os%sphere%mesh%np, block_size
            size = min(block_size, os%sphere%mesh%np - sp + 1)
            do ist = 1, psib%nst_linear
              idim = min(psib%linear_to_idim(ist), os%ndim)
              call blas_gemv('N', size, os%norbs, R_TOTYPE(M_ONE), os%eorb_mesh(sp, 1, idim, psib%ik), &
              os%sphere%mesh%np, weight(1, ist), 1, R_TOTYPE(M_ONE), psib%zff_linear(sp, ist), 1)
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist),os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do sp = 1, os%sphere%np, block_size
              size = min(block_size, os%sphere%np - sp + 1)
              ep = sp - 1 + size
              sorb(sp:ep) = R_TOTYPE(M_ZERO)
              do iorb = 1, os%norbs
                call blas_axpy(size, weight(iorb, ist), os%eorb_submesh(sp, idim, iorb, psib%ik), 1, sorb(sp), 1)
              end do
            end do
            do ip = 1,os%sphere%np
              psib%zff_linear(os%sphere%map(ip), ist) = &
                    psib%zff_linear(os%sphere%map(ip), ist) + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
        if(.not. os%submesh) then
          do sp = 1, os%sphere%mesh%np, block_size
            size = min(block_size, os%sphere%mesh%np - sp + 1)
            do ist = 1, psib%nst_linear
              idim = min(psib%linear_to_idim(ist), os%ndim)
              call blas_gemv('N', size, os%norbs, R_TOTYPE(M_ONE), os%X(orb)(sp, idim, 1), &
                    os%sphere%mesh%np * os%ndim, weight(1, ist), 1, R_TOTYPE(M_ONE),            &
                       psib%X(ff_linear)(sp, ist), 1)
            end do
          end do 
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb, ist)
              do ip = 1, os%sphere%np
                sorb(ip) = sorb(ip) + os%X(orb)(ip, idim, iorb) * tmp
              end do
            end do
            do ip = 1, os%sphere%np
              psib%X(ff_linear)(os%sphere%map(ip), ist) = &
                    psib%X(ff_linear)(os%sphere%map(ip), ist) + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
      end if

    case(BATCH_PACKED)
      !
      if(psib%has_phase) then
#ifdef R_TCOMPLEX
        if(.not. os%submesh) then
          do sp = 1, os%sphere%mesh%np, block_size
            ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
            do iorb = 1, os%norbs
              do ist = 1, psib%nst_linear - 4 + 1, 4
                idim1 = min(psib%linear_to_idim(ist),   os%ndim)
                idim2 = min(psib%linear_to_idim(ist+1), os%ndim)
                idim3 = min(psib%linear_to_idim(ist+2), os%ndim)
                idim4 = min(psib%linear_to_idim(ist+3), os%ndim)

                do ip = sp,ep
                  psib%zff_pack(ist, ip) = &
                     psib%zff_pack(ist, ip) + weight(iorb, ist)*os%eorb_mesh(ip, iorb, idim1, psib%ik)
                  psib%zff_pack(ist + 1, ip) = &
                     psib%zff_pack(ist + 1, ip) + weight(iorb, ist + 1)*os%eorb_mesh(ip, iorb, idim2, psib%ik)
                  psib%zff_pack(ist + 2, ip) = &
                     psib%zff_pack(ist + 2, ip) + weight(iorb, ist + 2)*os%eorb_mesh(ip, iorb, idim3, psib%ik)
                  psib%zff_pack(ist + 3, ip) = &
                    psib%zff_pack(ist + 3, ip) + weight(iorb, ist + 3)*os%eorb_mesh(ip, iorb, idim4, psib%ik)
                end do
              end do

              do ist = ist, psib%nst_linear
                idim = min(psib%linear_to_idim(ist), os%ndim)
                do ip = sp, ep
                  psib%zff_pack(ist, ip) = &
                     psib%zff_pack(ist, ip) + weight(iorb, ist)*os%eorb_mesh(ip, iorb, idim, psib%ik)
                end do
              end do
            end do
          end do

        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb, ist)
              do ip = 1, os%sphere%np
                sorb(ip) = sorb(ip) + os%eorb_submesh(ip, idim, iorb, psib%ik) * tmp
              end do
            end do
            do ip = 1, os%sphere%np
              psib%zff_pack(ist, os%sphere%map(ip)) = psib%zff_pack(ist, os%sphere%map(ip)) &
                                + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
#endif
      else
        if(.not. os%submesh) then 
          do iorb = 1, os%norbs
            do sp = 1, os%sphere%mesh%np, block_size
              ep = sp - 1 + min(block_size, os%sphere%mesh%np - sp + 1)
              do ist = 1, psib%nst_linear - 4 + 1, 4
                idim1 = min(psib%linear_to_idim(ist)  , os%ndim)
                idim2 = min(psib%linear_to_idim(ist+1), os%ndim)
                idim3 = min(psib%linear_to_idim(ist+2), os%ndim)
                idim4 = min(psib%linear_to_idim(ist+3), os%ndim)

                do ip = sp, ep
                  psib%X(ff_pack)(ist, ip) = &
                     psib%X(ff_pack)(ist, ip)   + weight(iorb, ist  ) * os%X(orb)(ip, idim1, iorb)
                  psib%X(ff_pack)(ist+1, ip) = &
                     psib%X(ff_pack)(ist+1, ip) + weight(iorb, ist+1) * os%X(orb)(ip, idim2, iorb)
                  psib%X(ff_pack)(ist+2, ip) = &
                     psib%X(ff_pack)(ist+2, ip) + weight(iorb, ist+2) * os%X(orb)(ip, idim3, iorb)
                  psib%X(ff_pack)(ist+3, ip) = &
                     psib%X(ff_pack)(ist+3, ip) + weight(iorb, ist+3) * os%X(orb)(ip, idim4, iorb)
                end do
              end do

              do ist = ist, psib%nst_linear
                idim = min(psib%linear_to_idim(ist), os%ndim)
                do ip = sp, ep
                  psib%X(ff_pack)(ist, ip) = &
                     psib%X(ff_pack)(ist, ip) + weight(iorb, ist) * os%X(orb)(ip, idim, iorb)
                end do
              end do
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            idim = min(psib%linear_to_idim(ist), os%ndim)
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb, ist)
              do ip = 1, os%sphere%np
                sorb(ip) = sorb(ip) + os%X(orb)(ip, idim, iorb) * tmp
              end do
            end do
            do ip = 1, os%sphere%np
              psib%X(ff_pack)(ist,os%sphere%map(ip)) = psib%X(ff_pack)(ist,os%sphere%map(ip)) &
                                + sorb(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(sorb)
        end if
      end if
    end select
  end if

  POP_SUB(X(orbitalset_add_to_batch))
  call profiling_out(prof)
end subroutine X(orbitalset_add_to_batch)


