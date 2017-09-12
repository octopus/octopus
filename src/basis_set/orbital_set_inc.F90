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
!!
!! $Id$

subroutine X(orbital_set_get_coefficients)(os, st_d, psi, ik, has_phase, dot)
  type(orbital_set_t),  intent(in) :: os
  type(states_dim_t),   intent(in) :: st_d
  R_TYPE,               intent(in) :: psi(:,:)
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:)

  integer :: im, ip, idim
  type(profile_t), save :: prof
  R_TYPE, allocatable :: spsi(:)

  call profiling_in(prof, "ORBSET_GET_COEFFICIENTS")

  PUSH_SUB(X(orbital_set_get_coefficients))

  do im = 1, os%norbs
    !If we need to add the phase, we explicitly do the operation using the sphere
    if(has_phase) then
#ifdef R_TCOMPLEX
      if(simul_box_is_periodic(os%sphere%mesh%sb)) then
        do idim = 1,st_d%dim
          dot(idim,im) = zmf_dotp(os%sphere%mesh, os%eorb_mesh(1:os%sphere%mesh%np,im,ik),&
                              psi(1:os%sphere%mesh%np,idim))
        end do
      else
        do idim = 1,st_d%dim
          dot(idim, im) = submesh_to_mesh_dotp(os%sphere, os%eorb_submesh(1:os%sphere%np,im,ik),&
                              psi(1:os%sphere%mesh%np, idim))
        end do
      endif 
#endif
    else
      do idim = 1,st_d%dim
        dot(idim,im) = submesh_to_mesh_dotp(os%sphere, os%X(orb)(1:os%sphere%np,im),&
                               psi(1:os%sphere%mesh%np, idim))
      end do
    end if
  end do

  POP_SUB(X(orbital_set_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbital_set_get_coefficients)


subroutine X(orbital_set_get_coeff_batch)(os, st_d, psib, ik, has_phase, dot)
  type(orbital_set_t),  intent(in) :: os
  type(states_dim_t),   intent(in) :: st_d
  type(batch_t),        intent(in) :: psib
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:,:)

  integer :: ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: spsi(:,:), psi(:,:)

  call profiling_in(prof, "ORBSET_GET_COEFF_BATCH")

  PUSH_SUB(X(orbital_set_get_coeff_batch))

  SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:st_d%dim))
  do ist = 1, psib%nst
    call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
    call X(orbital_set_get_coefficients)(os, st_d, psi, ik, has_phase, dot(1:st_d%dim,1:os%norbs,ist))
  end do
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(orbital_set_get_coeff_batch))
  call profiling_out(prof)
end subroutine X(orbital_set_get_coeff_batch)

subroutine X(orbital_set_add_to_psi)(os, st_d, psi, ik, has_phase, weight)
  type(orbital_set_t),  intent(in) :: os
  type(states_dim_t),   intent(in) :: st_d
  R_TYPE,            intent(inout) :: psi(:)
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,               intent(in) :: weight(:)

  integer :: im, ip
  type(profile_t), save :: prof

  call profiling_in(prof, "ORBSET_ADD_TO_PSI")

  PUSH_SUB(X(orbital_set_add_to_psi))

  do im = 1, os%norbs
    !In case of phase, we have to apply the conjugate of the phase here
    if(has_phase) then
#ifdef R_TCOMPLEX
      if(simul_box_is_periodic(os%sphere%mesh%sb)) then
        call lalg_axpy(os%sphere%mesh%np, weight(im), os%eorb_mesh(1:os%sphere%mesh%np,im,ik), &
                                psi(1:os%sphere%mesh%np))
      else
        call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,im,ik), &
                                psi(1:os%sphere%mesh%np), weight(im))
      endif
#endif
    else
      call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, im), &
                                psi(1:os%sphere%mesh%np), weight(im))
    end if
  end do

  POP_SUB(X(orbital_set_add_to_psi))
  call profiling_out(prof)
end subroutine X(orbital_set_add_to_psi)


subroutine X(orbital_set_add_to_batch)(os, st_d, psib, ik, has_phase, weight)
  type(orbital_set_t),  intent(in) :: os
  type(states_dim_t),   intent(in) :: st_d
  type(batch_t),     intent(inout) :: psib
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,               intent(in) :: weight(:,:)

  integer :: ip, iorb, ii, ist, idim, bind
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:), sorb(:)
  R_TYPE :: tmp

  call profiling_in(prof, "ORBSET_ADD_TO_BATCH")

  PUSH_SUB(X(orbital_set_add_to_batch))

  call batch_pack_was_modified(psib)

  if(os%sphere%mesh%use_curvilinear .or. batch_status(psib) == BATCH_CL_PACKED) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:st_d%dim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      do iorb = 1, os%norbs
        do idim = 1, st_d%dim
         bind = batch_ist_idim_to_linear(psib, (/ist, idim/))
          !In case of phase, we have to apply the conjugate of the phase here
          if(has_phase) then
#ifdef R_TCOMPLEX
            if(simul_box_is_periodic(os%sphere%mesh%sb)) then
              call lalg_axpy(os%sphere%mesh%np, weight(iorb,bind),os%eorb_mesh(1:os%sphere%mesh%np,iorb,ik), &
                                   psi(1:os%sphere%mesh%np,idim))
            else
              call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,iorb,ik), &
                                  psi(1:os%sphere%mesh%np,idim), weight(iorb,bind))
            endif
#endif
          else
            call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, iorb), &
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
        if(simul_box_is_periodic(os%sphere%mesh%sb)) then
          do ist = 1, psib%nst_linear
            do iorb = 1, os%norbs
              tmp =  weight(iorb,ist)
              forall(ip = 1:os%sphere%mesh%np)
                psib%states_linear(ist)%zpsi(ip) = &
                   psib%states_linear(ist)%zpsi(ip) + tmp*os%eorb_mesh(ip,iorb,ik)
              end forall
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb,ist)
              forall (ip = 1:os%sphere%np)
                sorb(ip) = sorb(ip) + os%eorb_submesh(ip,iorb,ik)*tmp
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
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%X(orb)(ip,iorb)*tmp
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
      if(has_phase) then
#ifdef R_TCOMPLEX
        if(simul_box_is_periodic(os%sphere%mesh%sb)) then
          do ist = 1, psib%nst_linear
            do iorb = 1, os%norbs
              tmp = weight(iorb,ist)
              forall(ip = 1:os%sphere%mesh%np)
                psib%pack%zpsi(ist,ip) = &
                   psib%pack%zpsi(ist,ip) + tmp*os%eorb_mesh(ip,iorb,ik)
              end forall 
            end do
          end do
        else
          SAFE_ALLOCATE(sorb(1:os%sphere%np))
          do ist = 1, psib%nst_linear
            sorb(:) = R_TOTYPE(M_ZERO)
            do iorb = 1, os%norbs
              tmp = weight(iorb,ist)
              forall (ip = 1:os%sphere%np)
                sorb(ip) = sorb(ip) + os%eorb_submesh(ip,iorb,ik)*tmp
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
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%X(orb)(ip,iorb)*tmp
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

  POP_SUB(X(orbital_set_add_to_batch))
  call profiling_out(prof)
end subroutine X(orbital_set_add_to_batch)
