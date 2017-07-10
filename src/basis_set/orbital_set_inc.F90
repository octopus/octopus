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
  R_TYPE,            intent(inout) :: dot(:)

  integer :: im, ip, idim
  type(profile_t), save :: prof
  R_TYPE, allocatable :: spsi(:)

  call profiling_in(prof, "ORBSET_GET_COEFFICIENTS")

  PUSH_SUB(X(orbital_set_get_coefficients))

  if(st_d%ispin == SPINORS) then
   message(1) = "orbital_set_get_coefficients is not implemented with spinors."
   call messages_fatal(1)
  end if

  if(os%sphere%mesh%use_curvilinear) then
    do im = 1, os%norbs
      !If we need to add the phase, we explicitly do the operation using the sphere
      if(has_phase) then
#ifdef R_TCOMPLEX
        dot(im) = submesh_to_mesh_dotp(os%sphere, st_d%dim, os%eorb(1:os%sphere%np,im,ik),&
                              psi(1:os%sphere%mesh%np, 1:st_d%dim))
#endif
      else
        dot(im) = submesh_to_mesh_dotp(os%sphere, st_d%dim, os%X(orb)(1:os%sphere%np,im),&
                               psi(1:os%sphere%mesh%np, 1:st_d%dim))
      end if
    end do
    !
  else
    !
    SAFE_ALLOCATE(spsi(1:os%sphere%np))
    !
    call X(submesh_copy_from_mesh)(os%sphere, psi(1:os%sphere%mesh%np,1), spsi(1:os%sphere%np))
    !
    !If we need to add the phase, we explicitly do the operation using the sphere
    if(has_phase) then
#ifdef R_TCOMPLEX
      call blas_gemv('C', os%sphere%np, os%norbs, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
        os%eorb(1,1,ik),  os%sphere%np, spsi(1), 1, R_TOTYPE(M_ZERO), dot(1), 1)
#endif
    else
      call blas_gemv('T', os%sphere%np, os%norbs, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
        os%X(orb)(1,1),  os%sphere%np, spsi(1), 1, R_TOTYPE(M_ZERO), dot(1), 1)
    end if
    !
    SAFE_DEALLOCATE_A(spsi)
  end if

  POP_SUB(X(orbital_set_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbital_set_get_coefficients)


subroutine X(orbital_set_get_coeff_batch)(os, st_d, psib, ik, has_phase, dot)
  type(orbital_set_t),  intent(in) :: os
  type(states_dim_t),   intent(in) :: st_d
  type(batch_t),        intent(in) :: psib
  integer,              intent(in) :: ik
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:,:)

  integer :: ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: spsi(:,:), psi(:,:)

  call profiling_in(prof, "ORBSET_GET_COEFF_BATCH")

  PUSH_SUB(X(orbital_set_get_coeff_batch))

  if(st_d%ispin == SPINORS) then
   message(1) = "orbital_set_get_coeff_batch is not implemented with spinors."
   call messages_fatal(1)
  end if

  if(os%sphere%mesh%use_curvilinear .or. batch_status(psib) == BATCH_CL_PACKED) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:st_d%dim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      call X(orbital_set_get_coefficients)(os, st_d, psi, ik, has_phase, dot(1:os%norbs,ist))
    end do
    SAFE_DEALLOCATE_A(psi)
    !
  else
    !
    select case(batch_status(psib))
    case(BATCH_NOT_PACKED)
      SAFE_ALLOCATE(spsi(1:os%sphere%np, 1:psib%nst))

      call X(submesh_copy_from_mesh_batch)(os%sphere, psib, spsi)
      !
      !If we need to add the phase, we explicitly do the operation using the sphere
      if(has_phase) then
#ifdef R_TCOMPLEX
        call blas_gemm('C', 'N', os%norbs, psib%nst, os%sphere%np, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
          os%eorb(1,1,ik),  os%sphere%np, spsi(1,1), os%sphere%np, R_TOTYPE(M_ZERO), dot(1,1), os%norbs)
#endif
      else
        call blas_gemm('T', 'N', os%norbs, psib%nst, os%sphere%np, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
          os%X(orb)(1,1),  os%sphere%np, spsi(1,1), os%sphere%np, R_TOTYPE(M_ZERO), dot(1,1), os%norbs)
      end if

    case(BATCH_PACKED)
      SAFE_ALLOCATE(spsi(1:psib%nst, 1:os%sphere%np))
      !
      call X(submesh_copy_from_mesh_batch)(os%sphere, psib, spsi)
      !
      !If we need to add the phase, we explicitly do the operation using the sphere
      if(has_phase) then
#ifdef R_TCOMPLEX
        call blas_gemm('C', 'T', os%norbs, psib%nst, os%sphere%np, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
          os%eorb(1,1,ik),  os%sphere%np, spsi(1,1), psib%nst, R_TOTYPE(M_ZERO), dot(1,1), os%norbs)
#endif
      else
        call blas_gemm('T', 'T', os%norbs, psib%nst, os%sphere%np, R_TOTYPE(os%sphere%mesh%vol_pp(1)), &
          os%X(orb)(1,1),  os%sphere%np, spsi(1,1), psib%nst, R_TOTYPE(M_ZERO), dot(1,1), os%norbs)
      end if
    end select
    !
    SAFE_DEALLOCATE_A(spsi)
  end if

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

  if(st_d%ispin == SPINORS) then
   message(1) = "orbital_set_add_to_batch is not implemented with spinors."
   call messages_fatal(1)
  end if

  !
  do im = 1, os%norbs
    !In case of phase, we have to apply the conjugate of the phase here
    if(has_phase) then
#ifdef R_TCOMPLEX
      call submesh_add_to_mesh(os%sphere, os%eorb(1:os%sphere%np,im,ik), &
                                psi(1:os%sphere%mesh%np), weight(im))
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

  integer :: im, ip, iorb, ii, ist
  type(profile_t), save :: prof
  R_TYPE, allocatable :: psi(:,:), sorb(:)
  R_TYPE :: tmp

  call profiling_in(prof, "ORBSET_ADD_TO_BATCH")

  PUSH_SUB(X(orbital_set_add_to_batch))

  if(st_d%ispin == SPINORS) then
   message(1) = "orbital_set_add_to_batch is not implemented with spinors."
   call messages_fatal(1)
  end if

  if(batch_is_packed(psib)) then
    call batch_pack_was_modified(psib)
  end if


  if(os%sphere%mesh%use_curvilinear .or. batch_status(psib) == BATCH_CL_PACKED) then
    !
    SAFE_ALLOCATE(psi(1:os%sphere%mesh%np, 1:st_d%dim))
    do ist = 1, psib%nst
      call batch_get_state(psib, ist, os%sphere%mesh%np, psi)
      do im = 1, os%norbs
        !In case of phase, we have to apply the conjugate of the phase here
        if(has_phase) then
#ifdef R_TCOMPLEX
          call submesh_add_to_mesh(os%sphere, os%eorb(1:os%sphere%np,im,ik), &
                                    psi(1:os%sphere%mesh%np, 1), weight(im,ist))
#endif
        else
          call submesh_add_to_mesh(os%sphere, os%X(orb)(1:os%sphere%np, im), &
                                    psi(1:os%sphere%mesh%np, 1), weight(im,ist))
        end if
      end do
      call batch_set_state(psib, ist, os%sphere%mesh%np, psi)
    end do
    SAFE_DEALLOCATE_A(psi)
    !
  else
    !
    SAFE_ALLOCATE(sorb(1:os%sphere%np))
    select case(batch_status(psib))
    case(BATCH_NOT_PACKED)
      !
      if(has_phase) then
#ifdef R_TCOMPLEX
        do ist = 1, psib%nst_linear
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%eorb(ip,iorb,ik)*tmp
            end forall
          end do
          do ip = 1,os%sphere%np
            psib%states_linear(ist)%X(psi)(os%sphere%map(ip)) = &
                  psib%states_linear(ist)%X(psi)(os%sphere%map(ip)) + sorb(ip)
          end do
        end do
#endif
      else
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
      end if

    case(BATCH_PACKED)
      if(has_phase) then
#ifdef R_TCOMPLEX
        do ist = 1, psib%nst_linear
          sorb(:) = R_TOTYPE(M_ZERO)
          do iorb = 1, os%norbs
            tmp = weight(iorb,ist)
            forall (ip = 1:os%sphere%np)
              sorb(ip) = sorb(ip) + os%eorb(ip,iorb,ik)*tmp
            end forall
          end do
          do ip = 1,os%sphere%np
            psib%pack%X(psi)(ist,os%sphere%map(ip)) = psib%pack%X(psi)(ist,os%sphere%map(ip)) &
                              + sorb(ip)
          end do
        end do
#endif
      else
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
      end if
    end select
    SAFE_DEALLOCATE_A(sorb) 
  end if

  POP_SUB(X(orbital_set_add_to_batch))
  call profiling_out(prof)
end subroutine X(orbital_set_add_to_batch)
