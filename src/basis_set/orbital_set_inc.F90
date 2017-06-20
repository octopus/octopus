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
   message(1) = "states_calc_projections is not implemented with spinors."
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

