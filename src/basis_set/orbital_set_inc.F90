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

subroutine X(orbital_set_get_coefficients)(os, psi, ik, d_dim, has_phase, dot)
  type(orbital_set_t),  intent(in) :: os
  R_TYPE,               intent(in) :: psi(:,:)
  integer,              intent(in) :: ik, d_dim
  logical,              intent(in) :: has_phase !True if the wavefunction has an associated phase
  R_TYPE,            intent(inout) :: dot(:)

  integer :: im
  type(profile_t), save :: prof

  call profiling_in(prof, "ORBSET_GET_COEFFICIENTS")

  PUSH_SUB(X(orbital_set_get_coefficients))

  do im = 1, os%norbs
    !If we need to add the phase, we explicitly do the operation using the sphere
    if(has_phase) then
#ifdef R_TCOMPLEX
      dot(im) = submesh_to_mesh_dotp(os%sphere, d_dim, os%orbitals(im)%eorb(1:os%sphere%np,ik),&
                              psi(1:os%sphere%mesh%np, 1:d_dim))
#endif
    else
      dot(im) = submesh_to_mesh_dotp(os%sphere, d_dim, os%orbitals(im)%X(orb)(1:os%sphere%np),&
                               psi(1:os%sphere%mesh%np, 1:d_dim))
    end if
  end do

  POP_SUB(X(orbital_set_get_coefficients))
  call profiling_out(prof)
end subroutine X(orbital_set_get_coefficients)

