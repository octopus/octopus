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

! ---------------------------------------------------------
!> This routine returns the atomic orbital basis -- provided
!! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(get_atomic_orbital) (geo, mesh, sm, iatom, ii, ll, jj, os, orbind, radius, d_dim)
  type(mesh_t),             intent(in)    :: mesh
  type(geometry_t), target, intent(in)    :: geo
  type(submesh_t),          intent(inout) :: sm
  integer,                  intent(in)    :: iatom, ii, ll
  FLOAT,                    intent(in)    :: jj
  type(orbital_set_t),      intent(inout) :: os
  integer,                  intent(in)    :: orbind
  FLOAT,                    intent(in)    :: radius
  integer,                  intent(in)    :: d_dim

  type(species_t), pointer :: spec
  #ifdef R_TCOMPLEX
  FLOAT, allocatable :: tmp(:)
  #endif
  integer :: is, mm, kappa
  FLOAT :: mu, coeff

  PUSH_SUB(X(get_atomic_orbital))

  spec => geo%atom(iatom)%species

  if(sm%np == -1) then
    
    if(mesh%sb%box_shape == MINIMUM .and. radius > mesh%sb%rsize) then
      message(1) = "The radius of an orbital set is bigger than the radius of the simulatio box."
      message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
      write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
      call messages_fatal(3)
    end if
 
    if(mesh%sb%box_shape == SPHERE .or. mesh%sb%box_shape == CYLINDER) then
      if(radius > mesh%sb%rsize) then
       message(1) = "The radius of an orbital set is bigger than the radius of the simulatio box."
       message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
       write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
       call messages_fatal(3) 
      end if
      if(mesh%sb%box_shape == CYLINDER .and. radius > mesh%sb%xsize) then
       message(1) = "The radius of an orbital set is bigger than the length of the cylinder box."
       message(2) = "Increase the value of XLength or decrease the value of OrbitalsThreshold_LDAU."
       write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
       call messages_fatal(3) 
      end if
    end if 

    if(mesh%sb%box_shape == SPHERE ) then
      if(sqrt(sum(geo%atom(iatom)%x(1:mesh%sb%dim)**2)) + radius > mesh%sb%rsize) then
       message(1) = "An orbital set has points outside of the simulatio box."
       message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
       write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
       call messages_fatal(3)
      end if
    end if

    if(mesh%sb%box_shape == CYLINDER ) then
      if(sqrt(sum(geo%atom(iatom)%x(2:mesh%sb%dim)**2)) + radius > mesh%sb%rsize) then
       message(1) = "An orbital set has points outside of the simulatio box."
       message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
       write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
       call messages_fatal(3)
      end if
      if(abs(geo%atom(iatom)%x(1)) + radius > mesh%sb%xsize) then
       message(1) = "An orbital set has points outside of the simulatio box."
       message(2) = "Increase the value of Xlength or decrease the value of OrbitalsThreshold_LDAU."
       write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
       call messages_fatal(3)
      end if
    end if

 
    !We initialise the submesh corresponding to the orbital 
    call submesh_init(sm, mesh%sb, mesh, geo%atom(iatom)%x, radius)

  end if

  if(.not.associated(os%X(orb))) then
    SAFE_ALLOCATE(os%X(orb)(1:sm%np,1:os%ndim,1:os%norbs))
    os%X(orb)(:,:,:) = R_TOTYPE(M_ZERO)
  end if

  if(d_dim == 1) then

    mm = orbind-1-ll

    !We get the orbital from the pseudopotential
  #ifdef R_TCOMPLEX
    !In this case we want to get a real orbital and to store it in complex array
    SAFE_ALLOCATE(tmp(1:sm%np))
    call dspecies_get_orbital_submesh(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x, tmp)
    os%X(orb)(1:sm%np,1,orbind) = tmp(1:sm%np)
    SAFE_DEALLOCATE_A(tmp)
  #else
      call X(species_get_orbital_submesh)(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,1,orbind))
  #endif
  else
    !see for instance https://arxiv.org/pdf/1011.3433.pdf
    kappa = (ll-jj)*(M_TWO*jj+M_ONE)
    mu = orbind-1-abs(kappa)+M_HALF

    mm = int(mu-M_HALF)
    if(abs(mm) <= ll) then
      call X(species_get_orbital_submesh)(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,1,orbind))
      coeff = sqrt((kappa-mu+M_HALF)/(M_TWO*kappa+M_ONE)) 
      do is=1,sm%np
        os%X(orb)(is,1,orbind) = coeff*os%X(orb)(is,1,orbind)
      end do
    else
       os%X(orb)(1:sm%np,1,orbind) = M_ZERO
    end if

    mm = int(mu+M_HALF)
    if(abs(mm) <= ll) then
      call X(species_get_orbital_submesh)(spec, sm, ii, ll, mm, 2, geo%atom(iatom)%x,&
                                        os%X(orb)(1:sm%np,2,orbind))
      coeff = (-kappa/abs(kappa))*sqrt((kappa+mu+M_HALF)/(M_TWO*kappa+M_ONE))
      do is=1,sm%np
        os%X(orb)(is,2,orbind) = coeff*os%X(orb)(is,2,orbind)
      end do
    else
      os%X(orb)(1:sm%np,2,orbind) = M_ZERO
    end if

  end if

  POP_SUB(X(get_atomic_orbital))

end subroutine X(get_atomic_orbital)

