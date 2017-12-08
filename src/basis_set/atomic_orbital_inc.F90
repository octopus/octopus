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
  type(orbitalset_t),       intent(inout) :: os
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
    call datomic_orbital_get_submesh(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x, tmp)
    os%X(orb)(1:sm%np,1,orbind) = tmp(1:sm%np)
    SAFE_DEALLOCATE_A(tmp)
  #else
    call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,1,orbind))
  #endif
  else
    if(jj == ll+M_HALF) then
      mm = orbind - 2 - ll
      if(mm >= -ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,1,orbind))
        coeff = sqrt((ll+mm+M_ONE)/(M_TWO*ll+M_ONE)) 
        do is=1,sm%np
          os%X(orb)(is,1,orbind) = coeff*os%X(orb)(is,1,orbind)
        end do
      else
        os%X(orb)(1:sm%np,1,orbind) = M_ZERO
      end if
      if(mm < ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm+1, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,2,orbind))
        coeff = sqrt((ll-mm)/(M_TWO*ll+M_ONE))                           
        do is=1,sm%np
          os%X(orb)(is,2,orbind) = coeff*os%X(orb)(is,2,orbind)
        end do
      else
       os%X(orb)(1:sm%np,2,orbind) = M_ZERO
      end if
    else
      mm = orbind - ll
      call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm, 1, geo%atom(iatom)%x,&
                                        os%X(orb)(1:sm%np,2,orbind))
      coeff = -sqrt((ll+mm)/(M_TWO*ll+M_ONE))                           
      do is=1,sm%np
        os%X(orb)(is,2,orbind) = coeff*os%X(orb)(is,2,orbind)
      end do
      if(mm > -ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm-1, 1, geo%atom(iatom)%x,&
                                         os%X(orb)(1:sm%np,1,orbind))
        coeff = sqrt((ll-mm+M_ONE)/(M_TWO*ll+M_ONE))      
        do is=1,sm%np
          os%X(orb)(is,1,orbind) = coeff*os%X(orb)(is,1,orbind)
        end do
      else
       os%X(orb)(1:sm%np,1,orbind) = M_ZERO
      end if
    end if

  end if

  POP_SUB(X(get_atomic_orbital))

end subroutine X(get_atomic_orbital)



  ! ---------------------------------------------------------
  subroutine X(atomic_orbital_get_submesh)(species, submesh, ii, ll, mm, ispin, pos, phi, derivative)
    type(species_t), target, intent(in)  :: species       !< The species.
    type(submesh_t),         intent(in)  :: submesh    !< The submesh descriptor where the orbital will be calculated.
    integer,                 intent(in)  :: ii
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    integer,                 intent(in)  :: ispin      !< The spin index.
    FLOAT,                   intent(in)  :: pos(:)     !< The position of the atom.
    R_TYPE,                  intent(out) :: phi(:)     !< The function defined in the mesh where the orbitals is returned.
    logical,       optional, intent(in)  :: derivative !< If present and .true. returns the derivative of the orbital.

    integer :: ip, nn(3), idir
    FLOAT :: sqrtw, ww
    R_TYPE, allocatable :: ylm(:)
    type(ps_t), pointer :: ps
    type(spline_t) :: dur
    logical :: derivative_
    
    if(submesh%np == 0) return

    PUSH_SUB(X(atomic_orbital_get_submesh))

    derivative_ = optional_default(derivative, .false.)

    ASSERT(ubound(phi, dim = 1) >= submesh%np)

    if(species_represents_real_atom(species) .and. submesh%mesh%sb%dim == 3) then
      ps => species_ps(species)
      
      forall(ip = 1:submesh%np) phi(ip) = submesh%x(ip, 0)

      if(species_is_ps(species)) then
        if(.not. derivative_) then
          call spline_eval_vec(ps%ur(ii, ispin), submesh%np, phi)
        else
          call spline_init(dur)
          call spline_der(ps%ur(ii, ispin), dur)
          call spline_eval_vec(dur, submesh%np, phi)
          call spline_end(dur)
        end if
      else
        ! FIXME: cache result somewhat. e.g. re-use result for each m. and use recursion relation.
        do ip = 1, submesh%np
          ww = species_zval(species)*submesh%x(ip, 0)/ii
          phi(ip) = sqrt( (2*species_zval(species)/ii)**3 * factorial(ii - ll - 1) / (2*ii*factorial(ii+ll)) ) * &
            exp(-ww) * (2 * ww)**ll * loct_sf_laguerre_n(ii-ll-1, real(2*ll + 1, REAL_PRECISION), 2*ww)
        end do
      end if

      SAFE_ALLOCATE(ylm(1:submesh%np))

#ifdef R_TCOMPLEX
      ! complex spherical harmonics. FIXME: vectorize
      do ip = 1, submesh%np
        call ylmr(submesh%x(ip, 1), submesh%x(ip, 2), submesh%x(ip, 3), ll, mm, ylm(ip))
      end do
#else
      ! real spherical harmonics
      call loct_ylm(submesh%np, submesh%x(1, 1), submesh%x(1, 2), submesh%x(1, 3), ll, mm, ylm(1))
#endif      

      do ip = 1, submesh%np
        phi(ip) = phi(ip)*ylm(ip)
      end do

      SAFE_DEALLOCATE_A(ylm)

      nullify(ps)

    else
      
      ASSERT(.not. derivative_)
      ! Question: why not implemented derivatives here?
      ! Answer: because they are linearly dependent with lower-order Hermite polynomials.

      ww = species_omega(species)
      sqrtw = sqrt(ww)

      ! FIXME: this is a pretty dubious way to handle l and m quantum numbers. Why not use ylm?
      nn = (/ii, ll, mm/)

      do ip = 1, submesh%np
        phi(ip) = exp(-ww*submesh%x(ip, 0)**2/M_TWO)
        do idir = 1, submesh%mesh%sb%dim
          phi(ip) = phi(ip) * hermite(nn(idir) - 1, submesh%x(ip, idir)*sqrtw)
        end do
      end do
      
    end if

    POP_SUB(X(atomic_orbital_get_submesh))
  end subroutine X(atomic_orbital_get_submesh)

