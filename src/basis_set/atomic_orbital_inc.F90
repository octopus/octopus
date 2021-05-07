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
subroutine X(get_atomic_orbital) (ions, mesh, sm, iatom, ii, ll, jj, os, orbind, radius, d_dim, &
                                    use_mesh, normalize)
  type(mesh_t),             intent(in)    :: mesh
  type(ions_t),     target, intent(in)    :: ions
  type(submesh_t),          intent(inout) :: sm
  integer,                  intent(in)    :: iatom, ii, ll
  FLOAT,                    intent(in)    :: jj
  type(orbitalset_t),       intent(inout) :: os
  integer,                  intent(in)    :: orbind
  FLOAT,                    intent(in)    :: radius
  integer,                  intent(in)    :: d_dim
  logical,                  intent(in)    :: use_mesh
  logical,                  intent(in)    :: normalize

  type(species_t), pointer :: spec
  FLOAT, allocatable :: tmp(:)
  R_TYPE, allocatable :: ztmp(:,:)
  integer :: mm
  FLOAT :: coeff, norm

  PUSH_SUB(X(get_atomic_orbital))

  spec => ions%atom(iatom)%species

  if(sm%np == -1) then
    
    select type (box => mesh%sb%box)
    type is (box_minimum_t)
      if (radius > box%radius) then
        message(1) = "The radius of an orbital set is bigger than the radius of the simulatio box."
        message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if

    type is (box_sphere_t)
      if (radius > box%radius) then
        message(1) = "The radius of an orbital set is bigger than the radius of the simulatio box."
        message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if
      if (norm2(ions%pos(:, iatom) - box%center) + radius > box%radius) then
        message(1) = "An orbital set has points outside of the simulatio box."
        message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if

    type is (box_cylinder_t)
      if(radius > box%radius) then
        message(1) = "The radius of an orbital set is bigger than the radius of the simulatio box."
        message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if
      if (radius > box%half_length) then
        message(1) = "The radius of an orbital set is bigger than the length of the cylinder box."
        message(2) = "Increase the value of XLength or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if

      if (norm2(ions%pos(2:mesh%sb%dim,iatom) - box%center(2:mesh%sb%dim)) + radius > box%radius) then
        message(1) = "An orbital set has points outside of the simulatio box."
        message(2) = "Increase the value of Radius or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if
      if (abs(ions%pos(1, iatom) - box%center(1)) + radius > box%half_length) then
        message(1) = "An orbital set has points outside of the simulatio box."
        message(2) = "Increase the value of Xlength or decrease the value of OrbitalsThreshold_LDAU."
        write(message(3),'(a,f8.5,a,i5,a)') 'The value of the radius is ', radius, ' Bohr.'
        call messages_fatal(3)
      end if

    end select

    !We initialise the submesh corresponding to the orbital 
    call submesh_init(sm, ions%space, mesh, ions%latt, ions%pos(:, iatom), radius)

  end if

  if(.not. allocated(os%X(orb))) then
    if(use_mesh) then
      SAFE_ALLOCATE(os%X(orb)(1:mesh%np,1:os%ndim,1:os%norbs))
    else
      SAFE_ALLOCATE(os%X(orb)(1:sm%np,1:os%ndim,1:os%norbs))
    end if
    os%X(orb)(:,:,:) = R_TOTYPE(M_ZERO)
  end if

  if(d_dim == 1) then

    mm = orbind-1-ll

    !We get the orbital from the pseudopotential
    !In this case we want to get a real orbital and to store it in complex array
    SAFE_ALLOCATE(tmp(1:sm%np))
    call datomic_orbital_get_submesh(spec, sm, ii, ll, mm, 1, tmp)
    if(normalize) then
      norm = dsm_nrm2(os%sphere, tmp)
      call lalg_scal(os%sphere%np, M_ONE/norm, tmp)
    end if

    if(use_mesh) then
      call submesh_add_to_mesh(sm, tmp, os%X(orb)(1:mesh%np, 1, orbind))
    else
      os%X(orb)(1:sm%np, 1, orbind) = tmp(1:sm%np)
    end if
    SAFE_DEALLOCATE_A(tmp)

  else
    SAFE_ALLOCATE(ztmp(1:sm%np, 1:2))  

    if(jj == ll+M_HALF) then
      mm = orbind - 2 - ll
      if(mm >= -ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm, 1, ztmp(:, 1))
        coeff = sqrt((ll+mm+M_ONE)/(M_TWO*ll+M_ONE)) 
        call lalg_scal(sm%np, coeff, ztmp(:, 1))
      else
        ztmp(1:sm%np, 1) = M_ZERO
      end if
      if(mm < ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm+1, 1, ztmp(:,2))
        coeff = sqrt((ll-mm)/(M_TWO*ll+M_ONE))                           
        call lalg_scal(sm%np, coeff, ztmp(:, 2))
      else
       ztmp(1:sm%np, 2) = M_ZERO
      end if
    else
      mm = orbind - ll
      call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm, 1, ztmp(:,2))
      coeff = -sqrt((ll+mm)/(M_TWO*ll+M_ONE))  
      call lalg_scal(sm%np, coeff, ztmp(:, 2)) 
      if(mm > -ll) then
        call X(atomic_orbital_get_submesh)(spec, sm, ii, ll, mm-1, 1, ztmp(:,1))
        coeff = sqrt((ll-mm+M_ONE)/(M_TWO*ll+M_ONE))      
        call lalg_scal(sm%np, coeff, ztmp(:, 1))
      else
       ztmp(1:sm%np, 1) = M_ZERO
      end if
    end if

    if(normalize) then  
      norm = X(sm_nrm2)(os%sphere, ztmp(:,1))**2 
      norm = norm + X(sm_nrm2)(os%sphere, ztmp(:,2))**2
      norm = sqrt(norm)
      call lalg_scal(os%sphere%np, M_ONE/norm, ztmp(:,1))
      call lalg_scal(os%sphere%np, M_ONE/norm, ztmp(:,2))
    end if


    if(use_mesh) then
      call submesh_add_to_mesh(sm, ztmp(:, 1), os%X(orb)(1:mesh%np, 1, orbind))
      call submesh_add_to_mesh(sm, ztmp(:, 2), os%X(orb)(1:mesh%np, 2, orbind))
    else
      os%X(orb)(1:sm%np, 1, orbind) = ztmp(1:sm%np, 1)
      os%X(orb)(1:sm%np, 2, orbind) = ztmp(1:sm%np, 2)
    end if
    SAFE_DEALLOCATE_A(tmp)

  end if

  POP_SUB(X(get_atomic_orbital))

end subroutine X(get_atomic_orbital)



  ! ---------------------------------------------------------
  subroutine X(atomic_orbital_get_submesh)(species, submesh, ii, ll, mm, ispin, phi, derivative)
    type(species_t), target, intent(in)  :: species       !< The species.
    type(submesh_t),         intent(in)  :: submesh    !< The submesh descriptor where the orbital will be calculated.
    integer,                 intent(in)  :: ii
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    integer,                 intent(in)  :: ispin      !< The spin index.
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

      do ip = 1, submesh%np
        phi(ip) = submesh%r(ip)
      end do

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
          ww = species_zval(species)*submesh%r(ip)/ii
          phi(ip) = sqrt( (2*species_zval(species)/ii)**3 * factorial(ii - ll - 1) / (2*ii*factorial(ii+ll)) ) * &
            exp(-ww) * (2 * ww)**ll * loct_sf_laguerre_n(ii-ll-1, TOFLOAT(2*ll + 1), 2*ww)
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
        phi(ip) = exp(-ww*submesh%r(ip)**2/M_TWO)
        do idir = 1, submesh%mesh%sb%dim
          phi(ip) = phi(ip) * hermite(nn(idir) - 1, submesh%x(ip, idir)*sqrtw)
        end do
      end do
      
    end if

    POP_SUB(X(atomic_orbital_get_submesh))
  end subroutine X(atomic_orbital_get_submesh)


  ! ---------------------------------------------------------
  ! Does the same job as atomic_orbital_get_submesh, but with an extra check that
  ! all points can be evaluated first.
  ! In case it cannot, it creates a temporary submesh on which the points can be evaluated,
  ! call atomic_orbital_get_submesh for this one, and copies back the points on the original one
  subroutine X(atomic_orbital_get_submesh_safe)(species, submesh, ii, ll, mm, ispin, phi, derivative)
    type(species_t), target, intent(in)  :: species       !< The species.
    type(submesh_t),         intent(in)  :: submesh    !< The submesh descriptor where the orbital will be calculated.
    integer,                 intent(in)  :: ii
    integer,                 intent(in)  :: ll
    integer,                 intent(in)  :: mm
    integer,                 intent(in)  :: ispin      !< The spin index.
    R_TYPE,                  intent(out) :: phi(:)     !< The function defined in the mesh where the orbitals is returned.
    logical,       optional, intent(in)  :: derivative !< If present and .true. returns the derivative of the orbital.

    integer :: ip, is
    logical :: safe
    integer, allocatable :: map(:)
    type(submesh_t) :: tmp_sm
    R_TYPE, allocatable :: phi_tmp(:)
    FLOAT :: threshold
    type(ps_t), pointer :: ps

    if(submesh%np == 0) return

    PUSH_SUB(X(atomic_orbital_get_submesh_safe))

    safe = .true.
    if(species_is_ps(species)) then
      ps => species_ps(species)
      threshold = spline_range_max(ps%ur(ii, ispin))
      if(any(submesh%r(1:submesh%np) > threshold)) safe = .false.
    end if

    if(safe) then
 
      call X(atomic_orbital_get_submesh)(species, submesh, ii, ll, mm, ispin, phi, derivative)

    else
      ASSERT(species_is_ps(species))
      ps => species_ps(species)
      threshold = spline_range_max(ps%ur(ii, ispin))

      is = 0
      do ip = 1, submesh%np
        if(submesh%r(ip) <= threshold) then
          is = is + 1
        end if
      end do

      SAFE_ALLOCATE(map(1:is))
      tmp_sm%mesh => submesh%mesh
      tmp_sm%np = is
      tmp_sm%np_part = is
      SAFE_ALLOCATE(tmp_sm%x(1:tmp_sm%np_part, 1:submesh%mesh%sb%dim))
      SAFE_ALLOCATE(tmp_sm%r(1:tmp_sm%np_part))
      SAFE_ALLOCATE(phi_tmp(1:tmp_sm%np))
      is = 0
      do ip = 1, submesh%np
        if(submesh%r(ip) <= threshold) then
          is = is + 1
          map(is) = ip
          tmp_sm%x(is, :) = submesh%x(ip, :)
          tmp_sm%r(is) = submesh%r(ip)
        end if
      end do

      call X(atomic_orbital_get_submesh)(species, tmp_sm, ii, ll, mm, ispin, phi_tmp, derivative)

      phi = R_TOTYPE(M_ZERO)
      do ip = 1, tmp_sm%np
        phi(map(ip)) = phi_tmp(ip)
      end do

      call submesh_end(tmp_sm)
      SAFE_DEALLOCATE_A(map)

    end if

    POP_SUB(X(atomic_orbital_get_submesh_safe))
  end subroutine X(atomic_orbital_get_submesh_safe)

