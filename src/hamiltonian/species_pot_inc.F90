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
!! $Id$

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Places, in the function phi (defined in each point of the mesh), the
  !! atomic orbital. The orbitals are obtained from the species data
  !! type, and are numbered from one to species_niwfs(spec). It may happen
  !! that there are different orbitals for each spin-polarization direction,
  !! and therefore the orbital is also characterized by the label "is".
  !!
  !! In order to put the orbital in the mesh, it is necessary to know where
  !! the species is, and this is given by the vector "pos".
  !!
  !! \todo Most of this work should be done inside the species
  !! module, and we should get rid of species_iwf_i, species_ifw_l, etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine X(species_get_orbital)(spec, mesh, ii, ll, mm, ispin, pos, orb, scale)
    type(species_t), target, intent(in)     :: spec
    type(mesh_t),            intent(in)     :: mesh
    integer,                 intent(in)     :: ii
    integer,                 intent(in)     :: ll
    integer,                 intent(in)     :: mm
    integer,                 intent(in)     :: ispin   !< The spin index.
    FLOAT,                   intent(in)     :: pos(:)  !< The position of the atom.
    R_TYPE,                  intent(out)    :: orb(:)  !< The function defined in the mesh where the orbitals is returned.
    FLOAT, optional,         intent(in)     :: scale

    integer :: ip, icell, nn(3), idir
    FLOAT :: r2, x(1:MAX_DIM), radius, xfactor
    FLOAT, allocatable :: xf(:, :), lorb(:)
    R_TYPE, allocatable :: ylm(:)
    type(ps_t), pointer :: ps
    type(periodic_copy_t) :: pc

    PUSH_SUB(X(species_get_orbital))

    xfactor = CNST(1.0)/optional_default(scale, CNST(1.0))

    radius = min(species_get_iwf_radius(spec, ii, ispin)/xfactor, maxval(mesh%sb%lsize))

    call periodic_copy_init(pc, mesh%sb, pos, range = radius)

    orb = M_ZERO

    SAFE_ALLOCATE(lorb(1:mesh%np))

    if(species_represents_real_atom(spec)) then
      if(species_is_ps(spec)) then
        ps => species_ps(spec)
        ASSERT(ii <= ps%conf%p)
      endif
      SAFE_ALLOCATE(xf(1:mesh%np, 1:mesh%sb%dim))
      SAFE_ALLOCATE(ylm(1:mesh%np))
    endif

    do icell = 1, periodic_copy_num(pc)

      if(species_represents_real_atom(spec) .and. mesh%sb%dim == 3) then

        do ip = 1, mesh%np
          x(1:mesh%sb%dim) = (mesh%x(ip, 1:mesh%sb%dim) - periodic_copy_position(pc, mesh%sb, icell))*xfactor
          r2 = sum(x(1:mesh%sb%dim)**2)
          xf(ip, 1:mesh%sb%dim) = x(1:mesh%sb%dim)
          
          if(species_is_ps(spec)) then
            if(r2 < spline_range_max(ps%ur_sq(ii, ispin))) then
              lorb(ip) = spline_eval(ps%ur_sq(ii, ispin), r2)
            else
              lorb(ip) = M_ZERO
            end if
          else
            ! FIXME: cache result somewhat. e.g. re-use result for each m. and use recursion relation.
            lorb(ip) = sqrt( (2*species_zval(spec)/ii)**3 * factorial(ii - ll - 1) / (2*ii*factorial(ii+ll)) ) * &
              exp(-species_zval(spec) * sqrt(r2) / ii) * (2*species_zval(spec)*sqrt(r2)/ii)**ll * &
              loct_sf_laguerre_n(ii-ll-1, real(2*ll + 1, REAL_PRECISION), 2*species_zval(spec)*sqrt(r2)/ii)
          endif
        end do

#ifdef R_TCOMPLEX
        ! complex spherical harmonics. FIXME: vectorize
        do ip = 1, mesh%np
          call ylmr(xf(ip, 1), xf(ip, 2), xf(ip, 3), ll, mm, ylm(ip))
        enddo
#else
        ! real spherical harmonics
        call loct_ylm(mesh%np, xf(1, 1), xf(1, 2), xf(1, 3), ll, mm, ylm(1))
#endif

        do ip = 1, mesh%np
          orb(ip) = orb(ip) + lorb(ip)*ylm(ip)
        end do

      else

        ! FIXME: this is a pretty dubious way to handle l and m quantum numbers. Why not use ylm?
        ! also, these are far from normalized.
        nn = (/ii, ll, mm/)

        do ip = 1, mesh%np
          x(1:mesh%sb%dim) = (mesh%x(ip, 1:mesh%sb%dim) - periodic_copy_position(pc, mesh%sb, icell))*xfactor
          r2 = sum(x(1:mesh%sb%dim)**2)
          lorb = exp(-species_omega(spec)*r2/M_TWO)
          do idir = 1, mesh%sb%dim
            lorb(ip) = lorb(ip) * hermite(nn(idir) - 1, x(idir)*sqrt(species_omega(spec)))
          enddo
          orb(ip) = orb(ip) + lorb(ip)
        end do
      end if

    end do

    SAFE_DEALLOCATE_A(lorb)

    if(species_is_ps(spec)) then
      SAFE_DEALLOCATE_A(xf)
      SAFE_DEALLOCATE_A(ylm)
      nullify(ps)
    endif

    call periodic_copy_end(pc)

    POP_SUB(X(species_get_orbital))
  end subroutine X(species_get_orbital)


  ! ---------------------------------------------------------
  subroutine X(species_get_orbital_submesh)(spec, submesh, ii, ll, mm, ispin, pos, phi, derivative)
    type(species_t), target, intent(in)  :: spec       !< The species.
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

    PUSH_SUB(X(species_get_orbital_submesh))

    derivative_ = optional_default(derivative, .false.)

    ASSERT(ubound(phi, dim = 1) >= submesh%np)

    if(species_represents_real_atom(spec)) then
      ps => species_ps(spec)
      
      forall(ip = 1:submesh%np) phi(ip) = submesh%x(ip, 0)

      if(species_is_ps(spec)) then
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
          ww = species_zval(spec)*submesh%x(ip, 0)/ii
          phi(ip) = sqrt( (2*species_zval(spec)/ii)**3 * factorial(ii - ll - 1) / (2*ii*factorial(ii+ll)) ) * &
            exp(-ww) * (2 * ww)**ll * loct_sf_laguerre_n(ii-ll-1, real(2*ll + 1, REAL_PRECISION), 2*ww)
        enddo
      endif

      SAFE_ALLOCATE(ylm(1:submesh%np))

#ifdef R_TCOMPLEX
      ! complex spherical harmonics. FIXME: vectorize
      do ip = 1, submesh%np
        call ylmr(submesh%x(ip, 1), submesh%x(ip, 2), submesh%x(ip, 3), ll, mm, ylm(ip))
      enddo
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

      ww = species_omega(spec)
      sqrtw = sqrt(ww)

      ! FIXME: this is a pretty dubious way to handle l and m quantum numbers. Why not use ylm?
      nn = (/ii, ll, mm/)

      do ip = 1, submesh%np
        phi(ip) = exp(-ww*submesh%x(ip, 0)**2/M_TWO)
        do idir = 1, submesh%mesh%sb%dim
          phi(ip) = phi(ip) * hermite(nn(idir) - 1, submesh%x(ip, idir)*sqrtw)
        enddo
      end do
      
    end if

    POP_SUB(X(species_get_orbital_submesh))
  end subroutine X(species_get_orbital_submesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
