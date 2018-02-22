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
  ! ---------------------------------------------------------
  subroutine X(species_get_orbital_submesh)(species, submesh, ii, ll, mm, ispin, pos, phi, derivative)
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

    PUSH_SUB(X(species_get_orbital_submesh))

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

    POP_SUB(X(species_get_orbital_submesh))
  end subroutine X(species_get_orbital_submesh)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
