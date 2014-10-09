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
  !! iorb-th atomic orbital. The orbitals are obtained from the species data
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
  subroutine X(species_get_orbital)(spec, mesh, iorb, ispin, pos, orb, scale)
    type(species_t), target, intent(in)     :: spec
    type(mesh_t),            intent(in)     :: mesh
    integer,                 intent(in)     :: iorb
    integer,                 intent(in)     :: ispin   !< The spin index.
    FLOAT,                   intent(in)     :: pos(:)  !< The position of the atom.
    R_TYPE,                  intent(out)    :: orb(:)  !< The function defined in the mesh where the orbitals is returned.
    FLOAT, optional,         intent(in)     :: scale

    integer :: i, l, m, ip, icell, nn(3), idir
    FLOAT :: r2, x(1:MAX_DIM), radius, xfactor
    FLOAT, allocatable :: xf(:, :), lorb(:)
    R_TYPE, allocatable :: ylm(:)
    type(ps_t), pointer :: ps
    type(periodic_copy_t) :: pc

    PUSH_SUB(X(species_get_orbital))

    xfactor = CNST(1.0)/optional_default(scale, CNST(1.0))

    call species_iwf_ilm(spec, iorb, ispin, i, l, m)

    radius = min(species_get_iwf_radius(spec, iorb, ispin)/xfactor, maxval(mesh%sb%lsize))

    call periodic_copy_init(pc, mesh%sb, pos, range = radius)

    orb = M_ZERO

    SAFE_ALLOCATE(lorb(1:mesh%np))

    if(species_is_ps(spec)) then
      ps => species_ps(spec)
      SAFE_ALLOCATE(xf(1:mesh%np, 1:mesh%sb%dim))
      SAFE_ALLOCATE(ylm(1:mesh%np))
    endif

    do icell = 1, periodic_copy_num(pc)

      if(species_is_ps(spec)) then

        do ip = 1, mesh%np
          x(1:mesh%sb%dim) = (mesh%x(ip, 1:mesh%sb%dim) - periodic_copy_position(pc, mesh%sb, icell))*xfactor
          r2 = sum(x(1:mesh%sb%dim)**2)
          xf(ip, 1:mesh%sb%dim) = x(1:mesh%sb%dim)
          
          if(r2 < spline_range_max(ps%ur_sq(i, ispin))) then
            lorb(ip) = spline_eval(ps%ur_sq(i, ispin), r2)
          else
            lorb(ip) = M_ZERO
          end if

        end do

#ifdef R_TCOMPLEX
        ! complex spherical harmonics. FIXME: vectorize
        do ip = 1, mesh%np
          call ylmr(xf(ip, 1), xf(ip, 2), xf(ip, 3), l, m, ylm(ip))
        enddo
#else
        ! real spherical harmonics
        call loct_ylm(mesh%np, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))
#endif

        do ip = 1, mesh%np
          orb(ip) = orb(ip) + lorb(ip)*ylm(ip)
        end do

      else

        ! FIXME: this is a pretty dubious way to handle l and m quantum numbers. Why not use ylm?
        nn = (/i, l, m/)

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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
