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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module ion_interaction_m
  use datasets_m
  use global_m
  use simul_box_m
  use ps_m
  use specie_m
  use solids_m
  use geometry_m
  use splines_m
  use profiling_m

  implicit none

  private
  public ::                    &
    ion_interaction_calculate

contains
  ! ---------------------------------------------------------
  !
  ! This function calculates the energy and the forces between ions.
  !
  ! ---------------------------------------------------------
  subroutine ion_interaction_calculate(geo, sb, energy, force)
    type(geometry_t),  target, intent(in)    :: geo
    type(simul_box_t),         intent(in)    :: sb
    FLOAT,                     intent(out)   :: energy
    FLOAT,                     intent(out)   :: force(:, :)

    type(specie_t), pointer :: s
    FLOAT :: r, rc, xi(1:MAX_DIM), dd, zi, zj
    integer :: iatom, jatom, icopy
    type(periodic_copy_t) :: pc

    type(profile_t), save :: ion_ion_prof

    call profiling_in(ion_ion_prof, "ION_ION_INTERACTION")

    ! see
    ! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:Ion-Ion_interaction
    ! for details about this routine.

    energy = M_ZERO
    force = M_ZERO

    ! interaction inside the cell, calculated directly
    do iatom = 1, geo%natoms
      s => geo%atom(iatom)%spec
      zi = geo%atom(iatom)%spec%z_val

      if(s%type .eq. SPEC_JELLI) then
        energy = energy + (M_THREE/M_FIVE)*s%z_val**2/s%jradius
      end if

      do jatom = 1, geo%natoms

        if(iatom == jatom) cycle

        zj = geo%atom(jatom)%spec%z_val
        r = sqrt(sum((geo%atom(iatom)%x - geo%atom(jatom)%x)**2))

        !the force
        dd = zi*zj/r**3
        force(1:MAX_DIM, iatom) = force(1:MAX_DIM, iatom) + dd*(geo%atom(iatom)%x(1:MAX_DIM) - geo%atom(jatom)%x(1:MAX_DIM))

        !energy
        if(jatom > iatom) cycle
        energy = energy + zi*zj/r
        
      end do !jatom
      
    end do !iatom

    ! if the system is periodic we have to add the energy of the
    ! interaction with the copies
    if(simul_box_is_periodic(sb)) then

      ! the short range part is calculated directly
      do iatom = 1, geo%natoms
        s => geo%atom(iatom)%spec
        if (.not. specie_is_ps(s)) cycle

        rc = spline_cutoff_radius(s%ps%vion, s%ps%projectors_sphere_threshold)
        rc = max(rc, spline_cutoff_radius(s%ps%dvion, s%ps%projectors_sphere_threshold))

        call periodic_copy_init(pc, sb, geo%atom(iatom)%x, rc)

        do icopy = 1, periodic_copy_num(pc)

          xi = periodic_copy_position(pc, sb, icopy)

          ! do not consider atoms inside the cell
          if ( maxval(abs(xi(1:MAX_DIM) - geo%atom(iatom)%x(1:MAX_DIM))) < M_TEN*M_EPSILON ) cycle

          do jatom = 1, geo%natoms
            zj = -geo%atom(jatom)%spec%z_val
            r = sqrt( sum( (xi - geo%atom(jatom)%x)**2 ) )
            ! energy
            energy = energy + M_HALF*zj*spline_eval(s%ps%vion, r)
            ! force
            force(1:MAX_DIM, jatom) = force(1:MAX_DIM, jatom) + &
              M_HALF*zj*spline_eval(s%ps%dvion, r)/r*(geo%atom(jatom)%x(1:MAX_DIM) - xi(1:MAX_DIM))
          end do
          
        end do
        
        call periodic_copy_end(pc)
        
      end do

      ! And the long range part, using the potential that was already
      ! calculated via poisson equation

    end if

    call profiling_out(ion_ion_prof)
  end subroutine ion_interaction_calculate

end module ion_interaction_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
