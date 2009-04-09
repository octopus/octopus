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
  use loct_math_m
  use simul_box_m
  use ps_m
  use species_m
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

    type(species_t), pointer :: s
    FLOAT :: r, dd, zi, zj
    integer :: iatom, jatom

    FLOAT, parameter :: alpha = CNST(1.1313708)

    type(profile_t), save :: ion_ion_prof

    call profiling_in(ion_ion_prof, "ION_ION_INTERACTION")

    ! see
    ! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:Ion-Ion_interaction
    ! for details about this routine.

    energy = M_ZERO
    force = M_ZERO

    if(simul_box_is_periodic(sb)) then
      call ion_interaction_periodic(geo, sb, energy, force)
    else
      ! only interaction inside the cell
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
          force(1:sb%dim, iatom) = force(1:sb%dim, iatom) + &
            dd*(geo%atom(iatom)%x(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))

          !energy
          if(jatom > iatom) cycle
          energy = energy + zi*zj/r

        end do !jatom
      end do !iatom

    end if
    call profiling_out(ion_ion_prof)

  end subroutine ion_interaction_calculate

  subroutine ion_interaction_periodic(geo, sb, energy, force)
    type(geometry_t),  target, intent(in)    :: geo
    type(simul_box_t),         intent(in)    :: sb
    FLOAT,                     intent(out)   :: energy
    FLOAT,                     intent(out)   :: force(:, :)

    type(species_t), pointer :: s
    FLOAT :: r, xi(1:MAX_DIM), zi, zj
    integer :: iatom, jatom, icopy
    type(periodic_copy_t) :: pc
    integer :: ix, iy, iz, isph, ss
    FLOAT   :: gg(1:MAX_DIM), gg2
    FLOAT   :: factor, charge
    CMPLX   :: sumatoms
    FLOAT, parameter :: alpha = CNST(1.1313708)

    ! see
    ! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:Ion-Ion_interaction
    ! for details about this routine.

    energy = M_ZERO
    force = M_ZERO

    ! if the system is periodic we have to add the energy of the
    ! interaction with the copies
    
    ! the short-range part is calculated directly
    do iatom = 1, geo%natoms
      s => geo%atom(iatom)%spec
      if (.not. species_is_ps(s)) cycle
      zi = geo%atom(iatom)%spec%z_val

      call periodic_copy_init(pc, sb, geo%atom(iatom)%x, CNST(5.0))
      
      do icopy = 1, periodic_copy_num(pc)
        xi = periodic_copy_position(pc, sb, icopy)
        
        do jatom = 1, geo%natoms
          zj = -geo%atom(jatom)%spec%z_val
          r = sqrt( sum( (xi(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))**2 ) )
          
          if(r < CNST(1e-5)) cycle
          
          ! energy
          energy = energy + M_HALF*zj*zi*(M_ONE - loct_erf(alpha*r))
          
          ! force
          force(1:sb%dim, jatom) = force(1:sb%dim, jatom) + &
            M_HALF*zj*zi*(M_ONE - loct_erf(alpha*r))/r*(geo%atom(jatom)%x(1:sb%dim) - xi(1:sb%dim))
        end do
        
      end do
      
      call periodic_copy_end(pc)
    end do

    ! And the long-range part, using an Ewald sum
    
    isph = 100
    do ix = -isph, isph
      do iy = -isph, isph
        do iz = -isph, isph
          
          ss = ix**2 + iy**2 + iz**2
          
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:sb%dim) = ix*sb%klattice(1:sb%dim, 1) + iy*sb%klattice(1:sb%dim, 2) + iz*sb%klattice(1:sb%dim, 3)
          gg2 = sum(gg(1:sb%dim)**2)
          
          ! k=0 must be removed from the sum
          if(gg2 == M_ZERO) cycle

          factor = M_TWO*M_PI/sb%rcell_volume*exp(-CNST(0.25)*gg2/alpha**2)/gg2
          
          sumatoms = M_Z0
          do iatom = 1, geo%natoms
            zi = geo%atom(iatom)%spec%z_val
            xi(1:sb%dim) = geo%atom(iatom)%x(1:sb%dim)
            sumatoms = sumatoms + zi*exp(-M_ZI*sum(gg(1:sb%dim)*xi(1:sb%dim)))
          end do
          energy = energy + factor*sumatoms*conjg(sumatoms)
          
          do iatom = 1, geo%natoms
            zi = geo%atom(iatom)%spec%z_val
            force(1:sb%dim, iatom) = -M_TWO*zi*factor*sumatoms
          end do
          
        end do
      end do
    end do
    
    ! remove self-interaction
    charge = M_ZERO
    do iatom = 1, geo%natoms
      zi = geo%atom(iatom)%spec%z_val
      charge = charge + zi
      energy = energy - alpha*zi**2/sqrt(M_PI) 
    end do
    
    ! This term is added in abinit, I am not sure where it comes
    ! from and whether we should add it.
    !
    ! energy = energy - M_PI*charge**2/(M_TWO*alpha**2*sb%rcell_volume)
    
  end subroutine ion_interaction_periodic
  
end module ion_interaction_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
