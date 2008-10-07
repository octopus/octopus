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

module energy_m
  use calc_mode_m
  use datasets_m
  use external_pot_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use smear_m
  use states_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                      &
    total_energy,                &
    dcalculate_eigenvalues,      &
    zcalculate_eigenvalues,      &
    delectronic_kinetic_energy,  &
    zelectronic_kinetic_energy,  &
    delectronic_external_energy, &
    zelectronic_external_energy

contains

  ! ---------------------------------------------------------
  ! This subroutine calculates the total energy of the system. Basically, it
  ! adds up the KS eigenvalues, and then it subtracts the whatever double
  ! counts exist (see TDDFT theory for details).
  subroutine total_energy(h, gr, st, iunit, full)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: iunit
    logical, optional,   intent(in)    :: full

    logical :: full_

    call push_sub('energy.energy_calculate')

    full_ = .false.
    if(present(full)) full_ = full

    select case(h%theory_level)
    case(INDEPENDENT_PARTICLES)
      h%eeigen = states_eigenvalues_sum(st)
      h%etot   = h%ep%eii + h%eeigen

    case(HARTREE)
      if(st%wfs_type == M_REAL) then
        h%t0     = delectronic_kinetic_energy(h, gr, st)
        h%eext   = delectronic_external_energy(h, gr, st)
      else
        h%t0     = zelectronic_kinetic_energy(h, gr, st)
        h%eext   = zelectronic_external_energy(h, gr, st)
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot = h%ep%eii + M_HALF*(h%eeigen + h%t0 + h%eext)

    case(HARTREE_FOCK)
      if(st%wfs_type == M_REAL) then
        h%t0     = delectronic_kinetic_energy(h, gr, st)
        h%eext   = delectronic_external_energy(h, gr, st)
      else
        h%t0     = zelectronic_kinetic_energy(h, gr, st)
        h%eext   = zelectronic_external_energy(h, gr, st)
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot = h%ep%eii + M_HALF*(h%eeigen + h%t0 + h%eext - h%epot) + h%ec

    case(KOHN_SHAM_DFT)
      if(full_) then
        if(st%wfs_type == M_REAL) then
          h%t0     = delectronic_kinetic_energy(h, gr, st)
          h%eext   = delectronic_external_energy(h, gr, st)
        else
          h%t0     = zelectronic_kinetic_energy(h, gr, st)
          h%eext   = zelectronic_external_energy(h, gr, st)
        end if
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot   = h%ep%eii + h%eeigen - h%ehartree + h%ex + h%ec - h%epot

    end select
    
    h%entropy = smear_calc_entropy(st%smear, st%eigenval, st%d%nik, st%nst, st%d%kweights)

    if(gauge_field_is_applied(h%ep%gfield)) then
      h%etot = h%etot + gauge_field_get_energy(h%ep%gfield, gr%sb)
    end if

    if (iunit > 0) then
      write(message(1), '(6x,a, f18.8)')'Total       = ', h%etot     / units_out%energy%factor
      write(message(2), '(6x,a, f18.8)')'Free        = ', (h%etot+h%entropy) / units_out%energy%factor
      write(message(3), '(6x,a)') '-----------'
      call write_info(3, iunit)

      write(message(1), '(6x,a, f18.8)')'Ion-ion     = ', h%ep%eii   / units_out%energy%factor
      write(message(2), '(6x,a, f18.8)')'Eigenvalues = ', h%eeigen   / units_out%energy%factor
      write(message(3), '(6x,a, f18.8)')'Hartree     = ', h%ehartree / units_out%energy%factor
      write(message(4), '(6x,a, f18.8)')'Int[n*v_xc] = ', h%epot     / units_out%energy%factor
      write(message(5), '(6x,a, f18.8)')'Exchange    = ', h%ex       / units_out%energy%factor
      write(message(6), '(6x,a, f18.8)')'Correlation = ', h%ec       / units_out%energy%factor
      write(message(7), '(6x,a, f18.8)')'-TS         = ', h%entropy  / units_out%energy%factor
      call write_info(7, iunit)
      if(full_) then
        write(message(1), '(6x,a, f18.8)')'Kinetic     = ', h%t0 / units_out%energy%factor
        write(message(2), '(6x,a, f18.8)')'External    = ', h%eext / units_out%energy%factor
        call write_info(2, iunit)
      end if
    end if

    call pop_sub()
  end subroutine total_energy

#include "undef.F90"
#include "real.F90"
#include "energy_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "energy_inc.F90"

end module energy_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
