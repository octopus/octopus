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
  use batch_m
  use berry_m
  use datasets_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use io_m
  use lalg_basic_m
  use mesh_m
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                      &
    total_energy,                &
    delectronic_kinetic_energy,  &
    zelectronic_kinetic_energy,  &
    delectronic_external_energy, &
    zelectronic_external_energy, &
    energy_calculate_eigenvalues

contains

  ! ---------------------------------------------------------
  ! This subroutine calculates the total energy of the system. Basically, it
  ! adds up the KS eigenvalues, and then it subtracts whatever double
  ! counts exist (see TDDFT theory for details).
  subroutine total_energy(hm, gr, st, iunit, full)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: iunit
    logical, optional,   intent(in)    :: full

    logical :: full_

    PUSH_SUB(total_energy)

    full_ = .false.
    if(present(full)) full_ = full

    hm%energy%eigenvalues = states_eigenvalues_sum(st)

    select case(hm%theory_level)
    case(INDEPENDENT_PARTICLES)
      hm%energy%total   = hm%ep%eii + hm%energy%eigenvalues

    case(HARTREE)
      if(states_are_real(st)) then
        hm%energy%kinetic     = delectronic_kinetic_energy(hm, gr, st)
        hm%energy%extern   = delectronic_external_energy(hm, gr, st)
      else
        hm%energy%kinetic     = zelectronic_kinetic_energy(hm, gr, st)
        hm%energy%extern   = zelectronic_external_energy(hm, gr, st)
      end if
      hm%energy%total = hm%ep%eii + M_HALF * (hm%energy%eigenvalues + hm%energy%kinetic + hm%energy%extern)

    case(HARTREE_FOCK)
      if(states_are_real(st)) then
        hm%energy%kinetic     = delectronic_kinetic_energy(hm, gr, st)
        hm%energy%extern   = delectronic_external_energy(hm, gr, st)
      else
        hm%energy%kinetic     = zelectronic_kinetic_energy(hm, gr, st)
        hm%energy%extern   = zelectronic_external_energy(hm, gr, st)
      end if
      hm%energy%total = hm%ep%eii + &
        M_HALF*(hm%energy%eigenvalues + hm%energy%kinetic + hm%energy%extern - hm%energy%intnvxc) + hm%energy%correlation

    case(KOHN_SHAM_DFT)
      if(full_) then
        if(states_are_real(st)) then
          hm%energy%kinetic     = delectronic_kinetic_energy(hm, gr, st)
          hm%energy%extern   = delectronic_external_energy(hm, gr, st)
        else
          hm%energy%kinetic     = zelectronic_kinetic_energy(hm, gr, st)
          hm%energy%extern   = zelectronic_external_energy(hm, gr, st)
        end if
      end if
      hm%energy%total = hm%ep%eii + hm%energy%eigenvalues &
        - hm%energy%hartree + hm%energy%exchange + hm%energy%correlation - hm%energy%intnvxc

    case(CLASSICAL)
      st%eigenval = M_ZERO
      hm%energy%eigenvalues = M_ZERO
      hm%energy%total = hm%ep%eii
    end select
    
    hm%energy%entropy = smear_calc_entropy(st%smear, st%eigenval, st%d%nik, st%nst, st%d%kweights, st%occ)
    if(st%smear%method == SMEAR_FIXED_OCC) then ! no temperature available
      hm%energy%TS = M_ZERO
    else
      hm%energy%TS = st%smear%dsmear * hm%energy%entropy
    endif

    if(gauge_field_is_applied(hm%ep%gfield)) then
      hm%energy%total = hm%energy%total + gauge_field_get_energy(hm%ep%gfield, gr%sb)
    end if

    if(associated(hm%vberry)) then
      hm%energy%berry = berry_energy_correction(st, gr%mesh, &
        hm%ep%E_field(1:gr%sb%periodic_dim), hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
      hm%energy%total = hm%energy%total + hm%energy%berry
    else
      hm%energy%berry = M_ZERO
    endif

    if (iunit > 0) then
      write(message(1), '(6x,a, f18.8)')'Total       = ', units_from_atomic(units_out%energy, hm%energy%total)
      write(message(2), '(6x,a, f18.8)')'Free        = ', units_from_atomic(units_out%energy, hm%energy%total - hm%energy%TS)
      write(message(3), '(6x,a)') '-----------'
      call messages_info(3, iunit)

      write(message(1), '(6x,a, f18.8)')'Ion-ion     = ', units_from_atomic(units_out%energy, hm%ep%eii)
      write(message(2), '(6x,a, f18.8)')'Eigenvalues = ', units_from_atomic(units_out%energy, hm%energy%eigenvalues)
      write(message(3), '(6x,a, f18.8)')'Hartree     = ', units_from_atomic(units_out%energy, hm%energy%hartree)
      write(message(4), '(6x,a, f18.8)')'Int[n*v_xc] = ', units_from_atomic(units_out%energy, hm%energy%intnvxc)
      write(message(5), '(6x,a, f18.8)')'Exchange    = ', units_from_atomic(units_out%energy, hm%energy%exchange)
      write(message(6), '(6x,a, f18.8)')'Correlation = ', units_from_atomic(units_out%energy, hm%energy%correlation)
      write(message(7), '(6x,a, f18.8)')'Entropy     = ', hm%energy%entropy ! the dimensionless sigma of Kittel&Kroemer
      write(message(8), '(6x,a, f18.8)')'-TS         = ', -units_from_atomic(units_out%energy, hm%energy%TS)
      call messages_info(8, iunit)
      if(full_) then
        write(message(1), '(6x,a, f18.8)')'Kinetic     = ', units_from_atomic(units_out%energy, hm%energy%kinetic)
        write(message(2), '(6x,a, f18.8)')'External    = ', units_from_atomic(units_out%energy, hm%energy%extern)
        call messages_info(2, iunit)
      end if
      if(associated(hm%ep%E_field) .and. simul_box_is_periodic(gr%sb)) then
        write(message(1), '(6x,a, f18.8)')'Berry       = ', units_from_atomic(units_out%energy, hm%energy%berry)
        call messages_info(1, iunit)
      endif  
    end if

    POP_SUB(total_energy)
  end subroutine total_energy

  ! --------------------------------------------------------------------
  
  subroutine energy_calculate_eigenvalues(hm, der, st, time)
    type(hamiltonian_t), intent(inout) :: hm
    type(derivatives_t), intent(inout) :: der
    type(states_t),      intent(inout) :: st
    FLOAT,   optional,   intent(in)    :: time
    
    PUSH_SUB(energy_calculate_eigenvalues)

    if(states_are_real(st)) then
      call dcalculate_eigenvalues(hm, der, st, time)
    else
      call zcalculate_eigenvalues(hm, der, st, time)
    end if

    POP_SUB(energy_calculate_eigenvalues)
  end subroutine energy_calculate_eigenvalues

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
