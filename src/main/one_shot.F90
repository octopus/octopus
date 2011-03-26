!! Copyright (C) 2002-2007 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: geom_opt.F90 3970 2008-03-29 11:38:27Z acastro $

#include "global.h"

module one_shot_m
  use datasets_m
  use density_m
  use energy_m
  use epot_m
  use geometry_m
  use global_m
  use output_m
  use hamiltonian_m
  use parser_m
  use loct_m
  use loct_math_m
  use mesh_m
  use messages_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use varinfo_m
  use xc_m
  use xc_OEP_m

  implicit none

  private
  public :: one_shot_run

contains

  ! ---------------------------------------------------------
  subroutine one_shot_run(sys, hm)
    type(system_t),              intent(inout) :: sys
    type(hamiltonian_t),         intent(inout) :: hm

    integer :: ierr
    FLOAT :: e_tot, e_t, e_ext

    PUSH_SUB(one_shot_run)

    ! load wavefunctions
    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call restart_read(trim(tmpdir)//GS_DIR, sys%st, sys%gr, sys%geo, ierr, exact = .true.)

    ! generate density
    call density_calc(sys%st, sys%gr, sys%st%rho)

    ! calculate the KS potential
    call v_ks_calc(sys%ks, hm, sys%st, calc_eigenval = .true.)

    e_t = M_ZERO
    e_ext = M_ZERO

    ! Probably this is particular to DFT
    select case(sys%ks%theory_level)
    case(KOHN_SHAM_DFT)
      ! kinetic energy + local potential + Hartree + xc
      if(states_are_real(sys%st)) then
        e_t     = delectronic_kinetic_energy(hm, sys%gr, sys%st)
        e_ext   = delectronic_external_energy(hm, sys%gr, sys%st)
      else
        e_t     = zelectronic_kinetic_energy(hm, sys%gr, sys%st)
        e_ext   = zelectronic_external_energy(hm, sys%gr, sys%st)
      end if

    case(INDEPENDENT_PARTICLES)
      ! there is nothing to do
    case default
      message(1) = "Fatal: the one shot calculation mode has not been implemented for"
      message(2) = "       this theory level."
      call messages_fatal(2)
    end select

    e_tot = e_t + e_ext + hm%energy%hartree + hm%energy%exchange + hm%energy%correlation + hm%ep%eii

    call messages_print_stress(stdout, "Energy")
    write(message(1), '(6x,a, f18.8)') 'Total       = ', units_to_atomic(units_out%energy, e_tot)
    write(message(2), '(6x,a, f18.8)') 'Ion-ion     = ', units_to_atomic(units_out%energy, hm%ep%eii)
    write(message(3), '(6x,a, f18.8)') 'Kinetic     = ', units_to_atomic(units_out%energy, e_t)
    write(message(4), '(6x,a, f18.8)') 'External    = ', units_to_atomic(units_out%energy, e_ext)
    write(message(5), '(6x,a, f18.8)') 'Hartree     = ', units_to_atomic(units_out%energy, hm%energy%hartree)
    write(message(6), '(6x,a, f18.8)') 'Exchange    = ', units_to_atomic(units_out%energy, hm%energy%exchange)
    write(message(7), '(6x,a, f18.8)') 'Correlation = ', units_to_atomic(units_out%energy, hm%energy%correlation)
    call messages_info(7, stdout)
    call messages_print_stress(stdout)

    call output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, STATIC_DIR)

    call states_deallocate_wfns(sys%st)

    POP_SUB(one_shot_run)
  end subroutine one_shot_run

end module one_shot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
