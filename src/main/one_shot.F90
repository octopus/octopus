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
  use external_pot_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use loct_parser_m
  use loct_m
  use loct_math_m
  use mesh_m
  use messages_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m
  use v_ks_m
  use varinfo_m
  use specie_pot_m
  use xc_m

  implicit none

  private
  public :: one_shot_run

contains

  ! ---------------------------------------------------------
  subroutine one_shot_run(sys, h)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: h

    integer :: ierr
    real(8) :: E_tot, E_t, E_ext, E_Hartree, E_x, E_c

    ! load wave-functions
    call states_allocate_wfns(sys%st, sys%gr%m)
    call restart_read(trim(tmpdir)//'gs', sys%st, sys%gr, sys%geo, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions"
      call write_fatal(1)
    end if

    ! generate density
    call states_calc_dens(sys%st, sys%gr%m%np, sys%st%rho)

    ! kinetic energy + local potential + Hartree + xc
    if(sys%st%wfs_type == M_REAL) then
      E_t     = delectronic_kinetic_energy(h, sys%gr, sys%st)
      E_ext   = delectronic_external_energy(h, sys%gr, sys%st)
    else
      E_t     = zelectronic_kinetic_energy(h, sys%gr, sys%st)
      E_ext   = zelectronic_external_energy(h, sys%gr, sys%st)
    end if

    ! Get the Hartree energy
    call v_ks_hartree(sys%gr, sys%st, h)
    E_Hartree = h%ehartree

    ! Get exchange-correlation energies
    call xc_get_vxc(sys%gr, sys%ks%xc, sys%st%rho, sys%st%d%ispin, E_x, E_c, &
      M_ZERO, sys%st%qtot)

    E_tot = E_t + E_ext + E_Hartree + E_x + E_c + h%ep%eii

    call messages_print_stress(stdout, "Energy")
    write(message(1), '(6x,a, f15.8)')'Total       = ', E_tot     / units_out%energy%factor
    write(message(2), '(6x,a, f15.8)')'Ion-ion     = ', h%ep%eii  / units_out%energy%factor
    write(message(3), '(6x,a, f15.8)')'Kinetic     = ', E_t       / units_out%energy%factor
    write(message(4), '(6x,a, f15.8)')'External    = ', E_ext     / units_out%energy%factor
    write(message(5), '(6x,a, f15.8)')'Hartree     = ', E_Hartree / units_out%energy%factor
    write(message(6), '(6x,a, f15.8)')'Exchange    = ', E_x       / units_out%energy%factor
    write(message(7), '(6x,a, f15.8)')'Correlation = ', E_c       / units_out%energy%factor
    call write_info(7, stdout)
    call messages_print_stress(stdout)

    call states_deallocate_wfns(sys%st)

  end subroutine one_shot_run

end module one_shot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
