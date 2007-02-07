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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module phonons_m
  use datasets_m
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lib_adv_alg_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use output_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m
  use v_ks_m

  implicit none

  private
  public :: &
       phonons_t, &
       phonons_init, &
       phonons_end,  &
       phonons_diagonalize_dm

  type phonons_t
    integer :: dim
    FLOAT, pointer :: dm(:,:), freq(:)

    FLOAT :: disp
  end type phonons_t

contains

  ! ---------------------------------------------------------
  subroutine phonons_init(ph, sys)
    type(phonons_t),     intent(out) :: ph
    type(system_t),      intent(inout) :: sys

    ph%dim = sys%geo%natoms*sys%gr%sb%dim
    ALLOCATE(ph%dm(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%freq(ph%dim), ph%dim)

  end subroutine phonons_init


  ! ---------------------------------------------------------
  subroutine phonons_end(ph)
    type(phonons_t),     intent(inout) :: ph

    deallocate(ph%dm)
    deallocate(ph%freq)

  end subroutine phonons_end

  ! ---------------------------------------------------------
  subroutine phonons_diagonalize_dm(ph)
    type(phonons_t),      intent(inout) :: ph
    
    FLOAT, allocatable :: tmpdm(:,:)

    !we need a temporary copy of DM, to avoid passing the same array twice
    ALLOCATE(tmpdm(ph%dim, ph%dim), ph%dim*ph%dim)
    
    tmpdm(1:ph%dim,1:ph%dim)=ph%dm(1:ph%dim,1:ph%dim)

    ! diagonalize DM
    call lalg_eigensolve(ph%dim, tmpdm, ph%dm, ph%freq)

    deallocate(tmpdm)

  end subroutine phonons_diagonalize_dm

end module phonons_m
