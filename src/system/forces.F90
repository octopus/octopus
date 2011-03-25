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

module forces_m
  use born_charges_m
  use comm_m
  use datasets_m
  use derivatives_m
  use epot_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use index_m
  use io_m
  use kpoints_m
  use lalg_adv_m
  use lalg_basic_m
  use lasers_m
  use linear_response_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multigrid_m
  use parser_m
  use poisson_m
  use poisson_cutoff_m
  use poisson_sete_m
  use profiling_m
  use projector_m
  use ps_m
  use simul_box_m
  use smear_m
  use solids_m
  use species_m
  use species_pot_m
  use splines_m
  use spline_filter_m
  use states_m
  use states_dim_m
  use submesh_m
  use varinfo_m

  implicit none

  private
  public ::                    &
    forces_calculate,          &
    dforces_from_potential,    &
    zforces_from_potential,    &
    dforces_born_charges,      &
    zforces_born_charges

  type(profile_t), save :: prof_comm

contains

  ! ---------------------------------------------------------
  subroutine forces_calculate(gr, geo, ep, st, t)
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(inout) :: geo
    type(epot_t),     intent(inout) :: ep
    type(states_t),   intent(inout) :: st
    FLOAT,     optional, intent(in) :: t

    integer :: i, j
    FLOAT :: x(MAX_DIM), time
    
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(forces_calculate)

    x(:) = M_ZERO
    time = M_ZERO
    if(present(t)) time = t

    ! the ion-ion term is already calculated
    do i = 1, geo%natoms
      geo%atom(i)%f(1:gr%sb%dim) = ep%fii(1:gr%sb%dim, i)
    end do
    
    if (states_are_real(st) ) then 
      call dforces_from_potential(gr, geo, ep, st, time)
    else
      call zforces_from_potential(gr, geo, ep, st, time)
    end if
    
    !\todo forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, ep%no_lasers
        select case(laser_kind(ep%lasers(j)))
        case(E_FIELD_ELECTRIC)
          call laser_field(ep%lasers(j), gr%sb, x, t)
          do i = 1, geo%natoms
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) + &
              P_PROTON_CHARGE * species_zval(geo%atom(i)%spec) * x(1:gr%mesh%sb%dim)
          end do

        case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL, E_FIELD_SCALAR_POTENTIAL)
          write(message(1),'(a)') 'The forces are currently not properly calculated if time-dependent'
          write(message(2),'(a)') 'magnetic fields are present.'
          call messages_fatal(2)
        end select
      end do
    end if

    if(associated(ep%E_field)) then
      do i = 1, geo%natoms
        geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) + &
          P_PROTON_CHARGE * species_zval(geo%atom(i)%spec) * ep%E_field(1:gr%mesh%sb%dim)
      end do
    end if
    
    POP_SUB(forces_calculate)
    call profiling_out(forces_prof)

  end subroutine forces_calculate

#include "undef.F90"
#include "real.F90"
#include "forces_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "forces_inc.F90"

end module forces_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
