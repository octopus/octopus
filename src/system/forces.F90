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

#include "global.h"

module forces_m
  use batch_m
  use batch_ops_m
  use born_charges_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use comm_m
  use derivatives_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use lasers_m
  use linear_response_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use projector_m
  use octcl_kernel_m
  use opencl_m
  use simul_box_m
  use species_m
  use species_pot_m
  use states_m
  use states_dim_m
  use symm_op_m
  use symmetrizer_m
  use types_m

  implicit none

  private
  public ::                    &
    forces_calculate,          &
    dforces_from_potential,    &
    zforces_from_potential,    &
    dforces_derivative,        &
    zforces_derivative,        &
    dforces_born_charges,      &
    zforces_born_charges,      &
    total_force_calculate

  type(profile_t), save :: prof_comm

contains

  ! ---------------------------------------------------------
  !> This computes the total forces on the ions created by the electrons
  !! (it excludes the force due to possible time-dependent external fields).
  subroutine total_force_calculate(gr, geo, ep, st, x)
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(in)    :: geo
    type(epot_t),     intent(inout) :: ep
    type(states_t),   intent(inout) :: st
    FLOAT, intent(inout)            :: x(MAX_DIM)

    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(total_force_calculate)

    x = M_ZERO
    if (states_are_real(st) ) then 
      call dtotal_force_from_potential(gr, geo, ep, st, x)
    else
      call ztotal_force_from_potential(gr, geo, ep, st, x)
    end if

    POP_SUB(total_force_calculate)
    call profiling_out(forces_prof)
  end subroutine total_force_calculate


  ! ---------------------------------------------------------
  subroutine forces_calculate(gr, geo, hm, st, t, dt)
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: st
    FLOAT,     optional, intent(in)    :: t
    FLOAT,     optional, intent(in)    :: dt

    integer :: i, j, iatom, idir
    FLOAT :: x(MAX_DIM), time
    FLOAT, allocatable :: force(:, :)
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(forces_calculate)

    x(:) = M_ZERO
    time = M_ZERO
    if(present(t)) time = t

    ! the ion-ion term is already calculated
    do i = 1, geo%natoms
      geo%atom(i)%f(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, i)
    end do

    SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
    
    if (states_are_real(st) ) then 
      call dforces_from_potential(gr, geo, hm, st, force)
    else
      call zforces_from_potential(gr, geo, hm, st, force)
    end if

    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + force(idir, iatom)
      end do
    end do

    SAFE_DEALLOCATE_A(force)
    
    !\todo forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(j)))
        case(E_FIELD_ELECTRIC)
          x(1:gr%sb%dim) = M_ZERO
          call laser_field(hm%ep%lasers(j), x(1:gr%sb%dim), t)
          do i = 1, geo%natoms
            ! Here the proton charge is +1, since the electric field has the usual sign.
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
             + species_zval(geo%atom(i)%spec)*x(1:gr%mesh%sb%dim)
          end do
    
        case(E_FIELD_VECTOR_POTENTIAL)
          ! Forces are correctly calculated only if the time-dependent
          ! vector potential has no spatial dependence.
          ! The full force taking account of the spatial dependence of A should be:
          ! F = q [- dA/dt + v x \nabla x A]

          x(1:gr%sb%dim) = M_ZERO
          call laser_electric_field(hm%ep%lasers(j), x(1:gr%sb%dim), t, dt) !convert in E field (E = -dA/ c dt)
          do i = 1, geo%natoms
            ! Also here the proton charge is +1
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
             + species_zval(geo%atom(i)%spec)*x(1:gr%mesh%sb%dim)
          end do

        case(E_FIELD_MAGNETIC, E_FIELD_SCALAR_POTENTIAL)
          write(message(1),'(a)') 'The forces are currently not properly calculated if time-dependent'
          write(message(2),'(a)') 'magnetic fields are present.'
          call messages_fatal(2)
        end select
      end do
    end if

    if(associated(hm%ep%E_field)) then
      do i = 1, geo%natoms
        ! Here the proton charge is +1, since the electric field has the usual sign.
        geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
          + species_zval(geo%atom(i)%spec)*hm%ep%E_field(1:gr%mesh%sb%dim)
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
