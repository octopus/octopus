!! Copyright (C) 2015 X. Andrade
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

module hirshfeld_m
  use messages_m
  use geometry_m
  use global_m
  use mesh_m
  use mesh_function_m
  use profiling_m
  use species_pot_m
  use states_m
  
  implicit none

  private
  public ::                 &
    hirshfeld_t,            &
    hirshfeld_init,         &
    hirshfeld_end,          &
    hirshfeld_partition

  type hirshfeld_t
    private
    type(mesh_t),     pointer     :: mesh
    type(geometry_t), pointer     :: geo
    type(states_t),   pointer     :: st
    FLOAT,            pointer     :: total_density(:)
    FLOAT,            pointer     :: free_charge(:)
  end type hirshfeld_t
    
contains

  subroutine hirshfeld_init(this, mesh, geo, st)
    type(hirshfeld_t),         intent(out)   :: this
    type(mesh_t),      target, intent(in)    :: mesh
    type(geometry_t),  target, intent(in)    :: geo
    type(states_t),    target, intent(in)    :: st
    
    integer :: iatom, ip, ispin
    FLOAT, allocatable :: atom_density(:, :)
    
    PUSH_SUB(hirshfeld_init)

    this%mesh => mesh
    this%geo  => geo
    this%st   => st

    SAFE_ALLOCATE(this%total_density(1:this%mesh%np))
    SAFE_ALLOCATE(this%free_charge(1:geo%natoms))
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))

    this%total_density = CNST(0.0)
    
    do iatom = 1, geo%natoms
      call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)
      this%free_charge(iatom) = CNST(0.0)
      do ispin = 1, st%d%nspin
        this%free_charge(iatom) = this%free_charge(iatom) + dmf_integrate(this%mesh, atom_density(:, ispin))
      end do
      forall(ip = 1:this%mesh%np) this%total_density(ip) = this%total_density(ip) + sum(atom_density(ip, 1:st%d%nspin))
    end do

    SAFE_DEALLOCATE_A(atom_density)
    
    POP_SUB(hirshfeld_init)    
  end subroutine hirshfeld_init

  ! ------------------------------------------------
  
  subroutine hirshfeld_end(this)
    type(hirshfeld_t), intent(inout) :: this

    PUSH_SUB(hirshfeld_end)

    SAFE_DEALLOCATE_P(this%total_density)
    SAFE_DEALLOCATE_P(this%free_charge)

    nullify(this%mesh)
    nullify(this%geo)
    nullify(this%st)
    
    POP_SUB(hirshfeld_end)    
  end subroutine hirshfeld_end

  ! -----------------------------------------------

  subroutine hirshfeld_partition(this, iatom, density, hirshfeld_charge, hirshfeld_volume)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: hirshfeld_charge
    FLOAT, optional,           intent(out)   :: hirshfeld_volume

    integer :: ip
    FLOAT :: dens_ip
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)
    
    PUSH_SUB(hirshfeld_partition)

    ASSERT(associated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    
    call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

    do ip = 1, this%mesh%np
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))
      if(abs(dens_ip) > CNST(1e-12)) then
        hirshfeld_density(ip) = sum(density(ip, 1:this%st%d%nspin))*dens_ip/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do

    hirshfeld_charge = dmf_integrate(this%mesh, hirshfeld_density)

    if(present(hirshfeld_volume)) then
      hirshfeld_volume = hirshfeld_charge/this%free_charge(iatom)
    end if

    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)
    
    POP_SUB(hirshfeld_partition)
  end subroutine hirshfeld_partition
  
end module hirshfeld_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
