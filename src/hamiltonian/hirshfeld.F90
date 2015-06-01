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
  use derivatives_m
  use messages_m
  use geometry_m
  use global_m
  use mesh_m
  use mesh_function_m
  use profiling_m
  use ps_m
  use species_pot_m
  use states_m
  use species_m
  
  implicit none

  private
  public ::                       &
    hirshfeld_t,                  &
    hirshfeld_init,               &
    hirshfeld_end,                &
    hirshfeld_charge,             &
    hirshfeld_volume_ratio,       &
    hirshfeld_density_derivative, &
    hirshfeld_position_derivative
  
  type hirshfeld_t
    private
    type(mesh_t),     pointer     :: mesh
    type(geometry_t), pointer     :: geo
    type(states_t),   pointer     :: st
    FLOAT,            pointer     :: total_density(:)
    FLOAT,            pointer     :: free_volume(:)
  end type hirshfeld_t
    
contains

  subroutine hirshfeld_init(this, mesh, geo, st)
    type(hirshfeld_t),         intent(out)   :: this
    type(mesh_t),      target, intent(in)    :: mesh
    type(geometry_t),  target, intent(in)    :: geo
    type(states_t),    target, intent(in)    :: st
    
    integer :: iatom, ip
    FLOAT :: rr
    FLOAT, allocatable :: atom_density(:, :)
    
    PUSH_SUB(hirshfeld_init)

    this%mesh => mesh
    this%geo  => geo
    this%st   => st

    SAFE_ALLOCATE(this%total_density(1:this%mesh%np))
    SAFE_ALLOCATE(this%free_volume(1:geo%natoms))
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))

    this%total_density = CNST(0.0)
    
    do iatom = 1, geo%natoms
      call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

      forall(ip = 1:this%mesh%np) this%total_density(ip) = this%total_density(ip) + sum(atom_density(ip, 1:st%d%nspin))
      
      do ip = 1, this%mesh%np
        rr = sqrt(sum((this%mesh%x(ip, 1:this%mesh%sb%dim) - this%geo%atom(iatom)%x(1:this%mesh%sb%dim))**2))
        atom_density(ip, 1) = sum(atom_density(ip, 1:this%st%d%nspin))*rr**3
      end do

      this%free_volume(iatom) = dmf_integrate(this%mesh, atom_density(:, 1))
    end do

    SAFE_DEALLOCATE_A(atom_density)
    
    POP_SUB(hirshfeld_init)    
  end subroutine hirshfeld_init

  ! ------------------------------------------------
  
  subroutine hirshfeld_end(this)
    type(hirshfeld_t), intent(inout) :: this

    PUSH_SUB(hirshfeld_end)

    SAFE_DEALLOCATE_P(this%total_density)
    SAFE_DEALLOCATE_P(this%free_volume)

    nullify(this%mesh)
    nullify(this%geo)
    nullify(this%st)
    
    POP_SUB(hirshfeld_end)    
  end subroutine hirshfeld_end

  ! -----------------------------------------------
  
  subroutine hirshfeld_charge(this, iatom, density, charge)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: charge

    integer :: ip
    FLOAT :: dens_ip
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)
    
    PUSH_SUB(hirshfeld_charge)

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

    charge = dmf_integrate(this%mesh, hirshfeld_density)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)
    
    POP_SUB(hirshfeld_charge)
  end subroutine hirshfeld_charge

  ! -----------------------------------------------    

  subroutine hirshfeld_volume_ratio(this, iatom, density, volume_ratio)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: volume_ratio

    integer :: ip
    FLOAT :: dens_ip, rr
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)
    
    PUSH_SUB(hirshfeld_volume_ratio)

    ASSERT(associated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    
    call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

    do ip = 1, this%mesh%np
      rr = sqrt(sum((this%mesh%x(ip, 1:this%mesh%sb%dim) - this%geo%atom(iatom)%x(1:this%mesh%sb%dim))**2))
            
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))
      if(abs(dens_ip) > CNST(1e-12)) then
        hirshfeld_density(ip) = rr**3*sum(density(ip, 1:this%st%d%nspin))*dens_ip/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do
    
    volume_ratio = dmf_integrate(this%mesh, hirshfeld_density)/this%free_volume(iatom)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)
    
    POP_SUB(hirshfeld_volume_ratio)
  end subroutine hirshfeld_volume_ratio
  
  ! -----------------------------------------------

  subroutine hirshfeld_density_derivative(this, iatom, ddensity)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(out)   :: ddensity(:)

    integer :: ip
    FLOAT :: dens_ip, rr
    FLOAT, allocatable :: atom_density(:, :)
    
    PUSH_SUB(hirshfeld_density_derivative)

    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    
    call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

    do ip = 1, this%mesh%np
      rr = sqrt(sum((this%mesh%x(ip, 1:this%mesh%sb%dim) - this%geo%atom(iatom)%x(1:this%mesh%sb%dim))**2))
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))

      if(abs(dens_ip) > CNST(1e-12)) then
        ddensity(ip) = rr**3*dens_ip/(this%total_density(ip)*this%free_volume(iatom))
      else
        ddensity(ip) = CNST(0.0)
      end if

    end do

    SAFE_DEALLOCATE_A(atom_density)
    
    POP_SUB(hirshfeld_density_derivative)
  end subroutine hirshfeld_density_derivative

  ! -----------------------------------------------

  subroutine hirshfeld_position_derivative(this, der, iatom, density, dposition)
    type(hirshfeld_t),         intent(in)    :: this
    type(derivatives_t),       intent(in)    :: der
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: dposition(:)

    integer :: ip, idir
    FLOAT :: dens_ip, rr
    FLOAT, allocatable :: tdensity(:), grad(:, :), atom_density(:, :)
    
    PUSH_SUB(hirshfeld_position_derivative)

    SAFE_ALLOCATE(atom_density(1:this%mesh%np, 1:this%st%d%nspin))
    SAFE_ALLOCATE(tdensity(1:this%mesh%np_part))
    SAFE_ALLOCATE(grad(1:this%mesh%np, 1:this%mesh%sb%dim))
    
    call species_atom_density(this%mesh, this%mesh%sb, this%geo%atom(iatom), this%st%d%nspin, atom_density)

    ! we use the same trick as for the forces:
    !  df(r-R)/dR = f(r-R)*d/dr
    
    do ip = 1, this%mesh%np
      tdensity(ip) = sum(density(ip, 1:this%st%d%nspin))
    end do

    call dderivatives_grad(der, tdensity, grad)
    
    do ip = 1, this%mesh%np
      rr = sqrt(sum((this%mesh%x(ip, 1:this%mesh%sb%dim) - this%geo%atom(iatom)%x(1:this%mesh%sb%dim))**2))
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))

      do idir = 1, this%mesh%sb%dim
        if(abs(dens_ip) > CNST(1e-12)) then
          grad(ip, idir) = grad(ip, idir)*rr**3*dens_ip/this%total_density(ip)
        else
          grad(ip, idir) = CNST(0.0)
        end if
      end do
      
    end do

    do idir = 1, this%mesh%sb%dim
      dposition(idir) = -dmf_integrate(this%mesh, grad(:, idir))/this%free_volume(iatom)
    end do
    
    SAFE_DEALLOCATE_A(atom_density)
    
    POP_SUB(hirshfeld_position_derivative)
  end subroutine hirshfeld_position_derivative
  
end module hirshfeld_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
