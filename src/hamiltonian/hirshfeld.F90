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

#include "global.h"

module hirshfeld_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use ions_oct_m
  use lattice_vectors_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ps_oct_m
  use species_pot_oct_m
  use states_elec_oct_m
  use species_oct_m
  use splines_oct_m
 
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
    type(mesh_t),        pointer     :: mesh
    type(ions_t),        pointer     :: ions
    type(states_elec_t), pointer     :: st
    FLOAT,               allocatable :: total_density(:)  !< (mesh%np)
    FLOAT,               allocatable :: free_volume(:)    !< (natoms)
    FLOAT,               allocatable :: free_vol_r3(:,:)  !< (natoms,mesh%np)
  end type hirshfeld_t

  FLOAT, parameter, public :: TOL_HIRSHFELD = CNST(1e-9)
    
contains

  subroutine hirshfeld_init(this, mesh, ions, st)
    type(hirshfeld_t),           intent(out)   :: this
    type(mesh_t),        target, intent(in)    :: mesh
    type(ions_t),        target, intent(in)    :: ions
    type(states_elec_t), target, intent(in)    :: st
    
    integer :: iatom, ip, isp
    FLOAT :: rr, pos(ions%space%dim), rmax
    FLOAT, allocatable :: atom_density(:, :), atom_density_acc(:)
    type(ps_t), pointer :: ps
    type(lattice_iterator_t) :: latt_iter
    integer :: icell
    type(profile_t), save :: prof
    PUSH_SUB(hirshfeld_init)

    call profiling_in(prof, "HIRSHFELD_INIT")

    this%mesh => mesh
    this%ions  => ions
    this%st   => st

    SAFE_ALLOCATE(this%total_density(1:mesh%np))
    SAFE_ALLOCATE(this%free_volume(1:ions%natoms))
    SAFE_ALLOCATE(this%free_vol_r3(1:ions%natoms,1:mesh%np))
    SAFE_ALLOCATE(atom_density(1:mesh%np, st%d%nspin))
    SAFE_ALLOCATE(atom_density_acc(1:mesh%np))

    this%total_density = CNST(0.0)
  
    do iatom = 1, ions%natoms
      ps => species_ps(ions%atom(iatom)%species)
      atom_density_acc(1:mesh%np) = M_ZERO

      rmax = CNST(0.0)
      do isp = 1, st%d%nspin
        rmax = max(rmax, spline_cutoff_radius(ps%density(isp), ps%projectors_sphere_threshold))
      end do

      latt_iter = lattice_iterator_t(ions%latt, rmax)
      do icell = 1, latt_iter%n_cells
        pos = this%ions%pos(:, iatom) + latt_iter%get(icell)
        !We get the non periodized density
        !We need to do it to have the r^3 correctly computed for periodic systems
        call species_atom_density_np(ions%atom(iatom)%species, ions%namespace, pos, mesh, st%d%nspin, atom_density)

        do ip = 1, mesh%np
          this%total_density(ip) = this%total_density(ip) + sum(atom_density(ip, 1:st%d%nspin))
        end do

        do ip = 1, mesh%np
          rr = norm2(mesh%x(ip, :) - pos)
          atom_density_acc(ip) = atom_density_acc(ip) + sum(atom_density(ip, 1:this%st%d%nspin))*rr**3  
        end do
      end do
      this%free_volume(iatom) = dmf_integrate(this%mesh, atom_density_acc(:), reduce = .false.)
      this%free_vol_r3(iatom,:) = atom_density_acc(:)
    end do

    if(this%mesh%parallel_in_domains) then
      call this%mesh%allreduce(this%free_volume)
    end if

    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(atom_density_acc)

    call profiling_out(prof)

    POP_SUB(hirshfeld_init)    
  end subroutine hirshfeld_init

  ! ------------------------------------------------
  
  subroutine hirshfeld_end(this)
    type(hirshfeld_t), intent(inout) :: this

    PUSH_SUB(hirshfeld_end)

    SAFE_DEALLOCATE_A(this%total_density)
    SAFE_DEALLOCATE_A(this%free_volume)
    SAFE_DEALLOCATE_A(this%free_vol_r3)

    nullify(this%mesh)
    nullify(this%ions)
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
    type(profile_t), save :: prof
    
    PUSH_SUB(hirshfeld_charge)

    call profiling_in(prof, "HIRSHFELD_CHARGE")

    ASSERT(allocated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    
    call species_atom_density(this%ions%atom(iatom)%species, this%ions%namespace, this%ions%space, this%ions%latt, &
      this%ions%pos(:, iatom), this%mesh, this%st%d%nspin, atom_density)

    do ip = 1, this%mesh%np
      dens_ip = sum(atom_density(ip, 1:this%st%d%nspin))
      if(abs(dens_ip) > TOL_HIRSHFELD) then
        hirshfeld_density(ip) = sum(density(ip, 1:this%st%d%nspin))*dens_ip/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do

    charge = dmf_integrate(this%mesh, hirshfeld_density)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)

    call profiling_out(prof)
    
    POP_SUB(hirshfeld_charge)
  end subroutine hirshfeld_charge

  ! -----------------------------------------------    

  subroutine hirshfeld_volume_ratio(this, iatom, density, volume_ratio)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: volume_ratio

    integer :: ip
    FLOAT, allocatable :: atom_density(:, :), hirshfeld_density(:)

    type(profile_t), save :: prof
    
    PUSH_SUB(hirshfeld_volume_ratio)

    call profiling_in(prof, "HIRSHFELD_VOLUME_RATIO")

    ASSERT(allocated(this%total_density))
    
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    SAFE_ALLOCATE(hirshfeld_density(1:this%mesh%np))
    

    do ip = 1, this%mesh%np
      if( this%total_density(ip) > TOL_HIRSHFELD) then
        hirshfeld_density(ip) = this%free_vol_r3(iatom,ip)*sum(density(ip, 1:this%st%d%nspin))/this%total_density(ip)
      else
        hirshfeld_density(ip) = CNST(0.0)
      end if
    end do
    
    volume_ratio = dmf_integrate(this%mesh, hirshfeld_density)/this%free_volume(iatom)
    
    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(hirshfeld_density)

    call profiling_out(prof)
    
    POP_SUB(hirshfeld_volume_ratio)
  end subroutine hirshfeld_volume_ratio
  
  ! -----------------------------------------------

  subroutine hirshfeld_density_derivative(this, iatom, ddensity)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    FLOAT,                     intent(out)   :: ddensity(:)

    integer :: ip
    FLOAT, allocatable :: atom_density(:, :)
    type(profile_t), save :: prof   
 
    PUSH_SUB(hirshfeld_density_derivative)

    call profiling_in(prof, "HIRSHFELD_DENSITY_DER")

    SAFE_ALLOCATE(atom_density(1:this%mesh%np, this%st%d%nspin))
    
    do ip = 1, this%mesh%np

      if(abs(this%total_density(ip)) > TOL_HIRSHFELD) then
        ddensity(ip) = this%free_vol_r3(iatom,ip)/(this%total_density(ip)*this%free_volume(iatom))
      else
        ddensity(ip) = CNST(0.0)
      end if
    end do

    SAFE_DEALLOCATE_A(atom_density)
    
    call profiling_out(prof)
    
    POP_SUB(hirshfeld_density_derivative)
  end subroutine hirshfeld_density_derivative

  ! -----------------------------------------------
  !dvadrr_ij = \frac{\delta V_i}{\delta \vec{x_j}}
  subroutine hirshfeld_position_derivative(this, iatom, jatom, density, dposition)
    type(hirshfeld_t),         intent(in)    :: this
    integer,                   intent(in)    :: iatom
    integer,                   intent(in)    :: jatom
    FLOAT,                     intent(in)    :: density(:, :)
    FLOAT,                     intent(out)   :: dposition(:)

    integer :: ip, idir, icell, jcell, isp
    FLOAT :: atom_dens, atom_der,rri, rrj, tdensity, pos_i(this%ions%space%dim), pos_j(this%ions%space%dim), rmax_i, rmax_j, &
             rij, rmax_isqu, rmax_jsqu
    FLOAT, allocatable :: grad(:, :), atom_density(:, :), atom_derivative(:, :)
    type(lattice_iterator_t) :: latt_iter_i, latt_iter_j
    type(ps_t), pointer :: ps_i, ps_j
    type(profile_t), save :: prof
    FLOAT :: tmp, xxi(this%ions%space%dim), xxj(this%ions%space%dim)

    FLOAT :: TOL_SPACING

    PUSH_SUB(hirshfeld_position_derivative)

    TOL_SPACING = maxval(this%mesh%spacing(1:this%ions%space%dim))

    call profiling_in(prof, "HIRSHFELD_POSITION_DER")


    SAFE_ALLOCATE(grad(1:this%mesh%np, 1:this%ions%space%dim))
    SAFE_ALLOCATE(atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))
    SAFE_ALLOCATE(atom_density(1:this%mesh%np, 1:this%st%d%nspin))

    dposition(1:this%ions%space%dim) = M_ZERO
    grad(1:this%mesh%np, 1:this%ions%space%dim) = M_ZERO

    ps_i => species_ps(this%ions%atom(iatom)%species)
    ps_j => species_ps(this%ions%atom(jatom)%species)

    rmax_i = CNST(0.0)
    rmax_j = CNST(0.0)
    do isp = 1, this%st%d%nspin
      rmax_i = max(rmax_i, spline_cutoff_radius(ps_i%density(isp), ps_i%projectors_sphere_threshold))
      rmax_j = max(rmax_j, spline_cutoff_radius(ps_j%density_der(isp), ps_j%projectors_sphere_threshold))
    end do

    rmax_isqu = rmax_i**2
    rmax_jsqu = rmax_j**2

    latt_iter_j = lattice_iterator_t(this%ions%latt, rmax_j)
    do jcell = 1, latt_iter_j%n_cells

      pos_j = this%ions%pos(:, jatom) + latt_iter_j%get(jcell)
      atom_derivative(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO
      call species_atom_density_derivative_np(this%ions%atom(jatom)%species, this%ions%namespace, pos_j, this%mesh, &
        this%st%d%spin_channels, atom_derivative(1:this%mesh%np, 1:this%st%d%nspin))

      latt_iter_i = lattice_iterator_t(this%ions%latt, (rmax_j+rmax_i)) ! jcells further away from this distance cannot respect the following 'if' condition with respect to the i atom in this icell
      do icell = 1, latt_iter_i%n_cells

        pos_i = pos_j + latt_iter_i%get(icell) + (this%ions%pos(:, iatom) - this%ions%pos(:, jatom))
        rij =  norm2(pos_i - pos_j)
          
        if(rij - (rmax_j+rmax_i) < TOL_SPACING) then 
 

          atom_density(1:this%mesh%np, 1:this%st%d%nspin) = M_ZERO

          !We get the non periodized density
          !We need to do it to have the r^3 correctly computed for periodic systems
          call species_atom_density_np(this%ions%atom(iatom)%species, this%ions%namespace, pos_i, this%mesh, this%st%d%nspin, &
            atom_density(1:this%mesh%np, 1:this%st%d%nspin))

          do ip = 1, this%mesh%np
            if(this%total_density(ip)< TOL_HIRSHFELD) cycle
            
            xxi = this%mesh%x(ip, :) - pos_i
            rri = sum(xxi**2)
            if(rri - rmax_isqu > TOL_SPACING) cycle ! In this case atom_dens = 0

            xxj = this%mesh%x(ip, :) - pos_j
            rrj = sum(xxj**2)
            if(rrj - rmax_jsqu > TOL_SPACING) cycle ! In this case atom_der = 0

            rri = sqrt(rri)
            rrj = sqrt(rrj)
              
            tdensity = sum(density(ip, 1:this%st%d%nspin))
            atom_dens = sum(atom_density(ip, 1:this%st%d%nspin))

            tmp = rri**3*atom_dens*tdensity/this%total_density(ip)**2

            atom_der = sum(atom_derivative(ip, 1:this%st%d%nspin))

            if(rrj > TOL_HIRSHFELD) then
              do idir = 1, this%ions%space%dim
                grad(ip, idir) = grad(ip, idir) - tmp*atom_der*xxj(idir)/rrj
              end do
            end if

            !Only if we really have the same atoms
            if(iatom == jatom .and. rij < TOL_HIRSHFELD) then
              grad(ip, :) = grad(ip, :) + (CNST(3.0)*rri*atom_dens + rri**2*atom_der)*tdensity/this%total_density(ip)*xxi
            end if

          end do
            
        end if
      end do
    end do

    do idir = 1, this%ions%space%dim
      dposition(idir) = dmf_integrate(this%mesh, grad(1:this%mesh%np, idir), reduce = .false.) &
                             /this%free_volume(iatom)
    end do

    if(this%mesh%parallel_in_domains) then
      call this%mesh%allreduce(dposition, dim = this%ions%space%dim)
    end if

    SAFE_DEALLOCATE_A(atom_density)
    SAFE_DEALLOCATE_A(atom_derivative)
    SAFE_DEALLOCATE_A(grad)

    call profiling_out(prof)

    POP_SUB(hirshfeld_position_derivative)
  end subroutine hirshfeld_position_derivative
  
end module hirshfeld_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
