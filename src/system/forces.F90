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

#include "global.h"

module forces_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use density_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use index_oct_m
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use species_oct_m
  use species_pot_oct_m
  use states_oct_m
  use states_dim_oct_m
  use symm_op_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                    &
    forces_calculate,          &
    dforces_from_potential,    &
    zforces_from_potential,    &
    total_force_calculate,     &
    forces_write_info

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
  subroutine forces_calculate(gr, geo, hm, st, ks, vhxc_old, t, dt)
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: st
    type(v_ks_t),        intent(in)      :: ks
    FLOAT,     optional, intent(in)    :: vhxc_old(:,:)
    FLOAT,     optional, intent(in)    :: t
    FLOAT,     optional, intent(in)    :: dt

    integer :: j, iatom, idir
    FLOAT :: x(MAX_DIM), time, global_force(1:MAX_DIM)
    FLOAT, allocatable :: force(:, :), force_loc(:, :), force_nl(:, :), force_u(:, :)
    FLOAT, allocatable :: force_scf(:, :)
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(forces_calculate)

    x(:) = M_ZERO
    time = M_ZERO
    if(present(t)) time = t

    !We initialize the different components of the force to zero
    do iatom = 1, geo%natoms
      geo%atom(iatom)%f_ii(1:gr%sb%dim) = M_ZERO
      geo%atom(iatom)%f_loc(1:gr%sb%dim) = M_ZERO
      geo%atom(iatom)%f_nl(1:gr%sb%dim) = M_ZERO
      geo%atom(iatom)%f_u(1:gr%sb%dim) = M_ZERO
      geo%atom(iatom)%f_fields(1:gr%sb%dim) = M_ZERO
      geo%atom(iatom)%f_scf(1:gr%sb%dim) = M_ZERO
    end do

    do iatom = 1, geo%natoms
      geo%atom(iatom)%f(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, iatom)
      geo%atom(iatom)%f_ii(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, iatom)
    end do

    if(present(t)) then
      call epot_global_force(hm%ep, geo, time, global_force)

      ! the ion-ion term is already calculated
      do iatom = 1, geo%natoms
        geo%atom(iatom)%f(1:gr%sb%dim) = geo%atom(iatom)%f(1:gr%sb%dim) + global_force(1:gr%sb%dim)
        geo%atom(iatom)%f_ii(1:gr%sb%dim) = geo%atom(iatom)%f_ii(1:gr%sb%dim) + global_force(1:gr%sb%dim)
      end do
    end if

    SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
    SAFE_ALLOCATE(force_loc(1:gr%mesh%sb%dim, 1:geo%natoms))
    SAFE_ALLOCATE(force_nl(1:gr%mesh%sb%dim, 1:geo%natoms))
    SAFE_ALLOCATE(force_u(1:gr%mesh%sb%dim, 1:geo%natoms)) 
    SAFE_ALLOCATE(force_scf(1:gr%mesh%sb%dim, 1:geo%natoms))

    if (states_are_real(st) ) then 
      call dforces_from_potential(gr, geo, hm, st, force, force_loc, force_nl, force_u)
    else
      call zforces_from_potential(gr, geo, hm, st, force, force_loc, force_nl, force_u)
    end if

    if(present(vhxc_old)) then
      call forces_from_scf(gr, geo, hm, st, force_scf, vhxc_old)
    else
      force_scf = M_ZERO
    end if

    if(hm%ep%force_total_enforce) then
      call forces_set_total_to_zero(geo, force)
      call forces_set_total_to_zero(geo, force_loc)
      call forces_set_total_to_zero(geo, force_nl)
      call forces_set_total_to_zero(geo, force_u)
      call forces_set_total_to_zero(geo, force_scf)
    end if

    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + force(idir, iatom) + force_scf(idir, iatom)
        geo%atom(iatom)%f_loc(idir) = force_loc(idir, iatom)
        geo%atom(iatom)%f_nl(idir) = force_nl(idir, iatom)
        geo%atom(iatom)%f_u(idir) = force_u(idir, iatom)
        geo%atom(iatom)%f_scf(idir) = force_scf(idir, iatom)
      end do
    end do

    SAFE_DEALLOCATE_A(force)
    SAFE_DEALLOCATE_A(force_loc)
    SAFE_DEALLOCATE_A(force_nl)
    SAFE_DEALLOCATE_A(force_u)
    SAFE_DEALLOCATE_A(force_scf)
 
    if(associated(hm%ep%E_field)) then
      do iatom = 1, geo%natoms
        ! Here the proton charge is +1, since the electric field has the usual sign.
        geo%atom(iatom)%f(1:gr%mesh%sb%dim) = geo%atom(iatom)%f(1:gr%mesh%sb%dim) &
          + species_zval(geo%atom(iatom)%species)*hm%ep%E_field(1:gr%mesh%sb%dim)
        geo%atom(iatom)%f_fields(1:gr%mesh%sb%dim) = geo%atom(iatom)%f_fields(1:gr%mesh%sb%dim) &
               + species_zval(geo%atom(iatom)%species)*hm%ep%E_field(1:gr%mesh%sb%dim)
      end do
    end if
    
    POP_SUB(forces_calculate)
    call profiling_out(forces_prof)

  end subroutine forces_calculate

  ! ----------------------------------------------------------------------

  subroutine forces_set_total_to_zero(geo, force)
    type(geometry_t),    intent(in)    :: geo
    FLOAT,               intent(inout) :: force(:, :)

    FLOAT, allocatable :: total_force(:)
    integer :: iatom

    PUSH_SUB(forces_set_total_to_zero)

    SAFE_ALLOCATE(total_force(1:geo%space%dim))

    total_force(1:geo%space%dim) = CNST(0.0)
    do iatom = 1, geo%natoms
      total_force(1:geo%space%dim) = total_force(1:geo%space%dim) + force(1:geo%space%dim, iatom)/geo%natoms
    end do

    do iatom = 1, geo%natoms
      force(1:geo%space%dim, iatom) = force(1:geo%space%dim, iatom) - total_force(1:geo%space%dim)
    end do

    SAFE_DEALLOCATE_A(total_force)
    POP_SUB(forces_set_total_to_zero)
  end subroutine forces_set_total_to_zero

 ! ----------------------------------------------------------------------
  subroutine forces_write_info(iunit, geo, sb, dir)
    integer,             intent(in)    :: iunit
    type(geometry_t),    intent(in)    :: geo
    type(simul_box_t),   intent(in)    :: sb
    character(len=*),    intent(in)    :: dir

    integer :: iatom, idir, ii, iunit2
    FLOAT:: rr(1:3), ff(1:3), torque(1:3)

    if(.not.mpi_grp_is_root(mpi_world)) return    

    PUSH_SUB(forces_write_info)

    write(iunit,'(3a)') 'Forces on the ions [', trim(units_abbrev(units_out%force)), "]"
    write(iunit,'(a,10x,99(14x,a))') ' Ion', (index2axis(idir), idir = 1, sb%dim)
    do iatom = 1, geo%natoms
      write(iunit,'(i4,a10,10f15.6)') iatom, trim(species_label(geo%atom(iatom)%species)), &
              (units_from_atomic(units_out%force, geo%atom(iatom)%f(idir)), idir=1, sb%dim)
    end do
    write(iunit,'(1x,100a1)') ("-", ii = 1, 13 + sb%dim * 15)
    write(iunit,'(a14, 10f15.6)') " Max abs force", &
            (units_from_atomic(units_out%force, maxval(abs(geo%atom(1:geo%natoms)%f(idir)))), idir=1, sb%dim)
    write(iunit,'(a14, 10f15.6)') " Total force", &
            (units_from_atomic(units_out%force, sum(geo%atom(1:geo%natoms)%f(idir))), idir=1, sb%dim)

    if(geo%space%dim == 2 .or. geo%space%dim == 3) then
      rr = M_ZERO
      ff = M_ZERO
      torque = M_ZERO
      do iatom = 1, geo%natoms
        rr(1:geo%space%dim) = geo%atom(iatom)%x(1:geo%space%dim)
        ff(1:geo%space%dim) = geo%atom(iatom)%f(1:geo%space%dim)
        torque(1:3) = torque(1:3) + dcross_product(rr, ff)
      end do
      write(iunit,'(a14, 10f15.6)') ' Total torque', &
              (units_from_atomic(units_out%force*units_out%length, torque(idir)), idir = 1, 3)
    end if


    iunit2 = io_open(trim(dir)//'/forces', action='write', position='asis')
    write(iunit2,'(a)') ' # Total force (x,y,z) Ion-Ion (x,y,z) Local (x,y,z) NL (x,y,z) Fields (x,y,z) Hubbard(x,y,z) SCF(x,y,z)'
    do iatom = 1, geo%natoms
       write(iunit2,'(i4,a10,24e15.6)') iatom, trim(species_label(geo%atom(iatom)%species)), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_ii(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_loc(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_nl(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_fields(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_u(idir)), idir=1, sb%dim), &
                 (units_from_atomic(units_out%force, geo%atom(iatom)%f_scf(idir)), idir=1, sb%dim)
    end do
    call io_close(iunit2) 

    POP_SUB(forces_write_info)

  end subroutine forces_write_info

  ! Implementation of the term from Chan et al.,  Phys. Rev. B 47, 4771 (1993).
  ! Here we make the approximation that the "atomic densities" are just the one 
  ! from the pseudopotential.  
  ! NTD : No idea if this is good or bad, but this is easy to implement 
  !       and works well in practice
  subroutine forces_from_scf(gr, geo, hm, st, force_scf, vhxc_old)
    type(grid_t),                   intent(inout) :: gr
    type(geometry_t),               intent(inout) :: geo
    type(hamiltonian_t),            intent(in)    :: hm
    type(states_t),                 intent(inout) :: st
    FLOAT,                          intent(out)   :: force_scf(:, :)
    FLOAT,                          intent(in)    :: vhxc_old(:,:)

    integer :: is, iatom, idir
    FLOAT, allocatable :: dvhxc(:,:), drho(:,:,:)

    PUSH_SUB(forces_from_scf)

    SAFE_ALLOCATE(dvhxc(1:gr%mesh%np, 1:hm%d%spin_channels))
    SAFE_ALLOCATE(drho(1:gr%mesh%np, 1:hm%d%spin_channels, 1:gr%mesh%sb%dim))

    !We average over spin channels
    do is = 1, hm%d%spin_channels
      dvhxc(1:gr%mesh%np, is) = hm%vhxc(1:gr%mesh%np, is) - vhxc_old(1:gr%mesh%np, is)
    end do

    force_scf = M_ZERO

    do iatom = geo%atoms_dist%start, geo%atoms_dist%end
      if(species_type(geo%atom(iatom)%species) == SPECIES_PSEUDO) then
        if(ps_has_density(species_ps(geo%atom(iatom)%species))) then

          call species_atom_density_grad(gr%mesh, gr%mesh%sb, geo%atom(iatom), &
            hm%d%spin_channels, drho)

          do idir = 1, gr%mesh%sb%dim
            do is = 1, hm%d%spin_channels
              force_scf(idir, iatom) = force_scf(idir, iatom) &
                -dmf_dotp(gr%mesh, drho(:,is,idir), dvhxc(:,is))
            end do
          end do
        end if
      end if
    end do

    SAFE_DEALLOCATE_A(dvhxc)
    SAFE_DEALLOCATE_A(drho)

    if(geo%atoms_dist%parallel) call dforces_gather(geo, force_scf) 

    POP_SUB(forces_from_scf)
  end subroutine forces_from_scf

#include "undef.F90"
#include "real.F90"
#include "forces_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "forces_inc.F90"

end module forces_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
