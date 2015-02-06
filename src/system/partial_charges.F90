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
!! $Id:$

#include "global.h"

! This module is an interface to use the routines of the Bader program,
! located in external_libs/bader/.
!
!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
! Version 0.28a (07/12/12)
!
! Authors:
!   Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Authors of the multipole code:
!   Sebastien Lebegue <Sebastien.Lebegue@crm2.uhp-nancy.fr>
!   Angyan Janos <Janos.Angyan@crm2.uhp-nancy.fr>
!   Emmanuel Aubert <emmanuel.aubert@crm2.uhp-nancy.fr>
!
! Contributers:
!   Johannes Voss (DTU), Erik McNellis (FHI), Matthew Dyer (Liverpool),
!   Sören Wohlthat (Sydney)
!
! Based on algorithms described in the following publications:
!
!   A fast and robust algorithm for Bader decomposition of charge density
!   G. Henkelman, A. Arnaldsson, and H. Jonsson
!   Comput. Mater. Sci. 36, 254-360 (2006).
!
!   An improved grid-based algorithm for Bader charge allocation
!   E. Sanville, S. Kenny, R. Smith, and G. Henkelman
!   J. Comput. Chem. 28, 899-908 (2007).
!
!   A grid-based Bader analysis algorithm without lattice bias
!   W. Tang, E. Sanville, and G. Henkelman
!   J. Phys.: Condens. Matter 21, 084204 (2009)
!
!-----------------------------------------------------------------------------------!

module partial_charges_m
  use batch_m
  use batch_ops_m
  use comm_m
  use cube_m
  use cube_function_m
  use derivatives_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
  use profiling_m
  use projector_m
  use ps_m
  use restart_m
  use simul_box_m
  use species_m
  use splines_m
  use states_m
  use states_dim_m
  use submesh_m
  use symmetries_m
  use symmetrizer_m
  use symm_op_m
  use unit_m
  use unit_system_m
  use varinfo_m

  ! these modules comes from the bader program

  use charge_mod
  use ions_mod
  use kind_mod
  use options_mod
  use matrix_mod
  use bader_mod
  use voronoi_mod

  implicit none

  private

  type partial_charges_t
    type(options_obj) :: options
  end type partial_charges_t

  public ::                             &
    partial_charges_t,                            &
    partial_charges_init,                         &
    partial_charges_end,                          &
    partial_charges_calculate

contains

  subroutine partial_charges_init(this)
    type(partial_charges_t), intent(out)   :: this

    PUSH_SUB(partial_charges_init)

    !default options
    this%options%out_opt = this%options%out_chgcar4
    this%options%in_opt = this%options%in_auto
    ! print options
    this%options%vac_flag = .false.
    this%options%weight_flag = .false.
    this%options%vacval = 1E-3
    this%options%print_all_atom = .false.
    this%options%print_all_bader = .false.
    this%options%print_sel_atom = .false.
    this%options%print_sel_bader = .false.
    this%options%print_sum_atom = .false.
    this%options%print_sum_bader = .false.
    this%options%print_bader_index = .false.
    this%options%print_atom_index = .false.
    ! end of print options
    this%options%bader_opt = this%options%bader_neargrid
    this%options%quit_opt = this%options%quit_known
    this%options%refine_edge_itrs = -1
    this%options%bader_flag = .true.
    this%options%voronoi_flag = .false.
    this%options%dipole_flag = .false.
    this%options%ldos_flag = .false.
    this%options%verbose_flag = .false.
    this%options%badertol = 1.0e-4_q2
    this%options%stepsize = 0.0_q2
    this%options%ref_flag = .false.
    

    POP_SUB(partial_charges_init)
  end subroutine partial_charges_init

  !----------------------------------------------

  subroutine partial_charges_calculate(this, mesh, st, geo, bader_charges, voronoi_charges)
    type(partial_charges_t), intent(in)    :: this
    type(mesh_t),            intent(in)    :: mesh
    type(states_t),          intent(in)    :: st
    type(geometry_t),        intent(in)    :: geo
    FLOAT, optional,         intent(out)   :: bader_charges(:)
    FLOAT, optional,         intent(out)   :: voronoi_charges(:)

    integer :: idir, iatom, ix, iy, iz, ip
    FLOAT :: vol
    type(ions_obj)   :: ions
    type(charge_obj) :: charge
    type(cube_t) :: cube
    type(cube_function_t) :: density_cube
    type(bader_obj) :: bdr
    type(voronoi_obj) :: vor
    real(q2) :: dlat(1:3), dcar(1:3)
    FLOAT    :: offset(MAX_DIM)
    FLOAT, allocatable :: total_density(:)
    type(profile_t), save :: prof

    PUSH_SUB(partial_charges_calculate)
    call profiling_in(prof, 'PARTIAL_CHARGES')

    ! put values in a nice cube
    call cube_init(cube, mesh%idx%ll, mesh%sb)
    call cube_function_null(density_cube)
    call dcube_function_alloc_RS(cube, density_cube)

    SAFE_ALLOCATE(total_density(1:mesh%np))
    
    do ip = 1, mesh%np
      total_density(ip) = sum(st%rho(ip, 1:st%d%nspin))
    end do

    call dmesh_to_cube(mesh, total_density, cube, density_cube, local = .true.)

    SAFE_DEALLOCATE_A(total_density)

    ions%nions = geo%natoms

    SAFE_ALLOCATE(ions%r_car(1:ions%nions, 3))
    SAFE_ALLOCATE(ions%r_dir(1:ions%nions, 3))
    SAFE_ALLOCATE(ions%ion_chg(1:ions%nions))
    SAFE_ALLOCATE(ions%atomic_num(1:ions%nions))

    offset = M_ZERO
    offset(1:3) = units_from_atomic(units_out%length, -matmul(mesh%sb%rlattice_primitive(1:3,1:3), mesh%sb%lsize(1:3)))

    do idir = mesh%sb%periodic_dim+1, 3
      offset(idir) = units_from_atomic(units_out%length, -(cube%rs_n_global(idir) - 1)/2*mesh%spacing(idir))
    end do

    charge%org_car(1:3) = offset(1:3)

    do idir = 1, 3
      ions%lattice(idir, 1:3) = mesh%spacing(idir)*mesh%sb%rlattice_primitive(1:3, idir)
    end do

    charge%npts(1:3) = cube%rs_n_global(1:3)

    charge%lat2car = transpose(ions%lattice)

    charge%i_npts = 1.0_q2/real(charge%npts, q2)

    do idir = 1, 3
      ions%lattice(idir, 1:3) = ions%lattice(idir, 1:3)*real(charge%npts(idir), q2)
    end do

    ions%dir2car = transpose(ions%lattice)
    call matrix_3x3_inverse(ions%dir2car, ions%car2dir)
    call matrix_3x3_inverse(charge%lat2car, charge%car2lat)

    vol = matrix_volume(ions%lattice)
    do iatom = 1, geo%natoms
      ions%atomic_num(iatom) = species_z(geo%atom(iatom)%spec)
      ions%ion_chg(iatom) = 0.0_q2
      ions%r_car(iatom, 1:3) = geo%atom(iatom)%x(1:3)
      call matrix_vector(ions%car2dir, ions%r_car(iatom, :) - charge%org_car(:), ions%r_dir(iatom, :))
    end do

    charge%org_lat = (/1.0_q2, 1.0_q2, 1.0_q2/)

    SAFE_ALLOCATE(ions%r_lat(1:ions%nions, 1:3))

    do iatom = 1, ions%nions
      call matrix_vector(charge%car2lat, ions%r_car(iatom, 1:3) - charge%org_car(1:3), ions%r_lat(iatom, 1:3))
      ions%r_lat(iatom, 1:3) = ions%r_lat(iatom, 1:3) + charge%org_lat(1:3)
      ! folds the position into the cell
      call pbc_r_lat(ions%r_lat(iatom, :), charge%npts)
    end do

    charge%nrho = product(charge%npts(1:3))

    SAFE_ALLOCATE(charge%rho(1:charge%npts(1), 1:charge%npts(2), 1:charge%npts(3)))

    do ix = 1, charge%npts(1)
      do iy = 1, charge%npts(2)
        do iz = 1, charge%npts(3)
          charge%rho(ix, iy, iz) = vol*density_cube%drs(ix, iy, iz)
        end do
      end do
    end do

    ! distance between neighboring points
    do ix = -1, 1
      dlat(1) = real(ix, q2)
      do iy = -1, 1
        dlat(2) = real(iy, q2)
        do iz = -1, 1
          dlat(3) = real(iz, q2)

          call matrix_vector(charge%lat2car, dlat, dcar)

          charge%lat_dist(ix, iy, iz) = sqrt(sum(dcar*dcar))

          if ((ix == 0) .and. (iy == 0) .and. (iz == 0)) then
            charge%lat_i_dist(ix, iy, iz)=0.0_q2
          else
            charge%lat_i_dist(ix, iy, iz)=1.0_q2/charge%lat_dist(ix, iy, iz)
          end if

        end do
      end do
    end do

    if(present(bader_charges)) then
      
      call bader_calc(bdr, ions, charge, this%options)
      
      do iatom = 1, geo%natoms
        bader_charges(iatom) = bdr%ionchg(iatom)
      end do
    end if

    if(present(voronoi_charges)) then
      
      call voronoi(vor, ions, charge)
      
      do iatom = 1, geo%natoms
        voronoi_charges(iatom) = vor%vorchg(iatom)
      end do
    end if


    SAFE_DEALLOCATE_A(ions%r_car)
    SAFE_DEALLOCATE_A(ions%r_dir)
    SAFE_DEALLOCATE_A(ions%ion_chg)
    SAFE_DEALLOCATE_A(ions%atomic_num)
    SAFE_DEALLOCATE_A(ions%r_lat)
    SAFE_DEALLOCATE_A(charge%rho)

    call dcube_function_free_RS(cube, density_cube)
    call cube_end(cube)

    call profiling_out(prof)
    POP_SUB(partial_charges_calculate)

  end subroutine partial_charges_calculate

  ! ---------------------------------------------------------

  subroutine partial_charges_end(this)
    type(partial_charges_t), intent(inout) :: this

    PUSH_SUB(partial_charges_end)

    POP_SUB(partial_charges_end)
  end subroutine partial_charges_end

end module partial_charges_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
