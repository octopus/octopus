!! Copyright (C) 2018 N. Tancogne-Dejean
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

module constrain_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hirshfeld_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::             &
    constrain_guess_moment

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 77

contains

  ! ---------------------------------------------------------
  !> transform the density into a density that matches the constrain
  !> Instead of atomic densities used in LCAO, we use the atomic density in the sens of
  !> Hirshfeld
  subroutine constrain_guess_moment(st, namespace, gr, sb, geo)
    type(states_elec_t), intent(in)    :: st
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in)    :: gr
    type(simul_box_t),   intent(in)    :: sb
    type(geometry_t),    intent(in)    :: geo

    integer :: ia, is, idir, gmd_opt
    integer, save :: iseed = 321
    type(block_t) :: blk
    FLOAT :: rr, rnd, phi, theta, mag(1:3), lmag, n1, n2
    FLOAT, allocatable :: atom_rho(:,:), rho(:,:)
    logical :: parallelized_in_atoms
    type(hirshfeld_t) :: hirshfeld
    FLOAT :: charge


    PUSH_SUB(constrain_guess_moment)

    parallelized_in_atoms = .false.

    call hirshfeld_init(hirshfeld, namespace, gr%mesh, geo, st)

    if (st%d%spin_channels == 1) then
      gmd_opt = INITRHO_PARAMAGNETIC
    else
      call parse_variable(namespace, 'GuessMagnetDensity', INITRHO_FERROMAGNETIC, gmd_opt)
      if(gmd_opt == INITRHO_RANDOM) then
        message(1) = "Constrain with random magnetization makes no sens."
        call messages_fatal(1, namespace=namespace)
      end if
    end if

    SAFE_ALLOCATE(rho(1:gr%fine%mesh%np, 1:st%d%nspin))
    rho = M_ZERO
    select case (gmd_opt)
    case (INITRHO_PARAMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:st%d%spin_channels))

      parallelized_in_atoms = .true.

      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call hirshfeld_charge(hirshfeld, namespace, ia, st%rho, charge, atom_rho)

        rho(1:gr%fine%mesh%np, 1:st%d%spin_channels) = rho(1:gr%fine%mesh%np, 1:st%d%spin_channels) + &
                                                  atom_rho(1:gr%fine%mesh%np, 1:st%d%spin_channels)
      end do

      if (st%d%spin_channels == 2) then
        rho(1:gr%fine%mesh%np, 1) = M_HALF*(sum(rho(1:gr%fine%mesh%np, 1:2), dim=2))
        rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 1)
        if(st%d%nspin > st%d%spin_channels) then
          rho(1:gr%fine%mesh%np, 3) = M_ZERO
          rho(1:gr%fine%mesh%np, 4) = M_ZERO
        end if
      end if

    case (INITRHO_FERROMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:2))

      parallelized_in_atoms = .true.

      atom_rho = M_ZERO
      rho = M_ZERO
      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call hirshfeld_charge(hirshfeld, namespace, ia, st%rho, charge, atom_rho)
        rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
      end do
      if(st%d%nspin > st%d%spin_channels) then
        rho(1:gr%fine%mesh%np, 3) = M_ZERO
        rho(1:gr%fine%mesh%np, 4) = M_ZERO
      end if

    case (INITRHO_USERDEF) ! User-defined
      
      if(parse_block(namespace, 'AtomsMagnetDirection', blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined."
        call messages_fatal(1, namespace=namespace)
      end if

      if (parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows."
        call messages_fatal(1, namespace=namespace)
      end if

      SAFE_ALLOCATE(atom_rho(1:gr%fine%mesh%np, 1:2))
      do ia = 1, geo%natoms
        !Read from AtomsMagnetDirection block 
        if (st%d%nspin == 2) then
          call parse_block_float(blk, ia-1, 0, mag(1))
          lmag = abs(mag(1))
        elseif (st%d%nspin == 4) then
          do idir = 1, 3
            call parse_block_float(blk, ia-1, idir-1, mag(idir))
            if (abs(mag(idir)) < CNST(1.0e-20)) mag(idir) = M_ZERO
          end do
          lmag = sqrt(dot_product(mag(1:3), mag(1:3)))
        end if

        !Get atomic density
        call hirshfeld_charge(hirshfeld, namespace, ia, st%rho, charge, atom_rho)

        !Scale magnetization density
        n1 = dmf_integrate(gr%fine%mesh, atom_rho(:, 1))
        n2 = dmf_integrate(gr%fine%mesh, atom_rho(:, 2))
        if (lmag > n1 + n2) then
          mag = mag*(n1 + n2)/lmag
          lmag = n1 + n2
        elseif (abs(lmag) <= M_EPSILON) then
          if (abs(n1 - n2) <= M_EPSILON) then
            rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
          else
            atom_rho(:, 1) = (atom_rho(:, 1) + atom_rho(:, 2))/M_TWO
            rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 1)
            rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + atom_rho(1:gr%fine%mesh%np, 1)
          end if
          cycle
        end if
        if (n1 - n2 /= lmag .and. n2 /= M_ZERO) then
          if (n1 - n2 < lmag) then
            atom_rho(:, 1) = atom_rho(:, 1) + (lmag - n1 + n2)/M_TWO/n2*atom_rho(:, 2)
            atom_rho(:, 2) = (n1 + n2 - lmag)/M_TWO/n2*atom_rho(:, 2)
          elseif (n1 - n2 > lmag) then
            atom_rho(:, 2) = atom_rho(:, 2) + (n1 - n2 - lmag)/M_TWO/n1*atom_rho(:, 1)
            atom_rho(:, 1) = (lmag + n1 + n2)/M_TWO/n1*atom_rho(:, 1)
          end if
        end if

        !Rotate magnetization density
        if (st%d%nspin == 2) then
          if (mag(1) > M_ZERO) then
            rho(1:gr%fine%mesh%np, 1:2) = rho(1:gr%fine%mesh%np, 1:2) + atom_rho(1:gr%fine%mesh%np, 1:2)
          else
            rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + atom_rho(1:gr%fine%mesh%np, 2)
            rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + atom_rho(1:gr%fine%mesh%np, 1)
          end if

        elseif (st%d%nspin == 4) then
          theta = acos(mag(3)/lmag)
          if (abs(mag(1)) <= M_EPSILON) then
            if (abs(mag(2)) <= M_EPSILON) then
              phi = M_ZERO
            elseif (mag(2) < M_ZERO) then
              phi = M_PI*CNST(3.0/2.0)
            elseif (mag(2) > M_ZERO) then
              phi = M_PI*M_HALF
            end if
          else
            if (mag(2) < M_ZERO) then
              phi = M_TWO*M_PI - acos(mag(1)/sin(theta)/lmag)
            elseif (mag(2) >= M_ZERO) then
              phi = acos(mag(1)/sin(theta)/lmag)
            end if
          end if

          rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 2) = rho(1:gr%fine%mesh%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 1) &
               + cos(theta/M_TWO)**2*atom_rho(1:gr%fine%mesh%np, 2)
          rho(1:gr%fine%mesh%np, 3) = rho(1:gr%fine%mesh%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
          rho(1:gr%fine%mesh%np, 4) = rho(1:gr%fine%mesh%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:gr%fine%mesh%np, 1) - atom_rho(1:gr%fine%mesh%np, 2))
        end if
      end do

      call parse_block_end(blk)

    end select


#ifdef HAVE_MPI
    if(geo%atoms_dist%parallel .and. parallelized_in_atoms) then
      ! NOTE: if random or user_defined are made parallelized in atoms, below should be st%d%nspin instead of spin_channels
      do is = 1, st%d%spin_channels
        atom_rho(1:gr%fine%mesh%np, 1) = rho(1:gr%fine%mesh%np, is)
        call MPI_Allreduce(atom_rho(1, 1), rho(1, is), gr%fine%mesh%np, &
          MPI_FLOAT, MPI_SUM, geo%atoms_dist%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    do is = 1, st%d%nspin
      call lalg_copy(gr%fine%mesh%np, rho(:,is), st%rho(:,is))
    end do

    SAFE_DEALLOCATE_A(atom_rho)

    call hirshfeld_end(hirshfeld)        

    POP_SUB(constrain_guess_moment)
  end subroutine constrain_guess_moment

end module constrain_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
