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
!! $Id: simul_box.F90 6396 2010-03-26 08:51:58Z mjv500 $

#include "global.h"

module ob_simul_box_m
  use c_pointer_m
  use calc_mode_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use loct_m
  use lookup_m
  use math_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use species_m
  use string_m
  use simul_box_m
  use symmetries_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    ob_simul_box_init

contains

  !--------------------------------------------------------------
  subroutine ob_simul_box_init(sb, lead_sb, ucells, dir, geo)
    type(simul_box_t), intent(inout) :: sb
    type(simul_box_t), intent(inout) :: lead_sb(:)
    integer,           intent(in)    :: ucells(:)
    character(len=*),  intent(in)    :: dir(:)
    type(geometry_t),  intent(inout) :: geo

    ! some local stuff
    FLOAT :: def_h, def_rsize
    integer :: idir, il

    PUSH_SUB(ob_simul_box_init)

    ! Open boundaries are only possible for rectangular simulation boxes.
    if(sb%box_shape.ne.PARALLELEPIPED) then
      message(1) = 'Open boundaries are only possible with a parallelepiped'
      message(2) = 'simulation box.'
      call write_fatal(2)
    end if
    ! Simulation box must not be periodic in transport direction.
    if(sb%periodic_dim .eq. 1) then
      message(1) = 'When using open boundaries, you cannot use periodic boundary'
      message(2) = 'conditions in the x-direction.'
      call write_fatal(2)
    end if

    call ob_read_lead_unit_cells(sb, lead_sb, dir)
    ! Adjust the size of the simulation box by adding the proper number
    ! of unit cells to the simulation region.
    do il = LEFT, RIGHT
      sb%lsize(TRANS_DIR) = sb%lsize(TRANS_DIR) + ucells(il)*lead_sb(il)%lsize(TRANS_DIR)
    end do
    call simul_box_add_lead_atoms(sb, geo) ! Add the atoms of the lead unit cells that are
                                           ! included in the simulation box to geo.
    POP_SUB(ob_simul_box_init)

  end subroutine ob_simul_box_init


  !--------------------------------------------------------------
  ! Read the simulation boxes of the leads
  subroutine ob_read_lead_unit_cells(sb, lead_sb, dir)
    type(simul_box_t), intent(inout) :: sb
    type(simul_box_t), intent(inout) :: lead_sb(:)
    character(len=*),  intent(in)    :: dir(:)

    integer :: iunit, il

    PUSH_SUB(ob_read_lead_unit_cells)

    do il = 1, NLEADS
      iunit = io_open(trim(dir(il))//'/'//GS_DIR//'mesh', action = 'read', is_tmp = .true., grp = mpi_world)
      call simul_box_init_from_file(lead_sb(il), iunit)
      call io_close(iunit)

      ! Check whether
      ! * simulation box is a parallelepiped,
      ! * the extensions in y-, z-directions fit the central box,
      ! * the central simulation box x-length is an integer multiple of
      !   the unit cell x-length,
      ! * periodic in one dimension, and
      ! * of the same dimensionality as the central system.

      if(lead_sb(il)%box_shape .ne. PARALLELEPIPED) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' lead is not a parallelepiped.'
        call write_fatal(1)
      end if

      if(any(sb%lsize(2:sb%dim) .ne. lead_sb(il)%lsize(2:sb%dim))) then
        message(1) = 'The size in non-transport-directions of the ' // LEAD_NAME(il) // ' lead'
        message(2) = 'does not fit the size of the non-transport-directions of the central system.'
        call write_fatal(2)
      end if

      if(.not. is_integer_multiple(sb%lsize(1), lead_sb(il)%lsize(1))) then
        message(1) = 'The length in x-direction of the central simulation'
        message(2) = 'box is not an integer multiple of the x-length of'
        message(3) = 'the ' // trim(LEAD_NAME(il)) // ' lead.'
        call write_fatal(3)
      end if

      if(lead_sb(il)%periodic_dim .ne. 1) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' lead is not periodic in x-direction.'
  !      call write_fatal(1)
        call write_warning(1)
      end if
      if(lead_sb(il)%dim .ne. sb%dim) then
        message(1) = 'Simulation box of ' // LEAD_NAME(il) // ' has a different dimension than'
        message(2) = 'the central system.'
        call write_fatal(2)
      end if
    end do

    POP_SUB(ob_read_lead_unit_cells)
  end subroutine ob_read_lead_unit_cells

  !--------------------------------------------------------------
  ! Read the coordinates of the leads atoms and add them to the
  ! simulation box
  subroutine ob_simul_box_add_lead_atoms(sb, lead_sb, geo)
    type(simul_box_t), intent(inout) :: sb
    type(simul_box_t), intent(inout) :: lead_sb(:)
    type(geometry_t),  intent(inout) :: geo

    type(geometry_t)  :: central_geo
    type(geometry_t), allocatable  :: lead_geo(:)
    character(len=32) :: label_bak
    integer           :: il, icell, iatom, jatom, icatom, dir

    PUSH_SUB(ob_simul_box_add_lead_atoms)

    SAFE_ALLOCATE(lead_geo(1:NLEADS))
    do il = 1, NLEADS
      ! We temporarily change the current label to read the
      ! coordinates of another dataset, namely the lead dataset.
      label_bak     = current_label
      current_label = sb%lead_dataset(il)
      call geometry_init(lead_geo(il), print_info=.false.)
      current_label = label_bak
      call simul_box_atoms_in_box(sb%sb_lead_unit_cell(il), lead_geo(il), .true.)
    end do

    ! Merge the geometries of the lead and of the central region.
    call geometry_copy(central_geo, geo)

    ! Set the number of atoms and classical atoms to the number
    ! of atoms coming from left and right lead and central part.
    if(geo%natoms .gt. 0) then
      SAFE_DEALLOCATE_P(geo%atom)
    end if
    geo%natoms = central_geo%natoms +                 &
      sb%add_unit_cells(LEFT) * lead_geo(LEFT)%natoms + &
      sb%add_unit_cells(RIGHT) * lead_geo(RIGHT)%natoms
    SAFE_ALLOCATE(geo%atom(1:geo%natoms))
    if(geo%ncatoms .gt. 0) then
      SAFE_DEALLOCATE_P(geo%catom)
    end if
    geo%ncatoms = central_geo%ncatoms +                &
      sb%add_unit_cells(LEFT) * lead_geo(LEFT)%ncatoms + &
      sb%add_unit_cells(RIGHT) * lead_geo(RIGHT)%ncatoms
    SAFE_ALLOCATE(geo%catom(1:geo%ncatoms))

    geo%only_user_def = central_geo%only_user_def .and. all(lead_geo(:)%only_user_def)
    geo%nlpp          = central_geo%nlpp .or. any(lead_geo(:)%nlpp)
    geo%nlcc          = central_geo%nlcc .or. any(lead_geo(:)%nlcc)
    geo%atoms%start   = 1
    geo%atoms%end     = geo%natoms
    geo%atoms%nlocal  = geo%natoms

    ! 1. Put the atoms of the central region into geo.
    geo%atom(1:central_geo%natoms)   = central_geo%atom
    geo%catom(1:central_geo%ncatoms) = central_geo%catom

    ! 2. Put the atoms of the leads into geo and adjust their x-coordinates.
    iatom  = central_geo%natoms + 1
    icatom = central_geo%ncatoms + 1

    do il = 1, NLEADS
      dir = (-1)**il
      ! We start from the "outer" unit cells of the lead.
      do icell = 1, sb%add_unit_cells(il)
        do jatom = 1, lead_geo(il)%natoms
          geo%atom(iatom) = lead_geo(il)%atom(jatom)
          geo%atom(iatom)%x(TRANS_DIR) = geo%atom(iatom)%x(TRANS_DIR) + &
            dir * (sb%lsize(TRANS_DIR) - (2*(icell - 1) + 1) * sb%sb_lead_unit_cell(il)%lsize(TRANS_DIR))
          iatom = iatom + 1
        end do

        do jatom = 1, lead_geo(il)%ncatoms
          geo%catom(icatom) = lead_geo(il)%catom(jatom)
          geo%catom(icatom)%x(TRANS_DIR) = geo%catom(icatom)%x(TRANS_DIR) + &
            dir * (sb%lsize(TRANS_DIR) - (2 * (icell - 1) + 1) * sb%sb_lead_unit_cell(il)%lsize(TRANS_DIR))
        end do
      end do
    end do

    ! Initialize the species of the "extended" central system.
    if(geo%nspecies.gt.0) then
      SAFE_DEALLOCATE_P(geo%species)
    end if
    call geometry_init_species(geo, print_info=.false.)

    do il = 1, NLEADS
      call geometry_end(lead_geo(il))
    end do

    call geometry_end(central_geo)
    SAFE_DEALLOCATE_A(lead_geo)

    POP_SUB(ob_simul_box_add_lead_atoms)
  end subroutine ob_simul_box_add_lead_atoms


end module ob_simul_box_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
