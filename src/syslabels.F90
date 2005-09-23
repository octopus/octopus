!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module syslabels
  use messages
  use lib_oct_parser

  implicit none

  public :: syslabels_init, syslabels_end
  public :: check_inp, read_system_labels

  ! variables to treat multi subsytems
  character(len=32), allocatable :: subsys_label(:)
  integer,           allocatable :: subsys_runmode(:), subsys_run_order(:)
  integer,           allocatable :: no_of_states(:)
  integer, parameter             :: multi_subsys_mode = 1000
  integer                        :: current_subsystem = 1
  integer                        :: no_subsystems
  integer                        :: no_syslabels, no_subsys_runmodes
  character(len=32)              :: current_label, tmpdir


contains


  ! ---------------------------------------------------------
  subroutine syslabels_init(calc_mode)
    integer, intent(in) :: calc_mode

    no_syslabels  = 1
    no_subsystems = 1
    allocate(subsys_label(no_syslabels), subsys_runmode(no_syslabels))
    allocate(subsys_run_order(no_syslabels))
    allocate(no_of_states(no_syslabels))
    current_subsystem = 1
    subsys_label(current_subsystem) = ""
    subsys_runmode(current_subsystem) = calc_mode
    subsys_run_order(current_subsystem) = 1

    ! An initial assignment to avoid having to do it in the subprograms.
    current_label = trim(subsys_label(subsys_run_order(1)))
    current_subsystem = subsys_run_order(1)

  end subroutine syslabels_init


  ! ---------------------------------------------------------
  subroutine syslabels_end()

    deallocate(subsys_label, subsys_runmode, subsys_run_order)
    deallocate(no_of_states)

  end subroutine syslabels_end


  ! ---------------------------------------------------------
  character(len=64) function check_inp(variable) result(var_name)
    character(len = * ), intent(in)  :: variable
    character(len = 64)              :: composite_name

    composite_name = trim(subsys_label(current_subsystem))//trim(variable)

    if(loct_parse_isdef(composite_name).ne.0) then
       ! composite name has been defined in the input file
       var_name = composite_name
    else
       ! could not find composite name in the input;
       ! will use bare variable name
       var_name = variable
    endif

  end function check_inp


  ! ---------------------------------------------------------
  subroutine read_system_labels()
    integer               :: i
    integer(POINTER_SIZE) :: blk

    ! first we read the required information from the input file
    ! and prompt the user for possible errors in the input

    ! find out how many subsystem we want to treat ...
    if(loct_parse_block('SystemLabels', blk) == 0) then
       no_syslabels = loct_parse_block_cols(blk,0)
    else
       message(1) = "Could not find block SystemLabels in the input file."
       message(2) = "This block is mandatory for run mode MultiSubsystem."
       call write_fatal(2)
    endif

    allocate(subsys_label(no_syslabels))

    ! ... and how the user would like to call them.
    do i = 1, no_syslabels
       call loct_parse_block_string(blk, 0, i-1, subsys_label(i))
    enddo
    call loct_parse_block_end(blk)

    ! now we check what we have to run in the respective subsystem
    if(loct_parse_block('SystemRunModes', blk) == 0) then
       no_subsys_runmodes = loct_parse_block_cols(blk,0)
    else
       message(1) = "Could not find block SystemRunModes in the input file."
       message(2) = "This block is mandatory for run mode MultiSubsystem."
       call write_fatal(2)
    endif

    if(no_subsys_runmodes/=no_syslabels) then
       message(1) = "The blocks SystemLabels and SystemRunModes do not have"
       message(2) = "the same size. Please correct your input file."
       call write_fatal(2)
    endif

    allocate(subsys_runmode(no_subsys_runmodes))

    do i = 1, no_subsys_runmodes
       call loct_parse_block_int(blk, 0, i-1, subsys_runmode(i))
    enddo
    call loct_parse_block_end(blk)

    ! ... and in which order
    if(loct_parse_block('SystemRunOrder', blk) == 0) then
       no_subsystems = loct_parse_block_cols(blk,0)
    else
       message(1) = "Could not find block SystemRunOrder in the input file."
       message(2) = "This block is mandatory for run mode MultiSubsystem."
       call write_fatal(2)
    endif

    allocate(subsys_run_order(no_subsystems))
    allocate(no_of_states(no_subsystems))

    do i = 1, no_subsystems
       call loct_parse_block_int(blk, 0, i-1, subsys_run_order(i))
       if (subsys_run_order(i).gt.no_subsys_runmodes) then
          message(1) = "Subsystem number too large. Please correct the block"
          message(2) = "SystemRunOrder in the input file."
          call write_fatal(2)
       endif
    enddo
    call loct_parse_block_end(blk)

  end subroutine read_system_labels

end module syslabels
