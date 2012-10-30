!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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
!! $Id: ob_grid.F90 6247 2009-12-20 13:20:41Z xavier $

#include "global.h"

module ob_grid_m
  use datasets_m
  use geometry_m
  use global_m
  use ob_interface_m
  use parser_m
  use mesh_m
  use messages_m
  use simul_box_m
  use profiling_m

  implicit none

  private
  public ::                      &
    ob_grid_lead_t,              &
    ob_grid_t,                   &
    ob_grid_init,                &
    ob_grid_end,                 &
    ob_grid_write_info

  type ob_grid_lead_t
    type(simul_box_ob_info_t) :: info
    type(simul_box_t)         :: sb          !< The simulation box of the lead.
    type(mesh_t)              :: mesh        !< The mesh of the lead.
    type(interface_t)         :: intf        !< The interface (reduced unit cell) of the lead.
    character(len=2000)       :: td_bias     !< Td-potential of lead.
  end type ob_grid_lead_t

  type ob_grid_t
    logical           :: open_boundaries     !< Use open boundaries?
    logical           :: transport_mode      !< transport mode switched on (during open boundaries)
    type(ob_grid_lead_t), pointer :: lead(:) !< leads(1:NLEADS)
  end type ob_grid_t

contains

  ! ---------------------------------------------------------
  !> reads the open boundary information and allocates accordingly
  subroutine ob_grid_init(ob_grid)
    type(ob_grid_t),               intent(out) :: ob_grid

    integer            :: nr, tag, nrows, ncols
    type(block_t)      :: blk

    integer, parameter ::   &
      LEAD_DATASET     = 1, &
      LEAD_RESTART_DIR = 2, &
      LEAD_STATIC_DIR  = 3, &
      ADD_UNIT_CELLS   = 4, &
      TD_POT_FORMULA   = 5, &
      NUMBER_LEADS     = 6, &
      TRANSPORT_MODE   = 7, &
      transport_on     = 1, &
      transport_off    = 0

    integer :: il, ol, tmode

    PUSH_SUB(ob_grid_init)

    ob_grid%open_boundaries = .false.

    !%Variable OpenBoundaries
    !%Type block
    !%Section Mesh
    !%Description
    !% This feature is experimental.
    !% In transport mode it enables open boundaries in the <i>x</i>-direction
    !% and defines the character of the leads attached to the left and right
    !% of the finite central system.
    !% The more general situation (non-transport) is that a given number
    !% of leads (number_leads) are attached to the central region.
    !%
    !% The format is as follows:
    !%
    !% <pre>
    !% %OpenBoundaries
    !%  lead_dataset     | "dataset"   | "dataset"
    !%  lead_restart_dir | "directory" | "directory"
    !%  lead_static_dir  | "directory" | "directory"
    !%  add_unit_cells   | nl          | nr
    !%  td_pot_formula   | "formula"   | "formula"
    !%  transport_mode   | transport_on
    !%  number_leads     | 2
    !% %
    !% </pre>
    !%
    !% The left column specifies characteristics of the left lead and
    !% and the right column characteristics of the right lead analogously.
    !% If only one column is given, the value specified is used for both leads.
    !%
    !% All entries except <tt>lead_dataset</tt> are optional.
    !%
    !%Option lead_dataset 1
    !% Gives the name of the dataset used for the periodic calculation of the
    !% ground states of the leads. It is used, <i>e.g.</i>, to read in the coordinates of the
    !% atoms of the lead. Both entries for left and right have to be equal.
    !%Option lead_restart_dir 2
    !% <tt>lead_restart_dir</tt> gives the name of restart directory of the periodic
    !% ground-state calculation for the lead unit cell. Though
    !% one may give different datasets for the left and right lead, they have to
    !% be identical due to the algorithm used to obtain extended eigenstates.
    !% The default is <tt>&lt;lead_dataset&gt;restart</tt>.
    !%Option lead_static_dir 3
    !% The same as <tt>lead_restart_dir</tt> for the <tt>static</tt> directory.
    !% <tt>Octopus</tt> needs the Kohn-Sham potential of the leads. Therefore, the periodic
    !% run must include <tt>Output = potential</tt> in the input file. The default
    !% of this entry is <tt>&lt;lead_dataset&gt;static</tt>.
    !%Option add_unit_cells 4
    !% <tt>add_unit_cells</tt> specifies how many lead unit cells should
    !% be included in the computational domain. Suitable values highly depend
    !% on the system under study but the numbers <tt>nl</tt> and <tt>nr</tt> should
    !% be taken large enough for the potential to equilibrate because we assume
    !% instaneous metallic screening in the leads. Furthermore, note that in
    !% a ground-state calculation, one additional unit cell is added automatically
    !% for each lead to the computational domain because the propagation
    !% algorithm needs the knowledge of the initial state for the first unit cell
    !% outside the simulation box. If omitted, no unit cells are included in the
    !% simulation region (apart from the one which is automatically added in
    !% ground-state calculations).
    !%Option td_pot_formula 5
    !% Defines a spatially local time-dependent potential in the leads as an
    !% analytic expression. This describes the time-dependent bias in the leads.
    !%Option transport_mode 6
    !% If set to on (transport_on) the normal transport calculation is done,
    !% otherwise (transport_off) an open system without the source term.
    !% The initial state is to be assumed to be completely localized in
    !% the central region. Default is transport_on.
    !%Option number_leads 7
    !% In the non-transport mode it defines the number of leads connected
    !% to the central region.
    !%
    !%Option transport_on 1
    !% Use transport (default).
    !%Option transport_off 0
    !% Just use open boundaries.
    !%End
    NLEADS = 0
    if(parse_block(datasets_check('OpenBoundaries'), blk).eq.0) then

      call messages_experimental("Open boundaries")
      ob_grid%open_boundaries = .true.
      NLEADS = 2
      ! first of all get the number of leads
      nrows = parse_block_n(blk)
      do nr = 0, nrows - 1
        call parse_block_integer(blk, nr, 0, tag)

        if(tag.eq.NUMBER_LEADS) then
          call parse_block_integer(blk, nr, 2, NLEADS)
          if ((NLEADS < 0) .or. (NLEADS > 2 * MAX_DIM)) &
            call input_error('OpenBoundaries:number_leads')
        end if
      end do
      ! now we can allocate the leads
      SAFE_ALLOCATE(ob_grid%lead(1:NLEADS))

      ! fill with default values
      ob_grid%lead(:)%info%dataset        = ''
      ob_grid%lead(:)%info%restart_dir    = ''
      ob_grid%lead(:)%info%static_dir     = ''
      ob_grid%lead(:)%info%ucells         = 0
      ob_grid%lead(:)%td_bias             = '0'
      ob_grid%transport_mode      = .true.

      ! now read the rest of the block
      nrows = parse_block_n(blk)
      do nr = 0, nrows - 1
        call parse_block_integer(blk, nr, 0, tag)
        ncols = parse_block_cols(blk, nr)

        select case(tag)
        case(LEAD_DATASET)
          call parse_block_string(blk, nr, 1, ob_grid%lead(LEFT)%info%dataset)
          if(ncols .eq. 3) then
            call parse_block_string(blk, nr, 2, ob_grid%lead(RIGHT)%info%dataset)
            if(trim(ob_grid%lead(LEFT)%info%dataset) .ne. trim(ob_grid%lead(RIGHT)%info%dataset)) then
              message(1) = 'For now datasets for left and right lead unit cells must'
              message(2) = 'be equal, i.e. only symmetric leads are possible.'
              call messages_fatal(2)
            end if
          else
            do ol = 2, NLEADS
              ob_grid%lead(ol)%info%dataset = ob_grid%lead(LEFT)%info%dataset
            end do
          end if
        case(LEAD_RESTART_DIR)
          call parse_block_string(blk, nr, 1, ob_grid%lead(LEFT)%info%restart_dir)
          if(ncols .eq. 3) then
            call parse_block_string(blk, nr, 2, ob_grid%lead(RIGHT)%info%restart_dir)
            if(trim(ob_grid%lead(LEFT)%info%restart_dir) .ne. trim(ob_grid%lead(RIGHT)%info%restart_dir)) then
              message(1) = 'Restart directories for left and right lead'
              message(2) = 'unit cells must be equal, i.e. only symmetric'
              message(3) = 'leads are possible.'
              call messages_fatal(3)
            end if
          else
            do ol = 2, NLEADS
              ob_grid%lead(ol)%info%restart_dir = ob_grid%lead(LEFT)%info%restart_dir
            end do
          end if
        case(LEAD_STATIC_DIR)
          call parse_block_string(blk, nr, 1, ob_grid%lead(LEFT)%info%static_dir)
          if(ncols .eq. 3) then
            call parse_block_string(blk, nr, 2, ob_grid%lead(RIGHT)%info%static_dir)
            if(trim(ob_grid%lead(LEFT)%info%static_dir) .ne. trim(ob_grid%lead(RIGHT)%info%static_dir)) then
              message(1) = 'Static directories for left and right lead'
              message(2) = 'unit cells must be equal, i.e. only symmetric'
              message(3) = 'leads are possible.'
              call messages_fatal(3)
            end if
          else
            do ol = 2, NLEADS
              ob_grid%lead(ol)%info%static_dir = ob_grid%lead(LEFT)%info%static_dir
            end do
          end if
        case(ADD_UNIT_CELLS)
          call parse_block_integer(blk, nr, 1, ob_grid%lead(LEFT)%info%ucells)
          if(ncols .eq. 3) then
            call parse_block_integer(blk, nr, 2, ob_grid%lead(RIGHT)%info%ucells)
          else
            do ol = 2, NLEADS
              ob_grid%lead(ol)%info%ucells = ob_grid%lead(LEFT)%info%ucells
            end do
          end if
          if(any(ob_grid%lead(1:NLEADS)%info%ucells .lt. 0)) then
            message(1) = 'add_unit_cells in the OpenBoundaries block must not be negative.'
            call messages_fatal(1)
          end if
        case(TD_POT_FORMULA)
          call parse_block_string(blk, nr, 1, ob_grid%lead(LEFT)%td_bias)
          if(ncols .eq. 3) then
            call parse_block_string(blk, nr, 2, ob_grid%lead(RIGHT)%td_bias)
          else
            do ol = 2, NLEADS
              ob_grid%lead(ol)%td_bias = ob_grid%lead(LEFT)%td_bias
            end do
          end if
        case(TRANSPORT_MODE)
          call parse_block_integer(blk, nr, 1, tmode)
          select case(tmode)
            case(transport_on)
              ob_grid%transport_mode = .true.
            case(transport_off)
              ob_grid%transport_mode = .false.
            case default
          end select
        case default
        end select
      end do
      ! Check if necessary lead_dataset line has been provided.
      if(all(ob_grid%lead(1:NLEADS)%info%dataset .eq. '')) then
        call input_error('OpenBoundaries')
      end if
      ! Set default restart directory.
      if(all(ob_grid%lead(1:NLEADS)%info%restart_dir .eq. '')) then
        do il = 1, NLEADS
          ob_grid%lead(il)%info%restart_dir = trim(ob_grid%lead(il)%info%dataset) // 'restart'
        end do
      end if
      ! Set default static directory.
      if(all(ob_grid%lead(1:NLEADS)%info%static_dir .eq. '')) then
        do il = 1, NLEADS
          ob_grid%lead(il)%info%static_dir = trim(ob_grid%lead(il)%info%dataset) // 'static'
        end do
      end if
    else
      ob_grid%open_boundaries  = .false.
      ob_grid%transport_mode   = .false.
      nullify(ob_grid%lead)
    end if

    POP_SUB(ob_grid_init)
  end subroutine ob_grid_init

  ! ---------------------------------------------------------
  subroutine ob_grid_end(ob_grid)
    type(ob_grid_t), intent(inout) :: ob_grid

    PUSH_SUB(ob_grid_end)

    if(associated(ob_grid%lead)) then
      SAFE_DEALLOCATE_P(ob_grid%lead)
    end if

    POP_SUB(ob_grid_end)
  end subroutine ob_grid_end

  !--------------------------------------------------------------
  subroutine ob_grid_write_info(ob_grid, iunit)
    type(ob_grid_t), intent(inout) :: ob_grid
    integer,           intent(in)  :: iunit

    PUSH_SUB(ob_grid_write_info)

    if(ob_grid%open_boundaries) then
      write(message(1), '(a)')       'Open boundaries in x-direction:'
      write(message(2), '(a,2i4)') '  Additional unit cells left:    ', &
        ob_grid%lead(LEFT)%info%ucells
      write(message(3), '(a,2i4)') '  Additional unit cells right:   ', &
        ob_grid%lead(RIGHT)%info%ucells
      write(message(4), '(a)')     '  Left lead read from directory:  ' // &
        trim(ob_grid%lead(LEFT)%info%restart_dir)
      write(message(5), '(a)')     '  Right lead read from directory: ' // &
        trim(ob_grid%lead(RIGHT)%info%restart_dir)
      call messages_info(5, iunit)
    end if

    POP_SUB(ob_grid_write_info)
  end subroutine ob_grid_write_info


  !--------------------------------------------------------------
  subroutine ob_grid_dump(ob_grid, iunit)
    type(ob_grid_t), intent(inout) :: ob_grid
    integer,           intent(in)  :: iunit

    PUSH_SUB(ob_grid_dump)

    if(ob_grid%open_boundaries) then
      write(iunit, '(a20,2i4)') 'add_unit_cells=     ', ob_grid%lead(1:2)%info%ucells ! FIXME for NLEADS>2
      write(iunit, '(a20,a32)') 'lead_restart_dir(L)=', ob_grid%lead(LEFT)%info%restart_dir
      write(iunit, '(a20,a32)') 'lead_restart_dir(R)=', ob_grid%lead(RIGHT)%info%restart_dir
    end if

    POP_SUB(ob_grid_dump)
  end subroutine ob_grid_dump


end module ob_grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
