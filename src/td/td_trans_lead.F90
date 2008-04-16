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
!! $Id: td_transport.F90 3030 2007-06-25 16:45:05Z marques $

! Initialize leads.

#include "global.h"

module td_trans_lead_m
  use datasets_m
  use global_m
  use loct_parser_m
  use messages_m
  use string_m
  use td_trans_intf_m
  use varinfo_m

  implicit none

  private

  public ::       &
    trans_lead_t, &
    lead_init,    &
    lead_end

  integer, public, parameter :: &
    flat_lead = 1

  character(len=4), public, parameter, dimension(1) :: &
    lead_type_names = (/'flat'/)

  type trans_lead_t
    FLOAT, pointer      :: td_pot(:, :)           ! td_pot(t, il), the local external td potential.
    character(len=1024) :: td_pot_formula(NLEADS)
    integer             :: lead_type(NLEADS)      ! Right now: just flat leads.
  end type trans_lead_t

contains

  ! ---------------------------------------------------------
  ! Read lead specification from input file and calculate td potential.
  subroutine lead_init(leads, n_steps, tstep)
    type(trans_lead_t), intent(out) :: leads
    integer,            intent(in)  :: n_steps
    FLOAT,              intent(in)  :: tstep

    integer             :: nrows, ncols, il, it
    FLOAT               :: t
    FLOAT, allocatable  :: pot_im(:)
    character(len=1024) :: tmp_c_string
    type(block_t)       :: blk

    call push_sub('td_trans_lead.lead_init')

    ALLOCATE(leads%td_pot(0:n_steps+1, NLEADS), (n_steps+2)*NLEADS)
    ALLOCATE(pot_im(0:n_steps+1), n_steps+2)

    !%Variable TDTransLeads
    !%Type block
    !%Section Transport
    !%Description
    !% this block describes the semi-infinite periodic leads attached
    !% to the central device region. It may have on or two columns per row:
    !% on column means symmetric leads, if two are given the left one
    !% describes the left lead and the right column the right lead.
    !% But the blocks has a fixed number of rows which is two at the moment.
    !%
    !% Example:
    !%
    !% <tt>%TDTransLeads<br>
    !% &nbsp;flat_lead<br>
    !% &nbsp;"-1"      | "1"<br>
    !% %</tt>
    !%
    !% The first row describes the lead type. Only flat leads are possible
    !% right now.
    !% The second row gives the time dependent local external potential in
    !% the leads which has to be given as a formula string and may include
    !% t as parameter but not any of the spatial coordinates.
    !% In the example above the potential is time independent but asymmetric.
    !%
    !%Option flat_lead 1
    !% A flat lead is attached to the central device region
    !%
    !%End

    if(loct_parse_block(check_inp('TDTransLeads'), blk) == 0) then
      nrows = loct_parse_block_n(blk)

      if(nrows.ne.2) then
        message(1) = 'The "TDTransLeads" block must have exactly two rows.'
        message(2) = 'See the variable description for an explanation of the format.'
        call write_fatal(2)
      end if

      ! Lead type.
      ncols = loct_parse_block_cols(blk, 0)
      do il = 1, ncols
        call loct_parse_block_int(blk, 0, il-1, leads%lead_type(il))
      end do
      if(ncols.eq.1) then
        leads%lead_type(RIGHT) = leads%lead_type(LEFT)
      end if
      do il = 1, NLEADS
        if(leads%lead_type(il).ne.flat_lead) then
          message(1) = 'Only flat leads are possible for now. Put "flat_leads" in'
          message(2) = '"TDTransLeads" block.'
          call write_fatal(2)
        end if
      end do

      ! TD potential.
      ncols = loct_parse_block_cols(blk, 1)
      do il = 1, ncols
        call loct_parse_block_string(blk, 1, il-1, leads%td_pot_formula(il))
      end do
      if(ncols.eq.1) then
        leads%td_pot_formula(RIGHT) = leads%td_pot_formula(LEFT)
      end if
    else
      message(1) = 'Input: No TDTransLeads block found.'
      message(2) = 'Input: Defaulting to flat leads with zero potential.'
      call write_info(2)

      leads%lead_type      = flat_lead
      leads%td_pot_formula = '0'
    end if

    ! Calculate td potential.
    do il = 1, NLEADS
      tmp_c_string = leads%td_pot_formula(il)
      call conv_to_c_string(tmp_c_string)

      leads%td_pot(0, il) = M_ZERO
      do it = 1, n_steps+1
        t = it*tstep
        call loct_parse_expression(leads%td_pot(it, il), pot_im(it), &
          M_ZERO, M_ZERO, M_ZERO, M_ZERO, t, tmp_c_string)
      end do
    end do

    deallocate(pot_im)
    call pop_sub()
  end subroutine lead_init


  ! ---------------------------------------------------------
  ! Free memory of lead information.
  subroutine lead_end(leads)
    type(trans_lead_t), intent(inout) :: leads

    call push_sub('td_trans_lead.lead_end')

    if(associated(leads%td_pot)) then
      deallocate(leads%td_pot)
      nullify(leads%td_pot)
    end if

    call pop_sub()
  end subroutine lead_end
end module td_trans_lead_m
