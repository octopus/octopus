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
!! $Id$

#include "global.h"

module lcao_m
  use batch_m
  use datasets_m
  use distributed_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use simul_box_m
  use solids_m
  use species_m
  use states_m
  use states_dim_m
  use states_lalg_m
  use states_block_m
  use h_sys_output_m

  implicit none

  private
  public ::            &
    lcao_t,            &
    lcao_init,         &
    lcao_wf,           &
    lcao_end,          &
    lcao_is_available, &
    lcao_num_orbitals

  integer, public, parameter ::     &
    LCAO_START_NONE    = 0, &
    LCAO_START_STATES  = 2, &
    LCAO_START_FULL    = 3

  type lcao_t
    private
    integer           :: state ! 0 => non-initialized;
                               ! 1 => initialized (k, s and v1 matrices filled)
    type(states_t)    :: st
    FLOAT,  pointer   :: ds     (:, :, :) ! s is the overlap matrix;
    CMPLX,  pointer   :: zs     (:, :, :) ! s is the overlap matrix;
  end type lcao_t

contains

  ! ---------------------------------------------------------
  subroutine lcao_init(this, gr, geo, st)
    type(lcao_t),         intent(out)   :: this
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(in)    :: st

    integer :: ia, norbs, n

    call profiling_in(C_PROFILING_LCAO_INIT)
    call push_sub('lcao.lcao_init')

    call states_null(this%st)

    !this is to avoid a bug whe deallocating in gfortran 4.2.0 20060520
    nullify(this%ds)
    nullify(this%zs)
        
    ! Fix the dimension of the LCAO problem (this%dim)
    norbs = 0
    do ia = 1, geo%natoms
      norbs = norbs + geo%atom(ia)%spec%niwfs
    end do
    if( (st%d%ispin.eq.SPINORS) ) norbs = norbs * 2

    if(norbs < st%nst) then
      this%state = 0
      write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
      call write_warning(1)
      call pop_sub()
      return
    end if

    !%Variable LCAODimension
    !%Type integer
    !%Default 0
    !%Section SCF
    !%Description
    !% Before starting the SCF cycle, an initial LCAO calculation can be performed
    !% in order to obtain reasonable initial guesses for spin-orbitals and densities.
    !% For this purpose, the code calculates a number of atomic orbitals -- this
    !% number depends on the given species. The default dimension for the LCAO basis
    !% set will be the sum of all these numbers, unless this dimension is larger than
    !% twice the number of required orbitals for the full calculation. 
    !%
    !% This dimension however can be reduced (never increased) by making use of the 
    !% variable LCAODimension. Note that LCAODimension cannot be smaller than the 
    !% number of orbitals needed in the full calculation -- if LCAODimension is smaller, 
    !% it will be changed silently increased to meet this requirement. In the same way, 
    !% if LCAODimension is larger than the available number of atomic orbitals, 
    !% it will be reduced. If you want to use the largest possible number, set
    !% LCAODimension to a negative number.
    !%End
    call loct_parse_int(check_inp('LCAODimension'), 0, n)
    if((n > 0) .and. (n <= st%nst)) then
      norbs = st%nst
    elseif( (n > st%nst) .and. (n < norbs) ) then
      norbs = n
    elseif( n.eq.0) then
      norbs = min(norbs, 2*st%nst)
    end if

    this%st%nst = norbs
    this%st%st_start = 1
    this%st%st_end = norbs
    this%st%d%dim = st%d%dim
    this%st%d%nik = st%d%nik
    this%st%d%ispin = st%d%ispin
    this%st%d%block_size = st%d%block_size
    this%st%parallel_in_states = .false.
    call distributed_copy(st%d%kpt, this%st%d%kpt)
    call states_allocate_wfns(this%st, gr%m, st%wfs_type)

    if (wfs_are_real(this%st)) then
      call dlcao_init(this, gr, geo, st, norbs)
    else
      call zlcao_init(this, gr, geo, st, norbs)
    end if

    this%state = 1

    call pop_sub()
    call profiling_out(C_PROFILING_LCAO_INIT)
  end subroutine lcao_init

  ! ---------------------------------------------------------
  subroutine lcao_end(this, nst)
    type(lcao_t), intent(inout) :: this
    integer,         intent(in) :: nst

    integer :: wfs_type

    call push_sub('lcao.lcao_end')

    call distributed_end(this%st%d%kpt)

    wfs_type = this%st%wfs_type

    if(this%st%nst >= nst) then
      if(wfs_type == M_REAL) then
        if(associated(this%ds     )) deallocate(this%ds)
      else
        if(associated(this%zs     )) deallocate(this%zs)
      end if
    endif

    if( associated(this%st%dpsi) .and. (wfs_type == M_REAL) ) then
      deallocate(this%st%dpsi)
      nullify(this%st%dpsi)
    end if
    if( associated(this%st%zpsi) .and. (wfs_type == M_CMPLX) ) then
      deallocate(this%st%zpsi)
      nullify(this%st%zpsi)
    end if

    this%state = 0
    call pop_sub()
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(this, st, gr, geo, h, start)
    type(lcao_t),        intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(hamiltonian_t), intent(in)    :: h
    integer, optional,   intent(in)    :: start

    integer :: start_

    ASSERT(this%state == 1)

    call profiling_in(C_PROFILING_LCAO)
    call push_sub('lcao.lcao_wf')

    start_ = 1
    if(present(start)) start_ = start

    if (wfs_are_real(this%st)) then
      call dlcao_wf(this, st, gr, geo, h, start_)
    else
      call zlcao_wf(this, st, gr, geo, h, start_)
    end if

    call pop_sub()
    call profiling_out(C_PROFILING_LCAO)
  end subroutine lcao_wf

  logical function lcao_is_available(this) result(available)
    type(lcao_t),        intent(in) :: this
    
    available = this%state == 1
  end function lcao_is_available

  integer function lcao_num_orbitals(this) result(norbs)
    type(lcao_t),        intent(in) :: this

    norbs = this%st%nst
  end function lcao_num_orbitals

#include "undef.F90"
#include "real.F90"
#include "lcao_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lcao_inc.F90"


end module lcao_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
