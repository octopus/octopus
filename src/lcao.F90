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

module lcao_m
  use global_m
  use messages_m
  use datasets_m
  use profiling_m
  use lib_oct_parser_m
  use lib_oct_m
  use lib_oct_gsl_spline_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use functions_m
  use mesh_m
  use simul_box_m
  use specie_m
  use geometry_m
  use states_m
  use system_m
  use hamiltonian_m
  use grid_m

  use output_m

  implicit none

  private
  public ::          &
    lcao_t,          &
    lcao_init,       &
    lcao_wf,         &
    lcao_end

  type lcao_t
    integer           :: state ! 0 => non-initialized;
                               ! 1 => initialized (k, s and v1 matrices filled)
    type(states_t) :: st

    FLOAT,  pointer  :: dhamilt(:, :, :) ! hamilt stores the Hamiltonian in the LCAO subspace;
    FLOAT,  pointer  :: ds     (:, :, :) ! s is the overlap matrix;
    FLOAT,  pointer  :: dk     (:, :, :) ! k is the kinetic + spin orbit operator matrix;
    FLOAT,  pointer  :: dv     (:, :, :) ! v is the potential.

    CMPLX,  pointer  :: zhamilt(:, :, :) ! hamilt stores the Hamiltonian in the LCAO subspace;
    CMPLX,  pointer  :: zs     (:, :, :) ! s is the overlap matrix;
    CMPLX,  pointer  :: zk     (:, :, :) ! k is the kinetic + spin orbit operator matrix;
    CMPLX,  pointer  :: zv     (:, :, :) ! v is the potential.
  end type lcao_t

contains

  ! ---------------------------------------------------------
  subroutine lcao_init(gr, lcao_data, st, h)
    type(lcao_t),         intent(out)   :: lcao_data
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(in)    :: st
    type(hamiltonian_t),  intent(in)    :: h

    integer :: ia, norbs, n

    if(lcao_data%state == 1) return

    call profiling_in(C_PROFILING_LCAO_INIT)
    call push_sub('lcao.lcao_init')

    call states_null(lcao_data%st)

    ! Fix the dimension of the LCAO problem (lcao_data%dim)
    norbs = 0
    do ia = 1, gr%geo%natoms
      norbs = norbs + gr%geo%atom(ia)%spec%niwfs
    end do
    if( (st%d%ispin.eq.SPINORS) ) norbs = norbs * 2

    if(norbs < st%nst) then
      lcao_data%state = 0
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
    !% twice the number of required orbitlas for the full calculation. 
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

    lcao_data%st%nst = norbs
    lcao_data%st%st_start = 1
    lcao_data%st%st_end = norbs
    lcao_data%st%d%dim = st%d%dim
    lcao_data%st%d%nik = st%d%nik
    lcao_data%st%d%ispin = st%d%ispin
    lcao_data%st%d%wfs_type = st%d%wfs_type
    call states_allocate_wfns(lcao_data%st, gr%m)

    if (lcao_data%st%d%wfs_type == M_REAL) then
      call dlcao_init(lcao_data, gr, h, norbs)
    else
      call zlcao_init(lcao_data, gr, h, norbs)
    end if

    lcao_data%state = 1

    call pop_sub()
    call profiling_out(C_PROFILING_LCAO_INIT)
  end subroutine lcao_init

  ! ---------------------------------------------------------
  subroutine lcao_end(lcao_data, nst)
    type(lcao_t), intent(inout) :: lcao_data
    integer,         intent(in) :: nst

    call push_sub('lcao.lcao_end')

    if(lcao_data%st%nst >= nst) then
      if(associated(lcao_data%dhamilt)) deallocate(lcao_data%dhamilt)
      if(associated(lcao_data%ds     )) deallocate(lcao_data%ds)
      if(associated(lcao_data%dk     )) deallocate(lcao_data%dk)
      if(associated(lcao_data%dv     )) deallocate(lcao_data%dv)

      if(associated(lcao_data%zhamilt)) deallocate(lcao_data%zhamilt)
      if(associated(lcao_data%zs     )) deallocate(lcao_data%zs)
      if(associated(lcao_data%zk     )) deallocate(lcao_data%zk)
      if(associated(lcao_data%zv     )) deallocate(lcao_data%zv)
    endif

    if(associated(lcao_data%st%dpsi)) then
      deallocate(lcao_data%st%dpsi); nullify(lcao_data%st%dpsi)
    end if
    if(associated(lcao_data%st%zpsi)) then
      deallocate(lcao_data%st%zpsi); nullify(lcao_data%st%zpsi)
    end if

    lcao_data%state = 0
    call pop_sub()
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(lcao_data, st, m, h, start)
    type(lcao_t),        intent(inout) :: lcao_data
    type(states_t),      intent(inout) :: st
    type(mesh_t),        intent(in)    :: m
    type(hamiltonian_t), intent(in)    :: h
    integer, optional,   intent(in)    :: start

    integer :: start_

    ASSERT(lcao_data%state == 1)

    call profiling_in(C_PROFILING_LCAO)
    call push_sub('lcao.lcao_wf')

    start_ = 1
    if(present(start)) start_ = start

    if (lcao_data%st%d%wfs_type == M_REAL) then
      call dlcao_wf(lcao_data, st, m, h, start_)
    else
      call zlcao_wf(lcao_data, st, m, h, start_)
    end if

    call pop_sub()
    call profiling_out(C_PROFILING_LCAO)
  end subroutine lcao_wf


#include "undef.F90"
#include "real.F90"
#include "lcao_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lcao_inc.F90"


end module lcao_m
