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
  use loct_m
  use loct_parser_m
  use mpi_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use solids_m
  use species_m
  use species_pot_m
  use states_m
  use states_dim_m
  use states_calc_m
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
    integer           :: norbs
    integer           :: maxorbs
    integer, pointer  :: atom(:)
    integer, pointer  :: level(:)
    integer, pointer  :: ddim(:)
  end type lcao_t

contains
  ! ---------------------------------------------------------
  subroutine lcao_init(this, gr, geo, st)
    type(lcao_t),         intent(out)   :: this
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(in)    :: st

    integer :: ia, n, ii, jj, maxj, idim

    call push_sub('lcao.lcao_init')

    ! nullify everything so we can check for associated pointers when deallocating
    nullify(this%atom)
    nullify(this%level)
    nullify(this%ddim)

    ! count the number of orbitals available
    maxj = 0
    this%maxorbs = 0
    do ia = 1, geo%natoms
      maxj = max(maxj, geo%atom(ia)%spec%niwfs)
      this%maxorbs = this%maxorbs + geo%atom(ia)%spec%niwfs
    end do

    this%maxorbs = this%maxorbs*st%d%dim

    if(this%maxorbs < st%nst) then
      this%state = 0
      write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
      call write_warning(1)
      call pop_sub(); return
    end if

    ! generate tables to know which indexes each atomic orbital has

    SAFE_ALLOCATE( this%atom(1:this%maxorbs))
    SAFE_ALLOCATE(this%level(1:this%maxorbs))
    SAFE_ALLOCATE( this%ddim(1:this%maxorbs))

    ! Each atom provides niwfs pseudo-orbitals (this number is given in
    ! geo%atom(ia)%spec%niwfs for atom number ia). This number is
    ! actually multiplied by two in case of spin-unrestricted or spinors
    ! calculations.
    !
    ! The pseudo-orbitals are placed in order in the following way (Natoms
    ! is the total number of atoms).
    !
    ! n = 1 => first orbital of atom 1,
    ! n = 2 => first orbital of atom 2.
    ! n = 3 => first orbital of atom 3.
    ! ....
    ! n = Natoms => first orbital of atom Natoms
    ! n = Natoms + 1 = > second orbital of atom 1
    ! ....
    !
    ! If at some point in this loop an atom pseudo cannot provide the corresponding
    ! orbital (because the niws orbitals have been exhausted), it moves on to the following
    ! atom.
    !
    ! In the spinors case, it changes a bit:
    !
    ! n = 1 => first spin-up orbital of atom 1, assigned to the spin-up component of the spinor.
    ! n = 2 => first spin-down orbital of atom 1, assigned to the spin-down component of the spinor.
    ! n = 3 => first spin-up orbital of atom 2, assigned to the spin-up component of the spinor.
    
    ii = 1
    do jj = 1, maxj
      do ia = 1, geo%natoms
        do idim = 1,st%d%dim
          if(jj > geo%atom(ia)%spec%niwfs) cycle

          this%atom(ii) = ia
          this%level(ii) = jj
          this%ddim(ii) = idim

          ii = ii + 1
        end do
      end do
    end do

    ASSERT(ii - 1 == this%maxorbs)

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
    !% This dimension however can be changed by making use of this
    !% variable. Note that LCAODimension cannot be smaller than the
    !% number of orbitals needed in the full calculation -- if
    !% LCAODimension is smaller, it will be silently increased to meet
    !% this requirement. In the same way, if LCAODimension is larger
    !% than the available number of atomic orbitals, it will be
    !% reduced. If you want to use the largest possible number, set
    !% LCAODimension to a negative number.
    !%End
    call loct_parse_int(datasets_check('LCAODimension'), 0, n)

    if(n > 0 .and. n <= st%nst) then
      this%norbs = st%nst
    else if(n > st%nst .and. n <= this%maxorbs) then
      this%norbs = n
    else if(n == 0) then
      this%norbs = min(this%maxorbs, 2*st%nst)
    else
      this%norbs = this%maxorbs
    end if
   
    ASSERT(this%norbs >= st%nst)
    ASSERT(this%norbs <= this%maxorbs)

    this%state = 1
       
    call pop_sub()
  end subroutine lcao_init

  ! ---------------------------------------------------------
  subroutine lcao_end(this)
    type(lcao_t), intent(inout) :: this

    call push_sub('lcao.lcao_end')

    SAFE_DEALLOCATE_P(this%atom)
    SAFE_DEALLOCATE_P(this%level)
    SAFE_DEALLOCATE_P(this%ddim)

    this%state = 0
    call pop_sub()
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(this, st, gr, geo, hm, start)
    type(lcao_t),        intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(hamiltonian_t), intent(in)    :: hm
    integer, optional,   intent(in)    :: start

    integer :: start_
    type(profile_t), save :: prof

    ASSERT(this%state == 1)

    call profiling_in(prof, "LCAO")
    call push_sub('lcao.lcao_wf')

    start_ = 1
    if(present(start)) start_ = start

    if (states_are_real(st)) then
      call dlcao_wf(this, st, gr, geo, hm, start_)
    else
      call zlcao_wf(this, st, gr, geo, hm, start_)
    end if

    call pop_sub()
    call profiling_out(prof)
  end subroutine lcao_wf

  logical function lcao_is_available(this) result(available)
    type(lcao_t),        intent(in) :: this
    
    available = this%state == 1
  end function lcao_is_available

  integer function lcao_num_orbitals(this) result(norbs)
    type(lcao_t),        intent(in) :: this

    norbs = this%norbs
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
