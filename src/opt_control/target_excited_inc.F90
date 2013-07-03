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
!! $Id: target_excited_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_excited(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    integer :: ierr, ip
    PUSH_SUB(target_init_excited)

    message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
    call messages_info(1)

    call states_look (trim(restart_dir)//GS_DIR, gr%mesh%mpi_grp, ip, ip, tg%st%nst, ierr)
    tg%st%st_start = 1
    tg%st%st_end   = tg%st%nst

    SAFE_DEALLOCATE_P(tg%st%occ)
    SAFE_DEALLOCATE_P(tg%st%eigenval)
    SAFE_DEALLOCATE_P(tg%st%node)

    SAFE_ALLOCATE(     tg%st%occ(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(tg%st%eigenval(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(    tg%st%node(1:tg%st%nst))
    if(tg%st%d%ispin == SPINORS) then
      SAFE_DEALLOCATE_P(tg%st%spin)
      SAFE_ALLOCATE(tg%st%spin(1:3, 1:tg%st%nst, 1:tg%st%d%nik))
    end if
    call states_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)
    tg%st%node(:)  = 0

    call restart_read(trim(restart_dir)//GS_DIR, tg%st, gr, ierr, exact = .true.)
    call excited_states_init(tg%est, tg%st, "oct-excited-state-target") 

    POP_SUB(target_init_excited)
  end subroutine target_init_excited


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_excited(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi

    PUSH_SUB(target_j1_excited)

    j1 = abs(zstates_mpdotp(gr%mesh, tg%est, psi))**2

    POP_SUB(target_j1_excited)
  end function target_j1_excited


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
