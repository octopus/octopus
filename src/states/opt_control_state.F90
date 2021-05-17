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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"


!> This module holds the "opt_control_state_t" datatype, which contains a quantum-classical
!! state.
module opt_control_state_oct_m
  use global_oct_m
  use ions_oct_m
  use messages_oct_m
  use profiling_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m

  implicit none

  private
  public :: opt_control_state_t,       &
            opt_control_state_init,    &
            opt_control_get_qs,        &
            opt_control_get_classical, &
            opt_control_set_classical, &
            opt_control_point_qs,      &
            opt_control_point_q,       &
            opt_control_point_p,       &
            opt_control_state_copy,    &
            opt_control_state_end,     &
            opt_control_state_null


  !> This is the datatype that contains the objects that are propagated: in principle this
  !! could be both the quantum and the classical subsystems, but for the moment it is only
  !! the quantum subsystem. So this data type is merely a wrapper around the states_elec_t data type.
  type opt_control_state_t
    private
    type(states_elec_t) :: psi
    FLOAT, allocatable :: q(:, :)
    FLOAT, allocatable :: p(:, :)
    integer :: natoms, ndim
  end type opt_control_state_t

contains


  subroutine opt_control_state_null(ocs)
    type(opt_control_state_t), intent(inout) :: ocs

    PUSH_SUB(opt_control_state_null)

    SAFE_DEALLOCATE_A(ocs%q)
    SAFE_DEALLOCATE_A(ocs%p)
    call ocs%psi%nullify()

    POP_SUB(opt_control_state_null)
  end subroutine opt_control_state_null


  function opt_control_point_qs(ocs)
    type(opt_control_state_t), target, intent(in) :: ocs
    type(states_elec_t), pointer :: opt_control_point_qs

    PUSH_SUB(opt_control_point_qs)

    opt_control_point_qs => ocs%psi

    POP_SUB(opt_control_point_qs)
  end function opt_control_point_qs

  function opt_control_point_q(ocs)
    type(opt_control_state_t), target, intent(in) :: ocs
    FLOAT, pointer :: opt_control_point_q(:, :)

    PUSH_SUB(opt_control_point_q)

    opt_control_point_q => ocs%q

    POP_SUB(opt_control_point_q)
  end function opt_control_point_q

  function opt_control_point_p(ocs)
    type(opt_control_state_t), target, intent(in) :: ocs
    FLOAT, pointer :: opt_control_point_p(:, :)

    PUSH_SUB(opt_control_point_p)

    opt_control_point_p => ocs%p

    POP_SUB(opt_control_point_p)
  end function opt_control_point_p

  subroutine opt_control_get_qs(qstate, ocs)
    type(states_elec_t),       intent(inout) :: qstate
    type(opt_control_state_t), intent(in)    :: ocs

    PUSH_SUB(opt_control_get_qs)

    call states_elec_copy(qstate, ocs%psi)

    POP_SUB(opt_control_get_qs)
  end subroutine opt_control_get_qs

  subroutine opt_control_get_classical(ions, ocs)
    type(ions_t),              intent(inout) :: ions
    type(opt_control_state_t), intent(in)    :: ocs

    integer :: idim, iatom

    PUSH_SUB(opt_control_get_classical)

    do idim = 1, ions%space%dim
      do iatom = 1, ions%natoms
        ions%pos(idim, iatom) = ocs%q(iatom, idim)
        ions%vel(idim, iatom) = ocs%p(iatom, idim) / ions%mass(iatom)
      end do
    end do

    POP_SUB(opt_control_get_classical)
  end subroutine opt_control_get_classical

  subroutine opt_control_set_classical(ions, ocs)
    type(ions_t),              intent(inout) :: ions
    type(opt_control_state_t), intent(inout) :: ocs

    integer :: idim, iatom

    PUSH_SUB(opt_control_set_classical)

    do idim = 1, ions%space%dim
      do iatom = 1, ions%natoms
        ocs%q(iatom, idim) = ions%pos(idim, iatom)
        ocs%p(iatom, idim) = ions%vel(idim, iatom) * ions%mass(iatom)
      end do
    end do

    POP_SUB(opt_control_set_classical)
  end subroutine opt_control_set_classical

  subroutine opt_control_state_init(ocs, qstate, ions)
    type(opt_control_state_t), intent(inout) :: ocs
    type(states_elec_t),       intent(in)    :: qstate
    type(ions_t),              intent(in)    :: ions

    integer :: iatom, idim

    PUSH_SUB(opt_control_state_init)

    call states_elec_copy(ocs%psi, qstate)

    SAFE_DEALLOCATE_A(ocs%q)
    SAFE_DEALLOCATE_A(ocs%p)

    ocs%ndim   = ions%space%dim
    ocs%natoms = ions%natoms

    SAFE_ALLOCATE(ocs%q(1:ocs%natoms, 1:ocs%ndim))
    SAFE_ALLOCATE(ocs%p(1:ocs%natoms, 1:ocs%ndim))

    do idim = 1, ions%space%dim
      do iatom = 1, ions%natoms
        ocs%q(iatom, idim) = ions%pos(idim, iatom)
        ocs%p(iatom, idim) = ions%mass(iatom) * ions%vel(idim, iatom)
      end do
    end do
    
    POP_SUB(opt_control_state_init)
  end subroutine opt_control_state_init

  subroutine opt_control_state_end(ocs)
    type(opt_control_state_t), intent(inout) :: ocs

    PUSH_SUB(opt_control_state_end)

    call states_elec_end(ocs%psi)

    SAFE_DEALLOCATE_A(ocs%q)
    SAFE_DEALLOCATE_A(ocs%p)

    POP_SUB(opt_control_state_end)
  end subroutine opt_control_state_end

  subroutine opt_control_state_copy(ocsout, ocsin)
    type(opt_control_state_t), intent(in)    :: ocsin
    type(opt_control_state_t), intent(inout) :: ocsout

    PUSH_SUB(opt_control_state_copy)

    call states_elec_end(ocsout%psi)
    call states_elec_copy(ocsout%psi, ocsin%psi)
    ocsout%ndim = ocsin%ndim
    ocsout%natoms = ocsin%natoms
    SAFE_DEALLOCATE_A(ocsout%q)
    SAFE_DEALLOCATE_A(ocsout%p)
    if(allocated(ocsin%q)) then
      SAFE_ALLOCATE(ocsout%q(1:ocsout%natoms, 1:ocsout%ndim))
      ocsout%q = ocsin%q
    end if
    if(allocated(ocsin%p)) then
      SAFE_ALLOCATE(ocsout%p(1:ocsout%natoms, 1:ocsout%ndim))
      ocsout%p = ocsin%p
    end if

    POP_SUB(opt_control_state_copy)
  end subroutine opt_control_state_copy

end module opt_control_state_oct_m
