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
!! $Id: opt_control_state.F90 6827 2010-07-19 16:13:27Z acastro $

#include "global.h"


!> This module holds the "oct_state_t" datatype, which contains all the information
!! about the "state" that is propagated, i.e. for the moment being just the electronic
!! state, that is a state_t datatype. However, it is intended to include the ionic positions
!! and velocities, also, so that one can propagate (and optimize) simultaneously the quantum
!! and classical subsystems.
module opt_control_state_m
  use geometry_m
  use global_m
  use loct_pointer_m
  use messages_m
  use profiling_m
  use species_m
  use states_m

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
            opt_control_state_end

  !> This is the datatype that contains the objects that are propagated: in principle this
  !! could be both the quantum and the classical subsystems, but for the moment it is only
  !! the quantum subsystem. So this data type is merely a wrapper around the states_t data type.
  type opt_control_state_t
    private
    type(states_t) :: psi
    FLOAT, allocatable :: q(:, :)
    FLOAT, allocatable :: p(:, :)
    integer :: natoms, ndim
  end type opt_control_state_t

contains

  function opt_control_point_qs(ocs)
    type(states_t), pointer :: opt_control_point_qs
    type(opt_control_state_t), target :: ocs
    opt_control_point_qs => ocs%psi
  end function opt_control_point_qs

  function opt_control_point_q(ocs)
    FLOAT, pointer :: opt_control_point_q(:, :)
    type(opt_control_state_t), target :: ocs
    opt_control_point_q => ocs%q
  end function opt_control_point_q

  function opt_control_point_p(ocs)
    FLOAT, pointer :: opt_control_point_p(:, :)
    type(opt_control_state_t), target :: ocs
    opt_control_point_p => ocs%p
  end function opt_control_point_p

  subroutine opt_control_get_qs(qstate, ocs)
    type(states_t), intent(inout)         :: qstate
    type(opt_control_state_t), intent(in) :: ocs
    PUSH_SUB(opt_control_get_qs)

    call states_copy(qstate, ocs%psi)

    POP_SUB(opt_control_get_qs)
  end subroutine opt_control_get_qs

  subroutine opt_control_get_classical(geo, ocs)
    type(geometry_t), intent(inout)       :: geo
    type(opt_control_state_t), intent(in) :: ocs

    integer :: idim, iatom
    PUSH_SUB(opt_control_get_classical)

    do idim = 1, geo%space%dim
      do iatom = 1, geo%natoms
        geo%atom(iatom)%x(idim) = ocs%q(iatom, idim)
        geo%atom(iatom)%v(idim) = ocs%p(iatom, idim) / species_weight(geo%atom(iatom)%spec)
      end do
    end do

    POP_SUB(opt_control_get_classical)
  end subroutine opt_control_get_classical

  subroutine opt_control_set_classical(geo, ocs)
    type(geometry_t),          intent(inout) :: geo
    type(opt_control_state_t), intent(inout) :: ocs

    integer :: idim, iatom
    PUSH_SUB(opt_control_set_classical)

    do idim = 1, geo%space%dim
      do iatom = 1, geo%natoms
        ocs%q(iatom, idim) = geo%atom(iatom)%x(idim)
        ocs%p(iatom, idim) = geo%atom(iatom)%v(idim) * species_weight(geo%atom(iatom)%spec)
      end do
    end do

    POP_SUB(opt_control_set_classical)
  end subroutine opt_control_set_classical

  subroutine opt_control_state_init(ocs, qstate, geo)
    type(opt_control_state_t), intent(inout) :: ocs
    type(states_t), intent(in)               :: qstate
    type(geometry_t), intent(in)             :: geo

    integer :: iatom, idim
    PUSH_SUB(opt_control_state_init)

    call states_copy(ocs%psi, qstate)

    SAFE_DEALLOCATE_A(ocs%q)
    SAFE_DEALLOCATE_A(ocs%p)

    ocs%ndim   = geo%space%dim
    ocs%natoms = geo%natoms

    SAFE_ALLOCATE(ocs%q(1:ocs%natoms, 1:ocs%ndim))
    SAFE_ALLOCATE(ocs%p(1:ocs%natoms, 1:ocs%ndim))

    do idim = 1, geo%space%dim
      do iatom = 1, geo%natoms
        ocs%q(iatom, idim) = geo%atom(iatom)%x(idim)
        ocs%p(iatom, idim) = species_weight(geo%atom(iatom)%spec) * geo%atom(iatom)%v(idim)
      end do
    end do
    
    POP_SUB(opt_control_state_init)
  end subroutine opt_control_state_init

  subroutine opt_control_state_end(ocs)
    type(opt_control_state_t), intent(inout) :: ocs

    PUSH_SUB(opt_control_state_end)

    call states_end(ocs%psi)

    SAFE_DEALLOCATE_A(ocs%q)
    SAFE_DEALLOCATE_A(ocs%p)

    POP_SUB(opt_control_state_end)
  end subroutine opt_control_state_end

  subroutine opt_control_state_copy(ocsout, ocsin)
    type(opt_control_state_t), intent(in) :: ocsin
    type(opt_control_state_t), intent(inout) :: ocsout

    call states_copy(ocsout%psi, ocsin%psi)
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

  end subroutine opt_control_state_copy

end module opt_control_state_m
