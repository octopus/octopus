!! Copyright (C) 2020 Heiko Appel
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

module propagator_etrs_oct_m
  use clock_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_oct_m

  implicit none

  private
  public ::                            &
    propagator_etrs_t

  type, extends(propagator_t) :: propagator_etrs_t
    private
  end type propagator_etrs_t

  interface propagator_etrs_t
    procedure propagator_etrs_constructor
  end interface propagator_etrs_t

contains

  ! ---------------------------------------------------------
  function propagator_etrs_constructor(dt, predictor_corrector) result(this)
    FLOAT,                     intent(in) :: dt
    logical,                   intent(in) :: predictor_corrector
    type(propagator_etrs_t), pointer   :: this

    PUSH_SUB(propagator_etrs_constructor)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = predictor_corrector
    this%start_step = ETRS_START
    this%final_step = ETRS_FINISH

    ! In the following, we use the convention t1 = t, t2 = t+dt/2, and t3 = t+dt
    if(predictor_corrector) then

      call this%add(STORE_CURRENT_STATUS)
      call this%add(ETRS_STORE_STATE_T1)       ! store current state: psi_t1 = psi(t)
      call this%add(ETRS_STORE_HAMILTONIAN_T1) ! store current Hamiltonian: H_t1 = H(t)
      call this%add(ETRS_PROPAGATE_T1_DT)      ! predict state: psi(t+dt) = U(dt)_H_t1 psi(t)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)        ! compute: H_t3 from psi(t+dt)
      call this%add(ETRS_RESET_STATE_T1)       ! reset state: psi(t) = psi_t1
      call this%add(ETRS_PROPAGATE_T1_DT_2)    ! predict: psi(t+dt/2) = U(dt/2)_H_t1 psi(t)
      call this%add(ETRS_STORE_STATE_T2)       ! store current state: psi_t2 = psi(t+dt/2)
      call this%add(START_SCF_LOOP)
      call this%add(ETRS_RESET_STATE_T2)       ! reset state: psi(t+dt/2) = psi_t2(t), only for scf iter > 1
      call this%add(ETRS_PROPAGATE_T3_DT_2)    ! predict: psi(t+dt) = U(dt/2)_H_t3 psi(t+dt/2)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)        ! update: H_t3 from psi(t+dt)
      call this%add(END_SCF_LOOP)
      call this%add(FINISHED)

      this%max_scf_count = 10
      this%scf_tol = CNST(1e-6) ! At the moment arbitrary. This is system specific and should be adapted.

    else

      call this%add(STORE_CURRENT_STATUS)
      call this%add(ETRS_STORE_STATE_T1)       ! store current state: psi_t1 = psi(t)
      call this%add(ETRS_STORE_HAMILTONIAN_T1) ! store current Hamiltonian: H_t1 = H(t)
      call this%add(ETRS_PROPAGATE_T1_DT)      ! predict state: psi(t+dt) = U(dt)_H_t1 psi(t)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)        ! compute: H_t3 from psi(t+dt)
      call this%add(ETRS_RESET_STATE_T1)       ! reset state: psi(t) = psi_t1
      call this%add(ETRS_PROPAGATE_T1_DT_2)    ! predict: psi(t+dt/2) = U(dt/2)_H_t1 psi(t)
      call this%add(ETRS_PROPAGATE_T3_DT_2)    ! predict: psi(t+dt) = U(dt/2)_H_t3 psi(t+dt/2)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)        ! update: H_t3 from psi(t+dt)
      call this%add(FINISHED)

    end if

    ! The ETRS propagator has two algorithmic steps
    this%algo_steps = 2

    this%dt = dt

    POP_SUB(propagator_etrs_constructor)
  end function propagator_etrs_constructor

end module propagator_etrs_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
