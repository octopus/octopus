!! Copyright (C) 2020 Nicolas Tancogne-Dejean, Sebastian Ohlmann, Heiko Appel
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

module propagator_exp_mid_oct_m
  use clock_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use propagator_abst_oct_m

  implicit none

  private
  public ::                            &
    propagator_exp_mid_t

  type, extends(propagator_abst_t) :: propagator_exp_mid_t
    private
  end type propagator_exp_mid_t

  interface propagator_exp_mid_t
    procedure propagator_exp_mid_init
  end interface propagator_exp_mid_t

contains

  ! ---------------------------------------------------------
  function propagator_exp_mid_init(namespace, predictor_corrector) result(this)
    type(namespace_t),         intent(in) :: namespace
    logical,                   intent(in) :: predictor_corrector
    type(propagator_exp_mid_t), pointer   :: this

    PUSH_SUB(propagator_exp_mid_init)

    SAFE_ALLOCATE(this)

    this%predictor_corrector = predictor_corrector
    this%start_step = EXPMID_START
    this%final_step = EXPMID_FINISH

    if(predictor_corrector) then

      call this%add(STORE_CURRENT_STATUS)
      call this%add(EXPMID_PREDICT_DT_2)  ! predict: psi(t+dt/2) = 0.5*(U_H(dt) psi(t) + psi(t)) or via extrapolation
      call this%add(START_SCF_LOOP)
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)   ! update: H(t+dt/2) from psi(t+dt/2)
      call this%add(EXPMID_PREDICT_DT)    ! predict: psi(t+dt) = U_H(t+dt/2) psi(t)
      call this%add(EXPMID_CORRECT_DT_2)  ! correct: psi(t+dt/2) = 0.5*(psi(t+dt) + psi(t))
      call this%add(END_SCF_LOOP) 
      call this%add(FINISHED)

      this%max_scf_count = 10
      this%scf_tol = CNST(1e-6) ! At the moment arbitrary. This is system specific and should be adapted.

    else

      call this%add(EXPMID_PREDICT_DT_2)  ! predict: psi(t+dt/2) = 0.5*(U_H(dt) psi(t) + psi(t)) or via extrapolation
      call this%add(UPDATE_INTERACTIONS)
      call this%add(UPDATE_HAMILTONIAN)   ! update: H(t+dt/2) from psi(t+dt/2)
      call this%add(EXPMID_PREDICT_DT)    ! predict: psi(t+dt) = U_H(t+dt/2) psi(t)
      call this%add(FINISHED)

    end if

    ! Exponential midpoint has only one algorithmic step
    this%algo_steps = 1

    call this%parse_td_variables(namespace)

    POP_SUB(propagator_exp_mid_init)
  end function propagator_exp_mid_init

end module propagator_exp_mid_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
