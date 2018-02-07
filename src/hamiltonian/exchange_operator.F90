!! Copyright (C) 2018 N. Tancogne-Dejean
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

module exchange_operator_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use kpoints_oct_m
  use mesh_oct_m
  use messages_oct_m
  use par_vec_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scdm_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_parallel_oct_m

  implicit none

  private
  public ::                          &
    exchange_operator_t,             &
    exchange_operator_nullify,       &
    exchange_operator_init,          &
    exchange_operator_end,           &
    dexchange_operator_single,       &
    zexchange_operator_single,       &
    dexchange_operator_apply,        &
    zexchange_operator_apply,        &
    dexchange_operator_hartree_apply,&
    zexchange_operator_hartree_apply,&
    dexchange_operator_scdm_apply,   &
    zexchange_operator_scdm_apply,   &
    dexchange_operator_rdmft_apply,  &
    zexchange_operator_rdmft_apply

  type exchange_operator_t
    type(states_t), pointer :: st
    FLOAT :: exx_coef !< how much of EXX to mix
    FLOAT :: cam_omega
    FLOAT :: cam_alpha
    FLOAT :: cam_beta

    type(scdm_t)  :: scdm
  end type exchange_operator_t

contains

  subroutine exchange_operator_nullify(this)
    type(exchange_operator_t), intent(out) :: this

    PUSH_SUB(exchange_operator_nullify)

    nullify(this%st)

    this%exx_coef = M_ZERO
    this%cam_omega = M_ZERO
    this%cam_alpha = M_ZERO
    this%cam_beta  = M_ZERO

    POP_SUB(exchange_operator_nullify)
  end subroutine exchange_operator_nullify
 
  subroutine exchange_operator_init(this, st, exxcoef, omega, alpha, beta)
    type(exchange_operator_t), intent(out) :: this
    type(states_t), target,    intent(in)  :: st
    FLOAT,                     intent(in)  :: exxcoef, omega, alpha, beta

    PUSH_SUB(exchange_operator_init)

    this%st => st

    this%exx_coef = exxcoef
    this%cam_omega = omega
    this%cam_alpha = alpha
    this%cam_beta  = beta
    
    POP_SUB(exchange_operator_init)
  end subroutine exchange_operator_init

  subroutine exchange_operator_end(this)
    type(exchange_operator_t), intent(out) :: this

    PUSH_SUB(exchange_operator_end)

   ! this is a bit ugly, hf_st is initialized in v_ks_calc but deallocated here.
    if(associated(this%st))  then
      if(this%st%parallel_in_states .or. this%st%d%kpt%parallel) &
        call states_parallel_remote_access_stop(this%st)
      call states_end(this%st)
      SAFE_DEALLOCATE_P(this%st)
    end if
    
    POP_SUB(exchange_operator_end)
  end subroutine exchange_operator_end

#include "undef.F90"
#include "real.F90"
#include "exchange_operator_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "exchange_operator_inc.F90"

end module exchange_operator_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
