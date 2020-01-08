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
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scdm_oct_m
  use simul_box_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use wfs_elec_oct_m

  implicit none

  private
  public ::                          &
    exchange_operator_t,             &
    exchange_operator_nullify,       &
    exchange_operator_init,          &
    exchange_operator_reinit,        &
    exchange_operator_end,           &
    dexchange_operator_single,       &
    zexchange_operator_single,       &
    dexchange_operator_apply,        &
    zexchange_operator_apply,        &
    dexchange_operator_hartree_apply,&
    zexchange_operator_hartree_apply,&
    dexchange_operator_scdm_apply,   &
    zexchange_operator_scdm_apply,   &
    exchange_operator_rdmft_occ_apply

  type exchange_operator_t
    type(states_elec_t), public, pointer :: st
    FLOAT :: cam_omega
    FLOAT :: cam_alpha
    FLOAT :: cam_beta

    logical :: user_defined_cam

    type(scdm_t)  :: scdm

    type(namespace_t), pointer :: namespace
  end type exchange_operator_t
 
contains

  subroutine exchange_operator_nullify(this)
    type(exchange_operator_t), intent(out) :: this

    PUSH_SUB(exchange_operator_nullify)

    nullify(this%st)

    this%cam_omega = M_ZERO
    this%cam_alpha = M_ZERO
    this%cam_beta  = M_ZERO

    this%user_defined_cam = .false.

    POP_SUB(exchange_operator_nullify)
  end subroutine exchange_operator_nullify
 
  subroutine exchange_operator_init(this, namespace, st_d, sb, mesh, omega, alpha, beta)
    type(exchange_operator_t), intent(inout) :: this
    type(namespace_t), target, intent(in)    :: namespace
    type(states_elec_dim_t),   intent(in)    :: st_d
    type(simul_box_t),         intent(in)    :: sb
    type(mesh_t),              intent(in)    :: mesh
    FLOAT,                     intent(in)    :: omega, alpha, beta

    type(block_t) :: blk

    PUSH_SUB(exchange_operator_init)

    this%namespace => namespace

    this%cam_omega = omega
    this%cam_alpha = alpha
    this%cam_beta  = beta

    !%Variable UserDefinedCamParameters
    !%Type block
    !%Section Hamiltonian::XC
    !%Description
    !% (Experimental) This variable allows to define the CAM parameters of an range-separated hybrid
    !% thus overwriting the values that are returned by libxc.
    !% The values are defined in order as omega, alpha, and beta.
    !% If the functional does not depends on these parameters, the value is simply ignored by the code.
    !% At the moment only the weight of the exact can be tuned, and only for the HSE06 or the PBE0 functionals.
    !%End
    if(parse_block(namespace, 'UserDefinedCamParameters', blk)==0) then
      call parse_block_float(blk, 0, 0, this%cam_omega)
      call parse_block_float(blk, 0, 1, this%cam_alpha)
      call parse_block_float(blk, 0, 2, this%cam_beta)
      call parse_block_end(blk) 
      this%user_defined_cam = .true.
      call messages_experimental('UserDefinedCamParameters')
    else
      this%user_defined_cam = .false. 
    end if

    POP_SUB(exchange_operator_init)
  end subroutine exchange_operator_init

  subroutine exchange_operator_reinit(this, st, omega, alpha, beta)
    type(exchange_operator_t), intent(inout) :: this
    type(states_elec_t), target, intent(in)  :: st
    FLOAT,                     intent(in)    :: omega, alpha, beta

    PUSH_SUB(exchange_operator_reinit)

    this%st => st

    if(.not. this%user_defined_cam) then
      this%cam_omega = omega
      this%cam_alpha = alpha
      this%cam_beta  = beta
    end if

    POP_SUB(exchange_operator_reinit)
  end subroutine exchange_operator_reinit

  subroutine exchange_operator_end(this)
    type(exchange_operator_t), intent(inout) :: this

    PUSH_SUB(exchange_operator_end)

    if(associated(this%st)) then
      if(this%st%parallel_in_states) call states_elec_parallel_remote_access_stop(this%st)
      call states_elec_end(this%st)
      SAFE_DEALLOCATE_P(this%st)
    end if
    nullify(this%st)

    POP_SUB(exchange_operator_end)
  end subroutine exchange_operator_end

  subroutine exchange_operator_rdmft_occ_apply(this, mesh, hpsib)
    type(exchange_operator_t), intent(in) :: this
    type(mesh_t),              intent(in) :: mesh
    class(wfs_elec_t),      intent(inout) :: hpsib

    PUSH_SUB(exchange_operator_rdmft_occ_apply)

    ! multiply linear terms in hamiltonian with occupation number
    ! nonlinear occupation number dependency occurs only in the exchange, which is treated there
    call batch_scal(mesh%np, this%st%occ(:, hpsib%ik), hpsib)

    POP_SUB(exchange_operator_rdmft_occ_apply)
  end subroutine exchange_operator_rdmft_occ_apply


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
