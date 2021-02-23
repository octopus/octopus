!! Copyright (C) 2002-2018 M. Marques, A. Castro, A. Rubio, G. Bertsch, N. Tancogne-Dejean
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
  use fourier_space_oct_m
  use global_oct_m
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use scdm_oct_m
  use simul_box_oct_m
  use singularity_oct_m
  use space_oct_m
  use symmetries_oct_m
  use symmetrizer_oct_m
  use states_abst_oct_m
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
    exchange_operator_rdmft_occ_apply,&
    dexchange_operator_compute_potentials, &
    zexchange_operator_compute_potentials, &
    dexchange_operator_commute_r, &
    zexchange_operator_commute_r, &
    dexchange_operator_ACE, &
    zexchange_operator_ACE


  type ACE_t
    integer :: nst
    FLOAT, allocatable :: dchi(:,:,:,:)
    CMPLX, allocatable :: zchi(:,:,:,:)
  end type ACE_t

  type exchange_operator_t
    type(states_elec_t), public, pointer :: st
    FLOAT :: cam_omega
    FLOAT :: cam_alpha
    FLOAT :: cam_beta

    type(poisson_t) :: psolver      !< Poisson solver

    type(scdm_t)  :: scdm

    type(singularity_t) :: singul !< Coulomb singularity

    logical       :: useACE
    type(ACE_t)   :: ace
    type(states_elec_t) :: xst !< The states after the application of the Fock operator
                               !! This is needed to construct the ACE operator
  end type exchange_operator_t
 
contains

  subroutine exchange_operator_nullify(this)
    type(exchange_operator_t), intent(out) :: this

    PUSH_SUB(exchange_operator_nullify)

    nullify(this%st)

    this%cam_omega = M_ZERO
    this%cam_alpha = M_ZERO
    this%cam_beta  = M_ZERO

    this%useACE = .false.

    POP_SUB(exchange_operator_nullify)
  end subroutine exchange_operator_nullify
 
  subroutine exchange_operator_init(this, namespace, space, st, sb, der, mc, kpoints, omega, alpha, beta)
    type(exchange_operator_t), intent(inout) :: this
    type(namespace_t),          intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    type(states_elec_t),       intent(in)    :: st
    type(simul_box_t),         intent(in)    :: sb
    type(derivatives_t),       intent(in)    :: der
    type(multicomm_t),         intent(in)    :: mc
    type(kpoints_t),           intent(in)    :: kpoints
    FLOAT,                     intent(in)    :: omega, alpha, beta

    PUSH_SUB(exchange_operator_init)

    this%cam_omega = omega
    this%cam_alpha = alpha
    this%cam_beta  = beta

    !%Variable AdaptivelyCompressedExchange
    !%Type logical
    !%Default false
    !%Section Hamiltonian
    !%Description
    !% (Experimental) If set to yes, Octopus will use the adaptively compressed exchange
    !% operator (ACE) for HF and hybrid calculations, as defined in
    !%  Lin, J. Chem. Theory Comput. 2016, 12, 2242.
    !%End
    call parse_variable(namespace, 'AdaptivelyCompressedExchange', .false., this%useACE)
    if(this%useACE) call messages_experimental('AdaptivelyCompressedExchange')

    if(this%useACE) then
      call this%xst%nullify()
    end if

    call singularity_init(this%singul, namespace, st, sb, kpoints)
    if(states_are_real(st)) then
      call poisson_init(this%psolver, namespace, space, der, mc, st%qtot, &
             force_serial = .true., verbose = .false.)
    else
      call poisson_init(this%psolver, namespace, space, der, mc, st%qtot, &
             force_serial = .true., verbose = .false., force_cmplx = .true.)
    end if

    POP_SUB(exchange_operator_init)
  end subroutine exchange_operator_init

  subroutine exchange_operator_reinit(this, st, omega, alpha, beta)
    type(exchange_operator_t), intent(inout) :: this
    type(states_elec_t), target, intent(in)  :: st
    FLOAT,                     intent(in)    :: omega, alpha, beta

    PUSH_SUB(exchange_operator_reinit)

    this%st => st

    this%cam_omega = omega
    this%cam_alpha = alpha
    this%cam_beta  = beta

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

    if(this%useACE) then
      call states_elec_end(this%xst)
      this%ace%nst = 0
      SAFE_DEALLOCATE_A(this%ace%dchi)
      SAFE_DEALLOCATE_A(this%ace%zchi)
    end if

    call singularity_end(this%singul)
    call poisson_end(this%psolver)

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
