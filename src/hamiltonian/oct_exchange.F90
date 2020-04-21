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

module oct_exchange_oct_m
  use derivatives_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                          &
    oct_exchange_t,                  &
    oct_exchange_nullify,            &
    oct_exchange_enabled,            &
    oct_exchange_prepare,            &
    oct_exchange_set,                &
    oct_exchange_remove,             &
    doct_exchange_operator,          &
    zoct_exchange_operator

  type oct_exchange_t
    private
    logical :: oct_exchange
    type(states_elec_t), pointer :: oct_st
    FLOAT, pointer :: oct_fxc(:, :, :)
    FLOAT, pointer :: oct_pot(:, :)
    FLOAT, pointer :: oct_rho(:, :)
  end type oct_exchange_t

contains

  subroutine oct_exchange_nullify(this)
    type(oct_exchange_t), intent(out) :: this

    PUSH_SUB(oct_exchange_nullify)

    this%oct_exchange = .false.
    nullify(this%oct_st)
    nullify(this%oct_fxc)
    nullify(this%oct_pot)
    nullify(this%oct_rho)

    POP_SUB(oct_exchange_nullify)
  end subroutine oct_exchange_nullify
  
  ! ---------------------------------------------------------
  logical function oct_exchange_enabled(this) result(oct_exchange)
    type(oct_exchange_t), intent(in) :: this

    PUSH_SUB(oct_exchange_enabled)
    oct_exchange = this%oct_exchange

    POP_SUB(oct_exchange_enabled)
  end function oct_exchange_enabled


  ! ---------------------------------------------------------
  subroutine oct_exchange_set(this, st, mesh)
    type(oct_exchange_t),     intent(inout) :: this
    type(states_elec_t), target, intent(in) :: st
    type(mesh_t),             intent(in)    :: mesh

    integer :: np, nspin

    PUSH_SUB(oct_exchange_set)

    ! In this release, no non-local part for the QOCT Hamiltonian.
    nullify(this%oct_st)

    this%oct_st => st
    this%oct_exchange = .true.
    np = mesh%np
    nspin = this%oct_st%d%nspin

    SAFE_ALLOCATE(this%oct_fxc(1:np, 1:nspin, 1:nspin))
    SAFE_ALLOCATE(this%oct_pot(1:np, 1:nspin))
    SAFE_ALLOCATE(this%oct_rho(1:np, 1:nspin))

    this%oct_fxc = M_ZERO
    this%oct_pot = M_ZERO
    this%oct_rho = M_ZERO

    POP_SUB(oct_exchange_set)
  end subroutine oct_exchange_set


  ! ---------------------------------------------------------
  subroutine oct_exchange_prepare(this, mesh, psi, xc, psolver, namespace)
    type(oct_exchange_t), intent(inout) :: this
    type(mesh_t),         intent(in)    :: mesh
    CMPLX,                intent(in)    :: psi(:, :, :, :)
    type(xc_t),           intent(in)    :: xc
    type(poisson_t),      intent(in)    :: psolver
    type(namespace_t),    intent(in)    :: namespace

    integer :: jst, ip, ik
    CMPLX, allocatable :: psi2(:, :)

    PUSH_SUB(oct_exchange_prepare)

    SAFE_ALLOCATE(psi2(1:mesh%np, 1:this%oct_st%d%dim))

    select case(this%oct_st%d%ispin)
    case(UNPOLARIZED)
      ASSERT(this%oct_st%d%nik  ==  1)

      this%oct_pot = M_ZERO
      this%oct_rho = M_ZERO
      do jst = 1, this%oct_st%nst
        call states_elec_get_state(this%oct_st, mesh, jst, 1, psi2)
        do ip = 1, mesh%np
          this%oct_rho(ip, 1) = this%oct_rho(ip, 1) + this%oct_st%occ(jst, 1)*aimag(conjg(psi2(ip, 1))*psi(ip, 1, jst, 1))
        end do
      end do
      call dpoisson_solve(psolver, this%oct_pot(:, 1), this%oct_rho(:, 1), all_nodes = .false.)

    case(SPIN_POLARIZED)
      ASSERT(this%oct_st%d%nik  ==  2)

      this%oct_pot = M_ZERO
      this%oct_rho = M_ZERO
      do ik = 1, 2
        do jst = 1, this%oct_st%nst
          call states_elec_get_state(this%oct_st, mesh, jst, ik, psi2)
          do ip = 1, mesh%np
            this%oct_rho(ip, ik) = this%oct_rho(ip, ik) + this%oct_st%occ(jst, ik) * aimag(conjg(psi2(ip, 1))*psi(ip, 1, jst, ik))
          end do
        end do
      end do

      do ik = 1, 2
        call dpoisson_solve(psolver, this%oct_pot(:, ik), this%oct_rho(:, ik), all_nodes = .false.)
      end do

    end select

    this%oct_fxc = M_ZERO
    call xc_get_fxc(xc, mesh, namespace, this%oct_st%rho, this%oct_st%d%ispin, this%oct_fxc)

    SAFE_DEALLOCATE_A(psi2)
    POP_SUB(oct_exchange_prepare)
  end subroutine oct_exchange_prepare


  ! ---------------------------------------------------------
  subroutine oct_exchange_remove(this)
    type(oct_exchange_t), intent(inout) :: this

    PUSH_SUB(oct_exchange_remove)

    nullify(this%oct_st)
    this%oct_exchange = .false.
    SAFE_DEALLOCATE_P(this%oct_fxc)
    SAFE_DEALLOCATE_P(this%oct_pot)
    SAFE_DEALLOCATE_P(this%oct_rho)

    POP_SUB(oct_exchange_remove)
  end subroutine oct_exchange_remove

#include "undef.F90"
#include "real.F90"
#include "oct_exchange_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "oct_exchange_inc.F90"

end module oct_exchange_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
