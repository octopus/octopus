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
!! $Id$


  ! ---------------------------------------------------------
  !> Magnus propagator
  subroutine td_magnus(ks, hm, gr, st, tr, time, dt)
    type(v_ks_t), target,            intent(inout) :: ks
    type(hamiltonian_t), target,     intent(inout) :: hm
    type(grid_t),        target,     intent(inout) :: gr
    type(states_t),      target,     intent(inout) :: st
    type(propagator_t),  target,     intent(inout) :: tr
    FLOAT,                           intent(in)    :: time
    FLOAT,                           intent(in)    :: dt

    integer :: j, is, ist, ik, i
    FLOAT :: atime(2)
    FLOAT, allocatable :: vaux(:, :, :), pot(:)

    PUSH_SUB(propagator_dt.td_magnus)

    SAFE_ALLOCATE(vaux(1:gr%mesh%np, 1:st%d%nspin, 1:2))

    atime(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
    atime(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then
      do j = 1, 2
        call vksinterp_interpolate(tr%vksold, 3, gr%mesh%np, st%d%nspin, time, dt, atime(j)-dt, hm%vhxc)
        call hamiltonian_update(hm, gr%mesh)
      end do
    else
      vaux = M_ZERO
    end if

    do j = 1, 2
      ! WARNING: This should be carefully tested, and extended to allow for velocity-gauge laser fields.
      do i = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(i)))
        case(E_FIELD_ELECTRIC)
          SAFE_ALLOCATE(pot(1:gr%mesh%np))
          pot = M_ZERO
          call laser_potential(hm%ep%lasers(i), gr%mesh, pot, time - dt + atime(j))
          do is = 1, st%d%nspin
            vaux(:, is, j) = vaux(:, is, j) + pot(:)
          end do
          SAFE_DEALLOCATE_A(pot)
        case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
          write(message(1),'(a)') 'The Magnus propagator cannot be used with magnetic fields, or'
          write(message(2),'(a)') 'with an electric field described in the velocity gauge.'
          call messages_fatal(2)
        end select
      end do
    end do

    tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
    tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/CNST(12.0))*dt*(vaux(:, :, 2) - vaux(:, :, 1))

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        call exponential_apply(tr%te, gr%der, hm, st%zpsi(:,:, ist, ik), ist, ik, dt, M_ZERO, &
          vmagnus = tr%vmagnus)
      end do
    end do

    if(.not. hm%cmplxscl%space) then
      call density_calc(st, gr, st%rho)
    else
      call density_calc(st, gr, st%zrho%Re, st%zrho%Im)
    end if

    SAFE_DEALLOCATE_A(vaux)
    POP_SUB(propagator_dt.td_magnus)
  end subroutine td_magnus
