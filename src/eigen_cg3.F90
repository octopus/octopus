!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "config_F90.h"

subroutine eigen_solver_cg3(st, sys, h, tol, niter, converged, errorflag, diff, reorder)
  use trl_info
  use trl_interface

  type(states_type), target, intent(inout)   :: st
  type(system_type), target, intent(IN)      :: sys
  type(hamiltonian_type), target, intent(IN) :: h
  real(r8), intent(in)               :: tol
  integer, intent(inout)             :: niter
  integer, intent(out)               :: errorflag, converged
  real(r8), intent(out), optional    :: diff(1:st%nst,1:st%nik)
  logical, intent(in), optional      :: reorder

  integer :: unit, i
  integer :: nrow, lohi, ned, maxlan, mev, lwrk
  type(trl_info_t) :: info
!!$  R_TYPE, allocatable :: evec(:, :, :)

  call push_sub('eigen_solver_cg3')

  nrow = (sys%m%np+1)*st%dim; ned = st%nst; maxlan = (ned + min(6, ned)); mev = st%nst

!!$  allocate(evec(sys%m%np, st%dim, mev))

  call trl_init_info(info, nrow, maxlan, lohi = -1, ned = ned, tol = 1.0e-5_r8, trestart = 3)

  if(conf%verbose > 999) then
    call io_assign(unit)
    call trl_set_debug(info, msglvl = 10, iou = unit, file = 'tmp/trlan_log_')
  endif

  call trl_set_iguess(info, nec = 0, iguess = 1)
  do ik_trlan = 1, st%nik
     st%R_FUNC(psi)(:, :, ik_trlan, 1) = R_TOTYPE(M_ZERO)
     do i = 1, st%nst
        st%R_FUNC(psi)(:, :, 1, ik_trlan) = st%R_FUNC(psi)(:, :, 1, ik_trlan) + &
                                            st%R_FUNC(psi)(:, :, i, ik_trlan)/st%nst
     enddo
  enddo


  h_trlan   => h
  m_trlan   => sys%m
  st_trlan  => st
  sys_trlan => sys

  do ik_trlan = 1, st%nik
!!$     do i = 1, st%nst
!!$        evec(:, :, i) = st%R_FUNC(psi)(:, :, i, ik_trlan)
!!$     enddo
     call trlan(op, info, nrow, mev, &
                st%eigenval(:, ik_trlan), st%R_FUNC(psi)(:, :, :, ik_trlan), nrow)
!!$     call trlan(op, info, nrow, mev, &
!!$                st%eigenval(:, ik_trlan), evec(:,:,:), nrow)
     if(info%stat.ne.0) then
        write(message(1),'(a,i5)') 'trlan call returned with error flag', info%stat
        call write_fatal(1)
     endif
     converged = info%nec
     niter     = info%matvec

!!$     do i = 1, st%nst
!!$        st%R_FUNC(psi)(:, :, i, ik_trlan) = evec(:, :, i)
!!$     enddo

  enddo

  call R_FUNC(scal)(sys%m%np*st%dim*st%nst*st%nik, &
                    R_TOTYPE(1.0/sqrt(sys%m%vol_pp)), st%R_FUNC(psi), 1)

  if(conf%verbose > 999) then
    call io_close(unit)
  endif

  nullify(h_trlan, m_trlan, st_trlan, sys_trlan)
  call pop_sub()
end subroutine eigen_solver_cg3

subroutine op(nrow, ncol, xin, ldx, yout, ldy)
  integer, intent(in) :: nrow, ncol, ldx, ldy
  real(r8), dimension(ldx, ncol), intent(in)  :: xin
  real(r8), dimension(ldy, ncol), intent(out) :: yout

  integer :: i

  do i = 1, ncol
     call R_FUNC(Hpsi) (h_trlan, m_trlan, st_trlan, sys_trlan, ik_trlan, &
          xin(1:, i), yout(2:,i))
  enddo

end subroutine op
