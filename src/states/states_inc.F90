!! Copyright (C) 2011 X. Andrade
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
!!
!! $Id: states_inc.F90 7770 2011-04-26 00:51:17Z dstrubbe $

subroutine X(states_get_state1)(st, mesh, ist, iqn, psi)
  type(states_t), intent(in)    :: st
  type(mesh_t),   intent(in)    :: mesh
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: iqn
  R_TYPE,         intent(out)   :: psi(:, :)
  
  integer :: idim, ip

  if(states_are_real(st)) then
    do idim = 1, st%d%dim
      forall(ip = 1:mesh%np) psi(ip, idim) = st%dpsi(ip, idim, ist, iqn)
    end do
  else
#ifdef R_TREAL
    ASSERT(.false.)
#endif
    do idim = 1, st%d%dim
      forall(ip = 1:mesh%np) psi(ip, idim) = st%zpsi(ip, idim, ist, iqn)
    end do
  end if

end subroutine X(states_get_state1)

! ------------------------------------------------------------

subroutine X(states_get_state2)(st, mesh, idim, ist, iqn, psi)
  type(states_t), intent(in)    :: st
  type(mesh_t),   intent(in)    :: mesh
  integer,        intent(in)    :: idim
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: iqn
  R_TYPE,         intent(out)   :: psi(:)
  
  integer :: ip

  if(states_are_real(st)) then
    forall(ip = 1:mesh%np) psi(ip) = st%dpsi(ip, idim, ist, iqn)
  else
#ifdef R_TREAL
    ASSERT(.false.)
#endif
    forall(ip = 1:mesh%np) psi(ip) = st%zpsi(ip, idim, ist, iqn)
  end if

end subroutine X(states_get_state2)

! ------------------------------------------------------------

subroutine X(states_set_state1)(st, mesh, ist, iqn, psi)
  type(states_t), intent(inout) :: st
  type(mesh_t),   intent(in)    :: mesh
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: iqn
  R_TYPE,         intent(in)    :: psi(:, :)
  
  integer :: idim, ip

  if(states_are_real(st)) then
#ifdef R_TCOMPLEX
    ASSERT(.false.)
#endif
    do idim = 1, st%d%dim
      forall(ip = 1:mesh%np) st%dpsi(ip, idim, ist, iqn) = psi(ip, idim)
    end do
  else
    do idim = 1, st%d%dim
      forall(ip = 1:mesh%np) st%zpsi(ip, idim, ist, iqn) = psi(ip, idim)
    end do
  end if

end subroutine X(states_set_state1)

! ------------------------------------------------------------

subroutine X(states_set_state2)(st, mesh, idim, ist, iqn, psi)
  type(states_t), intent(inout) :: st
  type(mesh_t),   intent(in)    :: mesh
  integer,        intent(in)    :: idim
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: iqn
  R_TYPE,         intent(in)    :: psi(:)
  
  integer :: ip

  if(states_are_real(st)) then
#ifdef R_TCOMPLEX
    ASSERT(.false.)
#endif
    forall(ip = 1:mesh%np) st%dpsi(ip, idim, ist, iqn) = psi(ip)
  else
    forall(ip = 1:mesh%np) st%zpsi(ip, idim, ist, iqn) = psi(ip)
  end if
  
end subroutine X(states_set_state2)

! ------------------------------------------------------------
