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
!!
!! $Id$

#include "global.h"

subroutine eigen_solver_jdqz(m, f_der, st, h, tol, niter, converged, errorflag, diff, reorder, verbose)
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(IN)    :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(out)   :: errorflag
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,      optional, intent(in)    :: reorder
  logical,      optional, intent(in)    :: verbose

  call jdqz


end subroutine eigen_solver_jdqz
