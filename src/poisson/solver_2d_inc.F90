!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
subroutine poisson2D_init(this)
  type(poisson_t), intent(inout) :: this

  PUSH_SUB(poisson2D_init)

  if(this%method == POISSON_FFT) then
    call poisson_fft_init(this%fft_solver, this%der%mesh, this%cube, this%kernel)
  end if

  POP_SUB(poisson2D_init)

end subroutine poisson2D_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
