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

#include "global.h"

module atom
  use global
  use lib_oct_gsl_spline
  use mesh
  use specie

  implicit none

  type atom_type
    character(len=10) :: label
    type(specie_type), pointer :: spec ! pointer to specie

    FLOAT :: x(3), v(3), f(3) ! position/velocity/force of atom in real space

    logical :: move              ! should I move this atom in the optimization mode
  end type atom_type

  type atom_classical_type
    FLOAT :: x(3), v(3), f(3)
    FLOAT :: charge

    character(len=4) :: label
  end type atom_classical_type

contains

  subroutine atom_get_wf(m, atom, l, lm, ispin, psi)
    type(mesh_type),        intent(in)    :: m
    type(atom_type), intent(in)    :: atom
    integer, intent(in)   :: l, lm, ispin
    R_TYPE, intent(out) :: psi(m%np)
    
    integer :: j, d2, ll
    FLOAT :: x(3), a(3), r, p, ylm, g(3)
    type(loct_spline_type), pointer :: s

    call push_sub('atom_get_wf')
    
    a = atom%x
    if(atom%spec%local) then
      ! add a couple of harmonic oscilator functions
    else
      s => atom%spec%ps%Ur(l, ispin)

      ll = atom%spec%ps%conf%l(l)
      do j = 1, m%np
        call mesh_r(m, j, r, x=x, a=a)
        p = loct_splint(s, r)
        ylm = loct_ylm(x(1), x(2), x(3), ll, lm)
        psi(j) = p * ylm
      end do
    end if
 
    call pop_sub()
  end subroutine atom_get_wf

end module atom
