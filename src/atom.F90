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
    
    integer :: j, ll
    FLOAT :: x(3), a(3), r, p, ylm
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

  function atom_density(m, atom, spin_channels) result(rho)
    type(mesh_type), intent(in) :: m
    type(atom_type), intent(in) :: atom
    integer,         intent(in) :: spin_channels
    FLOAT                       :: rho(m%np, spin_channels)

    integer :: opt, i, in_points, n
    integer, save :: j = 1
    FLOAT :: r
    R_TYPE :: psi1
    type(specie_type), pointer :: s

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    s => atom%spec

    ! select what kind of density we should build
    if(conf%dim==3) then
      select case(s%label(1:5))
      case('usdef')
        opt = 1
      case('jelli', 'point')
        opt = 2
      case default
        opt = 3
      end select
    else
      opt = 1      
    end if

    rho = M_ZERO
    ! build density...
    select case (opt)
    case (1) ! ...from userdef
      rho(1:m%np, 1:spin_channels) = real(s%Z_val, PRECISION) &
                                     /(m%np*m%vol_pp*real(spin_channels, PRECISION))
    case (2) ! ...from jellium
      in_points = 0
      do i = 1, m%np
        call mesh_r(m, i, r, a=atom%x)
        if(r <= s%jradius) then
          in_points = in_points + 1
        end if
      end do
    
      if(in_points > 0) then
        do i = 1, m%np
          call mesh_r(m, i, r, a=atom%x)
          if(r <= s%jradius) then
            rho(i, 1:spin_channels) = real(s%z_val, PRECISION)/ &
                                      (real(in_points*spin_channels, PRECISION)*m%vol_pp)
          end if
        end do
      end if

    case (3) ! ...from pseudopotential
      do i = 1, m%np
        call mesh_r(m, i, r, a=atom%x)
        do n = 1, s%ps%conf%p
          if(r >= r_small) then
            select case(spin_channels)
            case(1)
              psi1 = loct_splint(s%ps%Ur(n, 1), r)
              rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
            case(2)
              ! This is still a bit weird, but let us see how it works...
              psi1 = loct_splint(s%ps%Ur(n, 1), r)
              rho(i, mod(j,2)+1)   = rho(i, mod(j,2)+1)   + s%ps%conf%occ(n, 1)*psi1*psi1 / (M_FOUR*M_PI)
              rho(i, mod(j+1,2)+1) = rho(i, mod(j+1,2)+1) + s%ps%conf%occ(n, 2)*psi1*psi1 / (M_FOUR*M_PI)
              j = j + 1
            end select
          end if
        end do
      end do

    end select

  end function atom_density

end module atom
