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

program oct_test
  use global
  use lib_oct
  use functions
  use mesh
  use geometry
  use output
  use poisson

  implicit none

  FLOAT :: val_charge, def_h, def_rsize, alpha, beta, lambda, r, x(3), ylm, gamma, il, rr(3), s1, mult
  type(geometry_type), pointer :: geo
  type(mesh_type), pointer     :: m
  type(f_der_type)             :: f_der
  type(output_type)            :: outp
  FLOAT, allocatable :: f(:), g(:), vh_exact(:), vh_numerical(:)
  integer :: ll, mm, i, lldfac, ierr, l, ml, j


  call global_init() 

#ifdef HAVE_FFT
  call fft_all_init()
#endif

  call units_init()

  allocate(geo)
  call geometry_init_xyz(geo)
  call geometry_init_species(geo, val_charge, def_h, def_rsize)

  allocate(m)
  call mesh_init(m, geo, def_h, def_rsize)
  call f_der_init(m, f_der)
  call mesh_create_xyz(m, f_der%n_ghost(1))
  call f_der_build(f_der)

  call mesh_write_info(m, stdout)

  call output_init(outp)

  call poisson_init(m)

  allocate(f(m%np), g(m%np), vh_exact(m%np), vh_numerical(m%np))
  f = M_ZERO
  g = M_ZERO

  ll = 3
  lldfac = 1; do j = 1, 2*ll+1, 2; lldfac = lldfac * j; enddo
  mm = 2


!!$  ! This builds a charge distribution with unity (ll,mm) moment.
!!$  alpha = M_TWO
!!$  beta = (2**(ll+2))/( alpha**(2*ll+3) * sqrt(M_PI) * lldfac )
!!$  gamma = ( sqrt(M_PI)*2**(ll+3) ) / lldfac
!!$
!!$  do i = 1, m%np
!!$     call mesh_r(m, i, r, x = x)
!!$     ylm  = loct_ylm(x(1), x(2), x(3), ll, mm)
!!$     f(i) = beta * r**ll * exp(-(r/alpha)**2)*ylm
!!$  enddo
!!$
!!$  ! This checks that indeed the multipole is the one we wanted.
!!$  do l = 0, 4
!!$     do ml = -l, l
!!$        mult = M_ZERO
!!$        do i = 1, m%np
!!$           call mesh_r(m, i, r, x = x)
!!$           s1 = f(i)*m%vol_pp(i)*r**l
!!$           ylm = loct_ylm(x(1), x(2), x(3), l, ml)
!!$           mult = mult + ylm*s1
!!$        enddo
!!$        write(*, '(a2,i2,a1,i3,a4,f12.8)') 'M[', l, ',', ml, '] = ', mult
!!$     enddo
!!$  enddo
!!$
!!$  ! This calculates the exact potential of the (ll,mm) moment charge distribution.
!!$  call vhexact_nlm(alpha, ll, mm, m, vh_exact)

  ! This builds a normalized Gaussian charge
  alpha = M_TWO
  beta  = M_ONE / ( alpha**3 * sqrt(M_PI)**3 )
  rr(1:3) = (/ M_TWO*m%h(1), -CNST(5.0)*m%h(2), CNST(10.0)*m%h(3) /)
!!$  rr(1:3) = (/ M_ZERO, M_ZERO, M_ZERO /)
  do i = 1, m%np
     call mesh_r(m, i, r, x = x, a = rr)
     f(i) = beta*exp(-(r/alpha)**2)
  enddo
  write(*, '(a,f14.6)') 'Total charge of the Gaussian distribution', dmf_integrate(m, f)

  ! And this builds the corresponding *exact* potential due to the Gaussian charge
  do i = 1, m%np
     call mesh_r(m, i, r, x = x, a = rr)
     if(r > r_small) then
        vh_exact(i) = loct_erf(r/alpha)/r
     else
        vh_exact(i) = (M_TWO/sqrt(M_PI))/alpha
     endif
  enddo

  call poisson_solve(m, f_der, vh_numerical, f)  

  write(*, '(a,e14.4)') 'Error: ', dmf_nrm2(m, vh_numerical - vh_exact)

  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), ".", "rho", m, f, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), ".", "vh_exact", m, vh_exact, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), ".", "vh_numerical", m, vh_numerical, M_ONE, ierr)

  contains

  FLOAT function isubl(l, x)
    integer, intent(in) :: l
    FLOAT, intent(in) :: x
!!$    isubl = ( sqrt(M_PI)/M_TWO ) * loct_erf(x)
    isubl = M_HALF*loct_gamma(l + M_HALF)*(M_ONE - loct_incomplete_gamma(l+M_HALF, x**2) )
  end function isubl

  subroutine vhexact_nlm(alpha, ll, mm, m, vh)
    FLOAT, intent(in) :: alpha
    integer, intent(in) :: ll, mm
    type(mesh_type), intent(in) :: m
    FLOAT, intent(out) :: vh(:)

    integer :: i, j, lldfac
    FLOAT   :: r, x(3), ylm, beta, gamma


    lldfac = 1; do j = 1, 2*ll+1, 2; lldfac = lldfac * j; enddo

    beta = (2**(ll+2))/( alpha**(2*ll+3) * sqrt(M_PI) * lldfac )
    gamma = ( sqrt(M_PI)*2**(ll+3) ) / lldfac
    do i = 1, m%np
       call mesh_r(m, i, r, x = x)
       ylm  = loct_ylm(x(1), x(2), x(3), ll, mm)
       if(r .ne. M_ZERO) then
         vh(i) = gamma * isubl( ll, r/alpha) * ylm / r**(ll+1)
       else
         vh(i) = gamma * ylm / alpha
       endif
    enddo
  end subroutine vhexact_nlm

end program oct_test
