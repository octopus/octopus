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


!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tests a given Poisson solver by calculating the Hartree potential generated
! by a Gaussian distribution of charge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
program oct_test
  use global
  use units
  use lib_oct
  use fft
  use functions
  use mesh_function
  use mesh
  use derivatives
  use geometry
  use output
  use poisson

  implicit none

  type(geometry_type), pointer :: geo
  type(mesh_type), pointer     :: m
  type(f_der_type)             :: f_der
  type(output_type)            :: outp
  integer :: ierr, i
  FLOAT :: val_charge, def_h, def_rsize, alpha, beta, r, x(3), rr(3)
  FLOAT, allocatable :: f(:), g(:), vh_exact(:), vh_numerical(:), l2f(:), l2fn(:)

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

  allocate(f(m%np), g(m%np), vh_exact(m%np), vh_numerical(m%np), l2f(m%np), l2fn(m%np))
  f = M_ZERO
  g = M_ZERO
  vh_exact = M_ZERO
  vh_numerical = M_ZERO

  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), ".", "grid", m, f, M_ONE, ierr)

  ! This builds a normalized Gaussian charge
  alpha = M_ONE
  beta  = M_ONE / ( alpha**3 * sqrt(M_PI)**3 )
  rr(1:3) = (/ M_TWO*m%h(1), -CNST(1.0)*m%h(2), CNST(6.0)*m%h(3) /)
!!$  rr(1:3) = M_ZERO
  do i = 1, m%np
     call mesh_r(m, i, r, x = x, a = rr)
     f(i) = beta*exp(-(r/alpha)**2)
  enddo
  write(*, '(/,a,f14.6,/)') 'Total charge of the Gaussian distribution', dmf_integrate(m, f)

  ! This builds the Laplacian of the Gaussian charge distribution.
  do i = 1, m%np
     call mesh_r(m, i, r, x = x, a = rr)
     l2f(i) = (M_FOUR*r**2/alpha**4 - CNST(6.0)/alpha**2)*f(i)
  enddo

  ! This calculates the numerical Laplacian of the Gaussian Charge distribution.
  call dderivatives_lapl(f_der%der_discr, f, l2fn)

  call doutput_function(output_fill_how("AxisX"), ".", "rho", m, f, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX"), ".", "L2rho", m, l2f, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX"), ".", "L2rhon", m, l2fn, M_ONE, ierr)

  ! And this compares.
  write(*, '(/,a,f14.6,/)') 'Error in the Laplacian', dmf_nrm2(m, l2f-l2fn)

  ! And this builds the corresponding *exact* potential due to the Gaussian charge
  do i = 1, m%np
     call mesh_r(m, i, r, x = x, a = rr)
     if(r > r_small) then
        vh_exact(i) = loct_erf(r/alpha)/r
     else
        vh_exact(i) = (M_TWO/sqrt(M_PI))/alpha
     endif
  enddo

  write(*, *) 'Calling poisson_solve...', loct_clock()/1.0e6
  call dmf_random(m, vh_numerical)
  call poisson_solve(m, f_der, vh_numerical, f)
  write(*, *) 'Done.', loct_clock()/1.0e6

  write(*, *) 'Calling poisson_solve (again) ...', loct_clock()/1.0e6
  call poisson_solve(m, f_der, vh_numerical, f)
  write(*, *) 'Done.', loct_clock()/1.0e6

  write(*, '(/,a,e14.4,/)') 'Error: ', dmf_nrm2(m, vh_numerical - vh_exact)

  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), &
                        ".", "rho", m, f, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), &
                        ".", "vh_exact", m, vh_exact, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX.and.AxisZ.and.AxisY"), &
                        ".", "vh_numerical", m, vh_numerical, M_ONE, ierr)


  call global_end()
end program oct_test



