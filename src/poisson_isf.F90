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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: poisson.F90 2544 2006-11-03 17:41:04Z xavier $

#include "global.h"

module poisson_isf_m
  use cube_function_m
  use global_m
  use messages_m
  use mesh_m
  use functions_m

  implicit none

  private

  public :: &
       poisson_isf_init, &
       poisson_isf_solve, & 
       poisson_isf_end

  type(dcf_t) :: rho_cf
  FLOAT, pointer :: karray(:, :, :)
  integer, parameter :: itype_scf = 8 
  integer :: nfft1, nfft2, nfft3

contains

  subroutine poisson_isf_init(m)
    type(mesh_t), intent(in) :: m
    integer :: n1k, n2k, n3k

    call push_sub('poisson_isf.poisson_isf_init')

    call dcf_new(m%l, rho_cf)    ! allocate cube function

    call Dimensions_FFT(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
         nfft1, nfft2, nfft3)

    n1k = nfft1/2 + 1
    n2k = nfft2/2 + 1
    n3k = nfft3/2 + 1

    ALLOCATE(karray(n1k, n2k, n3k), n1k * n2k * n3k)

    call Build_Kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
         nfft1, nfft2, nfft3, & 
         m%h(1), itype_scf, karray)
    
    call pop_sub()

  end subroutine poisson_isf_init

  subroutine poisson_isf_solve(m, pot, rho)
    type(mesh_t), intent(in) :: m
    FLOAT, intent(out) :: pot(:) ! pot(m%np)
    FLOAT, intent(in) :: rho(:) ! rho(m%np)

    FLOAT :: ehartree

    call push_sub('poisson_isf.poisson_isf_solve')

    !allocate cubic grids
    call dcf_alloc_RS(rho_cf)

    !put the density in the cubic grid
    call dmf2cf(m, rho, rho_cf)

    call PSolver_Kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
         nfft1, nfft2, nfft3,&
         m%h(1), karray, rho_cf%RS, ehartree)

    !recover the potential to our mesh (it is returned in the same array)
    call dcf2mf(m, rho_cf, pot)
    
    !deallocate the cubic meshes
    call dcf_free_RS(rho_cf)
    
    call pop_sub()

  end subroutine poisson_isf_solve

  subroutine poisson_isf_end()

    call push_sub('poisson_isf.poisson_isf_init')

    deallocate(karray)
    call dcf_free(rho_cf)

    call pop_sub()

  end subroutine poisson_isf_end

end module poisson_isf_m
