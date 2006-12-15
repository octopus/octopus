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

#if defined(HAVE_MPI)
  use mpi_m
  use par_vec_m
#endif

  implicit none

  private

  public ::               &
       poisson_isf_init,  &
       poisson_isf_solve, & 
       poisson_isf_end

  type(dcf_t)        :: rho_cf
  FLOAT, pointer     :: kernel(:, :, :)
  integer, parameter :: order_scaling_function = 8 
  integer            :: nfft1, nfft2, nfft3

contains

  subroutine poisson_isf_init(m)
    type(mesh_t), intent(in) :: m

    integer :: n1, n2, n3
#if defined(HAVE_MPI)
    integer :: m1, m2, m3, md1, md2, md3
#endif

    call push_sub('poisson_isf.poisson_isf_init')

    call dcf_new(m%l, rho_cf)

#if defined(HAVE_MPI)
    call par_calculate_dimensions(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      m1, m2, m3, n1, n2, n3, md1, md2, md3, nfft1, nfft2, nfft3, mpi_world%size)

    ALLOCATE(kernel(nfft1, nfft2, nfft3/mpi_world%size), nfft1*nfft2*nfft3/mpi_world%size)

    call par_build_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3),       &
      n1, n2, n3, nfft1, nfft2, nfft3, m%h(1), order_scaling_function, &
      mpi_world%rank, mpi_world%size, mpi_world%comm, kernel)
#else
    call calculate_dimensions(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      nfft1, nfft2, nfft3)

    n1 = nfft1/2 + 1
    n2 = nfft2/2 + 1
    n3 = nfft3/2 + 1

    ALLOCATE(kernel(n1, n2, n3), n1*n2*n3)

    call build_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      nfft1, nfft2, nfft3, m%h(1), order_scaling_function, kernel)
#endif

    call pop_sub()
  end subroutine poisson_isf_init

  subroutine poisson_isf_solve(m, pot, rho)
    type(mesh_t), intent(in) :: m
    FLOAT, intent(out)       :: pot(:)
    FLOAT, intent(in)        :: rho(:)

    FLOAT, allocatable :: rho_global(:)
    FLOAT, allocatable :: pot_global(:)

    call push_sub('poisson_isf.poisson_isf_solve')

    call dcf_alloc_RS(rho_cf)

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      ALLOCATE(rho_global(m%np_global), m%np_global)
      ALLOCATE(pot_global(m%np_global), m%np_global)

      ! At this point, dvec_allgather is required because the ISF solver
      ! uses another data distribution algorithm than for the mesh functions
      ! for which every node requires the full data.
      call dvec_allgather(m%vp, rho_global, rho)
      call dmf2cf(m, rho_global, rho_cf)
#endif
    else
      call dmf2cf(m, rho, rho_cf)
    end if

#if defined(HAVE_MPI)
    call par_psolver_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      nfft1, nfft2, nfft3, m%h(1), kernel, rho_cf%RS,              &
      mpi_world%rank, mpi_world%size, mpi_world%comm)
#else   
    call psolver_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      nfft1, nfft2, nfft3, m%h(1), kernel, rho_cf%RS)
#endif

    if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
      call dcf2mf(m, rho_cf, pot_global)
      call dvec_scatter(m%vp, pot_global, pot)
      deallocate(rho_global, pot_global)
#endif
    else
      call dcf2mf(m, rho_cf, pot)
    end if

    call dcf_free_RS(rho_cf)
    
    call pop_sub()
  end subroutine poisson_isf_solve

  subroutine poisson_isf_end()
    call push_sub('poisson_isf.poisson_isf_end')

    deallocate(kernel)
    call dcf_free(rho_cf)

    call pop_sub()
  end subroutine poisson_isf_end

end module poisson_isf_m
