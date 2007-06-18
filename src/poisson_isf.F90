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
!! $Id: poisson.F90 2544 2006-11-03 17:41:04Z xavier $

#include "global.h"

module poisson_isf_m
  use cube_function_m
  use global_m
  use messages_m
  use mesh_m
  use functions_m
  use mpi_m
  use par_vec_m

  implicit none

  private

  public ::               &
       poisson_isf_init,  &
       poisson_isf_solve, & 
       poisson_isf_end

  ! Datatype to store kernel values to solve poisson equation
  ! on different communicators (configurations).
  type isf_cnf_t
    FLOAT, pointer  :: kernel(:, :, :)
    integer         :: nfft1, nfft2, nfft3
    type(mpi_grp_t) :: mpi_grp
  end type isf_cnf_t

  ! Indices for the cnf array
  integer, parameter :: serial = 1
  integer, parameter :: world = 2
  integer, parameter :: domain = 3
  integer, parameter :: n_cnf = 3

  type(dcf_t)                  :: rho_cf
  type(isf_cnf_t)              :: cnf(1:3)
  integer, parameter           :: order_scaling_function = 8 

contains
  ! ---------------------------------------------------------
  subroutine poisson_isf_init(m)
    type(mesh_t), intent(in) :: m

    integer :: n1, n2, n3
#if defined(HAVE_MPI)
    integer :: m1, m2, m3, md1, md2, md3
    integer :: n(3)
    integer :: i_cnf
#endif

    call push_sub('poisson_isf.poisson_isf_init')

    call dcf_new(m%l, rho_cf)

    ! The serial version is always needed (as used, e.g., in the casida runmode)
    call calculate_dimensions(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
      cnf(serial)%nfft1, cnf(serial)%nfft2, cnf(serial)%nfft3)

    n1 = cnf(serial)%nfft1/2 + 1
    n2 = cnf(serial)%nfft2/2 + 1
    n3 = cnf(serial)%nfft3/2 + 1

    ALLOCATE(cnf(serial)%kernel(n1, n2, n3), n1*n2*n3)

    call build_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3),   &
      cnf(serial)%nfft1, cnf(serial)%nfft2, cnf(serial)%nfft3, &
      m%h(1), order_scaling_function, cnf(serial)%kernel)

#if defined(HAVE_MPI)
    ! Allocate to configurations. The initialisation, especially the kernel,
    ! depends on the number of nodes used for the calculations. To avoid
    ! recalculating the kernel on each call of poisson_isf_solve depending on
    ! the all_nodes argument, both kernels are calculated.
    cnf(world)%mpi_grp = mpi_world
    cnf(domain)%mpi_grp = m%mpi_grp

    ! Build the kernel for all configurations. At the moment, this is
    ! solving the poisson equation with all nodes (i_cnf == world) and
    ! with the domain nodes only (i_cnf == domain).
    do i_cnf = 2, n_cnf
      call par_calculate_dimensions(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3),         &
        m1, m2, m3, n1, n2, n3, md1, md2, md3, cnf(i_cnf)%nfft1, cnf(i_cnf)%nfft2, &
        cnf(i_cnf)%nfft3, cnf(i_cnf)%mpi_grp%size)

      ! Shortcuts to avoid to "line too long" errors.
      n(1) = cnf(i_cnf)%nfft1
      n(2) = cnf(i_cnf)%nfft2
      n(3) = cnf(i_cnf)%nfft3

      ALLOCATE(cnf(i_cnf)%kernel(n(1), n(2), n(3)/cnf(i_cnf)%mpi_grp%size), n(1)*n(2)*n(3) / cnf(i_cnf)%mpi_grp%size)

    call par_build_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), n1, n2, n3,     &
      cnf(i_cnf)%nfft1, cnf(i_cnf)%nfft2, cnf(i_cnf)%nfft3,                      &
      m%h(1), order_scaling_function,                                            &
      cnf(i_cnf)%mpi_grp%rank, cnf(i_cnf)%mpi_grp%size, cnf(i_cnf)%mpi_grp%comm, &
      cnf(i_cnf)%kernel)
    end do
#endif

    call pop_sub()
  end subroutine poisson_isf_init

  ! ---------------------------------------------------------
  subroutine poisson_isf_solve(m, pot, rho, all_nodes)
    type(mesh_t), intent(in)  :: m
    FLOAT,        intent(out) :: pot(:)
    FLOAT,        intent(in)  :: rho(:)
    logical,      intent(in)  :: all_nodes

    integer :: i_cnf
#if defined(HAVE_MPI)
    FLOAT, allocatable :: rho_global(:)
    FLOAT, allocatable :: pot_global(:)
#endif

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

    ! Choose configuration.
    i_cnf = serial

#if defined(HAVE_MPI)
    if(all_nodes) then
      i_cnf = world
    else if(m%parallel_in_domains) then
      i_cnf = domain
    end if
#endif

#if !defined(HAVE_MPI)
    ASSERT(i_cnf == serial)
#endif

    if(i_cnf == serial) then
      call psolver_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3),    &
        cnf(serial)%nfft1, cnf(serial)%nfft2, cnf(serial)%nfft3, &
        m%h(1), cnf(serial)%kernel, rho_cf%RS)
 #if defined(HAVE_MPI)
    else
      call par_psolver_kernel(rho_cf%n(1), rho_cf%n(2), rho_cf%n(3), &
        cnf(i_cnf)%nfft1, cnf(i_cnf)%nfft2, cnf(i_cnf)%nfft3,  &
        m%h(1), cnf(i_cnf)%kernel, rho_cf%RS,                      &
        cnf(i_cnf)%mpi_grp%rank, cnf(i_cnf)%mpi_grp%size, cnf(i_cnf)%mpi_grp%comm)
#endif
    endif

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

  ! ---------------------------------------------------------
  subroutine poisson_isf_end()
#if defined(HAVE_MPI)
    integer :: i_cnf
#endif

    call push_sub('poisson_isf.poisson_isf_end')

#if defined(HAVE_MPI)
    do i_cnf = 1, n_cnf
      deallocate(cnf(i_cnf)%kernel)
    end do
#else
    deallocate(cnf(serial)%kernel)
#endif

    call dcf_free(rho_cf)

    call pop_sub()
  end subroutine poisson_isf_end

end module poisson_isf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
