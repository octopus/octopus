!! Copyright (C) 2002-2007 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: geom_opt.F90 3970 2008-03-29 11:38:27Z acastro $

#include "global.h"

module gcm_m
  use global_m
  use datasets_m
  use messages_m
  use loct_parser_m
  use units_m
  use system_m
  use hamiltonian_m
  use states_m
  use restart_m
  use poisson_m
  use mesh_function_m
  use lalg_adv_m
  use lalg_basic_m


  implicit none

  private
  public :: &
    gcm_run


  contains

  ! ---------------------------------------------------------
  ! Very preliminary implementation of the generator coodinates DFT scheme
  ! proposed by K. Capelle [K. Capelle, J. Chem. Phys. 119, 1285 (2003).
  ! ---------------------------------------------------------
  subroutine gcm_run(sys, h)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h

    type(states_t), allocatable :: phi(:)
    type(block_t) :: blk
    FLOAT, allocatable :: etot(:)
    FLOAT, allocatable :: hpsi(:, :)
    FLOAT, allocatable :: rho(:), vh(:)
    FLOAT, allocatable :: hmatrix(:, :), smatrix(:, :)
    character(len=100), allocatable :: slatdetnames(:)
    FLOAT :: uh, overlap, kij, gamma, emin
    integer :: ierr, i, ndeterminants, j, k

    call push_sub('gcm.gcm_run')

    !%Variable GCMSlaterDeterminants
    !%Type block
    !%Section Generator Coordinates
    !%Description
    !%
    !%End
    if(loct_parse_block(check_inp('GCMSlaterDeterminants'), blk) .ne. 0) then
      write(message(1),'(a)') 'If you run in "CalculationMode = gcm" mode, then you must'
      write(message(2),'(a)') 'supply also a "GCMSlaterDeterminants" block.'
      call write_fatal(2)
    else
     ndeterminants = loct_parse_block_cols(blk, 0)
     ALLOCATE(slatdetnames(ndeterminants), ndeterminants)
     do i = 0, ndeterminants - 1
       call loct_parse_block_string(blk, 0, i, slatdetnames(i+1))
     end do
     call loct_parse_block_end(blk)
    end if

    ALLOCATE(phi(ndeterminants), ndeterminants)
    ALLOCATE(etot(ndeterminants), ndeterminants)
    ALLOCATE(smatrix(ndeterminants, ndeterminants), ndeterminants**2)
    ALLOCATE(hmatrix(ndeterminants, ndeterminants), ndeterminants**2)

    ! Copy the basic structure in sys%st to all of the members of phi.
    do i = 1, ndeterminants
      call states_copy(phi(i), sys%st)
      call states_allocate_wfns(phi(i), sys%gr%m)
      ALLOCATE(phi(i)%eigenval(phi(i)%nst, phi(i)%d%nik), phi(i)%nst*phi(i)%d%nik)
      ALLOCATE(phi(i)%momentum(3, phi(i)%nst, phi(i)%d%nik), phi(i)%nst*phi(i)%d%nik)
      ALLOCATE(phi(i)%occ(phi(i)%nst, phi(i)%d%nik), phi(i)%nst*phi(i)%d%nik)
      if(phi(i)%d%ispin == SPINORS) then
        ALLOCATE(phi(i)%spin(3, phi(i)%nst, phi(i)%d%nik), phi(i)%nst*phi(i)%d%nik*2)
        phi(i)%spin = M_ZERO
      end if
      phi(i)%eigenval = huge(REAL_PRECISION)
      phi(i)%occ      = M_ZERO
    end do

    call messages_print_stress(stdout, 'Reading Slater determinants. ')
    ! Read each of the Slater determinants.
    do i = 1, ndeterminants
      call restart_read (trim(slatdetnames(i)), phi(i), sys%gr, sys%geo, ierr)
    end do
    call messages_print_stress(stdout)


    ALLOCATE(rho(sys%gr%m%np), sys%gr%m%np)
    ALLOCATE(vh(sys%gr%m%np), sys%gr%m%np)
    rho = M_ZERO
    vh  = M_ZERO

    ! Calculate the total energies for each of the Slater determinants
    do i = 1, ndeterminants
      call states_calc_dens(phi(i), sys%gr%m%np, phi(i)%rho)
      do j = 1, sys%gr%m%np
        rho(j) = phi(i)%rho(j, 1)
      end do
      call dpoisson_solve(sys%gr, vh, rho)
      uh = M_HALF*dmf_integrate(sys%gr%m, vh * rho)
      etot(i) = delectronic_kinetic_energy(h, sys%gr, phi(i)) + &
                delectronic_external_energy(h, sys%gr, phi(i))
      !write(0, '(a,i1,a,f18.10)') 'One body: (', i, ')', etot(i) / units_out%energy%factor
      etot(i) = etot(i) + M_HALF*uh + h%ep%eii
      hmatrix(i, i) = etot(i)
      !write(0, '(a,i1,a,f18.10)') 'etot(',i,') = ', etot(i) / units_out%energy%factor
    end do

    call messages_print_stress(stdout, 'Total energies of the Slater determinants')
    do i = 1, ndeterminants
      write(message(1),'(a,i2.2,a,f20.8,a)') &
        'Etot(',i,') = ', etot(i)/units_out%energy%factor, ' ['//trim(units_out%energy%abbrev)//']'
      call write_info(1)
    end do
    call messages_print_stress(stdout)


    ALLOCATE(hpsi(sys%gr%m%np_part, phi(1)%d%dim), sys%gr%m%np_part*phi(1)%d%dim)


    ! Calculate the cross terms
    do i = 1, ndeterminants
      do j = i + 1, ndeterminants

        overlap = dstates_dotp(sys%gr%m, 1, phi(i)%dpsi(:, :, 1, 1), phi(j)%dpsi(:, :, 1, 1))
        smatrix(i, j) = overlap**2

        hpsi = M_ZERO

        call dkinetic (h, sys%gr, phi(j)%dpsi(:, :, 1, 1), hpsi)
        call dvexternal (h, sys%gr, phi(j)%dpsi(:, :, 1, 1), hpsi, 1)
        kij = dstates_dotp(sys%gr%m, 1, phi(i)%dpsi(:, :, 1, 1), hpsi(:, :))
        !write(0, *) 'kij, overlap = ', kij, overlap

        rho = M_ZERO
        vh = M_ZERO
        do k = 1, sys%gr%m%np
          rho(k) = phi(i)%dpsi(k, 1, 1, 1) * phi(j)%dpsi(k, 1, 1, 1)
        end do
        call dpoisson_solve(sys%gr, vh, rho)
        uh =  dmf_integrate(sys%gr%m, vh(:) * rho(:))
        !write(0, *) 'One body: ', M_TWO * overlap * kij
        !write(0, *) 'Two body: ', uh

        gamma = M_TWO * overlap * kij + uh + h%ep%eii
        hmatrix(i, j) = gamma

      end do 
    end do

    do i = 1, ndeterminants
      smatrix(i, i) = M_ONE
      do j = 1, i-1
        smatrix(i, j) = smatrix(j, i)
        hmatrix(i, j) = hmatrix(j, i)
      end do
    end do

    call messages_print_stress(stdout, 'H matrix')
    do i = 1, ndeterminants
      write(message(1),'(10f18.8)') (hmatrix(i, j), j=1, ndeterminants)
      call write_info(1)
    end do
    call messages_print_stress(stdout)

    call messages_print_stress(stdout, 'S matrix')
    do i = 1, ndeterminants
      write(message(1),'(10f18.8)') (smatrix(i, j), j=1, ndeterminants)
      call write_info(1)
    end do
    call messages_print_stress(stdout)

    call  lalg_geneigensolve(ndeterminants, hmatrix, smatrix, etot)

    call messages_print_stress(stdout, 'Eigenvalues')
    do i = 1, ndeterminants
      write(message(1), '(a,i2.2,a,f20.8,a)') &
        'E(',i,') = ', etot(i) / units_out%energy%factor, ' ['//trim(units_out%energy%abbrev)//']'
      call write_info(1)
    end do
    call messages_print_stress(stdout)


    ! Clean up.
    do i = 1, ndeterminants
      call states_end(phi(i))
    end do
    deallocate(phi, etot, hpsi, rho, vh, smatrix, hmatrix, slatdetnames)
    call pop_sub()
  end subroutine gcm_run


end module gcm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

