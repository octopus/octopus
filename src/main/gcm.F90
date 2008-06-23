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
  use grid_m
  use system_m
  use hamiltonian_m
  use states_m
  use excited_states_m
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
    type(grid_t), pointer :: gr
    type(states_t) :: opst
    FLOAT, allocatable :: etot(:)
    FLOAT, allocatable :: hpsi(:, :)
    FLOAT, allocatable :: rho(:), vh(:)
    FLOAT, allocatable :: hmatrix(:, :), smatrix(:, :), overlap_matrix(:, :, :)
    character(len=100), allocatable :: slatdetnames(:)
    FLOAT :: uh, overlap, kij, gamma, emin, ex
    integer :: ierr, i, ndeterminants, j, k, n

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
     ndeterminants = loct_parse_block_n(blk)
     ALLOCATE(slatdetnames(ndeterminants), ndeterminants)
     do i = 0, ndeterminants - 1
       call loct_parse_block_string(blk, i, 0, slatdetnames(i+1))
     end do
     call loct_parse_block_end(blk)
    end if

    ALLOCATE(phi(ndeterminants), ndeterminants)
    ALLOCATE(etot(ndeterminants), ndeterminants)
    ALLOCATE(smatrix(ndeterminants, ndeterminants), ndeterminants**2)
    ALLOCATE(hmatrix(ndeterminants, ndeterminants), ndeterminants**2)
    gr => sys%gr

    ! Copy the basic structure in sys%st to all of the members of phi.
    do i = 1, ndeterminants
      call states_copy(phi(i), sys%st)
      call states_allocate_wfns(phi(i), gr%m)
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
      call restart_read (trim(slatdetnames(i)), phi(i), gr, sys%geo, ierr)
    end do
    call messages_print_stress(stdout)


    ALLOCATE(rho(NP), NP)
    ALLOCATE(vh(NP), NP)
    rho = M_ZERO
    vh  = M_ZERO

    ! Calculate the total energies for each of the Slater determinants
    do i = 1, ndeterminants
      ! The total density may be needed.
      call states_calc_dens(phi(i), NP, phi(i)%rho)

      ! First, the one-body part of the total energy:
      etot(i) = delectronic_kinetic_energy(h, gr, phi(i)) + &
                delectronic_external_energy(h, gr, phi(i))

      ! Coulomb contribution.
      do j = 1, NP
        rho(j) = phi(i)%rho(j, 1)
      end do
      call dpoisson_solve(gr, vh, rho)
      uh = M_HALF*dmf_integrate(gr%m, vh * rho)

      !Exchange contribution
      if(phi(i)%nst > 1) then
        ex = M_ZERO
        do j = 1, phi(i)%nst
          do k = 1, phi(i)%nst
            do n = 1, NP
              rho(n) = phi(i)%dpsi(n, 1, j, 1)*phi(i)%dpsi(n, 1, k, 1)
            end do
            call dpoisson_solve(gr, vh, rho)
            ex = ex - dmf_integrate(gr%m, vh*rho)
          end do
        end do
      else
        ex = - M_HALF * uh
      end if

      etot(i) = etot(i) + uh + ex + h%ep%eii
      hmatrix(i, i) = etot(i)
    end do

    call messages_print_stress(stdout, 'Total energies of the Slater determinants')
    do i = 1, ndeterminants
      write(message(1),'(a,i2.2,a,f20.8,a)') &
        'Etot(',i,') = ', etot(i)/units_out%energy%factor, ' ['//trim(units_out%energy%abbrev)//']'
      call write_info(1)
    end do
    call messages_print_stress(stdout)


    ALLOCATE(hpsi(NP_PART, phi(1)%d%dim), NP_PART*phi(1)%d%dim)


    ! Calculate the cross terms
    do i = 1, ndeterminants
      do j = i + 1, ndeterminants

        ALLOCATE(overlap_matrix(phi(i)%nst, phi(i)%nst, 1), phi(i)%nst*phi(i)%nst)

        call states_copy(opst, phi(j))
        do k = 1, phi(j)%nst
          opst%dpsi(:, :, k, 1) = M_ZERO
          call dkinetic (h, gr, phi(j)%dpsi(:, :, k, 1), opst%dpsi(:, :, k, 1))
          call dvexternal (h, gr, phi(j)%dpsi(:, :, k, 1), opst%dpsi(:, :, k, 1), 1)
        end do
        kij = dstates_mpmatrixelement(gr%m, phi(i), phi(j), opst)
        call states_end(opst)

        call dstates_matrix(gr%m, phi(i), phi(j), overlap_matrix)
        smatrix(i, j) = dstates_mpdotp(gr%m, phi(i), phi(j), overlap_matrix)

        rho = M_ZERO
        vh = M_ZERO
        do k = 1, NP
          rho(k) = phi(i)%dpsi(k, 1, 1, 1) * phi(j)%dpsi(k, 1, 1, 1)
        end do
        call dpoisson_solve(gr, vh, rho)
        uh =  dmf_integrate(gr%m, vh(:) * rho(:))

        gamma = kij + uh + h%ep%eii
        hmatrix(i, j) = gamma

        deallocate(overlap_matrix)
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

