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
  use density_m
  use energy_m
  use messages_m
  use parser_m
  use unit_m
  use unit_system_m
  use grid_m
  use system_m
  use hamiltonian_m
  use hamiltonian_base_m
  use states_m
  use states_dim_m
  use states_calc_m
  use excited_states_m
  use restart_m
  use poisson_m
  use profiling_m
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
  subroutine gcm_run(sys, hm)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm

    type(states_t), allocatable :: phi(:)
    type(block_t) :: blk
    type(grid_t), pointer :: gr
    type(states_t) :: opst
    FLOAT, allocatable :: etot(:)
    FLOAT, allocatable :: hpsi(:, :)
    FLOAT, allocatable :: rho(:), vh(:)
    FLOAT, allocatable :: hmatrix(:, :), smatrix(:, :), overlap_matrix(:, :, :)
    character(len=100), allocatable :: slatdetnames(:)
    FLOAT :: uh, kij, gamma, ex
    integer :: ierr, i, ndeterminants, j, k, n

    PUSH_SUB(gcm_run)

    call messages_experimental('Generator coordinates method')

    !%Variable GCMSlaterDeterminants
    !%Type block
    !%Section Calculation Modes::Generator Coordinates
    !%Description
    !%
    !%End
    if(parse_block(datasets_check('GCMSlaterDeterminants'), blk) .ne. 0) then
      write(message(1),'(a)') 'If you run in "CalculationMode = gcm" mode, then you must'
      write(message(2),'(a)') 'supply also a "GCMSlaterDeterminants" block.'
      call messages_fatal(2)
    else
     ndeterminants = parse_block_n(blk)
     SAFE_ALLOCATE(slatdetnames(1:ndeterminants))
     do i = 0, ndeterminants - 1
       call parse_block_string(blk, i, 0, slatdetnames(i+1))
     end do
     call parse_block_end(blk)
    end if

    SAFE_ALLOCATE(    phi(1:ndeterminants))
    SAFE_ALLOCATE(   etot(1:ndeterminants))
    SAFE_ALLOCATE(smatrix(1:ndeterminants, 1:ndeterminants))
    SAFE_ALLOCATE(hmatrix(1:ndeterminants, 1:ndeterminants))
    gr => sys%gr

    ! Copy the basic structure in sys%st to all of the members of phi.
    do i = 1, ndeterminants
      call states_copy(phi(i), sys%st)
      call states_allocate_wfns(phi(i), gr%mesh)
      SAFE_ALLOCATE(phi(i)%eigenval(1:phi(i)%nst, 1:phi(i)%d%nik))
      SAFE_ALLOCATE(phi(i)%occ(1:phi(i)%nst, 1:phi(i)%d%nik))
      if(phi(i)%d%ispin == SPINORS) then
        SAFE_ALLOCATE(phi(i)%spin(1:3, 1:phi(i)%nst, 1:phi(i)%d%nik))
        phi(i)%spin = M_ZERO
      end if
      phi(i)%eigenval = huge(phi(i)%eigenval)
      phi(i)%occ      = M_ZERO
    end do

    call messages_print_stress(stdout, 'Reading Slater determinants. ')
    ! Read each of the Slater determinants.
    do i = 1, ndeterminants
      call restart_read (trim(slatdetnames(i)), phi(i), gr, sys%geo, ierr)
    end do
    call messages_print_stress(stdout)


    SAFE_ALLOCATE(rho(1:gr%mesh%np))
    SAFE_ALLOCATE( vh(1:gr%mesh%np))
    rho = M_ZERO
    vh  = M_ZERO

    ! Calculate the total energies for each of the Slater determinants
    do i = 1, ndeterminants
      ! The total density may be needed.
      call density_calc(phi(i), gr, phi(i)%rho)

      ! First, the one-body part of the total energy:
      etot(i) = delectronic_kinetic_energy(hm, gr, phi(i)) + &
                delectronic_external_energy(hm, gr, phi(i))

      ! Coulomb contribution.
      do j = 1, gr%mesh%np
        rho(j) = phi(i)%rho(j, 1)
      end do
      call dpoisson_solve(psolver, vh, rho)
      uh = M_HALF*dmf_integrate(gr%mesh, vh * rho)

      !Exchange contribution
      if(phi(i)%nst > 1) then
        ex = M_ZERO
        do j = 1, phi(i)%nst
          do k = 1, phi(i)%nst
            do n = 1, gr%mesh%np
              rho(n) = phi(i)%dpsi(n, 1, j, 1)*phi(i)%dpsi(n, 1, k, 1)
            end do
            call dpoisson_solve(psolver, vh, rho)
            ex = ex - dmf_integrate(gr%mesh, vh*rho)
          end do
        end do
      else
        ex = - M_HALF * uh
      end if

      etot(i) = etot(i) + uh + ex + hm%ep%eii
      hmatrix(i, i) = etot(i)
    end do

    call messages_print_stress(stdout, 'Total energies of the Slater determinants')
    do i = 1, ndeterminants
      write(message(1),'(a,i2.2,a,f20.8,a)') &
        'Etot(',i,') = ', units_from_atomic(units_out%energy, etot(i)), ' ['//trim(units_abbrev(units_out%energy))//']'
      call messages_info(1)
    end do
    call messages_print_stress(stdout)

    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part, 1:phi(1)%d%dim))

    ! Calculate the cross terms
    do i = 1, ndeterminants
      do j = i + 1, ndeterminants

        SAFE_ALLOCATE(overlap_matrix(1:phi(i)%nst, 1:phi(i)%nst, 1:1))

        call states_copy(opst, phi(j))
        do k = 1, phi(j)%nst
          opst%dpsi(:, :, k, 1) = M_ZERO
          call dhamiltonian_apply(hm, gr%der, phi(j)%dpsi(:, :, k, 1), opst%dpsi(:, :, k, 1), ist = k, ik = 1, &
            terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        end do
        kij = dstates_mpmatrixelement(gr%mesh, phi(i), phi(j), opst)
        call states_end(opst)

        call dstates_matrix(gr%mesh, phi(i), phi(j), overlap_matrix)
        smatrix(i, j) = dstates_mpdotp(gr%mesh, phi(i), phi(j), overlap_matrix)

        uh =  mpdotp_twobody(phi(i), phi(j), gr)

        gamma = kij + uh + hm%ep%eii
        hmatrix(i, j) = gamma

        SAFE_DEALLOCATE_A(overlap_matrix)
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
      call messages_info(1)
    end do
    call messages_print_stress(stdout)

    call messages_print_stress(stdout, 'S matrix')
    do i = 1, ndeterminants
      write(message(1),'(10f18.8)') (smatrix(i, j), j=1, ndeterminants)
      call messages_info(1)
    end do
    call messages_print_stress(stdout)

    call  lalg_geneigensolve(ndeterminants, hmatrix, smatrix, etot)

    call messages_print_stress(stdout, 'Eigenvalues')
    do i = 1, ndeterminants
      write(message(1), '(a,i2.2,a,f20.8,a)') &
        'E(',i,') = ', units_from_atomic(units_out%energy, etot(i)), ' ['//trim(units_abbrev(units_out%energy))//']'
      call messages_info(1)
    end do
    call messages_print_stress(stdout)


    ! Clean up.
    do i = 1, ndeterminants
      call states_end(phi(i))
    end do
    SAFE_DEALLOCATE_A(phi)
    SAFE_DEALLOCATE_A(etot)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(vh)
    SAFE_DEALLOCATE_A(smatrix)
    SAFE_DEALLOCATE_A(hmatrix)
    SAFE_DEALLOCATE_A(slatdetnames)
    POP_SUB(gcm_run)
  end subroutine gcm_run
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function mpdotp_twobody(st1, st2, gr) result (st1opst2)
    type(states_t), intent(in)    :: st1
    type(states_t), intent(in)    :: st2
    type(grid_t),   intent(inout) :: gr

    integer :: nst, k, k1, k2, l1, l2, i
    FLOAT, allocatable :: rho(:), vh(:)
    FLOAT, allocatable :: mat(:, :, :, :)

    PUSH_SUB(mpdotp_twobody)

    nst = st1%nst

    select case(nst)

    case(1)

      SAFE_ALLOCATE(rho(1:gr%mesh%np))
      SAFE_ALLOCATE( vh(1:gr%mesh%np))

      rho = M_ZERO
      vh = M_ZERO
      do k = 1, gr%mesh%np
        rho(k) = st1%dpsi(k, 1, 1, 1) * st2%dpsi(k, 1, 1, 1)
      end do
      call dpoisson_solve(psolver, vh, rho)
      st1opst2 =  dmf_integrate(gr%mesh, vh(:) * rho(:))

      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(vh)

    case default
      write(message(1),'(a)') 'GCM mode does not handle yet systems with more than one orbital.'
      call messages_fatal(1)

      SAFE_ALLOCATE(mat(1:nst, 1:nst, 1:nst, 1:nst))
      SAFE_ALLOCATE(rho(1:gr%mesh%np))
      SAFE_ALLOCATE( vh(1:gr%mesh%np))

      ! Build the matrix <k1 k2 | 1/|r1-r2| | l1 l2>
      do k1 = 1, nst
        do l1 = 1, nst
          do i = 1, gr%mesh%np
            rho(i) = st1%dpsi(i, 1, k1, 1)*st2%dpsi(i, 1, l1, 1)
          end do
          call dpoisson_solve(psolver, vh, rho)
          do k2 = 1, nst
            do l2 = 1, nst
              do i = 1, gr%mesh%np
                rho(i) = st1%dpsi(i, 1, k2, 1)*st2%dpsi(i, 1, l2, 1)
              end do
              mat(k1, l1, k2, l2) = dmf_integrate(gr%mesh, vh(:)*rho(:))
            end do
          end do

        end do
      end do

      SAFE_DEALLOCATE_A(mat)
      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(vh)
    end select

    POP_SUB(mpdotp_twobody)
  end function mpdotp_twobody
  ! ---------------------------------------------------------



end module gcm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
