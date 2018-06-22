!!  Copyright (C) 2012 I. Theophilou, N. Helbig
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module rdmft_oct_m
  use density_oct_m
  use eigensolver_oct_m
  use energy_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use loct_math_oct_m 
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use minimizer_oct_m
  use output_oct_m
  use output_me_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_restart_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use xc_oep_oct_m
 
  implicit none

  private
  public ::                      & 
       assign_eigfunctions,      &
       rdm_t,                    &
       scf_occ,                  &
       scf_orb,                  &
       scf_rdmft

  type rdm_t
    type(states_t) :: psi
    integer  :: max_iter
    integer  :: iter
    integer  :: n_twoint !number of unique two electron integrals
    logical  :: do_basis
    FLOAT    :: mu, occsum, qtot, scale_f, toler, conv_ener, maxFO
    FLOAT, allocatable   :: eone(:), eone_int(:,:), twoint(:), hartree(:,:), exchange(:,:), evalues(:)   
    FLOAT, allocatable   :: vecnat(:,:), Coul(:,:,:), Exch(:,:,:) 
    integer, allocatable :: i_index(:), j_index(:), k_index(:), l_index(:) 

    !>shortcuts
    type(states_t),   pointer :: st
    type(grid_t),     pointer :: gr
  end type rdm_t
 
  type(rdm_t), save :: rdm
 
contains
   

  ! ----------------------------------------

  ! scf for the occupation numbers and the natural orbitals
  subroutine scf_rdmft(gr, geo, st, ks, hm, outp, max_iter, restart_dump)
    type(grid_t),  target, intent(inout) :: gr  !< grid
    type(geometry_t),      intent(inout) :: geo !< geometry
    type(states_t),target, intent(inout) :: st  !< States
    type(v_ks_t),          intent(inout) :: ks  !< Kohn-Sham
    type(hamiltonian_t),   intent(inout) :: hm  !< Hamiltonian
    type(output_t),        intent(in)    :: outp !< output
    integer,               intent(in)    :: max_iter
    type(restart_t),       intent(in)    :: restart_dump
    
    integer :: iter, icount, ip, ist, iatom, ierr, maxcount
    FLOAT :: energy, energy_dif, energy_old, energy_occ, xpos, xneg, sum_charge, rr, rel_ener
    FLOAT, allocatable :: species_charge_center(:), psi(:, :), stepsize(:)
    logical :: conv, gs_run_
    character(len=MAX_PATH_LEN) :: dirname    

    PUSH_SUB(scf_rdmft)

    gs_run_ = .true.

    if (hm%d%ispin /= 1) then
      call messages_not_implemented("RDMFT exchange function not yet implemented for spin_polarized or spinors")
    end if

    ! problem is about k-points for exchange
    if (simul_box_is_periodic(gr%sb)) then
      call messages_not_implemented("Periodic system calculations for RDMFT")
    end if

    ! exchange routine needs all states on each processor currently
    if(st%parallel_in_states) then
      call messages_not_implemented("RDMFT parallel in states")
    end if

    call rdmft_init() 
   

    !set initial values
    energy_old = CNST(1.0e20)
    xpos = M_ZERO 
    xneg = M_ZERO
    energy = M_ZERO 
    conv = .false.
    if(rdm%do_basis.eqv..false.) then
      !stepsize for steepest decent
      SAFE_ALLOCATE(stepsize(1:st%nst))
      stepsize = 0.1
      maxcount = 10
    else
      maxcount = 50
    endif
    
    !If using a basis set, localize the starting orbitals
    if(rdm%do_basis.eqv..true.) then
      !Find the charge center of they system  
      SAFE_ALLOCATE(species_charge_center(1:geo%space%dim))
      call geometry_dipole(geo,species_charge_center)
      sum_charge = M_ZERO
      do iatom = 1, geo%natoms
        sum_charge = sum_charge + species_zval(geo%atom(iatom)%species)
      end do
      species_charge_center = species_charge_center/(sum_charge*P_PROTON_CHARGE)
      !Localize orbitals
      SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
      do ist = 1, st%nst
        call states_get_state(st, gr%mesh, ist, 1, psi)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, species_charge_center)
          if (st%eigenval(ist, 1) < M_ZERO) then
            psi(ip, 1) = psi(ip, 1)&
              *exp((-(sqrt(-M_TWO*st%eigenval(int(st%qtot*M_HALF), 1))) + sqrt(-M_TWO*st%eigenval(ist, 1)))*rr)
          else
            psi(ip, 1) = psi(ip, 1)*exp(-(sqrt(-M_TWO*st%eigenval(int(st%qtot*M_HALF), 1)))*rr)
          end if
        end do
        call states_set_state(st, gr%mesh, ist, 1, psi)
      end do
      ! Orthogonalize the resulting orbitals
      call dstates_orthogonalization_full(st,gr%mesh,1)

      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(species_charge_center)
    endif
    
    ! Start the actual minimization, first step is minimization of occupation numbers
    ! Orbital minimization is according to Piris and Ugalde, Vol.13, No. 13, J. Comput. Chem. (scf_orb) or
    ! using steepest decent (scf_orb_direct)
    do iter = 1, max_iter
      write(message(1), '(a)') '**********************************************************************'
      write(message(2),'(a, i4)') 'RDM Iteration:', iter
      call messages_info(2)
      call scf_occ(rdm, gr, hm, st, energy_occ)
      ! Diagonalization of the generalized Fock matrix 
      write(message(1), '(a)') 'Optimization of natural orbitals'
      call messages_info(1)
      do icount = 1, maxcount !still under investigation how many iterations we need
        if (rdm%do_basis.eqv..false.) then
          call scf_orb_direct(rdm, gr, geo, st, ks, hm, stepsize, energy)
        else
          call scf_orb(rdm, gr, geo, st, ks, hm, energy)
        end if
        energy_dif = energy - energy_old
        energy_old = energy
        if (rdm%do_basis.eqv. .true.) then
          if (abs(energy_dif)/abs(energy).lt. rdm%conv_ener.and.rdm%maxFO < 1.d3*rdm%conv_ener)  exit
          if (energy_dif < M_ZERO) then 
            xneg = xneg + 1
          else
            xpos = xpos + 1
          end if
          if (xneg > CNST(1.5e0)*xpos) then
            rdm%scale_f = CNST(1.01)*rdm%scale_f
          elseif (xneg < CNST(1.1e0)*xpos) then 
            rdm%scale_f = CNST(0.95)* rdm%scale_f 
          end if
        endif !rdm%do_basis
        rdm%iter = rdm%iter + 1
      end do !icount
      xneg = M_ZERO
      xpos = M_ZERO
      
      rel_ener = abs(energy_occ-energy)/abs(energy)

      write(message(1),'(a,es15.5)') 'Total energy ', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
      call messages_info(1)
      if (rdm%do_basis.eqv..true.) then
        if ((rel_ener < rdm%conv_ener).and.rdm%maxFO < 1.e3*rdm%conv_ener) then
          conv = .true.
        end if
      else
        if (rel_ener < rdm%conv_ener) then
          conv = .true.
        end if
      endif
 
      if (rdm%toler > 1e-4) rdm%toler = rdm%toler*1e-1 !Is this still okay or does it restrict the possible convergence?

      ! save restart information
      if ( (conv .or. (modulo(iter, outp%restart_write_interval) == 0) &
        .or. iter == max_iter)) then
        call states_dump(restart_dump, st, gr, ierr, iter=iter) 
        if (ierr /= 0) then
          message(1) = 'Unable to write states wavefunctions.'
          call messages_warning(1)
        end if
      endif
      ! write output for iterations if requested
      if(outp%what/=0 .and. outp%duringscf .and. outp%output_interval /= 0 &
        .and. gs_run_ .and. mod(iter, outp%output_interval) == 0) then
        write(dirname,'(a,a,i4.4)') trim(outp%iter_dir),"scf.",iter
        call output_all(outp, gr, geo, st, hm, ks, dirname)
      end if

      if (conv) exit
    end do

    if(conv) then 
      ! output final information
      !if (rdm%do_basis.eqv..false.) then
        !call scf_orb(rdm, gr, geo, st, ks, hm, energy) !NH this call changes the energy again, needs to be replaced
        !write(message(1),'(a,es15.5)')  'The maximum non diagonal element of the matrix F formed by the Lagrange multiplyers is ', &
        !                                rdm%maxFO
        !call messages_info(1)
      !endif
      write(message(1),'(a,i3,a)')  'The calculation converged after ',iter,' iterations'
      write(message(2),'(a,es15.5)')  'The total energy is ', units_from_atomic(units_out%energy,energy + hm%ep%eii)
      call messages_info(2)
      if(gs_run_) then 
        !call scf_write_static(STATIC_DIR, "info") !NH can we turn this on?
        call output_all(outp, gr, geo, st, hm, ks, STATIC_DIR)
      end if
    else
      write(message(1),'(a,i3,a)')  'The calculation did not converge after ', iter, ' iterations '
      write(message(2),'(a,es15.5)') 'Relative energy difference between the last two iterations ', rel_ener
      write(message(3),'(a,es15.5)') 'The maximal non-diagonal element of the Hermitian matrix F is ', rdm%maxFO
      call messages_info(3)
    end if

    SAFE_DEALLOCATE_A(stepsize)

    call rdmft_end()
 
    POP_SUB(scf_rdmft) 

  contains

    ! ---------------------------------------------------------

    subroutine rdmft_init()

      PUSH_SUB(scf_rdmft.rdmft_init)  

      if(st%nst < st%qtot + 1) then   
        message(1) = "Too few states to run RDMFT calculation"
        message(2) = "Number of states should be at least the number of electrons plus one"
        call messages_fatal(2)
      end if
   
      if (states_are_complex(st)) then
        call messages_not_implemented("Complex states for RDMFT")
      end if

 

      !%Variable RDMTolerance
      !%Type float
      !%Default 1e-7 Ha
      !%Section SCF::RDMFT
      !%Description
      !% Convergence criterion for stopping the occupation numbers minimization. Minimization is
      !% stopped when all derivatives of the energy wrt. each occupation number 
      !% are smaller than this criterion. The bisection for finding the correct mu that is needed
      !% for the occupation number minimization also stops according to this criterion.
      !%End

      call parse_variable('RDMTolerance', CNST(1.0e-7), rdm%toler)

      !%Variable RDMConvEner
      !%Type float
      !%Default 1e-6 Ha
      !%Section SCF::RDMFT
      !%Description
      !% Convergence criterion for stopping the overall minimization of the energy with
      !% respect to occupation numbers and the orbitals. The minimization of the 
      !% energy stops when the total energy difference between two subsequent 
      !% minimizations of the energy with respect to the occupation numbers and the
      !% orbitals is smaller than this criterion. It is also used to exit the orbital minimization.
      !%End

      call parse_variable('RDMConvEner', CNST(1.0e-7), rdm%conv_ener)
      
      !%Variable RDMBasis
      !%Type logical
      !%Default yes 
      !%Section SCF::RDMFT
      !%Description
      !% If true, all the energy terms and corresponding derivatives involved in RDMFT will
      !% not be calculated on the grid but on the basis of the initial orbitals
      !%End

      call parse_variable('RDMBasis',.true., rdm%do_basis)

      ! shortcuts
      rdm%gr   => gr
      rdm%st   => st
      
      if (rdm%do_basis.eqv..true.) then
        rdm%n_twoint = st%nst*(st%nst+1)*(st%nst**2+st%nst+2)/8 
        SAFE_ALLOCATE(rdm%eone_int(1:st%nst, 1:st%nst))
        SAFE_ALLOCATE(rdm%twoint(1:rdm%n_twoint))
        SAFE_ALLOCATE(rdm%i_index(1:rdm%n_twoint))
        SAFE_ALLOCATE(rdm%j_index(1:rdm%n_twoint))
        SAFE_ALLOCATE(rdm%k_index(1:rdm%n_twoint))
        SAFE_ALLOCATE(rdm%l_index(1:rdm%n_twoint))
        SAFE_ALLOCATE(rdm%vecnat(1:st%nst, 1:st%nst))
        SAFE_ALLOCATE(rdm%Coul(1:rdm%st%nst, 1:rdm%st%nst, 1:rdm%st%nst))
        SAFE_ALLOCATE(rdm%Exch(1:rdm%st%nst, 1:rdm%st%nst, 1:rdm%st%nst))
        rdm%eone_int = M_ZERO
        rdm%twoint = M_ZERO
        rdm%vecnat = M_ZERO
        rdm%i_index = M_ZERO
        rdm%j_index = M_ZERO
        rdm%k_index = M_ZERO
        rdm%l_index = M_ZERO
        rdm%Coul = M_ZERO
        rdm%Exch = M_ZERO
        do ist = 1, st%nst
          rdm%vecnat(ist,ist)= M_ONE
        end do
      end if

      SAFE_ALLOCATE(rdm%eone(1:st%nst))
      SAFE_ALLOCATE(rdm%hartree(1:st%nst, 1:st%nst))
      SAFE_ALLOCATE(rdm%exchange(1:st%nst, 1:st%nst))
      SAFE_ALLOCATE(rdm%evalues(1:st%nst))

      rdm%eone = M_ZERO
      rdm%hartree = M_ZERO
      rdm%exchange = M_ZERO
      rdm%evalues = M_ZERO
      rdm%mu = M_TWO*st%eigenval(int(st%qtot*M_HALF), 1)
      rdm%qtot = st%qtot
      rdm%occsum = M_ZERO
      rdm%scale_f = CNST(1e-2)
      rdm%maxFO = M_ZERO
      rdm%iter = 1

      POP_SUB(scf_rdmft.rdmft_init)

    end subroutine rdmft_init

    ! ----------------------------------------

    subroutine rdmft_end()

      PUSH_SUB(scf_rdmft.rdmft_end)
    
      nullify(rdm%gr)
      nullify(rdm%st)

      SAFE_DEALLOCATE_A(rdm%evalues)
      SAFE_DEALLOCATE_A(rdm%eone)
      SAFE_DEALLOCATE_A(rdm%hartree)
      SAFE_DEALLOCATE_A(rdm%exchange)

      if (rdm%do_basis.eqv..true.) then
        SAFE_DEALLOCATE_A(rdm%eone_int)
        SAFE_DEALLOCATE_A(rdm%twoint)
        SAFE_DEALLOCATE_A(rdm%i_index)
        SAFE_DEALLOCATE_A(rdm%j_index)
        SAFE_DEALLOCATE_A(rdm%k_index)
        SAFE_DEALLOCATE_A(rdm%l_index)
        SAFE_DEALLOCATE_A(rdm%vecnat)
        SAFE_DEALLOCATE_A(rdm%Coul)
        SAFE_DEALLOCATE_A(rdm%Exch)
      end if

      POP_SUB(scf_rdmft.rdmft_end)

    end subroutine rdmft_end

  end subroutine scf_rdmft
  
  ! ---------------------------------------------------------
  
  ! scf for the occupation numbers 
  subroutine scf_occ(rdm, gr, hm, st, energy)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)   :: energy

    integer :: ist, icycle, ierr, ik
    FLOAT ::  sumgi1, sumgi2, sumgim, mu1, mu2, mum, dinterv
    FLOAT, allocatable ::  occin(:,:)
    FLOAT, parameter :: smallocc = CNST(0.00001) 
    REAL(8), allocatable ::   theta(:)
    REAL(8) :: objective

    PUSH_SUB(scf_occ)

    write(message(1),'(a)') 'Optimization of occupation numbers'
    call messages_info(1)
    
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:st%nst))

    occin = M_ZERO
    theta  = M_ZERO
    energy = M_ZERO
    
    
    !Initialize the occin. Smallocc is used for numerical stability
    
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:) = smallocc 
    where(occin(:,:) > st%smear%el_per_state-smallocc) occin(:,:) = st%smear%el_per_state - smallocc

    !Renormalize the occupation numbers 
    rdm%occsum = M_ZERO

    do ist = 1, st%nst
      do ik = 1, st%d%nik
        rdm%occsum = rdm%occsum + occin(ist,ik)
      end do
    end do

    do ist = 1, st%nst
      do ik = 1, st%d%nik
        occin(ist, ik) = occin(ist,ik)*st%qtot/rdm%occsum
      end do
    end do 

    rdm%occsum = st%qtot

    st%occ = occin
    
    if((rdm%iter == 1).and. (rdm%do_basis.eqv. .true.))  then 
      call dstates_me_two_body(gr, st, rdm%n_twoint, rdm%i_index, rdm%j_index, rdm%k_index, rdm%l_index, rdm%twoint)
      call rdm_integrals(rdm,hm,st,gr) 
      call sum_integrals(rdm)
    end if

    call rdm_derivatives(rdm, hm, st, gr)

    !finding the chemical potential mu such that the occupation numbers sum up to the number of electrons
    !bisection to find the root of rdm%occsum-st%qtot=M_ZERO
    mu1 = rdm%mu   !initial guess for mu 
    mu2 = -CNST(1.0e-6)
    dinterv = M_HALF

    !use n_j=sin^2(2pi*theta_j) to treat pinned states, minimize for both intial mu
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call minimize_multidim(MINMETHOD_BFGS2, st%nst, theta, real(0.05,8), real(0.01, 8), &
        real(CNST(1e-8), 8), real(CNST(1e-8),8), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
    sumgi1 = rdm%occsum - st%qtot
    rdm%mu = mu2
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call minimize_multidim(MINMETHOD_BFGS2, st%nst, theta, real(0.05,8), real(0.01, 8), &
        real(CNST(1e-8), 8), real(CNST(1e-8),8), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
    sumgi2 = rdm%occsum - st%qtot

    ! Adjust the interval between the initial mu to include the root of rdm%occsum-st%qtot=M_ZERO
    do while (sumgi1*sumgi2 > M_ZERO) 
      if (sumgi2 > M_ZERO) then
        mu2 = mu1
        sumgi2 = sumgi1
        mu1 = mu1 - dinterv
        rdm%mu = mu1
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call minimize_multidim(MINMETHOD_BFGS2, st%nst, theta, real(0.05,8), real(0.01, 8), &
            real(CNST(1e-8), 8), real(CNST(1e-8),8), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
        sumgi1 = rdm%occsum - st%qtot 
      else
        mu1 = mu2
        sumgi1 = sumgi2
        mu2 = mu2 + dinterv
        rdm%mu = mu2
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call minimize_multidim(MINMETHOD_BFGS2, st%nst, theta, real(0.05,8), real(0.01, 8), &
            real(CNST(1e-8), 8), real(CNST(1e-8),8), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
        sumgi2 = rdm%occsum - st%qtot 
      end if
    end do

    do icycle = 1, 50
      mum = (mu1 + mu2)*M_HALF
      rdm%mu = mum
      theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
      call minimize_multidim(MINMETHOD_BFGS2, st%nst, theta, real(0.05,8), real(0.01, 8), &
          real(CNST(1e-8), 8), real(CNST(1e-8),8), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
      sumgim = rdm%occsum - st%qtot
      if (sumgi1*sumgim < M_ZERO) then
        mu2 = mum
      else
        mu1 = mum
        sumgi1 = sumgim
      end if
      if(abs(sumgim) < rdm%toler .or. abs((mu1-mu2)*M_HALF) < rdm%toler)  exit
      cycle
    end do
    if (icycle >= 50) then
      write(message(1),'(a,1x,f11.4)') 'Bisection ended without finding mu, sum of occupation numbers:', rdm%occsum
      call messages_fatal(1)
    end if

    do ist = 1, st%nst
      st%occ(ist, 1) = st%smear%el_per_state*sin(theta(ist)*M_PI*M_TWO)**2
    end do
 
    energy = objective
    
    write(message(1),'(a4,5x,a12)')'#st','Occupation'
    call messages_info(1)   

    do ist = 1, st%nst
      write(message(1),'(i4,3x,f11.4)') ist, st%occ(ist, 1)
      call messages_info(1)
    end do

    write(message(1),'(a,1x,f11.4)') 'Sum of occupation numbers', rdm%occsum
    write(message(2),'(a,es15.5)') 'Total energy ', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
    call messages_info(2)   

    SAFE_DEALLOCATE_A(occin)
    SAFE_DEALLOCATE_A(theta)
    POP_SUB(scf_occ)

  end subroutine scf_occ

  

    subroutine objective_rdmft(size, theta, objective, getgrad, df)
      integer,     intent(in)    :: size
      REAL_DOUBLE, intent(in)    :: theta(size)
      REAL_DOUBLE, intent(inout) :: objective
      integer,     intent(in)    :: getgrad
      REAL_DOUBLE, intent(inout) :: df(size)

      integer :: ist
      FLOAT, allocatable :: dE_dn(:),occ(:)
 
      PUSH_SUB(objective_rdmft)

      ASSERT(size == rdm%st%nst)

      SAFE_ALLOCATE(dE_dn(1:size))
      SAFE_ALLOCATE(occ(1:size))

      occ = M_ZERO

      do ist = 1, size
        occ(ist) = M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
      end do
      
      rdm%occsum = M_ZERO
      do ist = 1, size
        rdm%occsum = rdm%occsum + occ(ist)
      end do
      
      !calculate the total energy without nuclei interaction and the energy
      !derivatives with respect to the occupation numbers

      call total_energy_rdm(rdm, occ, objective, dE_dn)
      do ist = 1, size
        df(ist) = M_FOUR*M_PI*sin(M_FOUR*theta(ist)*M_PI)*(dE_dn(ist)-rdm%mu)
      end do
      objective = objective - rdm%mu*(rdm%occsum - rdm%qtot)

      SAFE_DEALLOCATE_A(dE_dn)
      SAFE_DEALLOCATE_A(occ)


      POP_SUB(objective_rdmft)

    end subroutine objective_rdmft

    subroutine write_iter_info_rdmft(iter, size, energy, maxdr, maxdf, theta)
        implicit none
        integer, intent(in) :: iter
        integer, intent(in) :: size
        real(8), intent(in) :: energy, maxdr, maxdf
        real(8), intent(in) :: theta(size)

       PUSH_SUB(write_iter_info_rdmft)

       ASSERT(size == rdm%st%nst)

       POP_SUB(write_iter_info_rdmft)

    end subroutine write_iter_info_rdmft
   
  ! --------------------------------------------
    
  ! scf for the natural orbitals
  subroutine scf_orb(rdm, gr, geo, st, ks, hm, energy)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr !< grid
    type(geometry_t),     intent(inout) :: geo !< geometry
    type(states_t),       intent(inout) :: st !< States
    type(v_ks_t),         intent(inout) :: ks !< Kohn-Sham
    type(hamiltonian_t),  intent(inout) :: hm !< Hamiltonian
    FLOAT ,               intent(out)   :: energy    
    
    integer :: ist, jst
    FLOAT, allocatable ::  lambda(:,:), FO(:,:)
    
    PUSH_SUB(scf_orb)

    !matrix of Lagrange Multipliers from  Equation (8), Piris and Ugalde, Vol.13, No. 13, J. Comput. Chem. 
    SAFE_ALLOCATE(lambda(1:st%nst,1:st%nst)) 
    SAFE_ALLOCATE(FO(1:st%nst, 1:st%nst))    !Generalized Fockian Equation (11) 

    lambda = M_ZERO
    FO = M_ZERO
    if (rdm%do_basis.eqv..false.) then 
      call density_calc (st,gr,st%rho)
      call v_ks_calc(ks,hm,st,geo)
      call hamiltonian_update(hm, gr%mesh, gr%der%boundaries)
    end if

    call construct_f(hm,st,gr,lambda,rdm)
    
    !Set up FO matrix 
    if (rdm%iter==1) then
      do ist = 1, st%nst
        do jst = 1, ist
          FO(ist, jst) = M_HALF*(lambda(ist, jst) + lambda(jst, ist))
          FO(jst, ist) = FO(ist, jst)
        end do
      end do
    else
      do ist = 1, st%nst
        do jst = 1, ist - 1
          FO(jst, ist) = - ( lambda(jst, ist) - lambda(ist ,jst))
        end do
      end do
      rdm%maxFO = maxval(abs(FO))
      do ist = 1, st%nst
        FO(ist, ist) = rdm%evalues(ist)
        do jst = 1, ist-1
          if(abs(FO(jst, ist)) > rdm%scale_f) then
            FO(jst, ist) = rdm%scale_f*FO(jst,ist)/abs(FO(jst, ist))
          end if
          FO(ist, jst) = FO(jst, ist)
        end do
      end do
    end if
 
    call lalg_eigensolve(st%nst, FO, rdm%evalues)
    call assign_eigfunctions(st, gr, FO)
    if (rdm%do_basis.eqv..true.) call sum_integrals(rdm) ! to calculate rdm%Coul and rdm%Exch with the new rdm%vecnat 
    call rdm_derivatives(rdm, hm, st, gr)
    call total_energy_rdm(rdm, st%occ(:,1), energy)

    SAFE_DEALLOCATE_A(lambda) 
    SAFE_DEALLOCATE_A(FO) 

    POP_SUB(scf_orb)

  end subroutine scf_orb

 
  !---------------------------------------------------------------
  ! Minimize the total energy wrt. an orbital by steepest descent
  !---------------------------------------------------------------

  subroutine scf_orb_direct(rdm, gr, geo, st, ks, hm, stepsize, energy)
    
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr !< grid
    type(geometry_t),     intent(inout) :: geo !< geometry
    type(states_t),       intent(inout) :: st !< States
    type(v_ks_t),         intent(inout) :: ks !< Kohn-Sham
    type(hamiltonian_t),  intent(inout) :: hm !< Hamiltonian
    FLOAT,                intent(inout) :: stepsize(1:st%nst) !< step for steepest decent
    FLOAT,                intent(out)   :: energy    

    type(states_t)     :: states_old

    integer            :: ist, jst, lst, ip, itry, trymax
    FLOAT              :: scaleup, scaledown, norm, projection
    FLOAT              :: energy_old, energy_diff, smallstep, thresh
    FLOAT, allocatable :: E_deriv(:), dpsi(:,:), dpsi2(:,:)  

    PUSH_SUB(scf_orb_direct)

    SAFE_ALLOCATE(E_deriv(1:gr%mesh%np_part))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))

    E_deriv = M_ZERO

    !set parameters for steepest decent
    scaleup = 2.4d0
    scaledown = 0.5d0
    trymax = 5
    smallstep = 1d-10
    thresh = 1d-10

    call density_calc (st, gr, st%rho)
    call v_ks_calc(ks, hm, st, geo)
    call hamiltonian_update(hm, gr%mesh, gr%der%boundaries)

    call rdm_derivatives(rdm, hm, st, gr)
    call total_energy_rdm(rdm, st%occ(:,1), energy)
    !store current energy
    energy_old = energy

    do ist = 1, st%nst
      !do not try to improve states that are converged
      if (stepsize(ist).lt. smallstep) cycle

      !store current states
      call states_copy(states_old, st)

      !calculate the derivative with respect to state ist (also removes changes in direction of ist and normalizes)
      call E_deriv_calc(gr, st, hm, ist, E_deriv)    
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      
      !add a step along the gradient
      do itry = 1, trymax
        forall(ip=1:gr%mesh%np_part)
          dpsi(ip,1) = dpsi(ip,1) - E_deriv(ip)*stepsize(ist)
        end forall      
        !orthogonalize new orbital to all orbitals with smaller index and normalize
        do jst = 1, ist - 1
          call states_get_state(st, gr%mesh, jst, 1, dpsi2)
          projection = dmf_dotp(gr%mesh, dpsi(:,1), dpsi2(:,1))
          forall(ip = 1:gr%mesh%np_part)
            dpsi(ip, 1) = dpsi(ip, 1) - projection*dpsi2(ip, 1)
          end forall
        enddo !jst
        norm = sqrt(dmf_dotp(gr%mesh, dpsi(:,1), dpsi(:,1)))
        forall(ip=1:gr%mesh%np_part)
          dpsi(ip, 1) = dpsi(ip, 1)/norm
        end forall
        call states_set_state(st, gr%mesh, ist, 1, dpsi)
        !orthogonalize all orbitals with larger index to new orbital and normalize
        do jst = ist + 1, st%nst
          call states_get_state(st, gr%mesh, jst, 1, dpsi)
          do lst = ist, jst - 1
            call states_get_state(st, gr%mesh, lst, 1, dpsi2)
            projection = dmf_dotp(gr%mesh, dpsi(:,1), dpsi2(:,1))
            forall(ip = 1:gr%mesh%np_part)
              dpsi(ip,1) = dpsi(ip,1) - projection*dpsi2(ip,1)
            end forall
          enddo !lst
          norm = sqrt(dmf_dotp(gr%mesh, dpsi(:,1), dpsi(:,1)))
          forall(ip=1:gr%mesh%np_part)
            dpsi(ip, 1) = dpsi(ip, 1)/norm
          end forall
          call states_set_state(st, gr%mesh, jst, 1, dpsi)
        enddo !jst
        
        !calculate total energy
        call density_calc (st, gr, st%rho)
        call v_ks_calc(ks, hm, st, geo)
        call hamiltonian_update(hm, gr%mesh, gr%der%boundaries)
        call rdm_derivatives(rdm, hm, st, gr)
        call total_energy_rdm(rdm, st%occ(:,1), energy)

        !check if step lowers the energy
        energy_diff = energy - energy_old

        if (energy_diff .lt. - thresh) then 
          !sucessful step
          stepsize(ist) = stepsize(ist)*scaleup !increase stepsize for this state
          energy_old = energy
          exit
        else 
          !not sucessful step
          stepsize(ist) = stepsize(ist)*scaledown !decrease stepsize
          !undo changes in st and dpsi (cannot use states_copy due to memory leak) 
          do jst = st%nst, ist, -1
            call states_get_state(states_old, gr%mesh, jst, 1, dpsi)
            call states_set_state(st, gr%mesh, jst, 1, dpsi)
          enddo
          !if (itry == trymax) then
          !  write(message(1),'(a,i3,2x,a)') 'for state', ist, 'energy could not be improved'
          !  call messages_info(1)
          !endif
        endif !energy_diff  
      enddo !itry
      call states_end(states_old)
      !scale stepsize for next iteration, check convergence of states
      stepsize(ist) = scaledown*stepsize(ist)
      if(abs(stepsize(ist)) .lt. smallstep) then
        write(message(1),'(a,i3,2x,a)') 'for state', ist, 'the minimal step size is reached'
        write(message(2),'(a)') 'state is converged'
        call messages_info(2)
      endif
    enddo !ist

    SAFE_DEALLOCATE_A(E_deriv)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)

    POP_SUB(scf_orb_direct)

  end subroutine scf_orb_direct

  !--------------------------------------------------------------------
  ! calculate the derivative of the energy with respect to orbital ist
  !--------------------------------------------------------------------
  subroutine E_deriv_calc(gr, st, hm, ist, E_deriv)
    type(grid_t),         intent(inout) :: gr !< grid
    type(states_t),       intent(inout) :: st !< States
    type(hamiltonian_t),  intent(inout) :: hm !< Hamiltonian
    integer,              intent(in)    :: ist !number of state
    FLOAT,                intent(out)   :: E_deriv(1:gr%mesh%np_part)

    integer            :: jst, ii, ip
    FLOAT              :: E_deriv_corr, norm, projection
    FLOAT, allocatable :: rho_spin(:,:), rho(:), pot(:)
    FLOAT, allocatable :: hpsi1(:,:), hpsi2(:,:), dpsi(:,:), dpsi2(:,:)
  

    PUSH_SUB(E_deriv_calc)

    E_deriv = M_ZERO

    SAFE_ALLOCATE(rho_spin(1:gr%mesh%np_part, 1:hm%d%ispin)) 
    SAFE_ALLOCATE(rho(1:gr%mesh%np_part))
    SAFE_ALLOCATE(pot(1:gr%mesh%np_part))
    SAFE_ALLOCATE(hpsi1(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi2(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))

    !Initialize all local arrays to zero
    hpsi1 = M_ZERO
    hpsi2 = M_ZERO
    dpsi  = M_ZERO
    dpsi2  = M_ZERO
    rho_spin = M_ZERO    
    rho = M_ZERO    
    pot = M_ZERO
    E_deriv_corr = M_ZERO    

    call density_calc(st, gr, rho_spin)
    do ii = 1, hm%d%ispin
      rho(:) = rho_spin(:, ii)
    end do

    SAFE_DEALLOCATE_A(rho_spin) 

    call dpoisson_solve(psolver, pot, rho, all_nodes=.false.) !the Hartree potential
    
    call states_get_state(st, gr%mesh, ist, 1, dpsi)

    call dhamiltonian_apply(hm, gr%der, dpsi, hpsi1, ist, 1, &
                            & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
    call dhamiltonian_apply(hm, gr%der, dpsi, hpsi2, ist, 1, &
                            & terms = TERM_OTHERS)
      
    !bare derivative wrt. state ist   
    forall(ip=1:gr%mesh%np_part)
      E_deriv(ip) = st%occ(ist, 1)*(hpsi1(ip, 1) + pot(ip)*dpsi(ip, 1)) 
     
      !only for the Mueller functional
      E_deriv(ip) = E_deriv(ip) + sqrt(st%occ(ist, 1))*hpsi2(ip, 1)
    end forall

    norm = sqrt(dmf_dotp(gr%mesh, E_deriv, E_deriv))

    !remove gradient in direction of orbital ist itself and normalize
    projection = dmf_dotp(gr%mesh, E_deriv, dpsi(:,1))
    forall(ip = 1:gr%mesh%np_part)
      E_deriv(ip) = E_deriv(ip) - projection*dpsi(ip,1)
    end forall
    norm = sqrt(dmf_dotp(gr%mesh, E_deriv(:), E_deriv(:)))
    if(norm > M_ZERO) then
      forall(ip=1:gr%mesh%np_part)
        E_deriv(ip) = E_deriv(ip)/norm
      end forall
    endif      
 
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(pot)
    SAFE_DEALLOCATE_A(hpsi1)
    SAFE_DEALLOCATE_A(hpsi2)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)

    POP_SUB(E_deriv_calc)

  end subroutine E_deriv_calc


  ! ----------------------------------------
  ! constructs the Lagrange multiplyers needed for the orbital minimization
  subroutine construct_f(hm, st, gr, lambda, rdm)
    type(hamiltonian_t),  intent(in) :: hm
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in) :: gr
    FLOAT,                intent(out):: lambda(:,:) !< (1:st%nst, 1:st%nst)
    type(rdm_t),          intent(inout) :: rdm
      
    FLOAT, allocatable :: hpsi(:,:), hpsi1(:,:), dpsi(:,:), dpsi2(:,:), fvec(:) 
    FLOAT, allocatable :: g_x(:,:), g_h(:,:), rho(:,:), rho_tot(:), pot(:), fock(:,:,:)
    integer :: ist, ip, iorb, jorb, jst

    PUSH_SUB(construct_f)

    lambda = M_ZERO

    if (rdm%do_basis.eqv..false.) then 
      SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part,1:st%d%dim))
      SAFE_ALLOCATE(hpsi1(1:gr%mesh%np_part,1:st%d%dim))
      SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part ,1:st%d%dim))
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part ,1:st%d%dim))
      SAFE_ALLOCATE(rho(1:gr%mesh%np_part,1:hm%d%ispin)) 
      SAFE_ALLOCATE(rho_tot(1:gr%mesh%np_part))
      SAFE_ALLOCATE(pot(1:gr%mesh%np_part))
      SAFE_ALLOCATE(g_x(1:st%nst,1:st%nst))
      SAFE_ALLOCATE(g_h(1:st%nst,1:st%nst))

      hpsi = M_ZERO
      hpsi1 = M_ZERO
      dpsi2 = M_ZERO
      rho = M_ZERO    
      rho_tot = M_ZERO    
      pot = M_ZERO
      dpsi = M_ZERO
      dpsi2 = M_ZERO
      g_x = M_ZERO
      g_h = M_ZERO
    
      !calculate the Lagrange multiplyer lambda matrix on the grid, Eq. (9), Piris and Ugalde, Vol.13, No. 13, J. Comput. Chem.
      call density_calc(st, gr, rho)
      do ist =1, hm%d%ispin
        rho_tot(:) = rho(:, ist)
      end do
      call dpoisson_solve(psolver, pot, rho_tot, all_nodes=.false.) !the Hartree potential
    
      do iorb = 1, st%nst
        call states_get_state(st, gr%mesh, iorb, 1, dpsi)
        call dhamiltonian_apply(hm,gr%der, dpsi, hpsi, iorb, 1, &
                             terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        call dhamiltonian_apply(hm, gr%der, dpsi, hpsi1, iorb, 1, &
                              terms = TERM_OTHERS)
        forall (ip=1:gr%mesh%np_part)
          dpsi(ip,1) = st%occ(iorb, 1)*pot(ip)*dpsi(ip,1)
        end forall
        do jorb = 1, st%nst  
          call states_get_state(st, gr%mesh, jorb, 1, dpsi2)
          lambda(jorb, iorb) = dmf_dotp(gr%mesh, dpsi2(:,1), hpsi(:,1))
          lambda(iorb, jorb) = lambda(jorb, iorb)
          g_h(iorb, jorb) = dmf_dotp(gr%mesh, dpsi(:,1), dpsi2(:, 1))
          g_x(iorb, jorb) = dmf_dotp(gr%mesh, dpsi2(:,1), hpsi1(:,1))
          g_x(iorb, jorb) = sqrt(st%occ(iorb,1))*g_x(iorb, jorb)
        end do
      end do
      
 
      do jorb = 1,st%nst
        do iorb = 1,st%nst
	  lambda(jorb,iorb) = st%occ(iorb,1)*lambda(jorb,iorb) + g_h(iorb, jorb)+ g_x(iorb, jorb)
	end do
      end do

    else ! calculate the same lambda matrix on the basis
      !call sum_integrals(rdm)
      SAFE_ALLOCATE(fvec(1:st%nst))
      SAFE_ALLOCATE(Fock(1:st%nst, 1:st%nst, 1:st%nst))
      Fock = M_ZERO
      
      do iorb = 1, st%nst
        do ist = 1, st%nst
          do jst = 1, ist
            Fock(ist,jst,iorb) = st%occ(iorb, 1)*rdm%eone_int(ist,jst)
            do jorb = 1, st%nst
              !The coefficimnt of the Exchange term below is only for the Mueller functional
              Fock(ist ,jst, iorb) =  Fock(ist ,jst, iorb) + st%occ(iorb, 1)*st%occ(jorb,1)*rdm%Coul(ist, jst, jorb)  &
                                      - sqrt(st%occ(iorb,1))*sqrt(st%occ(jorb,1))*rdm%Exch(ist, jst,jorb)
            enddo
            Fock(jst, ist, iorb) = Fock(ist, jst, iorb)
          enddo
        enddo
      enddo

      do jorb = 1, st%nst
        do ist = 1, st%nst
          fvec(ist) = M_ZERO
          do jst =1, st%nst
            fvec(ist) = fvec(ist)+Fock(ist,jst,jorb)*rdm%vecnat(jst,jorb)
          enddo
        enddo
        do iorb= 1, st%nst
          lambda(iorb, jorb) = M_ZERO
          do ist = 1, st%nst
            lambda(iorb,jorb) = lambda(iorb,jorb) + rdm%vecnat(ist,iorb)*fvec(ist)
          enddo
        enddo
      enddo
    end if

    

    if (rdm%do_basis.eqv..false.) then 
      SAFE_DEALLOCATE_A(g_x)
      SAFE_DEALLOCATE_A(g_h)
      SAFE_DEALLOCATE_A(hpsi)
      SAFE_DEALLOCATE_A(hpsi1)
      SAFE_DEALLOCATE_A(dpsi2)
      SAFE_DEALLOCATE_A(dpsi)
      SAFE_DEALLOCATE_A(rho) 
      SAFE_DEALLOCATE_A(rho_tot) 
      SAFE_DEALLOCATE_A(pot) 
    else
      SAFE_DEALLOCATE_A(fvec) 
      SAFE_DEALLOCATE_A(Fock) 
    end if
   
    POP_SUB(construct_f)

  end subroutine construct_f
   
  ! ----------------------------------------
  
  ! finds the new states after the minimization of the orbitals
  subroutine assign_eigfunctions(st, gr, lambda)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in)    :: gr
    FLOAT,                intent(in)    :: lambda(:, :)
    
    integer :: iqn, iorb, jorb, ist
    FLOAT, allocatable :: vecnat_new(:,:)

    PUSH_SUB(assign_eigenfunctions)

    if (rdm%do_basis.eqv..false.) then 
      do iqn = st%d%kpt%start, st%d%kpt%end
        if(states_are_real(st)) then
          call states_rotate(gr%mesh, st, lambda, iqn)
        else
          call states_rotate(gr%mesh, st, M_z1*lambda, iqn)
        end if
      end do
    else
      SAFE_ALLOCATE(vecnat_new(1:st%nst, 1:st%nst))
      do iorb = 1, st%nst
        do ist = 1, st%nst
          vecnat_new(ist,iorb) = M_ZERO
          do jorb = 1, st%nst
            vecnat_new(ist , iorb) = vecnat_new(ist , iorb) + rdm%vecnat(ist, jorb)*lambda(jorb,iorb)
          enddo
        enddo
      enddo
      
      rdm%vecnat = vecnat_new


      SAFE_DEALLOCATE_A(vecnat_new)
    end if

    
    POP_SUB(assign_eigenfunctions)

  end subroutine assign_eigfunctions
   
  ! --------------------------------------------

  ! calculates the total energy when only the occupation numbers are updated
  subroutine total_energy_rdm(rdm, occ, energy, dE_dn)
    type(rdm_t),          intent(inout)  :: rdm
    FLOAT,                intent(in)     :: occ(:)
    FLOAT,                intent(out)    :: energy
    FLOAT, optional,      intent(out)    :: dE_dn(:) !< (1:st%nst)
     
    integer :: ist, jst
    FLOAT, allocatable :: V_h(:), V_x(:)
     
    PUSH_SUB(total_energy_rdm)
  
    SAFE_ALLOCATE(V_h(1:rdm%st%nst))
    SAFE_ALLOCATE(V_x(1:rdm%st%nst))

    energy = M_ZERO
    V_h = M_ZERO
    V_x = M_ZERO
     

    !Calculate hartree and exchange contribution 
    !This is only for the Mueller functional and has to be changed
    do ist = 1, rdm%st%nst
      do jst = 1, rdm%st%nst
        V_h(ist) = V_h(ist) + occ(jst)*rdm%hartree(ist, jst)
        V_x(ist) = V_x(ist) - sqrt(occ(jst))*rdm%exchange(ist, jst) 
      end do
      V_x(ist) = V_x(ist)*M_HALF/max(sqrt(occ(ist)), CNST(1.0e-16))
    end do


    !Calculate the energy derivative with respect to the occupation numbers
    if (present(dE_dn)) then
      dE_dn(:) = rdm%eone(:) + V_h(:) + V_x(:)
    end if

    !Total energy calculation without nuclei interaction  
    do ist = 1, rdm%st%nst
      energy = energy + occ(ist)*rdm%eone(ist) &
                      + M_HALF*occ(ist)*V_h(ist) & 
                      + occ(ist)*V_x(ist)
    end do

    SAFE_DEALLOCATE_A(V_h)
    SAFE_DEALLOCATE_A(V_x)
    
    POP_SUB(total_energy_rdm)
   
  end subroutine total_energy_rdm
  
  ! ----------------------------------------
  ! calculates the derivatives of the energy terms with respect to the occupation numbers
  subroutine rdm_derivatives(rdm, hm, st, gr)
    type(rdm_t),          intent(inout) :: rdm
    type(hamiltonian_t),  intent(in)    :: hm 
    type(states_t),       intent(in)    :: st 
    type(grid_t),         intent(inout) :: gr
    
    FLOAT, allocatable :: hpsi(:,:), rho1(:), rho(:), dpsi(:,:), dpsi2(:,:)
    FLOAT, allocatable :: v_ij(:,:,:)
    FLOAT, allocatable :: lxc(:, :, :) !required input variable for doep_x, not used otherwise, might get used
    FLOAT              :: ex !required input variable for doep_x, not used otherwise, might get used
    FLOAT              :: dd
    
    integer :: ist, jst, nspin_, is, jdm, iorb, jorb

    PUSH_SUB(rdm_derivatives) 


    nspin_ = min(st%d%nspin, 2)
   
    if (rdm%do_basis.eqv..false.) then 
      SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(rho1(1:gr%mesh%np))
      SAFE_ALLOCATE(rho(1:gr%mesh%np))
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
      SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))
      SAFE_ALLOCATE(v_ij(1:gr%der%mesh%np, 1:st%nst, 1:st%nst))
      SAFE_ALLOCATE(lxc(1:gr%mesh%np, st%st_start:st%st_end, 1:nspin_))

      lxc = M_ZERO
      v_ij = M_ZERO
      rdm%eone = M_ZERO
      rdm%hartree = M_ZERO
      rdm%exchange = M_ZERO

      !derivative of one-electron energy with respect to the natural orbitals occupation number
      do ist = 1, st%nst
        call states_get_state(st, gr%mesh, ist, 1, dpsi)
        call dhamiltonian_apply(hm,gr%der, dpsi, hpsi, ist, 1, &
                              terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        rdm%eone(ist) = dmf_dotp(gr%mesh, dpsi(:, 1), hpsi(:, 1))
      end do

      !integrals used for the hartree and exchange parts of the total energy and their derivatives
      do is = 1, nspin_
        do jdm = 1, st%d%dim
          call doep_x(gr%der, st, is, jdm, lxc, ex, 1.d0, v_ij)
        end do
      end do
      do ist = 1, st%nst
        call states_get_state(st, gr%mesh, ist, 1, dpsi)
        rho1(1:gr%mesh%np) = dpsi(1:gr%mesh%np, 1)**2
        do jst = ist, st%nst
          rdm%hartree(ist, jst) = dmf_dotp(gr%mesh, rho1, v_ij(:,jst, jst))
          rdm%hartree(jst, ist) = rdm%hartree(ist, jst)
          call states_get_state(st, gr%mesh, jst, 1, dpsi2)
          rho(1:gr%mesh%np) = dpsi2(1:gr%mesh%np, 1)*dpsi(1:gr%mesh%np, 1)
          rdm%exchange(ist, jst) = dmf_dotp(gr%mesh, rho, v_ij(:, ist, jst))
          rdm%exchange(jst, ist) = rdm%exchange(ist, jst)
        end do
      end do
    
      SAFE_DEALLOCATE_A(hpsi)
      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(rho1)
      SAFE_DEALLOCATE_A(dpsi)
      SAFE_DEALLOCATE_A(dpsi2)
      SAFE_DEALLOCATE_A(lxc)
      SAFE_DEALLOCATE_A(v_ij)

    else !if energy derivatives are expanded in a basis set

      do iorb = 1, st%nst
        rdm%eone(iorb) = M_ZERO
        do ist = 1, st%nst
          do jst = 1, st%nst
            dd = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)  
            rdm%eone(iorb) =rdm%eone(iorb) + dd*rdm%eone_int(ist, jst)
          end do
        end do
      end do

      do iorb = 1, st%nst
        do jorb =1 , iorb
          rdm%hartree(iorb ,jorb)= M_ZERO; rdm%exchange(iorb,jorb)=M_ZERO     
          do ist =1, st%nst
            do jst =1, st%nst
              dd = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)
              rdm%hartree(iorb ,jorb) = rdm%hartree(iorb ,jorb)+rdm%Coul(ist,jst, jorb)*dd 
              rdm%exchange(iorb ,jorb) = rdm%exchange(iorb ,jorb)+rdm%Exch(ist,jst, jorb)*dd 
            end do
          end do
          rdm%hartree(jorb, iorb) = rdm%hartree(iorb, jorb); rdm%exchange(jorb, iorb) = rdm%exchange(iorb, jorb)
        end do
      end do
    end if
     
    POP_SUB(rdm_derivatives) 

  end subroutine rdm_derivatives
  
  ! --------------------------------------------
  !calculates the one electron integrals in the basis of the initial orbitals
  subroutine rdm_integrals(rdm, hm, st, gr)
    type(rdm_t),          intent(inout) :: rdm
    type(hamiltonian_t),  intent(in)    :: hm 
    type(states_t),       intent(in)    :: st 
    type(grid_t),         intent(inout) :: gr
    
    FLOAT, allocatable :: hpsi(:,:)
    FLOAT, allocatable :: dpsi(:,:), dpsi2(:,:)
    integer :: ist, jst, nspin_ 

    PUSH_SUB(rdm_integrals)
 
    nspin_ = min(st%d%nspin, 2)
    
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part,1:st%d%dim))

    !calculate integrals of the one-electron energy term with respect to the initial orbital basis
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      do jst = ist, st%nst
        call states_get_state(st, gr%mesh, jst, 1, dpsi2)
        call dhamiltonian_apply(hm,gr%der, dpsi, hpsi, ist, 1, &
                              terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        rdm%eone_int(jst, ist) = dmf_dotp(gr%mesh, dpsi2(:, 1), hpsi(:, 1))
        rdm%eone_int(ist, jst) = rdm%eone_int(jst, ist)
      end do
    end do

    SAFE_DEALLOCATE_A(hpsi)    
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)
    
    POP_SUB(rdm_integrals) 

  end subroutine rdm_integrals

  ! --------------------------------------------
  ! constructs the Hartree and Exchange part of the RDMFT Fock matrix
  subroutine sum_integrals(rdm)
    type(rdm_t),        intent(inout) :: rdm

    integer :: ist, jst, kst, lst, iorb, icount
    logical :: inv_pairs
    FLOAT :: two_int, wij, wik, wil, wjk, wjl, wkl
    FLOAT, allocatable                :: DM(:,:,:)

    PUSH_SUB(sum_integrals)

    SAFE_ALLOCATE(DM(1:rdm%st%nst, 1:rdm%st%nst, 1:rdm%st%nst))
     
    rdm%Coul = M_ZERO  
    rdm%Exch = M_ZERO 
    DM = M_ZERO 
  
    do iorb = 1, rdm%st%nst
      do ist = 1, rdm%st%nst
        do jst = 1, ist
          DM(ist, jst, iorb) = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)
          DM(jst, ist, iorb) = DM(ist, jst, iorb)
        enddo
      enddo
    enddo

    do icount = 1, rdm%n_twoint

      ist = rdm%i_index(icount) 
      jst = rdm%j_index(icount) 
      kst = rdm%k_index(icount) 
      lst = rdm%l_index(icount) 

      two_int = rdm%twoint(icount)
         
      ! create weights of unique integrals 
      if(ist == jst) then
        wij = 1.d0
      else
        wij = 2.d0
      endif
      if(kst == lst) then
        wkl = 1.d0
      else
        wkl = 2.d0
      endif

      if(ist == kst .and. jst /= lst) then
        wik = M_TWO
      else
        wik = M_ONE
      endif
      if(ist == lst .and. jst /= kst) then
        wil = M_TWO
      else
        wil = M_ONE
      endif
      if(jst == kst .and. ist /= lst) then
        wjk = M_TWO
      else
        wjk = M_ONE
      endif
      if(jst == lst .and. ist /= kst) then
        wjl = M_TWO
      else
        wjl = M_ONE
      endif

      inv_pairs = (ist /= kst .or. jst /= lst)

      do iorb = 1, rdm%st%nst 

        !the Hartree terms
        rdm%Coul(ist, jst, iorb) = rdm%Coul(ist, jst, iorb) + DM(kst, lst, iorb)*wkl*two_int
        if(inv_pairs) rdm%Coul(kst, lst, iorb) = rdm%Coul(kst, lst, iorb) + DM(ist, jst, iorb)*wij*two_int

        !the exchange terms
        !weights are only included if they can differ from one
        rdm%Exch(ist, kst, iorb) = rdm%Exch(ist, kst, iorb) + two_int*DM(jst, lst, iorb)*wik
        if (kst /= lst) then 
          rdm%Exch(ist, lst, iorb) = rdm%Exch(ist, lst, iorb) + two_int*DM(jst, kst, iorb)*wil
        endif
        if(ist /= jst) then
          if(jst >= kst) then
            rdm%Exch(jst, kst, iorb) = rdm%Exch(jst, kst, iorb) + two_int*DM(ist, lst, iorb)*wjk
          else            
            if(inv_pairs) rdm%Exch(kst, jst, iorb) = rdm%Exch(kst, jst, iorb) + two_int*DM(ist, lst, iorb)
          endif
        endif        
        if(ist /=jst .and. kst /= lst) then
          if(jst >= lst) then
            rdm%Exch(jst, lst, iorb) = rdm%Exch(jst, lst, iorb) + two_int*DM(ist, kst, iorb)*wjl
          else
            if(inv_pairs) rdm%Exch(lst, jst, iorb) = rdm%Exch(lst, jst, iorb) + two_int*DM(ist, kst, iorb)
          endif
        endif

      end do !iorb
    end do !icount
   
    do iorb =1, rdm%st%nst
      do ist = 1, rdm%st%nst
        do jst = 1, ist-1
          rdm%Coul(jst, ist, iorb) = rdm%Coul(ist, jst, iorb)
          rdm%Exch(jst, ist, iorb) = rdm%Exch(ist, jst, iorb)
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(DM)

    POP_SUB(sum_integrals)
 
  end subroutine sum_integrals



end module rdmft_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:




