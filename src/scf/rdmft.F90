! Copyright (C) 2012 I. Theophilou, N. Helbig
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
!! $Id$

#include "global.h"

module rdmft_m
  use datasets_m
  use density_m
  use eigensolver_m
  use energy_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_m
  use loct_math_m 
  use messages_m
  use mesh_function_m
  use output_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use states_calc_m
  use unit_m
  use unit_system_m
  use v_ks_m
 
  implicit none

  private
  public ::                      & 
       assign_eigfunctions,      &
       rdm_t,                    &
       rdmft_init,               &
       rdmft_end,                &
       scf_occ,                  &
       scf_orb,                  &
       scf_rdmft

  type rdm_t
    type(states_t) :: psi
    integer  :: max_iter, iter
    FLOAT    :: mu, occsum, qtot, scale_f, step, toler, conv_ener
    FLOAT, allocatable   :: eone(:), hartree(:,:), exchange(:,:), evalues(:)   
  end type rdm_t
  
contains
   
  subroutine rdmft_init(rdm, st)
    type(rdm_t),           intent(inout)    :: rdm
    type(states_t),        intent(inout)    :: st


    PUSH_SUB(rdmft_init)  

    if(st%nst < st%qtot+5) then   
      message(1) = "Too few states to run RDMFT calculation"
      message(2) = "Number of states should be at least the number of electrons plus five"
      call messages_fatal(2)
    endif

    SAFE_ALLOCATE(rdm%eone(1:st%nst))
    SAFE_ALLOCATE(rdm%hartree(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(rdm%exchange(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(rdm%evalues(1:st%nst))

    rdm%eone = M_ZERO
    rdm%hartree = M_ZERO
    rdm%exchange = M_ZERO
    rdm%mu = M_TWO*st%eigenval(int(st%qtot*M_HALF), 1)
    rdm%qtot = st%qtot
    rdm%occsum = M_ZERO
    rdm%scale_f = CNST(1e-3)
    rdm%iter = 1
    
   

    !%Variable RDMMaxIter
    !%Type integer
    !%Default 100
    !%Section SCF::RDMFT
    !%Description
    !% Even if the convergence criterion is not satisfied, the minimization will stop
    !% after this number of iterations. The default is 100.
    !%End 
    call parse_integer(datasets_check('RDMMaxIter'), 100, rdm%max_iter)
    
    !%Variable RDMStep
    !%Type float
    !%Default 1e-2 
    !%Section SCF::RDMFT
    !%Description 
    !% Initial step for the occupation number optimizer. The default is 1.0e-2.
    !%End
    call parse_float(datasets_check('RDMStep'), CNST(1.0e-2), rdm%step)

    !%Variable RDMTolerance
    !%Type float
    !%Default 1e-6 Ha
    !%Section SCF::RDMFT
    !%Description
    !% Convergence criterion, for stopping the minimization. Minimization is
    !% stopped when all derivatives of the energy wrt. the occupation number 
    !% are smaller than this criterion. The default is 1.0e-6.
    !%End
    call parse_float(datasets_check('RDMTolerance'), CNST(1.0e-6), rdm%toler)

    !%Variable RDMConvEner
    !%Type float
    !%Default 1e-6 Ha
    !%Section SCF::RDMFT
    !% Convergence criterion, for stopping the minimization of the energy with
    !% respect to occupation numbers and the orbitals. The minimization of the 
    !% energy stops when the total energy difference between two subsequent 
    !% minimizations of the energy with respect to the occupation numbers and the
    !% orbitals is smaller than this criterion. The default is 1.0e-6
    !%End

    call parse_float(datasets_check('RDMConvEner'), CNST(1.0e-6), rdm%conv_ener)
    
    
    POP_SUB(rdmft_init)

  end subroutine rdmft_init

  ! ----------------------------------------

  subroutine rdmft_end(rdm)
    type(rdm_t), intent(inout) :: rdm

    PUSH_SUB(rdmft_end)


    SAFE_DEALLOCATE_A(rdm%evalues)
    SAFE_DEALLOCATE_A(rdm%eone)
    SAFE_DEALLOCATE_A(rdm%hartree)
    SAFE_DEALLOCATE_A(rdm%exchange)

    POP_SUB(rdmft_end)

  end subroutine rdmft_end

  ! ----------------------------------------

  ! scf for the occupation numbers and the natural orbitals
  subroutine scf_rdmft(rdm, gr, geo, st, ks, hm, outp)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr  !< grid
    type(geometry_t),     intent(inout) :: geo !< geometry
    type(states_t),       intent(inout) :: st  !< States
    type(v_ks_t),         intent(inout) :: ks  !< Kohn-Sham
    type(hamiltonian_t),  intent(inout) :: hm  !< Hamiltonian
    type(output_t),       intent(in)    :: outp !< output
    
    integer :: iter, counter, icount
    FLOAT :: energy, energy_dif, energy_old, energy_occ,  xpos, xneg, fr
    
    PUSH_SUB(scf_rdmft)
    
    energy_old = CNST(1.0e20)
    xpos = M_ZERO 
    xneg = M_ZERO
    counter = 0
    fr = M_ZERO
    energy = M_ZERO 

    write(message(1),'(a)') 'Initial minimization of occupation numbers'
    call messages_info(1)

    do iter = 1, rdm%max_iter
      write(message(1),'(a, 1x, i4)') 'RDM Iteration:', rdm%iter
      call messages_info(1)
      call scf_occ(rdm, gr, hm, st, energy_occ)
      do icount = 1, 15
        call scf_orb(rdm, gr, geo, st, ks, hm, energy)
        energy_dif = energy - energy_old
        energy_old = energy
        if (abs(energy_dif).lt. rdm%conv_ener)  exit
        if (energy_dif < M_ZERO) then
          xneg = xneg + 1
        else
          xpos = xpos + 1
        end if
        counter = counter + 1
        if (xneg > CNST(1.5e0)*xpos) then
          rdm%scale_f = CNST(1.01)*rdm%scale_f
        elseif (xneg > CNST(1.1e0)*xpos) then
          rdm%scale_f = rdm%scale_f
        else
          rdm%scale_f = CNST(0.95)* rdm%scale_f 
        endif
        counter = 0
        xneg = M_ZERO
        xpos = M_ZERO
      end do
      rdm%iter = rdm%iter + 1
      write(message(1),'(a,es15.8)') ' etot RDMFT after orbital minim = ', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
      call messages_info(1)
      if (abs(energy_occ-energy).lt. rdm%conv_ener)  exit
    end do
   

    call output_states(st,gr,geo,STATIC_DIR,outp)    

    POP_SUB(scf_rdmft)

  end subroutine scf_rdmft
  
  ! ---------------------------------------------------------
  
  ! scf for the occupation numbers 
  subroutine scf_occ(rdm, gr, hm, st, energy)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)   :: energy

    integer :: ist, icycle
    FLOAT ::  sumgi1, sumgi2, sumgim, mu1, mu2, mum, dinterv
    FLOAT, allocatable ::  occin(:,:), hpsi(:,:), pot(:), rho(:), dpsi(:,:), dpsi2(:,:), theta(:)
    FLOAT, parameter :: smallocc = CNST(0.0000001) 

    PUSH_SUB(scf_occ)

    if (hm%d%ispin /= 1) then
      call messages_not_implemented("RDMFT exchange function not yet implemented for spin_polarized or spinors")
    end if
    
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:st%nst))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(pot (1:gr%mesh%np))
    SAFE_ALLOCATE(rho (1:gr%mesh%np))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))

    occin = M_ZERO
    theta  = M_ZERO
    hpsi = M_ZERO
    energy = M_ZERO

    !Initialize the occin. Smallocc is used for numerical stability
    
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:) = smallocc 
    where(occin(:,:) > st%smear%el_per_state-smallocc) occin(:,:) = st%smear%el_per_state - smallocc

    st%occ = occin
    

    call rdm_derivatives(rdm, hm, st, gr)

    call total_energy_rdm(rdm, st,gr, st%occ(:,1), energy)

    !finding the chemical potential mu such that the occupation numbers sum up to the number of electrons
    !bisection to find the root of rdm%occsum-st%qtot=M_ZERO
    mu1 = rdm%mu   !initial guess for mu 
    mu2 = -CNST(1.0e-6)
    dinterv = M_HALF
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call  multid_minimize(st%nst, 1000, theta, energy) 
    sumgi1 = rdm%occsum - st%qtot
    rdm%mu = mu2
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call  multid_minimize(st%nst, 1000, theta, energy) 
    sumgi2 = rdm%occsum - st%qtot

    do while (sumgi1*sumgi2 > M_ZERO)
      if (sumgi2 > M_ZERO) then
        mu2 = mu1
        sumgi2 = sumgi1
        mu1 = mu1 - dinterv
        rdm%mu = mu1
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call  multid_minimize(st%nst, 1000, theta, energy) 
        sumgi1 = rdm%occsum - st%qtot 
      else
        mu1 = mu2
        sumgi1 = sumgi2
        mu2 = mu2 + dinterv
        rdm%mu = mu2
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call  multid_minimize(st%nst, 1000, theta, energy) 
        sumgi2 = rdm%occsum - st%qtot 
      end if
    end do

    do icycle = 1, 1000
      mum = (mu1 + mu2)*M_HALF
      rdm%mu = mum
      theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
      call  multid_minimize(st%nst, 1000, theta, energy) 
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
    if (icycle >= 1000) then
      write(message(1),'(a,1x,f11.6)') 'Bisection ended without finding mu, sum of occupation numbers:', rdm%occsum
      call messages_fatal(1)
    endif

    do ist = 1, st%nst
      st%occ(ist, 1) = st%smear%el_per_state*sin(theta(ist)*M_PI*M_TWO)**2
    end do
    

    write(message(1),'(a,1x,f11.6)') 'Occupations sum', rdm%occsum
    call messages_info(1)
    write(message(1),'(a,es15.8)') ' etot RDMFT after occ minim = ', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
    write(message(2),'(a4,1x,a12)')'#st','Occupation'
    call messages_info(1)   

    do ist = 1, st%nst
      write(message(1),'(i4,3x,f11.6)') ist, st%occ(ist, 1)
      call messages_info(1)  
    end do

    SAFE_DEALLOCATE_A(occin)
    SAFE_DEALLOCATE_A(theta)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(pot)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)
    POP_SUB(scf_occ)

  contains
  
    subroutine multid_minimize(nst, max_iter, theta, objective)!iris maybe do not pass from the arguments the number of iterations 
      integer, intent(in)        :: nst
      integer, intent(in)        :: max_iter
      FLOAT, intent(inout)       :: theta(:)
      FLOAT, intent(out)         :: objective

      integer :: icycle, ist, iexit
      FLOAT :: objective_new
      FLOAT, allocatable :: theta_new(:), df(:)
 
      PUSH_SUB(scf_occ.multid_minimize)

      SAFE_ALLOCATE(theta_new(1:nst))
      SAFE_ALLOCATE(df(1:nst))

      objective = M_ZERO
      df = M_ZERO
      objective_new = CNST(-1.0e-8)
      theta_new = theta

      do icycle = 1, max_iter
        if (icycle /= 1) then
          if (objective_new < objective) then
            rdm%step = CNST(1.3)*rdm%step
            objective = objective_new
            theta = theta_new
          else
            rdm%step = CNST(0.9)*rdm%step
          end if
        end if

        do ist = 1, nst
          theta_new(ist) = theta(ist) - rdm%step*df(ist)
        end do
           
        call calcul_objective(nst, theta_new, df, objective_new)

        iexit = 0
        do ist = 1, nst 
          if(abs(df(ist)) < rdm%toler)  iexit = iexit + 1 
        end do
        if (iexit == nst) exit
        cycle
      end do
   
      if (iexit /= nst) then
        write(message(1),'(a)') 'Did not manage to minimize the energy with respect to occupation numbers for this mu'
        call messages_info(1)
      end if

      objective = objective_new
      theta = theta_new
 
      SAFE_DEALLOCATE_A(theta_new)
      SAFE_DEALLOCATE_A(df)

      POP_SUB(scf_occ.multid_minimize)
  
    end subroutine multid_minimize
   
    ! --------------------------------------------
    
    subroutine calcul_objective(nst, theta_new, df, objective_new) 
      integer, intent(in)  :: nst
      FLOAT,   intent(in)  :: theta_new(:)
      FLOAT,   intent(out) :: df(:) !< (1:nst)
      FLOAT,   intent(out) :: objective_new

      integer :: ist
      FLOAT ::  energy
      FLOAT, allocatable :: V_h(:), V_x(:), dE_dn(:),occ(:) 
 
      PUSH_SUB(scf_occ.calcul_objective)

      SAFE_ALLOCATE(V_h(1:nst))
      SAFE_ALLOCATE(V_x(1:nst))
      SAFE_ALLOCATE(dE_dn(1:nst))
      SAFE_ALLOCATE(occ(1:nst))

      V_h = M_ZERO
      V_x = M_ZERO
      occ = M_ZERO

      do ist = 1, nst
        occ(ist) = M_TWO*sin(theta_new(ist)*M_PI*M_TWO)**2
      end do

      rdm%occsum = M_ZERO
      do ist = 1, nst
        rdm%occsum = rdm%occsum + occ(ist)
      end do
      
      !calculate the total energy without nuclei interaction and the energy
      !derivatives with respect to the occupation numbers

      call total_energy_rdm(rdm, st,gr, occ, energy, dE_dn)

      do ist = 1, nst
        df(ist) = M_FOUR*M_PI*sin(M_FOUR*theta_new(ist)*M_PI)*(dE_dn(ist)-rdm%mu)
      end do

      objective_new = energy - rdm%mu*(rdm%occsum - rdm%qtot)

      SAFE_DEALLOCATE_A(V_h)
      SAFE_DEALLOCATE_A(V_x)
      SAFE_DEALLOCATE_A(dE_dn)
      SAFE_DEALLOCATE_A(occ)

      POP_SUB(scf_occ.calcul_objective)

    end subroutine calcul_objective

  end subroutine scf_occ
   
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

    SAFE_ALLOCATE(lambda(1:st%nst,1:st%nst)) !matrix of Lagrange Multipliers from Piris paper Equation (8) 
    SAFE_ALLOCATE(FO(1:st%nst, 1:st%nst))

    lambda = M_ZERO
    FO = M_ZERO
    call density_calc (st,gr,st%rho)
    call v_ks_calc(ks, hm,st,geo)
    call hamiltonian_update(hm, gr%mesh)
    call construct_f(hm,st,gr,FO)
 
    if (rdm%iter==1) then
      do ist = 1, st%nst
        do jst = 1, ist
          lambda(ist, jst) = M_HALF*(FO(ist, jst) + FO(jst, ist))
	  lambda(jst, ist) = lambda(ist, jst)
        enddo
      end do
    else 
      do ist=1, st%nst
        lambda(ist, ist) = rdm%evalues(ist)
        do jst=1, ist-1
          lambda(jst, ist) =-( FO(jst, ist) - FO(ist ,jst))
          if(abs(lambda(jst,ist)) > rdm%scale_f) then
            lambda(jst, ist)= rdm%scale_f*lambda(jst,ist)/abs(lambda(jst, ist))
   	  endif
          lambda(ist, jst) = lambda(jst, ist) 
        enddo
      enddo
    endif

    call lalg_eigensolve(st%nst, lambda, rdm%evalues)
    call assign_eigfunctions(st, gr, lambda)
      

    call rdm_derivatives(rdm, hm, st, gr)
    call total_energy_rdm(rdm, st, gr, st%occ(:,1), energy)
    

    SAFE_DEALLOCATE_A(lambda) 
    SAFE_DEALLOCATE_A(FO) 

    POP_SUB(scf_orb)

  end subroutine scf_orb

  ! ----------------------------------------
  subroutine construct_f(hm, st, gr, lambda)
    type(hamiltonian_t),  intent(in) :: hm
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in) :: gr
    FLOAT,                intent(out):: lambda(1:st%nst, 1:st%nst)
      
    FLOAT, allocatable :: hpsi(:,:), hpsi1(:,:), dpsi(:,:), dpsi2(:,:) 
    FLOAT, allocatable :: g_x(:,:), g_h(:,:), rho(:), pot(:)
    integer :: ist, kst, ip

    PUSH_SUB(construct_f)

    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part,1:st%d%dim))
    SAFE_ALLOCATE(hpsi1(1:gr%mesh%np_part,1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part ,1:st%d%dim))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part ,1:st%d%dim))
    SAFE_ALLOCATE(g_x(1:st%nst,1:st%nst))
    SAFE_ALLOCATE(g_h(1:st%nst,1:st%nst))
    SAFE_ALLOCATE(rho(1:gr%mesh%np_part)) 
    SAFE_ALLOCATE(pot(1:gr%mesh%np_part))

    hpsi = M_ZERO
    hpsi1 = M_ZERO
    dpsi2 = M_ZERO
    lambda = M_ZERO
    g_x = M_ZERO
    g_h = M_ZERO
    rho = M_ZERO    
    pot = M_ZERO
    dpsi = M_ZERO
    dpsi2 = M_ZERO
        
    !calculate the single-particle part of the lambda matrix, Eq. (9)
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      call dhamiltonian_apply(hm,gr%der, dpsi, hpsi, ist, 1, &
                            & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
      do kst=1, ist
        call states_get_state(st, gr%mesh, kst, 1, dpsi2)
        lambda(kst, ist) = dmf_dotp(gr%mesh, dpsi2(:,1), hpsi(:,1))
        lambda(ist, kst) = lambda(kst, ist)
      end do
    end do
      
    do kst = 1, st%nst 
      do ist = 1, st%nst
        lambda(kst, ist) = lambda(kst, ist)*st%occ(ist,1)
      end do
    end do

    !calculate the Hartree contribution to lambda
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      forall (ip=1:gr%mesh%np_part)
        rho(ip) = rho(ip) + st%occ(ist,1)*dpsi(ip, 1)**2
      end forall
    end do

    call dpoisson_solve(psolver, pot , rho) !the Hartree potential
    
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      forall (ip=1:gr%mesh%np_part)
       dpsi(ip,1) = st%occ(ist, 1)*pot(ip)*dpsi(ip,1)
      end forall
      do kst = 1, st%nst  
        call states_get_state(st, gr%mesh, kst, 1, dpsi2)
        g_h(ist, kst) = dmf_dotp(gr%mesh, dpsi(:,1), dpsi2(:, 1))
      end do
    end do

    do kst = 1,st%nst
      do ist = 1,st%nst
        lambda(kst,ist) = lambda(kst,ist) + g_h(ist,kst)
      end do
    end do
    
    !calculate the exchange part of lambda
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      call dhamiltonian_apply(hm, gr%der, dpsi, hpsi, ist, 1, &
                            & terms = TERM_OTHERS)
      do kst = 1, st%nst
        call states_get_state(st, gr%mesh, kst, 1, dpsi2)
        g_x(ist,kst) = dmf_dotp(gr%mesh, dpsi2(:,1), hpsi(:,1))
        g_x(ist,kst) = sqrt(st%occ(ist,1))*g_x(ist,kst)
      end do
    end do
    
    do kst=1,st%nst
      do ist=1,st%nst
        lambda(kst,ist) = lambda(kst,ist) + g_x(ist,kst)
      end do
    end do

    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(hpsi1)
    SAFE_DEALLOCATE_A(dpsi2)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(g_x)
    SAFE_DEALLOCATE_A(g_h)
    SAFE_DEALLOCATE_A(rho) 
    SAFE_DEALLOCATE_A(pot) 
   
    POP_SUB(construct_f)

  end subroutine construct_f
   
  ! ----------------------------------------
  subroutine assign_eigfunctions(st, gr, lambda)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(in)    :: gr
    FLOAT,                intent(in)    :: lambda(:, :)
    
    type(states_t)   :: psi2
    FLOAT, allocatable :: occ(:,:)
    integer :: ist, jst

    PUSH_SUB(assign_eigenfunctions)
    
    SAFE_ALLOCATE(occ(1:st%nst,1:st%d%nik))

    occ = M_ZERO

    call states_copy(psi2, st)
    if(states_are_real(st)) then
      call dstates_rotate(gr%mesh, st, psi2, transpose(lambda))
    else
      call zstates_rotate(gr%mesh, st, psi2, M_z0 * transpose(lambda))
    endif
 
   ! reordering occupation numbers if needed
    occ = st%occ

    do ist = 1, st%nst
      if (abs(dmf_dotp(gr%mesh,st%dpsi(:,1,ist,1),psi2%dpsi(:,1,ist,1))) < M_HALF) then
        do jst = 1, st%nst 
          if (abs(dmf_dotp(gr%mesh,st%dpsi(:,1,ist,1),psi2%dpsi(:,1,jst,1))) >= M_HALF) then
            occ(ist,1) = st%occ(jst,1)
          end if 
        end do
      end if
    end do
  
  
    call states_end(psi2) 

    SAFE_DEALLOCATE_A(occ)
    
    POP_SUB(assign_eigenfunctions)

  end subroutine assign_eigfunctions
   
  ! ----------------------------------------
  subroutine total_energy_rdm(rdm, st, gr, occ, energy, dE_dn)
    type(rdm_t),          intent(inout)  :: rdm
    type(states_t),       intent(in)     :: st 
    type(grid_t),         intent(in)     :: gr
    FLOAT,                intent(in)     :: occ(:)
    FLOAT,                intent(out)    :: energy
    FLOAT, optional,      intent(out)    :: dE_dn(:) !< (1:st%nst)
     
    integer :: ist, jst
    FLOAT, allocatable :: hpsi(:,:), hpsi1(:,:), dpsi(:,:) 
    FLOAT, allocatable :: V_h(:), V_x(:)
     
    PUSH_SUB(total_energy_rdm)
  
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(hpsi1(1:gr%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part ,1:st%d%dim))
    SAFE_ALLOCATE(V_h(1:st%nst))
    SAFE_ALLOCATE(V_x(1:st%nst))

    hpsi = M_ZERO
    hpsi1 = M_ZERO
    dpsi = M_ZERO
    energy = M_ZERO
    V_h = M_ZERO
    V_x = M_ZERO
     
    !Calculate hartree contribution 
    do ist = 1, st%nst
      do jst = 1, st%nst
        V_h(ist) = V_h(ist) + occ(jst)*rdm%hartree(ist, jst)
      enddo
    end do

    !Calculate exchange contribution
    do ist = 1, st%nst 
      do jst = 1, st%nst
        V_x(ist) = V_x(ist) - sqrt(occ(jst))*rdm%exchange(ist, jst)
      end do
      V_x(ist) = V_x(ist)*M_HALF/max(sqrt(occ(ist)), CNST(1.0e-16))
    end do

    if (present(dE_dn)) then
      !Calculate the energy derivative with respect to the occupation numbers
      dE_dn(:) = rdm%eone(:) + V_h(:) + V_x(:)
    end if

    !Total energy calculation without nuclei interaction  
    do ist = 1, st%nst
      energy = energy + occ(ist)*rdm%eone(ist) &
                     & + M_HALF*occ(ist)*V_h(ist) & 
                     & + occ(ist)*V_x(ist)
    end do
    
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(hpsi1)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(V_h)
    SAFE_DEALLOCATE_A(V_x)
    
    POP_SUB(total_energy_rdm)
   
  end subroutine total_energy_rdm
  
  ! ----------------------------------------
  subroutine rdm_derivatives(rdm, hm, st, gr)
    type(rdm_t),          intent(inout) :: rdm
    type(hamiltonian_t),  intent(in)    :: hm 
    type(states_t),       intent(in)    :: st 
    type(grid_t),         intent(in)    :: gr
    
    FLOAT, allocatable :: hpsi(:,:), pot(:), rho(:), dpsi(:,:), dpsi2(:,:)
    integer :: ist, jst
    PUSH_SUB(rdm_derivatives) 

    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(pot (1:gr%mesh%np))
    SAFE_ALLOCATE(rho (1:gr%mesh%np))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np_part, 1:st%d%dim))

    !derivative of one-electron energy with respect to the natural orbitals occupation number
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      call dhamiltonian_apply(hm,gr%der, dpsi, hpsi, ist, 1, &
                            & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
      rdm%eone(ist) = dmf_dotp(gr%mesh, dpsi(:, 1), hpsi(:, 1))
    end do
    
    !calculates the integrals used for the hartree part of the total energy and its derivative
    do ist = 1, st%nst
      pot = M_ZERO
      rho = M_ZERO
      call states_get_state(st, gr%mesh, ist, 1, dpsi)
      rho(1:gr%mesh%np) = dpsi(1:gr%mesh%np, 1)**2
      ! FIXME: poisson solves here should probably be all_nodes = .false.
      call dpoisson_solve(psolver, pot , rho)
      
      do jst = 1, ist
        call states_get_state(st, gr%mesh, jst, 1, dpsi)
        rdm%hartree(jst, ist) = dmf_dotp(gr%mesh, dpsi(:, 1)**2, pot(:))
        rdm%hartree(ist, jst) = rdm%hartree(jst, ist)
      end do
    end do
   
    !calculates the integrals used for the exchange part of the total energy and its derivative
    do ist = 1, st%nst 
      pot = M_ZERO
      rho = M_ZERO
      dpsi2 = M_ZERO      
      call states_get_state(st, gr%mesh, ist, 1, dpsi2)
      
      do jst = 1, ist
        dpsi = M_ZERO      
        call states_get_state(st, gr%mesh, jst, 1, dpsi)
        rho(1:gr%mesh%np) = dpsi2(1:gr%mesh%np, 1)*dpsi(1:gr%mesh%np, 1)
        call dpoisson_solve(psolver, pot, rho)
        rdm%exchange(ist, jst) = dmf_dotp(gr%mesh,rho, pot)
        rdm%exchange(jst, ist) = rdm%exchange(ist, jst)
      end do
    end do

    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(pot)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)
    POP_SUB(rdm_derivatives) 

  end subroutine rdm_derivatives

end module rdmft_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

