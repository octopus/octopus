!! Copyright (C) 2012 I. Theophilou, N. Helbig
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
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use lalg_adv_m
  use loct_m
  use loct_math_m 
  use messages_m
  use mesh_function_m
  use output_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use unit_m
  use unit_system_m
  use v_ks_m
 
  implicit none

  private
  public ::                      & 
       rdm_t,                    &
       rdmft_init,               &
       rdmft_end,                &
       scf_occ,                  &
       scf_rdmft


  type rdm_t
    integer  :: max_iter
    FLOAT    :: mu, occsum, qtot
    FLOAT, allocatable   :: eone(:), hartree(:,:), exchange(:,:)   
    REAL(8) :: step, toler
  end type rdm_t
  
  type(rdm_t), save :: rdm
  
contains
   
  subroutine rdmft_init(rdm, st)
    type(rdm_t),           intent(inout)    :: rdm
    type(states_t),        intent(inout)    :: st

    PUSH_SUB(rdmft_init)  

    SAFE_ALLOCATE(rdm%eone(1:st%nst))
    SAFE_ALLOCATE(rdm%hartree(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(rdm%exchange(1:st%nst, 1:st%nst))

    rdm%eone = M_ZERO
    rdm%hartree = M_ZERO
    rdm%exchange = M_ZERO
    rdm%mu = M_TWO*st%eigenval(int(st%qtot*M_HALF), 1)
    rdm%qtot = st%qtot
    rdm%occsum = M_ZERO

    !%Variable RDMMaxIter
    !%Type integer
    !%Default 3000
    !%Section SCF::RDMFT
    !%Description
    !% Even if the convergence criterion is not satisfied, the minimization will stop
    !% after this number of iterations. The default is 3000.
    !%End 
    call parse_integer(datasets_check('RDMMaxIter'), 3000, rdm%max_iter)
    
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
    
    POP_SUB(rdmft_init)

  end subroutine rdmft_init

  ! ----------------------------------------

  subroutine rdmft_end(rdm)
    type(rdm_t), intent(inout) :: rdm

    PUSH_SUB(rdmft_end)

    SAFE_DEALLOCATE_A(rdm%eone)
    SAFE_DEALLOCATE_A(rdm%hartree)
    SAFE_DEALLOCATE_A(rdm%exchange)

    POP_SUB(rdmft_end)

  end subroutine rdmft_end

  ! ----------------------------------------

  ! scf for the occupation numbers and the natural orbitals
  subroutine scf_rdmft(rdm, gr, hm, st)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st

    PUSH_SUB(scf_rdmft)

    call scf_occ(rdm, gr, hm, st) 
    
    POP_SUB(scf_rdmft)

  end subroutine scf_rdmft
  
  ! ---------------------------------------------------------
  
  subroutine scf_occ(rdm, gr, hm, st)
    type(rdm_t),          intent(inout) :: rdm
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st

    
    REAL(8),allocatable :: theta(:),df(:)
    REAL(8)::energy
    integer :: ist, jst, icycle
    FLOAT ::  sumgi1, sumgi2, sumgim, mu1, mu2, mum, dinterv
    FLOAT, allocatable :: occout(:,:), occin(:,:), hpsi(:,:), pot(:), rho(:), dpsi(:,:), dpsi2(:,:)
    FLOAT :: el_per_state !maximum number of electrons per state
    FLOAT, parameter :: smallocc = CNST(0.0000001) 

    PUSH_SUB(scf_occ)
    
    SAFE_ALLOCATE(occout(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:st%nst))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(pot (1:gr%mesh%np))
    SAFE_ALLOCATE(rho (1:gr%mesh%np))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(df(1:st%nst))
    occout = M_ZERO
    occin = M_ZERO
    theta  = M_ZERO
    hpsi = M_ZERO
    energy = M_ZERO

    if(hm%d%ispin == 1) el_per_state = M_TWO ! unpolarized
    if(hm%d%ispin == 2) el_per_state = M_ONE

    !Initialize the occin. Smallocc is used for numerical stability
    
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:) = smallocc 
    where(occin(:,:) > el_per_state-smallocc) occin(:,:) = el_per_state - smallocc

    st%occ = occin
    theta(:) = asin(sqrt(occin(:, 1)/M_TWO))*(M_HALF/M_PI)
    
    if (hm%d%ispin /= 1) then
      call messages_not_implemented("RDMFT exchange function not yet implemented for spin_polarized or spinors")
    end if
  
    !derivative of one electron energy with respect to the natural orbitals occupation number
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
      
      do jst = ist, st%nst
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
      
      do jst = ist, st%nst
        call states_get_state(st, gr%mesh, jst, 1, dpsi)
        rho(1:gr%mesh%np) = dpsi2(1:gr%mesh%np, 1)*dpsi(1:gr%mesh%np, 1)
        call dpoisson_solve(psolver, pot, rho)
        rdm%exchange(ist, jst) = dmf_dotp(gr%mesh,rho, pot)
        rdm%exchange(jst, ist) = rdm%exchange(ist, jst)
      end do
    end do

    !finding the chemical potential mu such that the occupation numbers sum up to the number of electrons
    !bisection to find the root of rdm%occsum-st%qtot=M_ZERO
    mu1 = rdm%mu   !initial guess for mu 
    mu2 = -CNST(1.0e-6)
    dinterv = M_HALF
    call  multid_minimize(st%nst, rdm%max_iter, theta, energy)
    sumgi1 = rdm%occsum - st%qtot 
    rdm%mu = mu2
    call  multid_minimize(st%nst, rdm%max_iter, theta, energy) 
    sumgi2 = rdm%occsum - st%qtot

    do while (sumgi1*sumgi2 > M_ZERO)
      if (sumgi2 > M_ZERO) then
        mu2 = mu1
        sumgi2 = sumgi1
        mu1 = mu1 - dinterv
        rdm%mu = mu1
        call  multid_minimize(st%nst, rdm%max_iter, theta, energy) 
        sumgi1 = rdm%occsum - st%qtot 
      else
        mu1 = mu2
        sumgi1 = sumgi2
        mu2 = mu2 + dinterv
        rdm%mu = mu2
        call  multid_minimize(st%nst, rdm%max_iter, theta, energy) 
        sumgi2 = rdm%occsum - st%qtot 
      end if
    end do

    do icycle = 1, rdm%max_iter
      mum = (mu1 + mu2)*M_HALF
      rdm%mu = mum
      call  multid_minimize(st%nst, rdm%max_iter, theta, energy) 
      sumgim = rdm%occsum - st%qtot
      if (sumgi1*sumgim < M_ZERO) then
        mu2 = mum
      else
        mu1 = mum
        sumgi1 = sumgim
      end if
      if(abs(sumgim) < rdm%toler*rdm%step .or. abs((mu1-mu2)*M_HALF) < rdm%toler*rdm%step)  exit
      cycle
    end do
    if (icycle >= rdm%max_iter) then
      write(message(1),'(a,1x,f11.6)') 'Bisection ended without finding mu, sum of occupation numbers:', rdm%occsum
      call messages_fatal(1)
    endif

    occout = M_ZERO
    do ist = 1, st%nst
      occout(ist, 1) = M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
    end do

    write(message(1),'(a,1x,f11.6)') 'Occupations sum', rdm%occsum
    call messages_info(1)
    write(message(1),'(a,es15.8)') ' etot RDMFT= ', units_from_atomic(units_out%energy,energy+hm%ep%eii) 
    write(message(2),'(a4,1x,a12)')'#st','Occupation'
    call messages_info(2)   

    do ist = 1, st%nst
      write(message(1),'(i4,3x,f11.6)') ist, occout(ist, 1)
      call messages_info(1)  
    end do
   
    st%occ = occout

    SAFE_DEALLOCATE_A(occout)
    SAFE_DEALLOCATE_A(occin)
    SAFE_DEALLOCATE_A(theta)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(pot)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(dpsi2)
    
    POP_SUB(scf_occ)

    contains
  
    subroutine multid_minimize(nst, max_iter, theta, objective) 
      integer, intent(in)        :: nst
      integer, intent(in)        :: max_iter
      FLOAT, intent(inout)       :: theta(1:nst)
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
        if(objective_new < objective) then
          rdm%step = CNST(1.3)*rdm%step
          objective = objective_new
          theta = theta_new
        else
          rdm%step = CNST(0.9)*rdm%step
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
        write(message(1),'(a)') 'Did not manage to minimize '
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
      INTEGER, intent(in)                 :: nst
      FLOAT, intent(in)                   :: theta_new(1:nst)
      FLOAT, intent(out)                  :: df(1:nst)
      FLOAT, intent(out)                  :: objective_new

      integer :: ist, jst
      FLOAT ::  energy_tot
      FLOAT, allocatable :: V_h(:), V_x(:), dE_dn(:),occ(:) 
 
      PUSH_SUB(scf_occ.calcul_objective)

      SAFE_ALLOCATE(V_h(1:nst))
      SAFE_ALLOCATE(V_x(1:nst))
      SAFE_ALLOCATE(dE_dn(1:nst))
      SAFE_ALLOCATE(occ(1:nst))
      V_h = M_ZERO
      V_x = M_ZERO
      dE_dn = M_ZERO
      occ = M_ZERO

      do ist = 1, nst
        occ(ist) = M_TWO*sin(theta_new(ist)*M_PI*M_TWO)**2
      end do

      !Calculate hartree contribution 
      do ist = 1,nst
        do jst = 1, nst
          V_h(ist) = V_h(ist) + occ(jst)*rdm%hartree(jst, ist)
        enddo
      end do

      !Calculate exchange contribution
      do ist = 1, nst 
        do jst = 1, nst
          V_x(ist) = V_x(ist) - sqrt(occ(jst))*rdm%exchange(ist, jst)
        end do
        V_x(ist) = V_x(ist)*M_HALF/max(sqrt(occ(ist)), CNST(1.0e-16))
      end do

      rdm%occsum = M_ZERO
      do ist = 1, nst
        rdm%occsum = rdm%occsum + occ(ist)
      end do

      !Calculate the energy derivative with respect to the occupation numbers
      dE_dn(:) = rdm%eone(:) + V_h(:) + V_x(:)
      do ist = 1, nst
        df(ist) = M_FOUR*M_PI*sin(M_FOUR*theta_new(ist)*M_PI)*(dE_dn(ist)-rdm%mu)
      end do

      !Total energy calculation without nuclei interaction  
    
      energy_tot = M_ZERO
      do ist = 1, nst
        energy_tot = energy_tot + occ(ist)*rdm%eone(ist) + &
                                 & M_HALF*occ(ist)*V_h(ist) + &
                                 & occ(ist)*V_x(ist)
      end do
      objective_new = energy_tot - rdm%mu*(rdm%occsum - rdm%qtot)

      POP_SUB(scf_occ.calcul_objective)
    end subroutine calcul_objective

  end subroutine scf_occ
end module rdmft_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

