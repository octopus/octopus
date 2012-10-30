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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 
!! 02111-1307, USA. 
!! 
!! $Id: scf.F90 9160 2012-06-23 20:38:20Z xavier $ 

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
  use messages_m
  use mesh_function_m
  use loct_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use system_m
  use unit_m
  use unit_system_m
 
  implicit none

  private
  public ::                &
    rdmft_init,            &
    rdmft_end,             &
    scf_occ


  type rdm_t
    integer  :: max_iter
    FLOAT    :: mu
    FLOAT, ALLOCATABLE :: eone(:)
    FLOAT, ALLOCATABLE :: hartree(:,:)
    FLOAT, ALLOCATABLE :: exchange(:,:)
    FLOAT, ALLOCATABLE :: V_h(:)
    FLOAT, ALLOCATABLE :: V_x(:)
    
    type(eigensolver_t) :: eigens
    ! shortcuts
    type(hamiltonian_t), pointer :: hm
    type(grid_t),        pointer :: gr
    type(states_t),      pointer :: st
    integer                      :: dim
    integer                      :: size
  end type rdm_t

  type(rdm_t), save :: rdm  
  

contains
   
  subroutine rdmft_init(sys, hm)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm

    PUSH_SUB(rdmft_init)  

    rdm%gr     => sys%gr
    rdm%st     => sys%st
    rdm%hm     => hm
    rdm%dim    =  sys%gr%mesh%sb%dim
    rdm%size = rdm%st%nst
    
    SAFE_ALLOCATE(rdm%eone(1:rdm%st%nst))
    SAFE_ALLOCATE(rdm%hartree(1:rdm%st%nst, 1:rdm%st%nst))
    SAFE_ALLOCATE(rdm%exchange(1:rdm%st%nst, 1:rdm%st%nst))
    SAFE_ALLOCATE(rdm%V_h(1:rdm%st%nst))
    SAFE_ALLOCATE(rdm%V_x(1:rdm%st%nst))

    rdm%eone = M_ZERO
    rdm%hartree = M_ZERO
    rdm%exchange = M_ZERO
    
    POP_SUB(rdmft_init)

  end subroutine rdmft_init

  subroutine rdmft_end()

    PUSH_SUB(rdmft_end)

    SAFE_DEALLOCATE_A(rdm%eone)
    SAFE_DEALLOCATE_A(rdm%hartree)
    SAFE_DEALLOCATE_A(rdm%exchange)
    SAFE_DEALLOCATE_A(rdm%V_h)
    SAFE_DEALLOCATE_A(rdm%V_x)

    POP_SUB(rdmft_end)


  end subroutine rdmft_end

  !> optimization of occupation numbers according to arXiv:1208.4699
  ! ----------------------------------------------------------------
  subroutine scf_occ(gr, hm, st, sys) 
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st
    type(system_t),       intent(inout) :: sys
    
    integer :: ist, jst
    integer :: icycle, iter, n_unpinned
    FLOAT :: abs_occ, occsum,  sumgi1, sumgi2, sumgim, mu1, mu2, mum, muinitial, muhat, dinterv 
    FLOAT :: el_per_state ! maximum number of electrons in each state
    FLOAT, allocatable :: occin(:,:), occout(:,:)
    FLOAT, allocatable :: hpsi(:,:), hpsi2(:,:), pot(:)
    FLOAT, allocatable :: rho(:), dpsi2(:,:)
    FLOAT, allocatable :: dE_dn(:), dE_dn2(:), ksmear(:), sigma(:)

    FLOAT, parameter :: conv = CNST(0.00001)
    FLOAT, parameter :: smallocc = CNST(0.001)

    PUSH_SUB(scf_occ)

    call rdmft_init(sys,hm)
    
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(occout(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi2(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(pot (1:gr%mesh%np))
    SAFE_ALLOCATE(rho (1:gr%mesh%np))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np ,1:st%d%dim))
    SAFE_ALLOCATE(dE_dn(1:st%nst))
    SAFE_ALLOCATE(dE_dn2(1:st%nst))
    SAFE_ALLOCATE(ksmear(1:st%nst))
    SAFE_ALLOCATE(sigma(1:st%nst))

    occout = M_ZERO
    hpsi = M_ZERO
    sigma = M_ZERO
    dE_dn = M_ZERO 
    dE_dn2 = M_ZERO
    ksmear = 0.5
 
    if(hm%d%ispin == 1) & ! unpolarized
      el_per_state = M_TWO
    if (hm%d%ispin ==2) & 
     el_per_state = M_ONE
 
    occout(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occout(:,:) < smallocc) occout(:,:)=smallocc 
    where(occout(:,:) > el_per_state-smallocc) occout(:,:)=el_per_state-smallocc

    if (rdm%hm%d%ispin.ne.1) then
      call messages_not_implemented("RDMFT not yet implemented for spin_polarized or spinors")
    end if


    !derivative of one electron energy with respect to the natural orbitals occupation number
    do ist = 1, st%nst
      call dhamiltonian_apply(hm,gr%der,st%dpsi(:,:,ist, 1), hpsi, ist, 1, &
                            & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
      rdm%eone(ist) = dmf_dotp(gr%mesh,st%dpsi(:,1,ist,1), hpsi(:,1))
    enddo

    !calculates the integrals used for the hartree part of the total energy and its derivative
    do ist = 1, st%nst
      pot = M_ZERO
      rho = M_ZERO
      rho(:) = st%dpsi(:, 1, ist, 1)*st%dpsi(:, 1, ist, 1)
      call dpoisson_solve (psolver, pot, rho)
      do jst = ist, st%nst
        rdm%hartree(jst, ist) = dmf_dotp(gr%mesh,st%dpsi(:,1,jst,1)**2 , pot)
        rdm%hartree(ist, jst) = rdm%hartree(jst, ist)  
      enddo
    enddo
    !calculates the integrals used for the exchange part of the total energy and its derivative
    do ist= 1, st%nst 
      pot=M_ZERO
      rho=M_ZERO
      call states_get_state(st, gr%mesh, ist, 1, dpsi2)
      do jst = ist, st%nst
        rho(:) = dpsi2(:, 1)*st%dpsi(:, 1, jst, 1)
        call dpoisson_solve(psolver, pot, rho)
        rdm%exchange(ist, jst) = dmf_dotp(gr%mesh, rho , pot)
        rdm%exchange(jst, ist) = rdm%exchange(ist, jst)
      end do
    end do

    dinterv = M_HALF
    muinitial = rdm%st%eigenval(int(st%qtot*M_HALF), 1) !initial guess for mu
    rdm%max_iter = 1000
    do icycle = 1, rdm%max_iter  !finding the eigenvalues of the effective diagonal hamiltonian
      rdm%st%occ = occout
      rdm%mu = mum
      call energy_derivatives(dE_dn, dE_dn2)
      sigma(:) = (log(max(el_per_state - rdm%st%occ(:, 1), 1d-16)/max(rdm%st%occ(:, 1), 1d-16)))
      do ist = 1, st%nst
        rdm%st%eigenval(ist, 1) = dE_dn(ist) + sigma(ist)*ksmear(ist)
      end do
      mu1 = muinitial
      mu2 = -smallocc
      !check if initial interval contains the root, and broaden interval if necessary
      call fill_occupations(rdm%st%eigenval, mu1, ksmear, el_per_state, occout, occsum, sumgi1)
      call fill_occupations(rdm%st%eigenval, mu2, ksmear, el_per_state, occout, occsum, sumgi2)
      do while (sumgi1*sumgi2.gt.M_ZERO)
        if(sumgi2.gt.M_ZERO) then
          mu2 = mu1
          sumgi2 = sumgi1
          mu1 = mu1 - dinterv
          call fill_occupations(rdm%st%eigenval, mu1, ksmear, el_per_state, occout, occsum, sumgi1)
        else
          mu1 = mu2
          sumgi1 = sumgi2
          mu2 = mu2 + dinterv
          call fill_occupations(rdm%st%eigenval, mu2, ksmear, el_per_state, occout, occsum, sumgi2)
        endif
      enddo

      do iter = 1, rdm%max_iter !bisection to find the root of sumocc-st%qtot=M_ZERO, for fixed set of occupation numbers find mu
        mum = (mu1 + mu2)*M_HALF
        call fill_occupations(rdm%st%eigenval, mum, ksmear, el_per_state, occout, occsum, sumgim)        
        if (sumgi1*sumgim.lt.M_ZERO) then
          mu2 = mum
        else
          mu1 = mum
          sumgi1=sumgim
        end if
        if (abs(sumgim).lt.conv.or.abs((mu1-mu2)*M_HALF).lt.conv)  exit
        cycle
      end do
      abs_occ = M_ZERO
      n_unpinned = M_ZERO
      muhat = M_ZERO
      do ist = 1, rdm%st%nst
        if (occout(ist, 1).lt.el_per_state-conv .and. occout(ist, 1).gt.conv) then
          muhat = muhat + dE_dn(ist)
          n_unpinned = n_unpinned + 1
        end if
      end do
      muhat = muhat/n_unpinned

      do ist = 1, rdm%st%nst
        if (occout(ist,1).lt.el_per_state-conv .and. occout(ist,1).gt.conv) then
          abs_occ = abs_occ + (dE_dn(ist)-muhat)**2
        end if
      end do
      abs_occ=abs_occ/n_unpinned
      if (abs_occ.lt.conv) exit
      cycle
    end do !iteration to find optimal occupation numbers

    !total energy without nuclei interaction  
    rdm%hm%energy%total = M_ZERO
    do ist = 1, rdm%st%nst
      rdm%hm%energy%total = rdm%hm%energy%total + rdm%st%occ(ist,1)*rdm%eone(ist) &
                                 & + M_HALF*rdm%st%occ(ist, 1)*rdm%V_h(ist) &
                                 & + rdm%st%occ(ist,1)*rdm%V_x(ist)
    end do

    
    write(message(1),'(a,1x,f11.5)') 'Occupations sum', occsum
    call messages_info(1)
    write(message(1),'(a,es15.8)') ' etot RDMFT= ', units_from_atomic(units_out%energy,rdm%hm%energy%total+rdm%hm%ep%eii) 
    write(message(2),'(a4,1x,a12)') '#st','Occupation'
    call messages_info(2)   
    do ist = 1,st%nst
      write(message(1),'(i4,3x,f11.5)') ist, occout(ist,1)
      call messages_info(1)
    end do

   call rdmft_end()

   SAFE_DEALLOCATE_A(occout)
   SAFE_DEALLOCATE_A(occin)
   SAFE_DEALLOCATE_A(hpsi)
   SAFE_DEALLOCATE_A(hpsi2)
   SAFE_DEALLOCATE_A(pot)
   SAFE_DEALLOCATE_A(rho)
   SAFE_DEALLOCATE_A(dpsi2)
   SAFE_DEALLOCATE_A(dE_dn)
   SAFE_DEALLOCATE_A(dE_dn2)
   SAFE_DEALLOCATE_A(ksmear)
   SAFE_DEALLOCATE_A(sigma)
STOP
   POP_SUB(scf_occ)
  end subroutine scf_occ
  
  subroutine fill_occupations(eigenval, mu, ksmear, el_per_state, occ, occsum, sumgi)  
    FLOAT,                intent(in)    :: eigenval(rdm%st%nst,rdm%st%d%nik)
    FLOAT,                intent(in)    :: mu
    FLOAT,                intent(inout) :: ksmear(1:rdm%st%nst)
    FLOAT,                intent(in)    :: el_per_state    
    FLOAT,                intent(out)   :: occ(1:rdm%st%nst,1)
    FLOAT,                intent(out)   :: occsum
    FLOAT,                intent(out)   :: sumgi
   
    INTEGER :: ik, ist 
    FLOAT :: xx, stepf

    PUSH_SUB(fill_occupations)

    do ik = 1, rdm%st%d%nik 
      do ist = 1, rdm%st%nst
        xx = (eigenval(ist, ik)-mu)/ksmear(ist)
        stepf = M_ONE/(M_ONE + exp(xx))
        occ(ist, ik) = stepf * el_per_state
      end do
    end do

    occsum = M_ZERO
    do ist = 1, rdm%st%nst
      occsum = occsum + occ(ist, 1)
    end do
    sumgi = occsum - rdm%st%qtot

    POP_SUB(fill_occupations)
  end subroutine fill_occupations

  subroutine energy_derivatives(dE_dn, dE_dn2) 
    FLOAT, intent(inout)         :: dE_dn(1:rdm%st%nst)             
    FLOAT, intent(inout)         :: dE_dn2(1:rdm%st%nst)
  
    INTEGER :: ist, jst
    FLOAT :: dVx2
       
    PUSH_SUB(energy_derivatives)
             
    !Calculate hartree contribution 
    rdm%V_h = M_ZERO
    do ist = 1, rdm%st%nst
      do jst = 1, rdm%st%nst
        rdm%V_h(ist) = rdm%V_h(ist) + rdm%st%occ(jst,1)*rdm%hartree(jst,ist)
      enddo
    end do
  
    !Calcualate exchange contribution
    rdm%V_x = M_ZERO
    do ist = 1, rdm%st%nst 
      do jst=1, rdm%st%nst
        rdm%V_x(ist) = rdm%V_x(ist) - sqrt(rdm%st%occ(jst,1))*rdm%exchange(ist,jst)
      end do
      rdm%V_x(ist) = rdm%V_x(ist)*M_HALF/max(sqrt(rdm%st%occ(ist,1)),1d-16)
    end do
 
    dE_dn(:) = rdm%eone(:) + rdm%V_h(:) + rdm%V_x(:)

    do ist = 1, rdm%st%nst
      dE_dn2(ist) = rdm%hartree(ist, ist)
    end do

    do ist = 1, rdm%st%nst
      dVx2 = M_ZERO
      do jst = 1, ist-1
        dVx2 = dVx2 + sqrt(rdm%st%occ(jst, 1))*rdm%exchange(ist, jst)
      end do
      do jst = ist+1, rdm%st%nst
        dVx2 = dVx2 + sqrt(rdm%st%occ(jst, 1))*rdm%exchange(ist, jst)
      end do
      dE_dn2(ist) = dE_dn2(ist) + (M_ONE/max((M_FOUR*rdm%st%occ(ist, 1)**1.5), 1d-16))*dVx2
    end do

    POP_SUB(energy_derivatives)
  end subroutine energy_derivatives

end module rdmft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
