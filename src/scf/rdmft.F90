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
  use energy_m
  use geometry_m
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
  public ::                   &
       rdmft_init,            &
       calc_point_rdmft,      &
       scf_occ


  type rdmft_opt_t
    integer  :: max_iter
    FLOAT    :: mu
    FLOAT, ALLOCATABLE :: eone(:), hartree(:,:), exchange(:,:)   
    ! shortcuts
    type(hamiltonian_t), pointer :: hm
    type(grid_t),        pointer :: gr
    type(states_t),      pointer :: st
    integer                      :: dim
    integer                      :: size
  end type rdmft_opt_t

  type(rdmft_opt_t) :: rdmft_opt  
  

  contains
   
  subroutine rdmft_init(sys,hm)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm

    PUSH_SUB(rdmft_init)  

    rdmft_opt%gr     => sys%gr
    rdmft_opt%st     => sys%st
    rdmft_opt%hm     => hm
    rdmft_opt%dim    =  sys%gr%mesh%sb%dim
    rdmft_opt%size = rdmft_opt%st%nst
    
    POP_SUB(rdmft_init)

   end subroutine rdmft_init

! Calculation of the energy and the derivative wrt. the occupation numbers
   
   subroutine calc_point_rdmft(size, theta, objective, getgrad, df)
    integer,     intent(in)    :: size
    REAL_DOUBLE, intent(inout) :: theta(size)
    REAL_DOUBLE, intent(inout) :: objective
    integer,     intent(in)    :: getgrad
    REAL_DOUBLE, intent(inout) :: df(size)
     
    FLOAT, allocatable :: V_h(:), dpsi2(:,:), V_x(:), dE_dn(:)  
    FLOAT :: occsum,u, objective_new
    integer :: ist, jst, icycle, iexit
    FLOAT ::  theta_new(size)

    PUSH_SUB(calc_point_rdmft)

    ASSERT(size == rdmft_opt%size)
    
    SAFE_ALLOCATE(V_h(1:rdmft_opt%st%nst))
    SAFE_ALLOCATE(V_x(1:rdmft_opt%st%nst))
    SAFE_ALLOCATE(dE_dn(1:rdmft_opt%st%nst))

    objective_new=-1d-8
    u=0.01
    theta_new=theta
    df=M_ZERO

    do icycle=1,rdmft_opt%max_iter

      if(objective_new.lt.objective) then
         u = 1.3*u
         objective = objective_new
         theta = theta_new
       else
         u = 0.9*u
       end if
 
      do ist=1,rdmft_opt%st%nst
        theta_new(ist)=theta(ist)-u*df(ist)
      end do

      do ist=1,rdmft_opt%st%nst
        rdmft_opt%st%occ(ist,1)=M_TWO*sin(theta_new(ist)*M_PI*M_TWO)**2
      end do

   !Calculate hartree contribution 
      V_h = M_ZERO
      do ist= 1,rdmft_opt%st%nst
        do jst= 1,rdmft_opt%st%nst
          V_h(ist) = V_h(ist)+rdmft_opt%st%occ(jst,1)*rdmft_opt%hartree(jst,ist)
        enddo
      end do

   !Calcualate exchange contribution
     V_x=M_ZERO
     do ist= 1,rdmft_opt%st%nst 
       V_x(ist)=M_ZERO
       do jst=1,rdmft_opt%st%nst
         V_x(ist)=V_x(ist)-sqrt(rdmft_opt%st%occ(jst,1))*rdmft_opt%exchange(ist,jst)
       end do
       V_x(ist)=V_x(ist)*M_HALF/max(sqrt(rdmft_opt%st%occ(ist,1)),1d-16)
     end do

    occsum = M_ZERO
    do ist = 1, rdmft_opt%st%nst
      occsum = occsum + rdmft_opt%st%occ(ist, 1)
    end do

  !Calculate the energy derivative with respect to the occupation numbers
    dE_dn=M_ZERO
    if (getgrad.eq.1) then
       dE_dn(:)=rdmft_opt%eone(:)+V_h(:)+V_x(:)
       do ist=1,rdmft_opt%st%nst
         df(ist)=M_FOUR*M_PI*sin(M_FOUR*theta_new(ist)*M_PI)*(dE_dn(ist) - rdmft_opt%mu)
       end do
    end if

  !Total energy calculation without nuclei interaction  
    rdmft_opt%hm%energy%total=M_ZERO
    do ist=1,rdmft_opt%st%nst
      rdmft_opt%hm%energy%total=rdmft_opt%hm%energy%total+rdmft_opt%st%occ(ist,1)*rdmft_opt%eone(ist)+&
                                 &M_HALF*rdmft_opt%st%occ(ist,1)*V_h(ist)+&
                                 &rdmft_opt%st%occ(ist,1)*V_x(ist)
     end do

     objective_new = rdmft_opt%hm%energy%total -rdmft_opt%mu*(occsum - rdmft_opt%st%qtot) 
     
      iexit=M_ZERO
      do ist=1,rdmft_opt%st%nst 
         if  (abs(df(ist)).lt.0.5d-5)  iexit=iexit+1 
      end do
      if (iexit==rdmft_opt%st%nst)  exit
    cycle
  end do
  if (iexit.ne.rdmft_opt%st%nst) then
   write(message(1),'(a)'), 'did not manage to minimize the energy for this mu'
   call messages_info(1)
  end if
    SAFE_DEALLOCATE_A(V_h)
    SAFE_DEALLOCATE_A(V_x)
    SAFE_DEALLOCATE_A(de_dn)

   POP_SUB(calc_point_rdmft)

   end subroutine calc_point_rdmft
   

    ! ---------------------------------------------------------
  subroutine scf_occ(gr, geo, hm, st, sys ) 
    type(geometry_t),     intent(in)    :: geo
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st
    type(system_t),       intent(inout) :: sys

    
    FLOAT :: conv_abs_occ, conv_rel_occ, conv_mu, abs_occ, rel_occ, mu, mu_old, energy,objective
    FLOAT :: occsum, smallocc, sumgi1, sumgi2, sumgim, mu1, mu2, mum
    FLOAT, allocatable :: occout(:,:), occin(:,:),theta(:),hpsi(:,:),pot(:),rho(:),dpsi2(:,:)
    integer :: ist, jst, ik, ip, ierr, getgrad, icycle
    FLOAT :: df(st%nst) 
    logical :: finish

    PUSH_SUB(scf_occ)
    
    !%ConvAbsOcc
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the occupation numbers <math> n_j </math> in a RDMFT calculation:
    !%
    !% <math> \epsilon = \sum_{j=1}^{M} \vert n_j^{out} - n_j^{inp} \vert </math>
    !%
    !% where <math> M </math> is the number of natural orbitals in the calculation.
    !% A zero value (the default) means do not use this criterion.
    !%End
    call parse_float(datasets_check('ConvAbsOcc'), M_ZERO, conv_abs_occ)
    !%ConvRelOcc
    !%Type float
    !%Default 1e-5
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the occupation numbers <math> n_j </math> in a RDMFT calculation:
    !%
    !% <math> \epsilon = {1\over N} ConvAbsOcc</math>.
    !%
    !%End
    call parse_float(datasets_check('ConvRelOcc'), CNST(1e-5), conv_rel_occ)

    SAFE_ALLOCATE(occout(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:st%nst))
    SAFE_ALLOCATE(rdmft_opt%eone(1:st%nst))
    SAFE_ALLOCATE(rdmft_opt%hartree(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(rdmft_opt%exchange(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(pot (1:gr%mesh%np))
    SAFE_ALLOCATE(rho (1:gr%mesh%np))
    SAFE_ALLOCATE(dpsi2(1:gr%mesh%np ,1:st%d%dim))

    !Initialize the occin. Smallocc should be no less than 1d-7 for numerical stability
    smallocc = 1d-7 
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:)=smallocc !isws na allaksw to where se kati allo an den uparxei allou
    where(occin(:,:) > 2.0-smallocc) occin(:,:)=2.0-smallocc
    occout = M_ZERO
    theta  = M_ZERO
    rdmft_opt%eone = M_ZERO
    hpsi = M_ZERO
    rdmft_opt%hartree = M_ZERO

    st%occ=occin

    call rdmft_init(sys,hm)
    
    if (rdmft_opt%hm%d%ispin.ne.1) then
      call messages_not_implemented("RDMFT exchange function not yet implemented for spin_polarized or spinors")
    end if
  
   !derivative of one electron energy with respect to the natural orbitals occupation number
    do ist = 1, rdmft_opt%st%nst
      call dhamiltonian_apply(rdmft_opt%hm,rdmft_opt%gr%der,rdmft_opt%st%dpsi(:,:,ist, 1), hpsi, ist, 1, &
                            & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
      rdmft_opt%eone(ist) = dmf_dotp(rdmft_opt%gr%mesh,rdmft_opt%st%dpsi(:,1,ist,1), hpsi(:,1))
    enddo

    !calculates the integrals used for the hartree part of the total energy and ist derivative
    pot=M_ZERO
    rho=M_ZERO
    do ist=1,rdmft_opt%st%nst
      pot=M_ZERO
      rho(:)=rdmft_opt%st%dpsi(:,1,ist,1)**2
      call dpoisson_solve (psolver, pot , rho)
      do jst=ist,rdmft_opt%st%nst
          rdmft_opt%hartree(jst,ist) = dmf_dotp(rdmft_opt%gr%mesh,rdmft_opt%st%dpsi(:,1,jst,1)**2, pot(:))
          rdmft_opt%hartree(ist,jst) = rdmft_opt%hartree(jst,ist)  
      enddo
    enddo
 
  !calculates the integrals used for the exchange part of the total energy and ist derivative
    do ist= 1, rdmft_opt%st%nst 
      do jst = ist, rdmft_opt%st%nst
        pot = M_ZERO
        rho = M_ZERO
        dpsi2 = M_ZERO

        rdmft_opt%exchange(ist,jst) = M_ZERO
        call states_get_state(rdmft_opt%st,rdmft_opt%gr%mesh, jst, 1, dpsi2)

        forall (ip = 1:rdmft_opt%gr%mesh%np)
          rho(ip) = rho(ip) + dpsi2(ip, 1)*rdmft_opt%st%dpsi(ip, 1,ist,1)
        end forall

        call dpoisson_solve(psolver, pot, rho)
        rdmft_opt%exchange(ist,jst)=dmf_dotp(rdmft_opt%gr%mesh, rdmft_opt%st%dpsi(:,1,ist,1)*dpsi2(:,1),pot)
        rdmft_opt%exchange(jst,ist)=rdmft_opt%exchange(ist,jst)
      end do
    end do

  !finding the chemical potential mu such that the occupation numbers sum up to the number of electrons
  !bisection to find the root of sumocc-st%qtot=M_ZERO
   rdmft_opt%max_iter=2500
   getgrad=1
   mu1=2.0d0*st%eigenval(int(st%qtot*M_HALF),1)   !initial guess for mu in the neighbourhood of homo
   mu2=0.1d0*st%eigenval(int(st%qtot*M_HALF),1) 
   do icycle=1,rdmft_opt%max_iter
     objective=M_ZERO
     getgrad=1
     rdmft_opt%mu=mu1

     do ist=1,st%nst
       theta(ist)=asin(sqrt(occin(ist,1)*M_HALF))*M_HALF/M_PI
     end do
    
     call calc_point_rdmft(rdmft_opt%size, theta, objective, getgrad, df)
    
     occsum=M_ZERO
     do ist=1,st%nst 
       occout(ist,1)=M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
       occsum=occsum+occout(ist,1)
     end do
    
     sumgi1=occsum-st%qtot

     objective=M_ZERO
     rdmft_opt%mu=mu2
     do ist=1,st%nst
       theta(ist)=asin(sqrt(occin(ist,1)*M_HALF))*M_HALF/M_PI
     end do
    
     call calc_point_rdmft(rdmft_opt%size, theta, objective, getgrad, df)
   
     occsum=M_ZERO
     do ist=1,st%nst 
       occout(ist,1)=M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
       occsum=occsum+occout(ist,1)
     end do
  
     sumgi2=occsum-st%qtot

     if (icycle ==1.and. sumgi1*sumgi2.gt.M_ZERO) then
      !broaden interval if the root is not in the initial one
       mu1=1.5d0*mu1
       mu2=0.3d0*mu2
     end if      

     objective=M_ZERO
     mum=(mu1+mu2)*M_HALF
     rdmft_opt%mu=mum
    
     do ist=1,st%nst
       theta(ist)=asin(sqrt(occin(ist,1)*M_HALF))*M_HALF/M_PI
     end do
    
     call calc_point_rdmft(rdmft_opt%size, theta, objective, getgrad, df)
   
     occsum=M_ZERO
      do ist=1,st%nst 
        occout(ist,1)=M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
        occsum=occsum+occout(ist,1)
      end do

     sumgim=occsum-st%qtot

     if (sumgi1*sumgim.lt.M_ZERO) then
       mu2=mum
     else
       mu1=mum
     end if

     if (abs(sumgim).lt.1d-7.or.abs((mu1-mu2)*M_HALF).lt.1d-7)  exit
     cycle
   end do

   write(message(1),'(a,1x,f11.6)'), 'Occupations sum', occsum
   call messages_info(1)
   write(message(1),'(a,es15.8)') ' etot RDMFT= ',   units_from_atomic(units_out%energy,objective+rdmft_opt%hm%ep%eii) 
   write(message(2),'(a4,1x,a12)')'#st','Occupation'
   call messages_info(2)   

   do ist = 1, st%nst
     write(message(1),'(i4,3x,f11.6)'), ist, occout(ist, 1)
     call messages_info(1)  
   end do
   
STOP
   !compute convergence criteria
   abs_occ = M_ZERO
   do ist = 1, st%nst
     do ik =1 , st%d%nik
       abs_occ = abs_occ + abs( occout(ist,ik) - occin(ist, ik))
     end do
   end do

   rel_occ = abs_occ / st%qtot
   ! are we finished?
   finish = &
           (conv_abs_occ  <= M_ZERO .or. abs_occ  <= conv_abs_occ)  .and. &
           (conv_rel_occ  <= M_ZERO .or. rel_occ  <= conv_rel_occ)



   SAFE_DEALLOCATE_A(occout)
   SAFE_DEALLOCATE_A(occin)
   POP_SUB(scf_occ)
 end subroutine scf_occ
  
 end module rdmft_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

