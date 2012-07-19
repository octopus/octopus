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

#include "global.h"

module rdmft_m
  use messages_m
  use excited_states_m
  use exponential_m
  use geometry_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lasers_m
  use loct_m
  use loct_math_m
  use mesh_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_dim_m
  use system_m
  use datasets_m
  use density_m
  use energy_calc_m
  use epot_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use loct_m
  use loct_math_m
  use parser_m
  use mesh_m
  use messages_m
  use profiling_m
  use restart_m
  use simul_box_m
  use species_pot_m
  use states_m
  use states_calc_m
  use system_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use varinfo_m
  use xyz_adjust_m
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
  use loct_math_m
  use parser_m
  use poisson_m
  use profiling_m
  use states_m
  use system_m
  
  implicit none

  private
  public ::                   &
       rdmft_init,            &
       mu_optimize,           &
       calc_point_rdmft,      &
       write_iter_info_rdmft, &
       scf_occ


  type rdmft_opt_t
    integer  :: method
    FLOAT    :: step
    FLOAT    :: tolgrad
    FLOAT    :: toldr
    integer  :: max_iter
    integer  :: what2minimize
    FLOAT    :: mu

    ! shortcuts
    type(geometry_t),    pointer :: geo
    type(hamiltonian_t), pointer :: hm
    type(system_t),      pointer :: sys
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

    rdmft_opt%gr     => sys%gr
    rdmft_opt%geo    => sys%geo
    rdmft_opt%st     => sys%st
    rdmft_opt%hm     => hm
    rdmft_opt%sys    => sys
    rdmft_opt%dim    =  sys%gr%mesh%sb%dim
    rdmft_opt%size   =  rdmft_opt%st%nst*rdmft_opt%st%d%nik

  end subroutine rdmft_init

!--------------------------------------------------------------------------------
! This subroutine finds the optimum chemical potential for the occupation numbers to sum to the correct particle number
! occupation number optimization is called from here
! NOTE: this part of the code is under development, i.e. it does not work yet.
  subroutine mu_optimize(states, size, theta, mu)
    type(states_t),       intent(inout) :: states
    integer,              intent(in)    :: size
    REAL_DOUBLE,          intent(inout) :: theta(size)
    FLOAT,                intent(inout) :: mu
   
    integer :: ierr, ist
    REAL_DOUBLE :: energy
    FLOAT :: occsum
    
    rdmft_opt%mu = mu 
 
    ierr = loct_minimize(2, states%nst, theta(1), real(M_HALF, 8),&
           real(CNST(0.001), 8), real(CNST(0.001), 8), 50 , &
           calc_point_rdmft, write_iter_info_rdmft, energy)

  
   
    occsum = M_ZERO
    do ist = 1, rdmft_opt%st%nst
      occsum = occsum + rdmft_opt%st%occ(ist, 1)
    enddo

    mu = rdmft_opt%mu - (occsum - rdmft_opt%st%qtot)

  end subroutine mu_optimize

!----------------------------------------------------------------------------------
! Calculation of the energy and the derivative wrt. the occupation numbers

  subroutine calc_point_rdmft(size, theta, objective, getgrad, df)
    integer,     intent(in)    :: size
    REAL_DOUBLE, intent(in)    :: theta(size)
    REAL_DOUBLE, intent(inout) :: objective
    integer,     intent(in)    :: getgrad
    REAL_DOUBLE, intent(inout) :: df(size)
     
    FLOAT, ALLOCATABLE :: hpsi(:,:)  
    FLOAT, allocatable :: eone(:,:), pot(:), rdmhartree(:,:), rho(:),dpsi2(:,:), Vx(:),DE_Dn(:,:)  
    FLOAT :: occsum
    integer :: ist, jst, ik, ip
    
    PUSH_SUB(calc_point_rdmft)
    ASSERT(size == rdmft_opt%size)
    
    do ik = 1, rdmft_opt%st%d%nik
      do ist = 1, rdmft_opt%st%nst
        rdmft_opt%st%occ(ist,ik) = M_TWO*sin(theta(ist + (ik-1)*rdmft_opt%st%nst)*M_PI*M_HALF)**2
      enddo
    enddo
    
    call density_calc(rdmft_opt%st, rdmft_opt%gr, rdmft_opt%st%rho)
    call energy_calc_total(rdmft_opt%hm, rdmft_opt%gr,rdmft_opt%st) 
    
    occsum = M_ZERO
    do ik = 1, rdmft_opt%st%d%nik
      do ist = 1, rdmft_opt%st%nst
        occsum = occsum + rdmft_opt%st%occ(ist, ik)
      enddo
    enddo
    
    objective = rdmft_opt%hm%energy%total - rdmft_opt%mu*(occsum - rdmft_opt%st%qtot) 

    if (getgrad.eq.1) then

    SAFE_ALLOCATE(hpsi( 1:rdmft_opt%gr%mesh%np , 1:rdmft_opt%st%d%dim ))
    SAFE_ALLOCATE(eone( 1:rdmft_opt%st%nst ,1:rdmft_opt%st%d%nik ))
    SAFE_ALLOCATE(pot ( 1:rdmft_opt%gr%mesh%np  ))
    SAFE_ALLOCATE(rdmhartree ( 1:rdmft_opt%st%nst ,1:rdmft_opt%st%d%nik ))
    SAFE_ALLOCATE(rho ( 1:rdmft_opt%gr%mesh%np ))
    SAFE_ALLOCATE(dpsi2 (1:rdmft_opt%gr%mesh%np ,1:rdmft_opt%st%d%dim) )
    SAFE_ALLOCATE(Vx(1:rdmft_opt%st%nst))
    SAFE_ALLOCATE(De_Dn( 1:rdmft_opt%st%nst ,1:rdmft_opt%st%d%nik ))
    hpsi = M_ZERO
    eone = M_ZERO
    pot = M_ZERO
    rdmhartree = M_ZERO
    rho = M_ZERO
    DE_Dn=M_ZERO
 
    !calculate Hartree potential and initialize xc potential to zero
    call dpoisson_solve (psolver, pot ,rdmft_opt%st%rho(:,1))
    Vx=M_ZERO

    if(rdmft_opt%st%d%dim == 1) then
      do ist=1,rdmft_opt%st%nst
        do ik=1,rdmft_opt%st%d%nik
          !one electron part of the energy
          call dhamiltonian_apply(rdmft_opt%hm,rdmft_opt%gr%der,rdmft_opt%st%dpsi(:,:,ist, ik), hpsi, ist, ik, &
                                 & terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
          eone(ist,ik) = dmf_dotp(rdmft_opt%gr%mesh,rdmft_opt%st%dpsi(:,1,ist,ik), hpsi(:,1))
          
          !Hartree contribution
          rdmhartree(ist,ik) = dmf_dotp(rdmft_opt%gr%mesh,rdmft_opt%st%dpsi(:,1,ist,ik)**2, pot(:))
          
          !xc part
          do jst = 1, rdmft_opt%st%nst
            pot = M_ZERO
            rho = M_ZERO
            dpsi2 = M_ZERO

            call states_get_state(rdmft_opt%st,rdmft_opt%gr%mesh, jst, ik, dpsi2)

            forall (ip = 1:rdmft_opt%gr%mesh%np)
              rho(ip) = rho(ip) + sqrt(rdmft_opt%st%occ(jst,1))*dpsi2(ip, 1)*rdmft_opt%st%dpsi(ip, 1,ist,ik)
            end forall

            call dpoisson_solve(psolver, pot, rho)

            forall (ip = 1:rdmft_opt%gr%mesh%np) pot(ip)=pot(ip)*dpsi2(ip, 1)
            Vx(ist) = Vx(ist) - dmf_dotp(rdmft_opt%gr%mesh, rdmft_opt%st%dpsi(:,1,ist,ik),pot)
          end do
          Vx(ist) = Vx(ist)/sqrt(rdmft_opt%st%occ(ist, 1))
          dE_dn(ist, ik)=eone(ist,ik)+rdmhartree(ist,ik)+Vx(ist)
        end do !ik 
      end do !ist
    else
      call messages_write ("RDMFT is not implemented yet for spin_polarized or spinors")
      call messages_fatal ()
    endif

    do ik = 1, rdmft_opt%st%d%nik
      do ist = 1, rdmft_opt%st%nst
        df(ist + (ik-1)*rdmft_opt%st%nst) = M_PI*sin(theta(ist + (ik-1)*rdmft_opt%st%nst)*M_PI) & 
	                                   *(dE_dn(ist, ik) - rdmft_opt%mu)
      enddo
    enddo

    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(eone)
    SAFE_DEALLOCATE_A(pot )
    SAFE_DEALLOCATE_A(rdmhartree )
    SAFE_DEALLOCATE_A(rho )
    SAFE_DEALLOCATE_A(dpsi2  )
    SAFE_DEALLOCATE_A(Vx)
    SAFE_DEALLOCATE_A(De_Dn)
    end if

    POP_SUB(calc_point_rdmft)

  end subroutine calc_point_rdmft
   
    ! ---------------------------------------------------------
   subroutine write_iter_info_rdmft(geom_iter, size, energy, maxdx, maxdf, coords)
    integer,     intent(in) :: geom_iter
    integer,     intent(in) :: size ! must equal dim * natoms
    REAL_DOUBLE, intent(in) :: energy, maxdx, maxdf
    REAL_DOUBLE, intent(in) :: coords(size)


    PUSH_SUB(write_iter_info_rdmft)

    message(1) = ""
    message(2) = ""
    message(3) = "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(message(4),'("+++++++++++++++++++++ MINIMIZATION ITER #: ",I4," ++++++++++++++++++++++")') geom_iter
    message(8) = message(3)
    message(9) = message(3)
    message(10) = ""
    message(11) = ""
    call messages_info(11)

    POP_SUB(write_iter_info_rdmft)
  end subroutine write_iter_info_rdmft

    ! ---------------------------------------------------------
  subroutine scf_occ(gr, geo, hm, st, sys ) 
    type(geometry_t),     intent(in)    :: geo
    type(grid_t),         intent(inout) :: gr
    type(hamiltonian_t),  intent(inout) :: hm
    type(states_t),       intent(inout) :: st
    type(system_t),       intent(inout) :: sys

    FLOAT :: conv_abs_occ, conv_rel_occ, conv_mu, abs_occ, rel_occ, mu, mu_old
    FLOAT :: smallocc !the minimum value occ can have for numerical stability
    FLOAT, allocatable :: occout(:,:), occin(:,:), occnew(:,:), theta(:)
    integer :: ist, jst, ik, ip, size
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

    size = st%nst*st%d%nik

    SAFE_ALLOCATE(occout(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(occnew(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:size))
    
    !Initialize the occin. Smallocc should be no less than 1d-8 for numerical stability
    smallocc = 1d-8
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:)=smallocc 
    where(occin(:,:) > 2.0-smallocc) occin(:,:)=2.0-smallocc
    occout = M_ZERO
    occnew = M_ZERO
    theta  = M_ZERO

    st%occ=occin

    call rdmft_init(sys,hm)
    
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        theta(ist + (ik-1)*st%nst) = asin(sqrt(occin(ist, ik)*M_HALF))*M_TWO/M_PI
      enddo
    enddo
    
    mu = -CNST(1e-2)
   ! scf cycle for mu the chemical potential 
    do
      mu_old=mu
      call mu_optimize (st, size, theta,mu)
      conv_mu=abs(mu - mu_old)
      if ((conv_mu).lt.real(CNST(1e-10)))  exit
    end do
    
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        occout(ist, ik) = M_TWO*sin(theta(ist + (ik-1)*st%nst)*M_PI*M_HALF)**2
      enddo
    enddo
     
    !compute convergence criteria
    abs_occ = M_ZERO
    do ist = 1, st%nst
      do ik =1 , st%d%nik
        abs_occ = abs_occ + abs(occout(ist, ik) - occin(ist, ik))
      end do
    end do

    rel_occ = abs_occ / st%qtot
    ! are we finished?
    finish = &
           (conv_abs_occ  <= M_ZERO .or. abs_occ  <= conv_abs_occ)  .and. &
           (conv_rel_occ  <= M_ZERO .or. rel_occ  <= conv_rel_occ)



    SAFE_DEALLOCATE_A(occout)
    SAFE_DEALLOCATE_A(occin)
    SAFE_DEALLOCATE_A(occnew)

    POP_SUB(scf_occ)
  end subroutine scf_occ
 
end module rdmft_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

