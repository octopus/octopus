!! Copyright (C) 2004 Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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

module rsdfpt

  use geometry
  use system
  use states
  use restart
  use hamiltonian 
  use eigen_solver
  use mix
  use mesh
  use poisson

  implicit none

  type rsdfpt_type
    ! response convergence to e-field/max iterations
    FLOAT :: conv_abs_dens_e  
    FLOAT :: abs_dens_e
    integer :: max_iter_e   

    ! response convergence to ion displacement/max iterations   
    FLOAT :: conv_abs_dens_r 
    FLOAT :: abs_dens_r    
    integer :: max_iter_r

    type(mix_type) :: mixer

    ! energy hessian 
    FLOAT, pointer :: rshess(:) 

    ! polarizability tensor 	
    FLOAT         :: rspolt(9) 
    ! mean static polarizability 1/3 trace of pol. tensor 
    FLOAT         :: msp

     ! density response
    FLOAT, pointer :: rsrho(:,:,:) 
    FLOAT, pointer :: rsrhoin(:,:,:)
    FLOAT, pointer :: rsrhonew(:,:,:) 
    
    ! KS orbitals response rsdpsi(m%np, 1, st%nst, st%d%nspin)
    FLOAT, pointer :: rsdpsi(:,:,:,:)

    ! change in hartree potential 
    FLOAT, pointer :: rsVhar(:)
    ! change in xc potential
    ! rsVxc(:,1,1) = \partial V_{xc}^{\alpha}/\partial \rho_{\alpha}
    ! rsVxc(:,1,2) = \partial V_{xc}^{\alpha}/\partial \rho_{\beta} 
    ! rsVxc(:,2,1) = \partial V_{xc}^{\beta}/\partial \rho_{\alpha} 
    ! rsVxc(:,2,2) = \partial V_{xc}^{\beta}/\partial \rho_{\beta} 
    FLOAT, pointer :: rsVxc(:,:,:)

    ! some working arrays...
    FLOAT, pointer :: Y0(:,:), Y(:,:)
  end type rsdfpt_type

contains

  integer function rsdfpt_run(sys, h, fromScratch) result(ierr)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical, intent(inout) :: fromScratch

    type(rsdfpt_type) :: lr

    ierr = 0
    call init_()

    ! load wave-functions
    if(X(restart_read) ("tmp/restart_gs", sys%st, sys%m).ne.sys%st%nst) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
      
      ierr = 1
      call end_()
      return
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    !if(.not.fromScratch) then ! try to load delta_psi
    !  if(X(restart_read) ("tmp/restart_lr_static_pol", sys%st, sys%m).ne.sys%st%nst) then
    !  
    !end if
    fromScratch = .true.

    call rsdfpt_e_init(sys%m, sys%st, h, lr)
    call pol_tensor(sys%m, sys%st, h, lr, sys%f_der)
    call rsdfpt_e_end(lr)

    call end_()

  contains

    subroutine init_()
      call push_sub('rsdfpt_run')

      ! allocate wfs
      allocate(sys%st%X(psi)(sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

    end subroutine init_

    subroutine end_()
      deallocate(sys%st%X(psi))
      
      call pop_sub()
    end subroutine end_
    
  end function rsdfpt_run


  subroutine rsdfpt_init(m, st, h, p)
    type(mesh_type),    intent(in)    :: m
    type(states_type),  intent(in)    :: st
    type(hamiltonian_type), intent(in) :: h
    type(rsdfpt_type),  intent(inout) :: p

    call push_sub('rsdfpt_init')

    allocate(p%rsrho(1, m%np, st%d%nspin))
    allocate(p%rsrhoin(1, m%np, st%d%nspin))
    allocate(p%rsrhonew(1, m%np, st%d%nspin))
    
    allocate(p%rsdpsi(m%np, st%dim, st%nst, st%d%nspin))
    allocate(p%rsVhar(m%np))
    allocate(p%rsVxc(m%np, st%d%nspin, st%d%nspin))

    ! the kernel does not change during the linear response calculation
    if(.not.h%ip_app) call build_fxc_kernel(m, st, p)

    allocate(p%Y0(m%np,st%dim), p%Y(m%np,st%dim))
    call mix_init(p%mixer, 1, m%np, st%d%nspin)

    call pop_sub() 
  end subroutine rsdfpt_init


  subroutine rsdfpt_e_init(m, st, h, p)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(hamiltonian_type), intent(in) :: h
    type(rsdfpt_type),     intent(inout) :: p
  
    integer :: i
    
    call push_sub('rsdfpt_e_init')
    
    call loct_parse_float("ConvAbsDens", CNST(1e-5), p%conv_abs_dens_e)
    call loct_parse_int("MaximumIter", 200, p%max_iter_e)

    call rsdfpt_init(m, st, h, p)

    do i = 1, 9
      p%rspolt(i) = M_ZERO 
    end do

    call pop_sub() 
  end subroutine rsdfpt_e_init

    
  subroutine rsdfpt_e_end(p)   
    type(rsdfpt_type),      intent(inout) :: p 

    call push_sub('rsdfpt_e_end')

    call rsdfpt_end(p)

    call pop_sub()
  end subroutine rsdfpt_e_end


  subroutine rsdfpt_end(p)   
    type(rsdfpt_type),      intent(inout) :: p 

    call push_sub('rsdfpt_end')

    deallocate(p%rsrho, p%rsrhoin, p%rsrhonew)
    deallocate(p%rsdpsi)
    deallocate(p%rsVhar)
    deallocate(p%rsVxc)
    deallocate(p%Y0,p%Y)

    call mix_end(p%mixer)

    call pop_sub()
  end subroutine rsdfpt_end


  subroutine build_fxc_kernel(m, st, p)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(rsdfpt_type),     intent(inout) :: p
    
    FLOAT :: fx(st%d%nspin, st%d%nspin)
    FLOAT :: fc(st%d%nspin, st%d%nspin)        
    integer :: i, j, k

    call push_sub('build_fxc_kernel')

    ! warning: irel is not used in any way
    ! warning: spin polarized pw92 is not implemented  
    ! warning: we assume that st contains ground state density
    do i = 1, m%np 
      call fxc_lda_exchange(st%d%nspin, 0, st%rho(i,:), fx(:,:))
      call fxc_pw92_correlation(st%d%nspin, 0, st%rho(i,:), fc(:,:))
      do j = 1, st%d%nspin
        do k = 1, st%d%nspin
          p%rsVxc(i,k,j) = fx(k,j)+fc(k,j)
        end do
      end do
    end do

    call pop_sub()
  end subroutine build_fxc_kernel


  subroutine build_response_dns(m, st, p)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(rsdfpt_type),     intent(inout) :: p
    
    integer :: i, nik, nst
    
    call push_sub('build_response_dns')
    
    do i = 1, m%np
      p%rsrho(1,i,:) = M_ZERO
      do nik = 1, st%d%nspin
        do nst = 1, st%nst
          if(st%occ(nst, nik) > M_ZERO) then 
            p%rsrho(1,i,nik) = p%rsrho(1,i,nik) + &
               M_TWO*st%occ(nst, nik)*st%dpsi(i,1,nst,nik)*p%rsdpsi(i,1,nst,nik)
          endif
        end do
      end do
    end do

    call pop_sub()
  end subroutine build_response_dns

  ! this computes fxc for LDA exchange, 
  ! ULDAC = -(\frac{1}{9 M_PI})^{1/3} 
  ! PLDAC = -(\frac{2}{9 M_PI})^{1/3} 
  ! EXPM23 = -\frac{2}{3}
  ! in the spin-polarized case:
  ! fxc(1,1) = \partial V_x^{\alpha}/\partial \rho_{\alpha}
  ! fxc(1,2) = \partial V_x^{\alpha}/\partial \rho_{\beta} 
  ! fxc(2,1) = \partial V_x^{\beta}/\partial \rho_{\alpha} 
  ! fxc(2,2) = \partial V_x^{\beta}/\partial \rho_{\beta} 
  subroutine fxc_lda_exchange(nsp, irel, ds, fxc)
    integer, intent(in)  :: nsp, irel
    FLOAT,   intent(in)  :: ds(:)    ! ds(nsp)
    FLOAT,   intent(out) :: fxc(:,:) ! fxc(nsp, nsp)

    FLOAT,   parameter   :: MINDEN = CNST(1e-15),          &
       ULDAC  = CNST(-0.328248341),   &
       EXPM23 = CNST(-0.666666666),   &   
       PLDAC  = CNST(-0.413566994)
    FLOAT                :: dns1, dns2

    call push_sub('fxc_lda_exchange')

    select case(nsp) 
    case(1)
      dns1 = max(M_ZERO, ds(1))
      if(dns1 < MINDEN) then 
        fxc(1,1) = M_ZERO
      else
        fxc(1,1) =ULDAC*(dns1**EXPM23)
      end if
    case(2) 
      dns1 = max(M_ZERO, ds(1))
      dns2 = max(M_ZERO, ds(2))
      if(dns1 < MINDEN) then 
        fxc(1,1) = M_ZERO
        fxc(1,2) = M_ZERO
      else 
        fxc(1, 1) = PLDAC*(dns1**EXPM23)
        fxc(1, 2) = M_ZERO
      endif
      if(dns2 < MINDEN) then 
        fxc(2, 1) = M_ZERO
        fxc(2, 2) = M_ZERO 
      else 
        fxc(2, 1) = M_ZERO
        fxc(2, 2) = PLDAC*(dns2**EXPM23) 
      endif
    end select

    call pop_sub()
  end subroutine fxc_lda_exchange

  !  computes deriv. of cor.pot for Perdew-Wang lda correlation phys.rev B 45, 13244 (1992)
  !  EXP13 = \frac{1}{3}
  !  PARS - array of parameters        1  2   3   4,  5,  6,  7
  !  PARS(:,1) - spin spin-unpolarized P, A, a1, b1, b2, b3, b4
  subroutine fxc_pw92_correlation(nsp, irel, ds, fxc)
    integer, intent(in)  :: nsp, irel
    FLOAT,   intent(in)  :: ds(nsp) 
    FLOAT,   intent(out) :: fxc(nsp,nsp)

    FLOAT,   parameter   :: MINDEN = CNST(1e-15),      &
       EXP13 = CNST(0.333333333), &
       PS(7,1) =   reshape( (/ CNST(1.0),     CNST(0.031091), &
       CNST(0.21370), CNST(7.5957),   &
       CNST(3.5876),  CNST(1.6382),   &
       CNST(0.49294) /), (/7, 1/) )

    FLOAT                :: dns1, dns2, rs
    FLOAT                :: q0, q1, dq0, dq1, d2q1, QQ
    FLOAT                :: dcor, d2cor
    
    call push_sub('fxc_pw92_correlation')

    select case(nsp) 
    case(1)
      dns1 = max(M_ZERO, ds(1))
      if(dns1 < MINDEN) then 
        fxc(1,1) = M_ZERO
      else 
        rs   = (0.75/(M_PI*dns1))**EXP13
        q0   = -2.0*PS(2,1)*(1.0+PS(3,1)*rs)
        dq0  = -2.0*PS(2,1)*PS(3,1)
        q1   = 2.0*PS(2,1)*(PS(4,1)*rs**0.5+PS(5,1)*rs     &
           + PS(6,1)*rs**1.5+PS(7,1)*rs**(PS(1,1)+1.0))
        dq1  = PS(2,1)*(PS(4,1)*rs**(-0.5)+2.0*PS(5,1)     &
           + 3.0*PS(6,1)*rs**0.5                         &
           + 2.0*(PS(1,1)+1.0)*PS(7,1)*rs**PS(1,1))
        d2q1 = PS(2,1)*(-0.5*PS(4,1)*rs**(-1.5)            &
           + 1.5*PS(6,1)*rs**(-0.5)                      &
           + 2.0*PS(1,1)*(PS(1,1)+1.0)*PS(7,1)*rs**(PS(1,1)-1))   
        QQ   = q1*q1+q1  
        dcor = dq0*log(1.0 + 1.0/q1) - q0*dq1/QQ
        d2cor= -(2.0*dq0*dq1+q0*d2q1)/QQ                  &
           + q0*dq1*dq1*(2.0*q1+1.0)/(QQ*QQ)
        fxc(1,1) = -2.0*rs*dcor/(9.0*dns1)                &
           + rs*rs*d2cor/(9.0*dns1) 
      end if
    end select

    call pop_sub()
  end subroutine fxc_pw92_correlation

  subroutine get_response_e(m, st, h, p, alpha, f_der)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(hamiltonian_type), intent(inout) :: h
    type(rsdfpt_type),      intent(inout) :: p
    integer,                intent(in)    :: alpha
    type(f_der_type),       intent(inout) :: f_der
    
    integer :: iter, nik, nik2, nst, np 
    FLOAT, allocatable :: tmp(:)
    logical :: finish
    
    call push_sub('get_response_e')

    allocate(tmp(m%np))

    call init_response_e(m, st, p)
    call orth_response(m, st, p)
    call build_response_dns(m, st, p)

    do iter=1, p%max_iter_e
      p%rsrhoin(:,:,:) = p%rsrho(:,:,:)

      do np = 1, m%np
        tmp(np) = sum(p%rsrho(1,np,:))
      end do

      if(.not.h%ip_app) call poisson_solve(m, f_der, p%rsVhar, tmp) 

      do nik = 1, st%d%nspin
        do nst = 1, st%nst 
          if (st%occ(nst, nik) > M_ZERO) then  
            p%Y(:,1) = (p%rsVhar(:) + m%x(:,alpha))*tmp(:) 
            do nik2 = 1, st%d%nspin
              p%Y(:,1) = p%Y(:,1) + p%rsVxc(:,nik,nik2)*p%rsrho(1,:,nik2)  
            end do
            p%Y(:,:) = -p%Y(:,:)*p%rsdpsi(:,:, nst, nik)
            call orth_vector(m, st, p%Y(:,:), nik)
            call solve_AXeY(m, st, p, f_der, h, nst, nik, 500, CNST(1.0e-5))
          endif
        end do
      end do

      call orth_response(m, st, p)
      call build_response_dns(m, st, p)

      p%rsrhonew(:,:,:) = M_ZERO
      call mixing(p%mixer, iter, 1, m%np, st%d%nspin, &
         p%rsrhoin, p%rsrho, p%rsrhonew)  

      p%abs_dens_e = M_ZERO
      do nik = 1, st%d%nspin
        tmp(:) = (p%rsrhoin(1,:,nik)-p%rsrho(1,:,nik))**2
        p%abs_dens_e = p%abs_dens_e + dmf_integrate(m, tmp)
      end do
      p%abs_dens_e = sqrt(p%abs_dens_e) 
      finish = (p%abs_dens_e <= p%conv_abs_dens_e) 

      if(finish) then 
        write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, 'iterations'  
        exit
      else  
        p%rsrho(:,:,:) = p%rsrhonew(:,:,:)  
      end if
    end do
    deallocate(tmp)  
    
    call pop_sub()
  end subroutine get_response_e
  
  subroutine solve_AXeY(m, st, p, f_der, h, nst, nik, MAXITER, tol)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(rsdfpt_type),      intent(inout) :: p
    type(f_der_type),       intent(inout) :: f_der
    type(hamiltonian_type), intent(inout) :: h
    integer,                intent(in)    :: nst, nik
    integer,                intent(in)    :: MAXITER
    FLOAT,                  intent(in)    :: tol
    
    integer :: iter  
    FLOAT, allocatable :: z(:,:), g(:,:), pp(:,:), t(:,:)
    FLOAT  :: sm1, sm2, sm3, alpha, beta

    call push_sub('solve_AXeY')

    allocate(z(m%np,st%dim), g(m%np,st%dim), pp(m%np,st%dim), t(m%np,st%dim))

    pp(:,:) = p%rsdpsi(:,:,nst,nik)
    call X(Hpsi)(h, m, f_der, pp, z, nik)

    z(:,:) = -z(:,:) + st%eigenval(nst,nik)*pp(:,:) + p%Y(:,:)
    pp(:,:) = z(:,:)
    do iter = 1, MAXITER     
      call X(Hpsi)(h, m, f_der, pp, t, nik)
      t(:,:) = t(:,:) - st%eigenval(nst, nik)*pp(:,:)

      sm1 = X(states_dotp) (m, 1, z, z)
      sm2 = X(states_dotp) (m, 1, z, t)
      alpha = sm1/sm2

      p%rsdpsi(:,:,nst,nik) = p%rsdpsi(:,:,nst,nik) + alpha*pp(:,:)
      g(:,:) = z(:,:)
      g(:,:) = g(:,:) -alpha*t(:,:)
      sm3 = X(states_dotp) (m, 1, g, g) 
      beta = sm3/sm1

      if(sm3 <= tol) then 
        exit
      else 
        pp(:,:) = beta*pp(:,:) + g(:,:)
        z(:,:) = g(:,:) 
      endif
    end do

    if(sm3 > tol) then 
      message(1) = "Response using CG: Not converged!"
      write(message(2),*) "iter = ", iter, " sm3 = ", sm3
      call write_warning(2)
    end if
    deallocate(z, g, pp, t)

    call pop_sub()
  end subroutine solve_AXeY

  ! orthogonalizes response of \alpha KS orbital to all \alpha KS orbitals  
  subroutine orth_response(m, st, p)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(rsdfpt_type),      intent(inout) :: p
    
    integer :: nst, nik
    call push_sub('orth_response')
    do nik = 1, st%d%nspin
      do nst = 1, st%nst 
        if(st%occ(nst, nik) > M_ZERO) then  
          call orth_vector(m, st, p%rsdpsi(:,:,nst,nik), nik) 
        endif
      end do
    end do
    
    call pop_sub()
  end subroutine orth_response
  
  subroutine orth_vector(m, st, v, nsp)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    FLOAT,                  intent(inout) :: v(:,:)
    integer,                intent(in)    :: nsp     
    
    FLOAT                                :: scalp
    integer                              :: nst
    call push_sub('orth_vector')
    
    do nst = 1, st%nst
      if(st%occ(nst, nsp) > M_ZERO) then
        scalp = R_REAL(X(states_dotp)(m, 1, v, st%X(psi)(:,:,nst,nsp)))
        v(:,:) = v(:,:) - scalp* st%X(psi)(:,:,nst,nsp)
      end if
    end do

    call pop_sub()
  end subroutine orth_vector

  subroutine init_response_e(m, st, p)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(rsdfpt_type),      intent(inout) :: p

    integer :: nik, nst, np
    FLOAT :: rd

    call push_sub('init_response_e')

    do nik = 1, st%d%nspin
      do nst = 1, st%nst
        if (st%occ(nst, nik) > M_ZERO) then
          do np = 1, m%np 
            call mesh_r(m, np, rd)
            p%rsdpsi(np,1,nst,nik) = st%dpsi(np,1,nst,nik)*rd*exp(-rd)
          end do
        endif
      end do
    end do
    p%rsVhar(:) = M_ZERO

    call pop_sub()
  end subroutine init_response_e

  subroutine pol_tensor(m, st, h, p, f_der)
    type(mesh_type),        intent(in) :: m
    type(states_type),      intent(in) :: st
    type(hamiltonian_type), intent(inout) :: h
    type(rsdfpt_type),      intent(inout) :: p
    type(f_der_type),       intent(inout) :: f_der

    integer :: i, j
    FLOAT :: rhov
    call push_sub('pol_tensor')

    do i = 1, conf%dim
      call get_response_e(m, st, h, p, i, f_der)

      do j = 1, m%np 
        rhov = sum(p%rsrho(1,j,:))*m%vol_pp(j)
        p%rspolt((i-1)*3+1) = p%rspolt((i-1)*3+1) - m%x(j,1)*rhov
        p%rspolt((i-1)*3+2) = p%rspolt((i-1)*3+2) - m%x(j,2)*rhov
        p%rspolt((i-1)*3+3) = p%rspolt((i-1)*3+3) - m%x(j,3)*rhov 
      end do
    end do

    p%msp = (p%rspolt(1)+p%rspolt(5)+p%rspolt(9))/M_THREE

    write(message(1), '(a,f12.6,a)') 'Info: Mean static polarizability', &
       p%msp, ' [a.u.]' 
    call write_info(1)

    call pop_sub()
  end subroutine pol_tensor

end module rsdfpt
