!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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

module fxc
  use mesh
  use states

  implicit none

  private
  public :: build_fxc_kernel

contains

  subroutine build_fxc_kernel(dim, m, rho, fxc)
    type(states_dim_type), intent(in)    :: dim
    type(mesh_type),       intent(in)    :: m
    FLOAT,                 intent(in)    :: rho(:,:)    ! rho(m%np, dim%nspin)
    FLOAT,                 intent(inout) :: fxc(:,:,:)  ! fxc(m%np, dim%nspin, dim%nspin)
    
    FLOAT :: fx(dim%nspin, dim%nspin)
    FLOAT :: fc(dim%nspin, dim%nspin)        
    integer :: i

    call push_sub('build_fxc_kernel')

    ! warning: irel is not used in any way
    ! warning: spin polarized pw92 is not implemented  
    ! warning: we assume that st contains ground state density
    do i = 1, m%np 
      call fxc_lda_exchange(dim%nspin, 0, rho(i,:), fx(:,:))
      call fxc_pw92_correlation(dim%nspin, 0, rho(i,:), fc(:,:))
      fxc(i,:,:) = fx(:,:) + fc(:,:)
    end do
    
    call pop_sub()
  end subroutine build_fxc_kernel

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

end module fxc
