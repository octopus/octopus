!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!!! This module contains routines necessary to the split operator
!!! methods defined in td_exp
module td_exp_split
  use global
  use lib_basic_alg
  use mesh
  use states
  use hamiltonian

  implicit none

contains

  !!! Calculates psi = exp{factor*T} psi
  !!! where T is the kinetic energy operator
  subroutine zexp_kinetic (m, h, psi, ik, cf, factor)
    type(mesh_type),        intent(IN) :: m
    type(hamiltonian_type), intent(IN) :: h
    CMPLX,                  intent(inout) :: psi(m%np, h%d%dim)
    integer,                intent(in) :: ik
    type(zcf),              intent(inout) :: cf
    CMPLX,                  intent(in) :: factor
    
    integer :: ix, iy, iz, k(3), idim
    FLOAT :: cutoff, temp(3), g2

    call push_sub('exp_kinetic')
    
    if(conf%periodic_dim>0) then
      message(1) = 'Internal error in exp_kinetic'
      call write_fatal(1)
    endif
    
    if(h%cutoff > M_ZERO) then
      cutoff = h%cutoff
    else
      cutoff = CNST(1e10)
    endif

    temp = M_ZERO
    temp(1:conf%dim) = (M_TWO*M_Pi)/(cf%n(1:conf%dim)*m%h(1:conf%dim))

    call zcf_alloc_RS(cf)
    call zcf_alloc_FS(cf)

    do idim = 1, h%d%dim
      call zmf2cf(m, psi(:, idim), cf)
      call zcf_RS2FS(cf)

      do iz = 1, cf%n(3)
        k(3) = pad_feq(iz, cf%n(3), .true.)
        do iy = 1, cf%n(2)
          k(2) = pad_feq(iy, cf%n(2), .true.)
          do ix = 1, cf%n(1)
            k(1) = pad_feq(ix, cf%n(1), .true.)

            g2 = min(cutoff, sum((temp(1:conf%dim)*k(1:conf%dim))**2))
            cf%FS(ix, iy, iz) = exp(factor*g2/M_TWO)*cf%FS(ix, iy, iz)
          end do
        end do
      end do

      call zcf_FS2RS(cf)
      call zcf2mf(m, cf, psi(:, idim))
    enddo

    call zcf_free_RS(cf)
    call zcf_free_FS(cf)
    
    call pop_sub()
  end subroutine zexp_kinetic

  !!! Calculates psi = exp{factor*V_KS(t)} psi
  !!! where V_KS is the Kohn-Sham potential
  subroutine zexp_vlpsi (m, h, psi, ik, t, factor)
    type(mesh_type),        intent(IN)    :: m
    type(hamiltonian_type), intent(IN)    :: h
    CMPLX,                  intent(inout) :: psi(m%np, h%d%dim)
    integer,                intent(in)    :: ik
    FLOAT,                  intent(in)    :: t
    CMPLX,                  intent(in)    :: factor

    integer :: k
    FLOAT :: x(3), f(3)
    
    call push_sub('vlpsi')
    
    ! WARNING: spinors not yet supported.
    select case(h%d%ispin)
    case(UNPOLARIZED)
      psi(:, 1) = exp(factor*(h%vpsl(:)+h%vhxc(:, 1)))*psi(:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        psi(:, 1) = exp(factor*(h%vpsl(:)+h%vhxc(:, 1)))*psi(:, 1)
      else
        psi(:, 1) = exp(factor*(h%vpsl(:)+h%vhxc(:, 2)))*psi(:, 1)
      end if
    case(SPINORS)
      message(1) = 'Internal error in exp_vlpsi'
      call write_fatal(1)
    end select
    
    if(h%ep%no_lasers > 0) then
      if(h%gauge.ne.1) then  ! only length gauge is supported
        message(1) = "Only the length gauge is supported in exp_vlpsi"
        call write_fatal(1)
      end if

      call epot_laser_field(h%ep, t, f)
      do k = 1, m%np
        call mesh_xyz(m, k, x)
        psi(k,:) = exp(factor*sum(x*f)) * psi(k,:)
      end do
    end if
    
    call pop_sub()
  end subroutine zexp_vlpsi

  !!! calculates psi = exp{factor V_nlpp} psi
  !!! where V_nlpp is the non-local part of the pseudpotential
  subroutine zexp_vnlpsi (m, h, psi, ik, factor, order)
    type(mesh_type),        intent(IN) :: m
    type(hamiltonian_type), intent(IN) :: h
    CMPLX,                  intent(inout) :: psi(m%np, h%d%dim)
    integer,                intent(in) :: ik
    CMPLX,                  intent(in) :: factor
    logical,                intent(in) :: order

    integer :: idim, ikbc, jkbc, &
         ivnl_start, ivnl_end, step, kbc_start, kbc_end, ivnl
    CMPLX :: uvpsi, p2, ctemp
    CMPLX, allocatable :: lpsi(:), lHpsi(:), initzpsi(:, :)

    call push_sub('vnlpsi')

    allocate(initzpsi(m%np, 1:h%d%dim))
    initzpsi = psi

    dimension_loop: do idim = 1, h%d%dim

    if(order) then
      step = 1;  ivnl_start = 1; ivnl_end = h%ep%nvnl
    else
      step = -1; ivnl_start = h%ep%nvnl; ivnl_end  = 1
    endif

    do ivnl = ivnl_start, ivnl_end, step
       allocate(lpsi(h%ep%vnl(ivnl)%n), lhpsi(h%ep%vnl(ivnl)%n))
       lpsi(:) = initzpsi(h%ep%vnl(ivnl)%jxyz(:), idim)
       lhpsi(:) = M_z0

       if(order) then
          kbc_start = 1; kbc_end = h%ep%vnl(ivnl)%c
       else
          kbc_start = h%ep%vnl(ivnl)%c; kbc_end = 1
       end if

       do ikbc = kbc_start, kbc_end, step
          do jkbc = kbc_start, kbc_end, step
             p2 = sum(h%ep%vnl(ivnl)%uv(:, ikbc)*h%ep%vnl(ivnl)%uv(:, ikbc))*m%vol_pp
             ctemp = h%ep%vnl(ivnl)%uvu(ikbc, jkbc)*p2*factor
             uvpsi = sum(h%ep%vnl(ivnl)%uv(:, ikbc)*lpsi(:))*m%vol_pp*(exp(ctemp)-M_z1)/p2
             lhpsi(:) = lhpsi(:) + uvpsi*h%ep%vnl(ivnl)%uv(:, jkbc)
          end do
       end do

       deallocate(lpsi, lhpsi)
    enddo

    end do dimension_loop

    deallocate(initzpsi)

    call pop_sub()
  end subroutine zexp_vnlpsi

end module td_exp_split
