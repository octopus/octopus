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

subroutine xc_gga(func, nlcc, m, st, pot, energy)
  integer, intent(in) :: func
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: pot(m%np, st%nspin), energy
  
  real(r8) :: dedd(st%nspin), dummy_dedd(st%nspin), &
              dedgd(3, st%nspin), dummy_dedgd(3, st%nspin), &
              e, dummy_e, dvol, den(3)
  real(r8), allocatable :: d(:,:), gd(:,:,:)
  integer :: i, is, in, ic, ind(3)

  sub_name = 'xc_gga'; call push_sub()

  allocate(d(0:m%np, st%nspin), gd(3, m%np, st%nspin))
  do is = 1, st%nspin
    d(0, is) = 0._r8
    d(1:m%np, is) = abs(st%rho(1:m%np, is))
    if(nlcc) then ! non-linear core corrections
      d(1:m%np, is) = d(1:m%np, is) + st%rho_core(1:m%np)/st%nspin
    end if
    call dmesh_derivatives(m, d(:, is), grad=gd(:,:, is))
  end do

  energy = M_ZERO
  do i = 1, m%np
    select case(func)      
    case(X_FUNC_GGA_PBE)
      call pbexc(0, st%nspin, d(i, :), gd(1:3, i, 1:st%nspin), e, dummy_e, &
                 dedd, dummy_dedd, dedgd, dummy_dedgd) 
    case(X_FUNC_GGA_PBER)
      call pbexc(1, st%nspin, d(i, :), gd(1:3, i, 1:st%nspin), e, dummy_e, &
                 dedd, dummy_dedd, dedgd, dummy_dedgd) 
    case(X_FUNC_GGA_LB94)
      call xc_x_lb94(st%nspin, d(i, :), gd(1:3,i,1:st%nspin), e, dedd, dedgd) 
    case(C_FUNC_GGA_PBE)
      call pbexc(0, st%nspin, d(i, :), gd(1:3, i, 1:st%nspin), dummy_e, e, &
                 dummy_dedd, dedd, dummy_dedgd, dedgd) 
    end select
    energy = energy + sum(d(i, :)) * e * m%vol_pp

    do is = 1, st%nspin
      pot(i, is) = pot(i, is) + dedd(is)

      do in = -m%d%norder , m%d%norder
        ind(1) = m%Lxyz_inv(m%Lxyz(1,i)+in,m%Lxyz(2,i),m%Lxyz(3,i))
        ind(2) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i)+in,m%Lxyz(3,i))
        ind(3) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i),m%Lxyz(3,i)+in)

#ifndef BOUNDARIES_ZERO_DERIVATIVE
        den = 0.0_r8
#endif
        if(ind(1) > 0)den(1) = d(ind(1), is)
        if(ind(2) > 0)den(2) = d(ind(2), is)
        if(ind(3) > 0)den(3) = d(ind(3), is)
        
        do ic = 1, 3
          if(ind(ic) > 0) then
            pot(ind(ic), is) = pot(ind(ic), is) + &
                 dedgd(ic, is) * m%d%dgidfj(in)/ m%h(ic)
          end if
        end do
      end do
    end do
  end do

  if(func == X_FUNC_GGA_LB94) then ! we have to calculate the energy
    ! Levy-Perdew relation (PRA 32, 2010 (1985))
    energy = 0._r8
    do is = 1, st%nspin
      call dmesh_derivatives(m, pot(:, is), grad=gd(:, :, is))
      do i = 1, m%np
        energy = energy + d(i, is) * sum(m%Lxyz(:,i)*m%h(:)*gd(:, i, is))
      end do
    end do
    energy = - energy * m%vol_pp
  end if

  call pop_sub(); return
end subroutine xc_gga

subroutine xc_x_lb94(nspin, dens, gdens, ex, dexdd, dexdgd)
  integer, intent(in)  :: nspin
  real(r8), intent(IN) :: dens(nspin), gdens(3, nspin)
  real(r8), intent(out):: ex, dexdd(nspin), dexdgd(3, nspin)
  
  ! Internal variables
  integer(i4) :: is 
  real(r8)    :: d(nspin), gd(3, nspin), gdm, x, f

! Lower bounds of density and its gradient to avoid divisions by zero
! plus some numerical constants
  real(r8), parameter :: &
      DENMIN = 1.E-20_r8,    GDMIN  = 1.E-20_r8,          &
      THRD   = 1._r8/3._r8,        &
      FTHRD  = 4._r8/3._r8, BETA   = 0.05_r8

! first we add the LDA potential
  call exchng(0, nspin, dens, ex, dexdd)

! Translate density and its gradient to new variables
  if (nspin .eq. 1) then
    d(1)     = dens(1) * M_HALF
    gd(:, 1) = gdens(:, 1) * M_HALF
  else
    d  = dens
    gd = gdens
  endif

  do is = 1, nspin
    gdm   = sqrt( gd(1,is)**2 + gd(2,is)**2 + gd(3,is)**2 )
    
    if(d(is) >= DENMIN .and. gdm >=GDMIN) then
      x = gdm / d(is)**FTHRD
      f = x**2/(1._r8 + 3._r8*BETA*x*oct_asinh(x))
      dexdd(is) = dexdd(is) - BETA * d(is)**THRD * f
    end if
  end do

 ! this is to be calculated afterwards
  ex = M_ZERO; dexdgd = M_ZERO

  return
end subroutine xc_x_lb94  
