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
  use liboct
  integer, intent(in) :: func
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: pot(m%np, st%nspin), energy
  
  real(r8) :: e, e_x, dvol, den(3), dpol, dtot, vpol, glob(4, 3)!, loc(st%spin_channels)
  real(r8), parameter :: tiny = 1.0e-12_r8

  real(r8), allocatable :: d(:, :),    &
                           gd(:,:,:), &
                           lpot(:, :)
  real(r8), allocatable :: rhoplus(:), rhominus(:), mag(:) 
  real(r8), allocatable :: grhoplus(:, :), grhominus(:, :), gmag(:, :)

  real(r8) :: locald(st%spin_channels), localgd(3, st%spin_channels), &
              localdedd(st%spin_channels), localdedgd(3, st%spin_channels), &
              localdedd_x(st%spin_channels), localdedgd_x(3, st%spin_channels)
  !complex(r8), allocatable :: u(:, :, :)

  integer :: i, j, is, in, ic, ind(3), k

  call push_sub('xc_gga')

  allocate(d     (     m%np, st%nspin), &
           lpot  (    0:m%np, st%spin_channels))
  allocate(rhoplus(m%np), rhominus(m%np))
  allocate(grhoplus(3, m%np), grhominus(3, m%np))

  ! Store in local variables d the density matrix
  ! (in the global reference system).
  call dcopy(m%np*st%nspin, st%rho, 1, d, 1)

  ! If the pseudo has non-local core corrections, add the core charge
  ! (to the diagonal of the density matrix)
  if(nlcc) then
    do is = 1, st%spin_channels
       call daxpy(m%np, M_ONE/st%spin_channels, st%rho_core(1), 1, d(1, is), 1)
    enddo
  end if

  grhoplus = M_ZERO; grhominus = M_ZERO
  do i = 1, m%np
    select case(st%ispin)
    case(UNPOLARIZED)
      rhoplus(i) = max(d(i, 1), M_ZERO)
    case(SPIN_POLARIZED)
      rhoplus(i)  = max(d(i, 1), M_ZERO)
      rhominus(i) = max(d(i, 2), M_ZERO)
    case(SPINORS)
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      rhoplus(i)  = max(M_HALF*(dtot+dpol), M_ZERO)
      rhominus(i) = max(M_HALF*(dtot-dpol), M_ZERO)
    end select
  enddo

  call dmesh_derivatives(m, rhoplus(:),  grad=grhoplus(:, :))
  if(st%ispin > UNPOLARIZED) call dmesh_derivatives(m, rhominus(:), grad=grhominus(:, :))
  energy = M_ZERO
  lpot = M_ZERO

  space_loop: do i = 1, m%np

    locald(1) = rhoplus(i)
    localgd(1:3, 1) = grhoplus(1:3, i)
    if(st%ispin > UNPOLARIZED) then
      locald(2) = rhominus(i)
      localgd(1:3, 2) = grhominus(1:3, i)
    endif

    ! Calculate the potential/gradient density in local reference frame.
    select case(func)      
    case(X_FUNC_GGA_PBE)
      call pbex(0, st%spin_channels, locald, localgd, e, localdedd, localdedgd)
    case(X_FUNC_GGA_PBER)
      call pbex(1, st%spin_channels, locald, localgd, e, localdedd, localdedgd)
    case(X_FUNC_GGA_LB94)
      call xc_x_lb94(st%spin_channels, locald, localgd, e, &
                     localdedd, localdedgd) 
    case(C_FUNC_GGA_PBE)
      call pbec(st%spin_channels, locald, localgd, e, localdedd, localdedgd)
    case(C_FUNC_GGA_PBEX)
      call pbex(1, st%spin_channels, locald, localgd, e_x, localdedd_x, localdedgd_x)
      call pbec(st%spin_channels, locald, localgd, e, localdedd, localdedgd)
      e = e + e_x; localdedd = localdedd + localdedd_x; localdedgd = localdedgd + localdedgd_x
    end select

    energy = energy + sum(d(i, :)) * e * m%vol_pp

    lpot(i, :) = lpot(i, :) + localdedd(:)
    do in = -m%d%norder , m%d%norder
       if(m%lxyz(1,i)+in<m%nr(1,1).or.m%lxyz(1,i)+in>m%nr(2,1)) cycle
       ind(1) = m%Lxyz_inv(m%Lxyz(1,i)+in,m%Lxyz(2,i),m%Lxyz(3,i))
       if(m%lxyz(2,i)+in<m%nr(1,2).or.m%lxyz(2,i)+in>m%nr(2,2)) cycle
       ind(2) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i)+in,m%Lxyz(3,i))
       if(m%lxyz(3,i)+in<m%nr(1,3).or.m%lxyz(3,i)+in>m%nr(2,3)) cycle
       ind(3) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i),m%Lxyz(3,i)+in)
       do ic = 1, 3
          if(ind(ic) > 0) then
            !loc(:) = localdedgd(ic, :)
            lpot(ind(ic), :) = lpot(ind(ic), :) + localdedgd(ic, :) * m%d%dgidfj(in)/ m%h(ic)
          end if
       end do
    end do

  end do space_loop

  ! And now we rotate back (do not need the rotation matrix for this).
  if(st%ispin == SPINORS) then
    do i = 1, m%np
       dtot = d(i, 1) + d(i, 2)
       dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
       vpol = (lpot(i, 1) - lpot(i, 2))*(d(i, 1) - d(i, 2)) / (dpol + tiny)
       pot(i, 1) = M_HALF*(lpot(i, 1) + lpot(i, 2) + vpol)
       pot(i, 2) = M_HALF*(lpot(i, 1) + lpot(i, 2) - vpol)
       pot(i, 3) = (lpot(i, 1) - lpot(i, 2))*d(i, 3) / (dpol + tiny)
       pot(i, 4) = (lpot(i, 1) - lpot(i, 2))*d(i, 4) / (dpol + tiny)
    enddo
  else
    pot = lpot
  endif

  ! If LB94, we have to calculate the energy 
  ! Levy-Perdew relation (PRA 32, 2010 (1985))
  if(func == X_FUNC_GGA_LB94) then
    energy = 0._r8
    do is = 1, st%nspin
      call dmesh_derivatives(m, pot(:, is), grad=gd(:, :, is))
      do i = 1, m%np
        energy = energy + d(i, is) * sum(m%Lxyz(:,i)*m%h(:)*gd(:, i, is))
      end do
    end do
    energy = - energy * m%vol_pp
  end if

  deallocate(d, lpot, rhoplus, rhominus, grhoplus, grhominus)
  call pop_sub()
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
