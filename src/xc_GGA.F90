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
  
  real(r8) :: e, e_x, dvol, den(3), dpol, dtot, vpol, glob(4, 3), loc(2)
  real(r8), parameter :: tiny = 1.0e-12_r8

  real(r8), allocatable :: d(:, :),    &
                           gd(:,:,:), &
                           lpot(:, :)

  real(r8) :: locald(st%spin_channels), localgd(3, st%spin_channels), &
              localdedd(st%spin_channels), localdedgd(3, st%spin_channels), &
              localdedd_x(st%spin_channels), localdedgd_x(3, st%spin_channels)
  complex(r8), allocatable :: u(:, :, :)

  integer :: i, is, in, ic, ind(3), k

  sub_name = 'xc_gga'; call push_sub()

  allocate(d     (   0:m%np, st%nspin), &
           gd    (3,   m%np, st%nspin), &
           lpot  (     m%np, st%spin_channels), &
           u     (2, 2, m%np))

  ! Store in local variables d the density matrix
  ! (in the global reference system).
  do is = 1, st%nspin
    d(0, is) = M_ZERO
    d(1:m%np, is) = st%rho(1:m%np, is)
  end do

  ! If the pseudo has non-local core corrections, add the core charge
  ! (to the diagonal of the density matrix)
  if(nlcc) then
    do is = 1, st%spin_channels
       d(1:m%np, is) = d(1:m%np, is) + st%rho_core(1:m%np)/st%spin_channels
    enddo
  end if

  ! Store in local variable gd the denssity gradient.
  do is = 1, st%nspin
    call dmesh_derivatives(m, d(:, is), grad=gd(:,:, is))
  enddo

  ! Build the rotation matrix u: By doing u*g*transpose[u*], it takes the
  ! operator g to the reference frame where the density is diagonal (in
  ! particular, u*d*transpose[u*] is a diagonal matrix).
  if(st%ispin == SPINORS) then
    do i = 1, m%np
       call rot_matrix(u(:, :, i), d(i, :))
    enddo
  endif

  energy = M_ZERO
  lpot = M_ZERO
  space_loop: do i = 1, m%np

    ! Rotate density to the local reference frame.
    if(st%ispin == SPINORS) then
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      locald(1) = max(M_HALF*(dtot + dpol), M_ZERO)
      locald(2) = max(M_HALF*(dtot - dpol), M_ZERO)
    else
      do is = 1, st%spin_channels
         locald(is) = max(d(i, is), M_ZERO)
      enddo
    endif

    ! Now we rotate the gradient. For this purpose, we do need the rotation matrix.
    if(st%ispin == SPINORS) then
      do ic = 1, 3
         call to_local(gd(ic, i, :), localgd(ic, :), u(:, :, i))
      enddo
    else
      localgd(1:3, 1:st%spin_channels) = gd(1:3, i, 1:st%spin_channels)
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
    if(st%ispin==SPINORS) then
      call to_global(localdedgd(1, :), glob(:, 1), u(:, :, i))
      call to_global(localdedgd(2, :), glob(:, 2), u(:, :, i))
      call to_global(localdedgd(3, :), glob(:, 3), u(:, :, i))
    endif
    do in = -m%d%norder , m%d%norder
       ind(1) = m%Lxyz_inv(m%Lxyz(1,i)+in,m%Lxyz(2,i),m%Lxyz(3,i))
       ind(2) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i)+in,m%Lxyz(3,i))
       ind(3) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i),m%Lxyz(3,i)+in)
       do ic = 1, 3
          if(ind(ic) > 0) then
            if(st%ispin == SPINORS) then
               call to_local (glob(:, ic),             loc,  u(:, :, ind(ic)))
            else
               loc(:) = localdedgd(ic, :)
            endif
            do is = 1, st%spin_channels
               lpot(ind(ic), is) = lpot(ind(ic), is) + loc(is) * m%d%dgidfj(in)/ m%h(ic)
            enddo
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

  deallocate(d, gd, lpot, u)
  call pop_sub(); return
  contains


  subroutine rot_matrix(u, dens)
    complex(r8), intent(out) :: u(2, 2)
    real(r8), intent(in)     :: dens(4)

    real(r8)    ::phi, theta, absrho12

    absrho12 = sqrt(dens(3)**2+dens(4)**2)
    if(absrho12<tiny) then
      phi = M_ZERO
    else
      phi   = -atan(dens(4)/dens(3))
    endif
    if(absrho12<tiny) then
      theta = M_ZERO
    elseif(abs(real(dens(1)-dens(2)))<tiny) then
      theta = M_PI/2
    else
      theta =  atan(M_TWO*absrho12/(dens(1)-dens(2)))
    endif
    u(1, 1) =  cos(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(1, 2) =  sin(M_HALF*theta)*exp(-M_zI*M_HALF*phi)
    u(2, 1) = -sin(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(2, 2) =  cos(M_HALF*theta)*exp(-M_zI*M_HALF*phi)

  end subroutine rot_matrix

  ! local = u*global*transpose[u*]
  subroutine to_local(global, local, u)
    implicit none
    real(r8), intent(in)  :: global(4)
    real(r8), intent(out) :: local (2)
    complex(r8), intent(in) :: u(2,2)

    complex(r8) :: g12, g21

    g12 = global(3) + M_zI*global(4)
    g21 = conjg(g12)
    local(1) = (u(1,1)*global(1) + u(1,2)*g21   )*conjg(u(1,1)) + &
               (u(1,1)*g12    + u(1,2)*global(2))*conjg(u(1,2))
    local(2) = (u(2,1)*global(1) + u(2,2)*g21   )*conjg(u(2,1)) + &
               (u(2,1)*g12    + u(2,2)*global(2))*conjg(u(2,2))
  end subroutine to_local

  ! global = transpose[u*]*local*u
  subroutine to_global(local, global, u)
    implicit none
    real(r8), intent(out) :: global(4)
    real(r8), intent(in)  :: local (2)
    complex(r8), intent(in) :: u(2,2)

    complex(r8) :: g12

    global(1) = local(1)*abs(u(1,1))**2 + local(2)*abs(u(2,1))**2
    global(2) = local(1)*abs(u(1,2))**2 + local(2)*abs(u(2,2))**2
    g12 = local(1)*conjg(u(1,1))*u(1,2)
    global(3) = real(g12)
    global(4) = aimag(g12)
  end subroutine to_global

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
