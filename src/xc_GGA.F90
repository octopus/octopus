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
  
  real(r8) :: dummy_dedd(st%nspin), dummy_dedgd(3, st%nspin), e, dummy_e, dvol, den(3), &
              dpol, dtot, vpol, glob(4), loc(2)
  real(r8), parameter :: tiny = 1.0e-12_r8

  real(r8), allocatable :: d(:, :),    &
                           ld(:, :),   &
                           gd(:,:,:), &
                           lgd(:, :, :), &
                           dedd(:, :), &
                           ldedd(:, :), &
                           dedgd(:, :, :), &
                           ldedgd(:, :, :), &
                           lpot(:, :)

  integer :: i, is, in, ic, ind(3)

  sub_name = 'xc_gga'; call push_sub()

  allocate(d     (   0:m%np, st%nspin), &
           gd    (3,   m%np, st%nspin), &
           dedd  (     m%np, st%nspin), &
           dedgd (3,   m%np, st%nspin), &
           ld    (   0:m%np, st%spin_channels), &
           lgd   (3,   m%np, st%spin_channels), &
           ldedd (     m%np, st%spin_channels), &
           ldedgd(3,   m%np, st%spin_channels), &
           lpot  (     m%np, st%spin_channels))

  ! Store in local variables d and gd the density matrix, and its gradient
  ! (in the global reference system).
  do is = 1, st%nspin
    d(0, is) = M_ZERO
    d(1:m%np, is) = abs(st%rho(1:m%np, is))
    call dmesh_derivatives(m, d(:, is), grad=gd(:,:, is))
  end do

  ! If the pseudo has non-local core corrections, add the core charge
  ! (to the diagonal of the density matrix)
  if(nlcc) then
    do is = 1, st%spin_channels
       d(1:m%np, is) = d(1:m%np, is) + st%rho_core(1:m%np)/st%spin_channels
    enddo
  end if

  energy = M_ZERO
  lpot = M_ZERO
  space_loop: do i = 1, m%np

    ! Rotate density to the local reference frame.
    if(st%ispin == SPINORS) then
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      ld(i, 1) = M_HALF*(dtot + dpol)
      ld(i, 2) = M_HALF*(dtot - dpol)
    else
      ld(i, 1:st%spin_channels) = d(i, 1:st%spin_channels)
    endif

    ! Now we rotate the gradient. For this purpose, we do need the rotation matrix.
    if(st%ispin == SPINORS) then
      do ic = 1, 3
         call to_local(gd(ic, i, :), lgd(ic, i, :), d(i, :))
      enddo
    else
      lgd(1:3, i, 1:st%spin_channels) = gd(1:3, i, 1:st%spin_channels)
    endif

    ! Calculate the potential/gradient density in local reference frame.
    select case(func)      
    case(X_FUNC_GGA_PBE)
      call pbexc(0, st%spin_channels, ld(i, :), lgd(1:3, i, 1:st%nspin), e, dummy_e, &
                 ldedd(i, :), dummy_dedd, ldedgd(:, i, :), dummy_dedgd) 
    case(X_FUNC_GGA_PBER)
      call pbexc(1, st%spin_channels, ld(i, :), lgd(1:3, i, 1:st%nspin), e, dummy_e, &
                 ldedd(i, :), dummy_dedd, ldedgd(:, i, :), dummy_dedgd) 
    case(X_FUNC_GGA_LB94)
      call xc_x_lb94(st%spin_channels, ld(i, :), lgd(1:3,i,1:st%nspin), e, &
                     ldedd(i, :), ldedgd(:, i, :)) 
    case(C_FUNC_GGA_PBE)
      call pbexc(0, st%spin_channels, ld(i, :), lgd(1:3, i, 1:st%nspin), dummy_e, e, &
                 dummy_dedd, ldedd(i, :), dummy_dedgd, ldedgd(:, i, :))
    end select
    energy = energy + sum(d(i, :)) * e * m%vol_pp

      lpot(i, :) = lpot(i, :) + ldedd(i, :)
      do in = -m%d%norder , m%d%norder
        ind(1) = m%Lxyz_inv(m%Lxyz(1,i)+in,m%Lxyz(2,i),m%Lxyz(3,i))
        ind(2) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i)+in,m%Lxyz(3,i))
        ind(3) = m%Lxyz_inv(m%Lxyz(1,i),m%Lxyz(2,i),m%Lxyz(3,i)+in)
        do ic = 1, 3
          if(ind(ic) > 0) then
            call to_global(ldedgd(ic, i, :), glob, d(i, :))
            call to_local (glob,             loc,  d(ind(ic), :))
            do is = 1, st%spin_channels
               lpot(ind(ic), is) = lpot(ind(ic), is) + loc(is) * m%d%dgidfj(in)/ m%h(ic)
            enddo
          end if
        end do
      end do

  end do space_loop

  ! And now we rotate back.
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

  deallocate(d, gd, dedd, dedgd, ld, lgd, ldedd, ldedgd, lpot)
  call pop_sub(); return
  contains

  ! This subroutine projects a 2x2 hemitian matrix (density or potential in spin space
  ! to the reference frame where the density (entry "dens") is diagonal. Of course
  ! the input "global" need not be diagonal, but nondiagonal terms are ignored.
  ! A 2x2 hermitian matrix x is represented by a real 4-vector, xx(1) = x(1, 1),
  ! xx(2) = x(2, 2), xx(3) = Re(x(1,2)), xx(4) = Im(x(1,2)).
  subroutine to_local(global, local, dens)
    implicit none
    real(r8), intent(in)  :: global(4)
    real(r8), intent(out) :: local (2)
    real(r8), intent(in)  :: dens(4)

    real(r8), parameter :: tiny = 1.0e-12_r8
    real(r8) :: phi, theta
    complex(r8) :: rho(2, 2), u(2, 2), ut(2, 2), b(2, 2), aux(2, 2), g(2, 2)

    ! Builds the density matrix
    rho(1, 1) = dens(1)
    rho(2, 2) = dens(2)
    rho(1, 2) = dens(3) + M_zI*dens(4)
    rho(2, 1) = dens(3) - M_zI*dens(4)

    g(1, 1) = global(1)
    g(2, 2) = global(2)
    g(1, 2) = global(3) + M_zI*global(4)
    g(2, 1) = global(3) - M_zI*global(4)

    if(abs(rho(1, 2))<tiny) then
      phi = M_ZERO
    else
      phi   = -atan(aimag(rho(1, 2))/real(rho(1, 2)))
    endif
    if(abs(rho(1,2))<tiny) then
      theta = M_ZERO
    elseif(abs(real(dens(1)-dens(2)))<tiny) then
      theta = M_PI/2
    else
      theta =  atan(M_TWO*abs(rho(1, 2))/(real(rho(1,1))-real(rho(2,2))))
    endif

    u(1, 1) =  cos(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(1, 2) =  sin(M_HALF*theta)*exp(-M_zI*M_HALF*phi)
    u(2, 1) = -sin(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(2, 2) =  cos(M_HALF*theta)*exp(-M_zI*M_HALF*phi)

    aux = matmul(u, matmul(g, transpose(conjg(u))))

    local(1) = aux(1, 1)
    local(2) = aux(2, 2)
    
  end subroutine to_local

  subroutine to_global(local, global, dens)
    implicit none
    real(r8), intent(out) :: global(4)
    real(r8), intent(in)  :: local (2)
    real(r8), intent(in)  :: dens(4)

    real(r8), parameter :: tiny = 1.0e-12_r8
    real(r8) :: phi, theta
    complex(r8) :: rho(2, 2), u(2, 2), ut(2, 2), b(2, 2), aux(2, 2), g(2, 2)

    ! Builds the density matrix
    rho(1, 1) = dens(1)
    rho(2, 2) = dens(2)
    rho(1, 2) = dens(3) + M_zI*dens(4)
    rho(2, 1) = dens(3) - M_zI*dens(4)

    g(1, 1) = local(1)
    g(2, 2) = local(2)
    g(1, 2) = (M_ZERO, M_ZERO)
    g(2, 1) = (M_ZERO, M_ZERO)

    if(abs(rho(1, 2))<tiny) then
      phi = M_ZERO
    else
      phi   = -atan(aimag(rho(1, 2))/real(rho(1, 2)))
    endif
    if(abs(rho(1,2))<tiny) then
      theta = M_ZERO
    elseif(abs(real(dens(1)-dens(2)))<tiny) then
      theta = M_PI/2
    else
      theta =  atan(M_TWO*abs(rho(1, 2))/(real(rho(1,1))-real(rho(2,2))))
    endif

    u(1, 1) =  cos(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(1, 2) =  sin(M_HALF*theta)*exp(-M_zI*M_HALF*phi)
    u(2, 1) = -sin(M_HALF*theta)*exp( M_zI*M_HALF*phi)
    u(2, 2) =  cos(M_HALF*theta)*exp(-M_zI*M_HALF*phi)

    aux = matmul(transpose(conjg(u)), matmul(g, u))

    global(1) = aux(1, 1)
    global(2) = aux(2, 2)
    global(3) = real (aux(1, 2))
    global(4) = aimag(aux(1, 2))
    
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
