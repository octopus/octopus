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

subroutine xc_gga(xcs, m, st, vxc, ex, ec, ip, qtot)
  type(xc_type),     intent(in)  :: xcs
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(in)  :: st
  FLOAT,          intent(out) :: vxc(m%np, st%nspin), ex, ec
  FLOAT,          intent(in)  :: ip, qtot
  
  FLOAT :: e, dpol, dtot, vpol, r

  FLOAT, allocatable :: d(:, :), lpot(:, :)
  FLOAT, allocatable :: rhoplus(:), rhominus(:)
  FLOAT, allocatable :: grhoplus(:, :), grhominus(:, :)
  FLOAT, allocatable :: vlocaldedgd(:,:,:), vlocaldedgd1(:,:)

  FLOAT :: x(3), locald(st%spin_channels), localgd(3, st%spin_channels), &
       localdedd(st%spin_channels), &
       localdedd_x(st%spin_channels), localdedgd_x(3, st%spin_channels)

  integer :: i, j, is, in, ic, ind(3), k, n
  integer :: ixc

  call push_sub('xc_gga')

  allocate(d(m%np, st%nspin), lpot(m%np, st%spin_channels))
  allocate(rhoplus(m%np), rhominus(m%np))
  allocate(grhoplus(3, m%np), grhominus(3, m%np))
  allocate(vlocaldedgd(m%np, 3, st%spin_channels), vlocaldedgd1(3, st%spin_channels))

  ! Store in local variables d the density matrix
  ! (in the global reference system).
  call dcopy(m%np*st%nspin, st%rho(1, 1), 1, d(1, 1), 1)

  ! If the pseudo has non-local core corrections, add the core charge
  ! (to the diagonal of the density matrix)
  if(xcs%nlcc) then
    do is = 1, st%spin_channels
      call daxpy(m%np, M_ONE/st%spin_channels, st%rho_core(1), 1, d(1, is), 1)
    end do
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

  call df_gradient(m, rhoplus, grhoplus)
  if(st%ispin > UNPOLARIZED) call df_gradient(m, rhominus, grhominus)

  lpot        = M_ZERO
  vlocaldedgd = M_ZERO
  space_loop: do i = 1, m%np
        
    locald(1) = rhoplus(i)
    localgd(1:3, 1) = grhoplus(1:3, i)
    if(st%ispin > UNPOLARIZED) then
      locald(2) = rhominus(i)
      localgd(1:3, 2) = grhominus(1:3, i)
    endif

    ! Calculate the potential/gradient density in local reference frame.
    functl_loop: do ixc = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(.not.btest(xcs%functl, ixc)) cycle
        
      select case(ibset(0, ixc))
      case(X_FUNC_GGA_PBE)
        call pbex(0, st%spin_channels, locald, localgd, e, localdedd, vlocaldedgd1(:,:))
      case(X_FUNC_GGA_PBER)
        call pbex(1, st%spin_channels, locald, localgd, e, localdedd, vlocaldedgd1(:,:))
      case(X_FUNC_GGA_LB94)
        call mesh_r(m, i, r)
        call lb94(st%spin_channels, locald, localgd, localdedd, &
             r, ip, qtot, xcs%lb94_modified, xcs%lb94_beta, xcs%lb94_threshold)
        e = M_ZERO
        vlocaldedgd1(:,:) = M_ZERO
      case(C_FUNC_GGA_PBE)
        call pbec(st%spin_channels, locald, localgd, e, localdedd, vlocaldedgd1(:,:))
      end select
      
      if(ixc < N_X_FUNCTL) then
        ex = ex + sum(d(i,:)) * e * m%vol_pp
      else
        ec = ec + sum(d(i,:)) * e * m%vol_pp
      end if
      lpot(i, :) = lpot(i, :) + localdedd(:)
      vlocaldedgd(i,:,:) = vlocaldedgd(i,:,:) + vlocaldedgd1(:,:)
      
    end do functl_loop

  end do space_loop

  ! We now add substract the divergence of the functional derivative of fxc with respect to
  ! the gradient of the density.
  do is = 1, st%spin_channels
    call df_divergence(m, vlocaldedgd(:,:,is), rhoplus(:))
    call daxpy(m%np, -M_ONE, rhoplus(1), 1, lpot(1, is), 1)
  end do
      
  ! And now we rotate back (do not need the rotation matrix for this).
  if(st%ispin == SPINORS) then
    do i = 1, m%np
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      vpol = (lpot(i, 1) - lpot(i, 2))*(d(i, 1) - d(i, 2)) / (dpol + tiny)
      vxc(i, 1) = vxc(i, 1) + M_HALF*(lpot(i, 1) + lpot(i, 2) + vpol)
      vxc(i, 2) = vxc(i, 2) + M_HALF*(lpot(i, 1) + lpot(i, 2) - vpol)
      vxc(i, 3) = vxc(i, 3) + (lpot(i, 1) - lpot(i, 2))*d(i, 3) / (dpol + tiny)
      vxc(i, 4) = vxc(i, 4) + (lpot(i, 1) - lpot(i, 2))*d(i, 4) / (dpol + tiny)
    enddo
  else
    vxc = vxc + lpot
  endif
  
  ! If LB94, we have to calculate the energy 
  ! Levy-Perdew relation (PRA 32, 2010 (1985))
  if(iand(xcs%functl, X_FUNC_GGA_LB94).ne.0) then
    do is = 1, st%nspin
      call df_gradient(m, vxc(:, is), grhoplus)
      do i = 1, m%np
        ex = ex - d(i, is) * sum(m%Lxyz(:,i)*m%h(:)*grhoplus(:, i)) * m%vol_pp
      end do
    end do
  end if

  deallocate(d, lpot, rhoplus, rhominus, grhoplus, grhominus, vlocaldedgd, vlocaldedgd1)
  call pop_sub()
end subroutine xc_gga

