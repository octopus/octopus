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
  real(r8),          intent(out) :: vxc(m%np, st%nspin), ex, ec
  real(r8),          intent(in)  :: ip, qtot
  
  real(r8) :: e, dpol, dtot, vpol, vtot, vdif, r
  real(r8), allocatable :: d(:, :), gradd(:, :, :), div(:), dedgd(:, :, :)
  real(r8) :: dedd(st%spin_channels), dedgd1(3, st%spin_channels)!, &
  integer :: i, j, is, in, ic, ind(3), k, n, ixc

  call push_sub('xc_gga')

  allocate(d(m%np, st%nspin), gradd(3, m%np, st%nspin), dedgd(m%np, 3, st%spin_channels))
  
  ! Store in local variable d the density matrix
  ! (in the global reference system).
  ! If the pseudo has non-local core corrections, add the core charge
  ! (to the diagonal of the density matrix)
  call dcopy(m%np*st%nspin, st%rho, 1, d, 1)
  if(xcs%nlcc) then
    do is = 1, st%spin_channels
      call daxpy(m%np, M_ONE/st%spin_channels, st%rho_core(1), 1, d(1, is), 1)
    end do
  end if

  ! Make sure nothing is below zero, and rotate to the local reference frame in the case of spinors.
  do i = 1, m%np
    select case(st%ispin)
    case(UNPOLARIZED)
      d(i, 1) = max(d(i, 1), M_ZERO)
    case(SPIN_POLARIZED)
      d(i, 1)  = max(d(i, 1), M_ZERO)
      d(i, 2)  = max(d(i, 2), M_ZERO)
    case(SPINORS)
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      d(i, 1)  = max(M_HALF*(dtot+dpol), M_ZERO)
      d(i, 2) = max(M_HALF*(dtot-dpol), M_ZERO)
    end select
  enddo

  ! Calculate the gradient.
  call df_gradient(m, d(:, 1), gradd(:, :, 1))
  if(st%ispin > UNPOLARIZED) call df_gradient(m, d(:, 2), gradd(:, :, 2))

  vxc = M_ZERO
  dedgd = M_ZERO
  space_loop: do i = 1, m%np
        
    ! Calculate the potential in local reference frame.
    functl_loop: do ixc = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(.not.btest(xcs%functl, ixc)) cycle
        
      select case(ibset(0, ixc))
      case(X_FUNC_GGA_PBE)
        call pbex(0, st%spin_channels, d(i, :), gradd(:, i, :), e, dedd, dedgd1(:,:))
      case(X_FUNC_GGA_PBER)
        call pbex(1, st%spin_channels, d(i, :), gradd(:, i, :), e, dedd, dedgd1(:,:))
      case(X_FUNC_GGA_LB94)
        call mesh_r(m, i, r)
        call lb94(st%spin_channels, d(i, :), gradd(:, i, :), dedd, &
             r, ip, qtot, xcs%lb94_modified, xcs%lb94_beta, xcs%lb94_threshold)
        e = M_ZERO
        dedgd1(:,:) = M_ZERO
      case(C_FUNC_GGA_PBE)
        call pbec(st%spin_channels, d(i, :), gradd(:, i, :), e, dedd, dedgd1(:,:))
      end select
      
      if(ixc < N_X_FUNCTL) then
        ex = ex + sum(d(i,:)) * e * m%vol_pp
      else
        ec = ec + sum(d(i,:)) * e * m%vol_pp
      end if
      vxc(i, 1:st%spin_channels) = vxc(i, 1:st%spin_channels) + dedd(:)
      dedgd(i,:,:) = dedgd(i,:,:) + dedgd1(:,:)
      
    end do functl_loop

  end do space_loop

  ! We now substract the divergence of the functional derivative of fxc with respect to
  ! the gradient of the density. Still in local reference frame.
  allocate(div(m%np))
  do is = 1, st%spin_channels
    call df_divergence(m, dedgd(:,:,is), div(:))
    call daxpy(m%np, -M_ONE, div, 1, vxc(1, is), 1)
  end do
  deallocate(div)
      
  ! And now we rotate back
  if(st%ispin == SPINORS) then
    do i = 1, m%np
      dtot = d(i, 1) + d(i, 2)
      dpol = sqrt( (d(i, 1)-d(i, 2))**2 + M_FOUR*(d(i, 3)**2+d(i, 4)**2) )
      vtot = vxc(i,1)+vxc(i,2)
      vdif = vxc(i,1)-vxc(i,2)
      vpol = (vxc(i, 1) - vxc(i, 2))*(d(i, 1) - d(i, 2)) / (dpol + tiny)
      vxc(i, 1) = vxc(i, 1) + M_HALF*(vtot + vpol)
      vxc(i, 2) = vxc(i, 2) + M_HALF*(vtot - vpol)
      vxc(i, 3) = vxc(i, 3) + vdif*d(i, 3) / (dpol + tiny)
      vxc(i, 4) = vxc(i, 4) + vdif*d(i, 4) / (dpol + tiny)
    enddo
  endif

  
  ! If LB94, we have to calculate the energy 
  ! Levy-Perdew relation (PRA 32, 2010 (1985))
  ! (Note that gradd is used to store the gradient of vxc)
  if(iand(xcs%functl, X_FUNC_GGA_LB94).ne.0) then
    do is = 1, st%nspin
      call df_gradient(m, vxc(:, is), gradd(:, :, 1))
      do i = 1, m%np
        ex = ex - d(i, is) * sum(m%Lxyz(:,i)*m%h(:)*gradd(:, i, is)) * m%vol_pp
      end do
    end do
  end if

  deallocate(d, gradd, dedgd)
  call pop_sub()
end subroutine xc_gga

