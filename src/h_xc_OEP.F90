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

!!! This file handles the evaluation of the OEP potential, in the KLI or full OEP
!!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)

!!! This file has to be outside the module xc, for it requires the Hpsi

subroutine X(h_xc_oep)(xcs, m, st, vxc, ex, ec)
  type(xc_type),     intent(in)    :: xcs
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(inout) :: st
  real(r8),          intent(inout) :: vxc(m%np, st%nspin), ex, ec
  
  R_TYPE,   allocatable :: lxc(:,:)
  real(r8) :: socc, sfact, e
  integer :: is, ixc, ifunc

  call push_sub('h_xc_oep')

  ! this routine is only prepared for finite systems, and ispin = 1, 2
  if(st%ispin > 2 .or. st%nik>st%ispin) then
    message(1) = "OEP only works for finite systems and collinear spin!"
    call write_fatal(1)
  end if

  call getSpinFactor(st%nspin, socc, sfact)
  allocate(lxc(m%np, st%nst))

  spin: do is = 1, min(st%nspin, 2)
    lxc     = M_ZERO
    e       = M_ZERO

    ! get lxc
    functl_loop: do ixc = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(.not.btest(xcs%functl, ixc)) cycle
      ifunc = ibset(0, ixc)
      if(.not.( &
           ifunc == X_FUNC_OEP_X    .or. &
           ifunc == X_FUNC_OEP_SIC  .or. &
           ifunc == C_FUNC_OEP_SIC)) cycle
      
      select case(ifunc)
      case(X_FUNC_OEP_X)
        call X(oep_x) (m, st, is, socc, sfact, lxc, e)
        
      case(X_FUNC_OEP_SIC)
        call X(oep_x_sic) (xcs, m, st, is, socc, sfact, lxc, e)
      case(C_FUNC_OEP_SIC)
        call X(oep_c_sic) (xcs, m, st, is, socc, sfact, lxc, e)
      end select
      
      if(ixc < N_X_FUNCTL) then
        ex = ex + e
      else
        ec = ec + e
      end if
    end do functl_loop
  
    ! solve the KLI equation
    call X(xc_KLI_solve) (m, st, xcs%oep_level, is, socc, lxc, vxc)

  end do spin

  deallocate(lxc)

  call pop_sub()
end subroutine X(h_xc_oep)
