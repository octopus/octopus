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
!!
!! $Id$


! ---------------------------------------------------------
! Orthogonalizes v against all the occupied states
! ---------------------------------------------------------
subroutine X(lr_orth_vector) (m, st, v, ik, tol)
  type(mesh_t),        intent(in)    :: m
  type(states_t),      intent(in)    :: st
  R_TYPE,              intent(inout) :: v(:,:)
  integer,             intent(in)    :: ik
  FLOAT,   optional,   intent(in)    :: tol

  R_TYPE  :: scalp
  integer :: ist

  call push_sub('linear_response_inc.Xlr_orth_vector')

  do ist = 1, st%nst
    if(st%occ(ist, ik) > lr_min_occ) then
      scalp = X(states_dotp)(m, st%d%dim, st%X(psi)(:,:, ist, ik), v)
      if(present(tol)) then 
        if( ABS(scalp) < tol ) cycle
      end if
      v(1:m%np, 1:st%d%dim) = v(1:m%np, 1:st%d%dim) - scalp*st%X(psi)(1:m%np, 1:st%d%dim, ist, ik)
    end if
  end do

  call pop_sub()
end subroutine X(lr_orth_vector)


! -------------------------------------------------------------------
! calculates
!    lr%dl_rho += sum_{i occ} psi_i^0 (r) * psi_i^1*(r)   <=   type=1
!    lr%dl_rho += sum_{i occ} psi_i^0*(r) * psi_i^1 (r)   <=   type=2
!    type 3 => type 1 + type 2
! --------------------------------------------------------------------
subroutine X(lr_build_dl_rho) (m, st, lr, type)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr
  integer,        intent(in)    :: type

  integer :: i, p, ik, sp
  CMPLX   :: c
  R_TYPE  :: d(4)

  call push_sub('linear_response_inc.Xlr_build_dl_rho')

  sp = 1
  if(st%d%ispin == SPIN_POLARIZED) sp = 2

  do ik = 1, st%d%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, m%np
        d(1) = st%d%kweights(ik)*st%occ(p, ik) * &
          st%X(psi)(i, 1, p, ik)*R_CONJ(lr%X(dl_psi)(i, 1, p, ik))

        select case(st%d%ispin)
        case(SPIN_POLARIZED)
          d(2) = st%d%kweights(ik+1)*st%occ(p, ik+1) * &
            st%X(psi)(i, 1, p, ik+1)*R_CONJ(lr%X(dl_psi)(i, 1, p, ik+1))

        case(SPINORS)
          lr%X(dl_rho)(i, 2) = lr%X(dl_rho)(i, 2) + st%d%kweights(ik)*st%occ(p, ik) * &
            st%X(psi)(i, 2, p, ik)*R_CONJ(lr%X(dl_psi)(i, 2, p, ik))

          c = st%X(psi)(i, 1, p, ik) * R_CONJ(lr%X(dl_psi)(i, 2, p, ik))

          d(3) = st%d%kweights(ik)*st%occ(p, ik) * R_REAL(c)
          d(4) = st%d%kweights(ik)*st%occ(p, ik) * R_AIMAG(c)
        end select

        if(type == 2) d(1:st%d%nspin) = R_CONJ(d(1:st%d%nspin))
        if(type == 3) d(1:st%d%nspin) = R_CONJ(d(1:st%d%nspin)) + d(1:st%d%nspin)
        lr%X(dl_rho)(i, 1:st%d%nspin) = lr%X(dl_rho)(i, 1:st%d%nspin) + d(1:st%d%nspin)
      end do
    end do
  end do

  call pop_sub()
end subroutine X(lr_build_dl_rho)


! ---------------------------------------------------------
! orthogonalizes response of \alpha KS orbital to all occupied
! \alpha KS orbitals
! ---------------------------------------------------------
subroutine X(lr_orth_response)(m, st, lr)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr
  
  integer :: ist, ik
  call push_sub('linear_response_inc.Xlr_orth_response')
  
  do ik = 1, st%d%nspin
    do ist = 1, st%nst
      if(st%occ(ist, ik) > lr_min_occ) then
        call X(lr_orth_vector) (m, st, lr%X(dl_psi)(:,:, ist, ik), ik)
      end if
    end do
  end do
  
  call pop_sub()
end subroutine X(lr_orth_response)

