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
subroutine X(lr_orth_vector) (m, st, v, ik)
  type(mesh_t),        intent(in)    :: m
  type(states_t),      intent(in)    :: st
  R_TYPE,              intent(inout) :: v(:,:)
  integer,             intent(in)    :: ik

  call push_sub('linear_response_inc.Xlr_orth_vector')

  call X(states_gram_schmidt)(m, st%nst, st%d%dim, st%X(psi)(:, :, :, ik), v(:, :))

  call pop_sub()
end subroutine X(lr_orth_vector)


! --------------------------------------------------------------------
subroutine X(lr_build_dl_rho) (m, st, lr, nsigma)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  type(lr_t),     intent(inout) :: lr(:)
  integer,        intent(in)    :: nsigma

  integer :: i, ist, ik, ik2, sp
  CMPLX   :: c
  R_TYPE  :: d

  call push_sub('linear_response_inc.Xlr_build_dl_rho')

  ! initialize density
  do i = 1, nsigma
    lr(i)%X(dl_rho)(:, :) = M_ZERO
  end do

  sp = 1
  if(st%d%ispin == SPIN_POLARIZED) sp = 2

  ! calculate density
  do ik = 1, st%d%nik, sp
    do ist  = st%st_start, st%st_end
      do i = 1, m%np

        do ik2 = ik, ik+sp-1 ! this loop takes care of the SPIN_POLARIZED case
          d = st%d%kweights(ik2)*st%occ(ist, ik2)

          if(nsigma == 1) then  ! either omega purely real or purely imaginary
            d = d * st%X(psi)(i, 1, ist, ik2)*R_CONJ(lr(1)%X(dl_psi)(i, 1, ist, ik2))
            lr(1)%X(dl_rho)(i, 1) = lr(1)%X(dl_rho)(i, 1) + d + R_CONJ(d)
          else
            c = d*(                                                         &
              R_CONJ(st%X(psi)(i, 1, ist, ik2))*lr(1)%X(dl_psi)(i, 1, ist, ik2) + &
              st%X(psi)(i, 1, ist, ik2)*R_CONJ(lr(2)%X(dl_psi)(i, 1, ist, ik2)))
            lr(1)%X(dl_rho)(i, 1) = lr(1)%X(dl_rho)(i, 1) + c
            lr(2)%X(dl_rho)(i, 1) = lr(2)%X(dl_rho)(i, 1) + R_CONJ(c)
          end if
        end do

        if(st%d%ispin == SPINORS) then
          message(1) = "Not yet implemented - please fix me"
          call write_fatal(1)
        end if

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
      if(st%occ(ist, ik) > lr_min_occ) call X(lr_orth_vector) (m, st, lr%X(dl_psi)(:,:, ist, ik), ik)
    end do
  end do
  
  call pop_sub()
end subroutine X(lr_orth_response)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
