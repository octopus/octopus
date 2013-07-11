!! Copyright (C) 2013 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!


subroutine X(pnfft_forward)(pnfft, in, out)
  type(pnfft_t), intent(inout)  :: pnfft
  R_TYPE,           intent(in)  :: in(:,:,:)
  CMPLX,            intent(out) :: out(:,:,:)

  integer :: i1, i2, i3

  PUSH_SUB(X(pnfft_forward))

!   print *, mpi_world%rank, "---> pnfft%N_local       ", pnfft%N_local 
!   print *, mpi_world%rank, "---> pnfft%M             ", pnfft%M 
!   print *, mpi_world%rank, "---> size(in)            ", size(in,1), size(in,2), size(in,3) 
!   print *, mpi_world%rank, "---> size(pnfft%f_hat)   ", size(pnfft%f_hat,1), size(pnfft%f_hat,2), size(pnfft%f_hat, 3) 
!   print *, mpi_world%rank, "---> size(out)           ", size(out,1), size(out,2), size(out,3)  
!   print *, mpi_world%rank, "---> size(pnfft%f)       ", size(pnfft%f,1), size(pnfft%f,2), size(pnfft%f,3) 
  
!   do i1 = 1, pnfft%N_local(1)
!     do i2 = 1, pnfft%N_local(2)
!       do i3 = 1, pnfft%N_local(3)
!         pnfft%f_hat(i1,i3,i2) = in(i1,i3,i2)
!       end do
!     end do
!   end do

! #ifdef R_TCOMPLEX 
!   pnfft%f_hat(:,:,:) = cmplx(real(in(:,:,:)),aimag(in(:,:,:)), C_DOUBLE_COMPLEX)
! #else
!   pnfft%f_hat(:,:,:) = cmplx(in(:,:,:), C_DOUBLE_COMPLEX)
! #endif

  pnfft%f_hat(:,:,:) = in(:,:,:)

  call pnfft_trafo(pnfft%plan)

  out(:,:,:) = pnfft%f(:,:,:)


!   do i1 = 1, pnfft%M(1)
!     do i2 = 1, pnfft%M(2)
!       do i3 = 1, pnfft%M(3)
! !         out(i1,i2,i3) = pnfft%f_lin(pnfft_idx_3to1(pnfft,i1,i2,i3))
!         out(i3,i2,i1) = pnfft%f(i1,i2,i3)
!       end do
!     end do
!   end do

  POP_SUB(X(pnfft_forward))
    
end subroutine X(pnfft_forward)


! ---------------------------------------------------------
subroutine X(pnfft_backward)(pnfft, in, out)
  type(pnfft_t), intent(inout)  :: pnfft
  CMPLX,            intent(in)  :: in (:,:,:)
  R_TYPE,           intent(out) :: out(:,:,:)

  integer :: i1, i2, i3

  PUSH_SUB(X(pnfft_backward))

!   print *, mpi_world%rank, "<--- pnfft%N_local       ", pnfft%N_local 
!   print *, mpi_world%rank, "<--- pnfft%M             ", pnfft%M 
!   print *, mpi_world%rank, "<--- size(in)            ", size(in,1), size(in,2), size(in,3) 
!   print *, mpi_world%rank, "<--- size(pnfft%f_hat)   ", size(pnfft%f_hat,1), size(pnfft%f_hat,2), size(pnfft%f_hat, 3) 
!   print *, mpi_world%rank, "<--- size(out)           ", size(out,1), size(out,2), size(out,3)  
!   print *, mpi_world%rank, "<--- size(pnfft%f)       ", size(pnfft%f,1), size(pnfft%f,2), size(pnfft%f,3) 


!   do i1 = 1, pnfft%M(1)
!     do i2 = 1, pnfft%M(2)
!       do i3 = 1, pnfft%M(3)
!         pnfft%f_lin(pnfft_idx_3to1(pnfft,i3,i2,i1)) = cmplx(real(in(i1,i2,i3)), aimag(in(i1,i2,i3)),C_DOUBLE_COMPLEX)
!       end do
!     end do
!   end do

  pnfft%f(:,:,:) = in(:,:,:)


  call pnfft_adj(pnfft%plan)

  out(:,:,:) = pnfft%f_hat(:,:,:)

!   do i1 = 1,pnfft%N_local(1)
!     do i2 = 1, pnfft%N_local(2)
!       do i3 = 1, pnfft%N_local(3)
!         out(i2,i3,i1) = pnfft%f_hat(i1,i3,i2)
!       end do
!     end do
!   end do

  out = out/pnfft%norm

  POP_SUB(X(pnfft_backward))

end subroutine X(pnfft_backward)
