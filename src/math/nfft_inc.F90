!! Copyright (C) 2011 U. De Giovannini
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


subroutine X(nfft_forward)(nfft, in, out)
  type(nfft_t), intent(in)  :: nfft
  R_TYPE,        intent(in)  :: in(:,:,:)
  CMPLX,        intent(out) :: out(:,:,:)

  integer :: ix, iy, iz, MM(3)

  PUSH_SUB(X(nfft_forward))
    
  select case(nfft%dim)
    case (1)
      MM(1) = nfft%M
      MM(2) = 1
      MM(3) = 1

    case (2)
      MM(1) = nfft%M
      MM(2) = nfft%M
      MM(3) = 1

    case (3)
      MM(1) = nfft%M
      MM(2) = nfft%M
      MM(3) = nfft%M

  end select 

  do ix = 1, nfft%N(1)
    do iy = 1, nfft%N(2)
      do iz = 1, nfft%N(3)
        call X(oct_set_f_hat)(nfft%plan, nfft%N(1), nfft%dim, in(ix,iy,iz), ix, iy, iz)
      end do
    end do
  end do


  call oct_nfft_trafo(nfft%plan)


  do ix = 1, MM(1)
    do iy = 1, MM(2)
      do iz = 1, MM(3)
        call X(oct_get_f)(nfft%plan, nfft%M, nfft%dim, out(ix,iy,iz), ix, iy, iz)
      end do
    end do
  end do

  POP_SUB(X(nfft_forward))
    
end subroutine X(nfft_forward)


! ---------------------------------------------------------
subroutine X(nfft_backward)(nfft, in, out)
  type(nfft_t), intent(in)  :: nfft
  CMPLX,        intent(in)  :: in (:,:,:)
  R_TYPE,        intent(out) :: out(:,:,:)

  integer :: ix, iy, iz, MM(3)

  PUSH_SUB(X(nfft_backward))

  select case(nfft%dim)
    case (1)
      MM(1) = nfft%M
      MM(2) = 1
      MM(3) = 1

    case (2)
      MM(1) = nfft%M
      MM(2) = nfft%M
      MM(3) = 1

    case (3)
      MM(1) = nfft%M
      MM(2) = nfft%M
      MM(3) = nfft%M

  end select

  do ix = 1, MM(1)
    do iy = 1, MM(2)
      do iz = 1, MM(3)
        call X(oct_set_f)(nfft%plan, nfft%M, nfft%dim, in(ix,iy,iz), ix, iy, iz)
      end do
    end do
  end do

  call oct_nfft_adjoint(nfft%plan)

  do ix = 1,nfft%N(1)
    do iy = 1, nfft%N(2)
      do iz = 1, nfft%N(3)
        call X(oct_get_f_hat)(nfft%plan, nfft%M, nfft%dim, out(ix,iy,iz), ix, iy, iz)
      end do
    end do
  end do

  out = out/nfft%norm

  POP_SUB(X(nfft_backward))

end subroutine X(nfft_backward)
