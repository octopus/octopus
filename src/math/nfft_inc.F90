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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

  !--------------------------------------------
  subroutine X(nfft_forward)(nfft, in, out)
    type(nfft_t), intent(in)  :: nfft
    R_TYPE,       intent(in)  :: in(:,:,:)
    CMPLX,        intent(out) :: out(:,:,:)

    integer :: ix, iy, iz

    PUSH_SUB(X(nfft_forward))

    do ix = 1, nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call X(oct_set_f_hat)(nfft%plan, nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call oct_nfft_trafo(nfft%plan)

    do ix = 1, nfft%M(1)
      do iy = 1, nfft%M(2)
        do iz = 1, nfft%M(3)
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
    R_TYPE,       intent(out) :: out(:,:,:)

    integer :: ix, iy, iz

    PUSH_SUB(znfft_backward)

    do ix = 1, nfft%M(1)
      do iy = 1, nfft%M(2)
        do iz = 1, nfft%M(3)
          call X(oct_set_f)(nfft%plan, nfft%M, nfft%dim, in(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    call oct_nfft_adjoint(nfft%plan)

    do ix = 1,nfft%N(1)
      do iy = 1, nfft%N(2)
        do iz = 1, nfft%N(3)
          call X(oct_get_f_hat)(nfft%plan, nfft%dim, out(ix,iy,iz), ix, iy, iz)
        end do
      end do
    end do

    out = out/nfft%norm

    POP_SUB(X(nfft_backward))
  end subroutine X(nfft_backward)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
