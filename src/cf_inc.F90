!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!!! The following routines handle creation/destruction of the cube
subroutine X(cf_new)(n, cf)
  integer, intent(in) :: n(3)
  type(X(cf)), intent(out) :: cf

  call push_sub('cf_new')

  ASSERT(all(n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = n

#if defined(HAVE_FFT)
  nullify(cf%fft)
#endif
  call pop_sub()
end subroutine X(cf_new)

subroutine X(cf_new_from)(cf, cf_i)
  type(X(cf)), intent(out) :: cf
  type(X(cf)), intent( in) :: cf_i

  call push_sub('cf_new_from')
  ASSERT(all(cf_i%n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = cf_i%n

#if defined(HAVE_FFT)
  if(associated(cf_i%fft)) then
    allocate(cf%fft)
    call fft_copy(cf_i%fft, cf%fft)
    cf%nx = cf_i%nx
  else
    nullify(cf%fft)
  end if
#endif
  call pop_sub()
end subroutine X(cf_new_from)

subroutine X(cf_alloc_RS)(cf)
  type(X(cf)), intent(inout) :: cf

  call push_sub('cf_alloc_RS')

  ASSERT(.not.associated(cf%RS))
  allocate(cf%RS(cf%n(1), cf%n(2), cf%n(3)))

  call pop_sub()
end subroutine X(cf_alloc_RS)

subroutine X(cf_free_RS)(cf)
  type(X(cf)), intent(inout) :: cf

  call push_sub('cf_free_RS')

  ASSERT(associated(cf%RS))
  deallocate(cf%RS)
  nullify(cf%RS)

  call pop_sub()
end subroutine X(cf_free_RS)

#if defined(HAVE_FFT)
subroutine X(cf_alloc_FS)(cf)
  type(X(cf)), intent(inout) :: cf

  call push_sub('cf_alloc_FS')
  ASSERT(.not.associated(cf%FS))
  ASSERT(associated(cf%fft))

  allocate(cf%FS(cf%nx, cf%n(2), cf%n(3)))

  call pop_sub()
end subroutine X(cf_alloc_FS)

subroutine X(cf_free_FS)(cf)
  type(X(cf)), intent(inout) :: cf

  call push_sub('cf_free_FS')
  ASSERT(associated(cf%FS))
  deallocate(cf%FS)
  nullify(cf%FS)

  call pop_sub()
end subroutine X(cf_free_FS)
#endif

subroutine X(cf_free)(cf)
  type(X(cf)), intent(inout) :: cf

  call push_sub('cf_free')
  if(associated(cf%RS)) then
    deallocate(cf%RS)
    nullify(cf%RS)
  end if

  if(associated(cf%FS)) then
    deallocate(cf%FS)
    nullify(cf%FS)
  end if

#if defined(HAVE_FFT)
  if(associated(cf%fft)) then
    call fft_end(cf%fft)
    deallocate(cf%fft)
    nullify(cf%fft)
  end if
#endif

  call pop_sub()
end subroutine X(cf_free)

#if defined(HAVE_FFT)
!!! initializes the ffts. As the dimension of the fft may be adjusted, this
!!! routine has to be called before allocating anything
subroutine X(cf_fft_init)(cf)
  type(X(cf)), intent(inout)  :: cf

  call push_sub('cf_fft_init')
  ASSERT(.not.associated(cf%fft))
  ASSERT(.not.associated(cf%RS))
  ASSERT(.not.associated(cf%FS))
  
  allocate(cf%fft)
#ifdef R_TREAL
  call fft_init(cf%n, fft_real, cf%fft)
  cf%nx = cf%n(1)/2 + 1
#else
  call fft_init(cf%n, fft_complex, cf%fft)
  cf%nx = cf%n(1)
#endif

  call pop_sub()
end subroutine X(cf_fft_init)

!!! The next routines convert the function between real space and fourier space
!!! Note that the dimensions of the function in FS are different wether f 
!!! is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cf_RS2FS)(cf)
  type(X(cf)), intent(inout)  :: cf

  ASSERT(associated(cf%RS))
  if(.not.associated(cf%FS)) call X(cf_alloc_FS)(cf)

  call X(fft_forward)(cf%fft, cf%RS, cf%FS)

end subroutine X(cf_RS2FS)

subroutine X(cf_FS2RS)(cf)
  type(X(cf)), intent(inout)  :: cf

  ASSERT(associated(cf%FS))
  if(.not.associated(cf%RS)) call X(cf_alloc_RS)(cf)

  call X(fft_backward)(cf%fft, cf%FS, cf%RS)

end subroutine X(cf_FS2RS)

! This function calculates the laplacian of a function in Fourier space.
! optionally, one can apply a cutoff
! i.e. \nabla^2 f = min(cutoff, G^2) * f
subroutine X(cf_FS_lapl)(m, cf, cutoff_)
  type(mesh_type), intent(in) :: m
  type(X(cf)), intent(inout)  :: cf
  FLOAT, intent(in), optional :: cutoff_

  FLOAT :: cutoff, temp(3), g2
  integer :: k(3), ix, iy, iz

  ASSERT(associated(cf%FS))

  cutoff = CNST(1e10)
  if(present(cutoff_)) then
    if(cutoff_>M_ZERO) cutoff = cutoff_
  end if

  temp = M_ZERO
  temp(1:conf%dim) = (M_TWO*M_PI)/(cf%n(1:conf%dim)*m%h(1:conf%dim))

  do iz = 1, cf%n(3)
    k(3) = pad_feq(iz, cf%n(3), .true.)
    do iy = 1, cf%n(2)
      k(2) = pad_feq(iy, cf%n(2), .true.)
      do ix = 1, cf%nx
        k(1) = pad_feq(ix, cf%n(1), .true.)
        g2 = min(cutoff, sum((temp(1:conf%dim)*k(1:conf%dim))**2))
        cf%FS(ix, iy, iz) = -g2*cf%FS(ix, iy, iz)
      end do
    end do
  end do

end subroutine X(cf_FS_lapl)

subroutine X(cf_FS_grad)(m, cf, j)
  type(mesh_type), intent(in) :: m
  type(X(cf)), intent(inout)  :: cf
  integer, intent(in) :: j

  FLOAT :: temp(3), g
  integer :: k(3), ix, iy, iz

  call push_sub('cf_FS_grad')

  temp = M_ZERO
  temp(1:conf%dim) = (M_TWO*M_PI)/(cf%n(1:conf%dim)*m%h(1:conf%dim))

  k = 0
  do iz = 1, cf%n(3)
    if(j == 3) k(3) = pad_feq(iz, cf%n(3), .true.)
    do iy = 1, cf%n(2)
      if(j == 2) k(2) = pad_feq(iy, cf%n(2), .true.)
      do ix = 1, cf%nx
        if(j == 1) k(1) = pad_feq(ix, cf%n(1), .true.)

        g = sum(temp(1:conf%dim)*k(1:conf%dim))
        cf%FS(ix, iy, iz) = g * M_zI * cf%FS(ix, iy, iz)
      end do
    end do
  end do

  call pop_sub()
end subroutine X(cf_FS_grad)

#endif
