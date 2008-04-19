!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
! The following routines handle creation/destruction of the cube
subroutine X(cf_new)(n, cf)
  integer, intent(in) :: n(MAX_DIM)
  type(X(cf_t)), intent(out) :: cf

  call push_sub('cf_inc.Xcf_new')

  ASSERT(all(n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = n

  nullify(cf%fft)
  call pop_sub()
end subroutine X(cf_new)


! ---------------------------------------------------------
subroutine X(cf_new_from)(cf, cf_i)
  type(X(cf_t)), intent(out) :: cf
  type(X(cf_t)), intent( in) :: cf_i

  call push_sub('cf_inc.Xcf_new_from')
  ASSERT(all(cf_i%n>0))

  nullify(cf%RS)
  nullify(cf%FS)
  cf%n = cf_i%n

  if(associated(cf_i%fft)) then
    ALLOCATE(cf%fft, 1)
    call fft_copy(cf_i%fft, cf%fft)
    cf%nx = cf_i%nx
  else
    nullify(cf%fft)
  end if
  call pop_sub()
end subroutine X(cf_new_from)


! ---------------------------------------------------------
subroutine X(cf_alloc_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_alloc_RS')

  ASSERT(.not.associated(cf%RS))
  ALLOCATE(cf%RS(cf%n(1), cf%n(2), cf%n(3)), cf%n(1)*cf%n(2)*cf%n(3))

  call pop_sub()
end subroutine X(cf_alloc_RS)


! ---------------------------------------------------------
subroutine X(cf_free_RS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_free_RS')

  ASSERT(associated(cf%RS))
  deallocate(cf%RS)
  nullify(cf%RS)

  call pop_sub()
end subroutine X(cf_free_RS)


! ---------------------------------------------------------
subroutine X(cf_alloc_FS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_alloc_FS')
  ASSERT(.not.associated(cf%FS))
  ASSERT(associated(cf%fft))

  ALLOCATE(cf%FS(cf%nx, cf%n(2), cf%n(3)), cf%nx*cf%n(2)*cf%n(3))

  call pop_sub()
end subroutine X(cf_alloc_FS)


! ---------------------------------------------------------
subroutine X(cf_free_FS)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_free_FS')
  ASSERT(associated(cf%FS))
  deallocate(cf%FS)
  nullify(cf%FS)

  call pop_sub()
end subroutine X(cf_free_FS)


! ---------------------------------------------------------
subroutine X(cf_free)(cf)
  type(X(cf_t)), intent(inout) :: cf

  call push_sub('cf_inc.Xcf_free')

  if(associated(cf%RS)) then
    deallocate(cf%RS)
    nullify(cf%RS)
  end if

  if(associated(cf%FS)) then
    deallocate(cf%FS)
    nullify(cf%FS)
  end if

  if(associated(cf%fft)) then
    call fft_end(cf%fft)
    deallocate(cf%fft)
    nullify(cf%fft)
  end if

  call pop_sub()
end subroutine X(cf_free)


! ---------------------------------------------------------
! initializes the ffts. As the dimension of the fft may be adjusted, this
! routine has to be called before allocating anything
subroutine X(cf_fft_init)(cf, sb)
  type(X(cf_t)),     intent(inout) :: cf
  type(simul_box_t), intent(in)    :: sb

  call push_sub('cf_inc.Xcf_fft_init')
  ASSERT(.not.associated(cf%fft))
  ASSERT(.not.associated(cf%RS))
  ASSERT(.not.associated(cf%FS))

  ALLOCATE(cf%fft, 1)
#ifdef R_TREAL
  call fft_init(cf%n, fft_real, cf%fft, optimize = .not.simul_box_is_periodic(sb))
  cf%nx = cf%n(1)/2 + 1
#else
  call fft_init(cf%n, fft_complex, cf%fft, optimize = .not.simul_box_is_periodic(sb))
  cf%nx = cf%n(1)
#endif

  call pop_sub()
end subroutine X(cf_fft_init)


! ---------------------------------------------------------
! The next routines convert the function between real space and fourier space
! Note that the dimensions of the function in FS are different wether f
! is real or complex, because the FFT representation is different (FFTW scheme).
subroutine X(cf_RS2FS)(cf)
  type(X(cf_t)), intent(inout)  :: cf

  ASSERT(associated(cf%RS))
  if(.not.associated(cf%FS)) call X(cf_alloc_FS)(cf)

  call X(fft_forward)(cf%fft, cf%RS, cf%FS)

end subroutine X(cf_RS2FS)

! ---------------------------------------------------------
subroutine X(cf_FS2RS)(cf)
  type(X(cf_t)), intent(inout)  :: cf

  ASSERT(associated(cf%FS))
  if(.not.associated(cf%RS)) call X(cf_alloc_RS)(cf)

  call X(fft_backward)(cf%fft, cf%FS, cf%RS)

end subroutine X(cf_FS2RS)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
