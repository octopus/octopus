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

subroutine mesh_init(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)

  integer :: i, j, k, morder

  sub_name = 'mesh_init'; call push_sub()

  call mesh_init_derivatives_coeff(m)

  call mesh1D_create(m)

  m%fft_n(1) = 2*m%nr(1) + 2; m%fft_n(2) = 0; m%fft_n(3) = 0
  m%hfft_n = m%fft_n(1)/2 + 1
  m%dplanf = int(-1, POINTER_SIZE)

  call oct_parse_double(C_string('DoubleFFTParameter'), 2.0_r8, m%fft_alpha)
  if (m%fft_alpha < 1.0_r8 .or. m%fft_alpha > 3.0_r8 ) then
    write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
    message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
    call write_fatal(2)
  end if

  m%fft_n2(1) = 2*nint(m%fft_alpha*m%nr(1)) + 2; m%fft_n2(2) = 0; m%fft_n2(3) = 0
  m%hfft_n2 = m%fft_n2(1)/2 + 1
  m%dplanf2 = int(-1, POINTER_SIZE)

  call oct_parse_int(C_string('DerivativesSpace'), REAL_SPACE, m%d%space)
  if(m%d%space < 0 .or. m%d%space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", m%d%space, &
         "' is not a valid DerivativesSpace"
    message(2) = '(0 <= DerivativesSpace <=1)'
    call write_fatal(2)
  end if

  ! order is very important
  ! init rs -> create -> init fft
  if(m%d%space == REAL_SPACE) then 
    message(1) = 'Info: Derivatives calculated in real-space'
  else
    call mesh_alloc_ffts(m, 1)
    message(1) = 'Info: Derivatives calculated in reciprocal-space'
  end if
  call write_info(1)

  call pop_sub()
end subroutine mesh_init

subroutine mesh_alloc_ffts(m, i)
  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: i

  sub_name = 'mesh_alloc_ffts'; call push_sub()

  if(i == 1 .and. m%dplanf == int(-1, POINTER_SIZE)) then
    call rfftwnd_f77_create_plan(m%dplanf, 1, (/ m%fft_n(1) /), &
         fftw_real_to_complex + fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftwnd_f77_create_plan(m%dplanb, 1, (/ m%fft_n(1) /), &
         fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
    call  fftwnd_f77_create_plan(m%zplanf, 1, (/ m%fft_n(1) /), &
         fftw_forward,  fftw_measure + fftw_threadsafe)
    call  fftwnd_f77_create_plan(m%zplanb, 1, (/ m%fft_n(1) /), &
         fftw_backward, fftw_measure + fftw_threadsafe) 
  else if(i == 2 .and. m%dplanf2 == int(-1, POINTER_SIZE)) then
    message(1) = "Info: FFTs used in a double box (for local potential)"
    write(message(2), '(6x,a,i4,a,i4,a,i4,a)') &
          'box size = (', m%fft_n2(1), ',', m%fft_n2(2), ',',m%fft_n2(3),')'
    write(message(3), '(6x,a,f12.5)') 'alpha = ', m%fft_alpha
    call write_info(3)

    call rfftwnd_f77_create_plan(m%dplanf2, 1, (/ m%fft_n2(1) /), &
         fftw_real_to_complex + fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftwnd_f77_create_plan(m%dplanb2, 1, (/ m%fft_n2(1) /),  &
         fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
  end if

  call pop_sub()
end subroutine mesh_alloc_ffts

subroutine mesh_gradient_in_FS(m, nx, n, f, grad)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nx, n(3)
  complex(r8), intent(IN)  :: f   (nx)
  complex(r8), intent(out) :: grad(nx)

  real(r8) :: temp, g2
  integer :: kx,ix
  
  temp = (2.0_r8*M_Pi)/(n(1)*m%h(1))
  do ix=1, nx
     kx = pad_feq(ix, n(1), .true.)
     g2 = temp*kx
     grad(ix) = g2 * M_zI * f(ix)
  enddo

end subroutine mesh_gradient_in_FS

subroutine mesh_laplacian_in_FS(m, nx, n, f, lapl)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nx, n(3)
  complex(r8), intent(IN)  :: f   (nx)
  complex(r8), intent(out) :: lapl(nx)

  real(r8) :: temp, g2
  integer :: kx, ix
  
  temp = (2.0_r8*M_Pi)/(n(1)*m%h(1))
  do ix=1, nx
     kx = pad_feq(ix, n(1), .true.)
     g2 = temp**2*kx**2
     lapl(ix) = - g2 * f(ix)
  enddo  

end subroutine mesh_laplacian_in_FS

! this actually adds to outp
subroutine phase_factor(m, n, vec, inp, outp)
  implicit none
  type(mesh_type), intent(IN) :: m
  integer, intent(in)         :: n(3)
  real(r8), intent(IN)        :: vec(3)
  complex(r8), intent(IN)     :: inp (n(1)/2+1)
  complex(r8), intent(inout)  :: outp(n(1)/2+1)
  
  complex(r8) :: k
  integer     :: ix, ixx
  
  k = M_zI * ((2.0_r8*M_Pi)/(n(1)*m%h(1)))

  do ix = 1, n(1)/2 + 1
     ixx = pad_feq(ix, n(1), .true.)
     outp(ix) = outp(ix) + exp( -(k*vec(1)*ixx) ) * inp(ix)
  end do

end subroutine phase_factor

#include "mesh1D_create.F90"

