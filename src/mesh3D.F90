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
  
  call mesh_create(m, natoms, atom)

  ! we will probably need ffts in a lot of places
#ifdef POLYMERS
  j = 2
#else
  j = 3
#endif
  m%fft_n(:)  = 2*m%nr(:) + 1
  do i = 1, j
    call oct_fft_optimize(m%fft_n(i), 7, 1) ! always ask for an odd number
  end do
  m%hfft_n = m%fft_n(1)/2 + 1
  m%dplanf = int(-1, POINTER_SIZE)
  
  call oct_parse_double(C_string('DoubleFFTParameter'), 2.0_r8, m%fft_alpha)
  if (m%fft_alpha < 1.0_r8 .or. m%fft_alpha > 3.0_r8 ) then
    write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
    message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
    call write_fatal(2)
  end if

  do i = 1, j
    m%fft_n2(i) = 2*nint(m%fft_alpha*maxval(m%nr)) + 1
    call oct_fft_optimize(m%fft_n2(i), 7, 1) ! always ask for an odd number
  end do
  do i = j+1, 3 ! for periodic systems
    m%fft_n2(i) = m%fft_n(i)
  end do
  
  ! But I will leave if for the moment.
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
    call rfftw3d_f77_create_plan(m%dplanf, m%fft_n(1), m%fft_n(2), m%fft_n(3), &
         fftw_real_to_complex + fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftw3d_f77_create_plan(m%dplanb, m%fft_n(1), m%fft_n(2), m%fft_n(3), &
         fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
    call  fftw3d_f77_create_plan(m%zplanf, m%fft_n(1), m%fft_n(2), m%fft_n(3), &
         fftw_forward,  fftw_measure + fftw_threadsafe)
    call  fftw3d_f77_create_plan(m%zplanb, m%fft_n(1), m%fft_n(2), m%fft_n(3), &
         fftw_backward, fftw_measure + fftw_threadsafe) 
  else if(i == 2 .and. m%dplanf2 == int(-1, POINTER_SIZE)) then
    message(1) = "Info: FFTs used in a double box (for poisson | local potential)"
    write(message(2), '(6x,a,i4,a,i4,a,i4,a)') &
          'box size = (', m%fft_n2(1), ',', m%fft_n2(2), ',',m%fft_n2(3),')'
    write(message(3), '(6x,a,f12.5)') 'alpha = ', m%fft_alpha
    call write_info(3)

    call rfftw3d_f77_create_plan(m%dplanf2, m%fft_n2(1), m%fft_n2(2), m%fft_n2(3), &
         fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftw3d_f77_create_plan(m%dplanb2, m%fft_n2(1), m%fft_n2(2), m%fft_n2(3), &
         fftw_backward, fftw_measure + fftw_threadsafe)
  end if

  call pop_sub()
end subroutine mesh_alloc_ffts

subroutine mesh_gradient_in_FS(m, nx, n, f, grad, dir)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nx, n(3), dir
  complex(r8), intent(IN)  :: f   (nx, n(2), n(3))
  complex(r8), intent(out) :: grad(nx, n(2), n(3))

  real(r8) :: temp(3), g2
  integer :: kx, ky, kz, ix, iy, iz
  
  temp(:) = (2.0_r8*M_Pi)/(n(:)*m%h(:))
  kx = 0; ky = 0; kz = 0

  do iz = 1, n(3)
    if(dir == 3) kz = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      if(dir == 2) ky = pad_feq(iy, n(2), .true.)
      do ix = 1, nx
        if(dir == 1) kx = pad_feq(ix, n(1), .true.)
        
        g2 = temp(1)*kx + temp(2)*ky + temp(3)*kz
        grad(ix, iy, iz) = g2 * M_zI * f(ix, iy, iz)
      end do
    end do
  end do

end subroutine mesh_gradient_in_FS

subroutine mesh_laplacian_in_FS(m, nx, n, f, lapl)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: nx, n(3)
  complex(r8), intent(IN)  :: f   (nx, n(2), n(3))
  complex(r8), intent(out) :: lapl(nx, n(2), n(3))

  real(r8) :: temp(3), g2
  integer :: kx, ky, kz, ix, iy, iz
  
  temp(:) = (2.0_r8*M_Pi)/(n(:)*m%h(:))

  do iz = 1, n(3)
    kz = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      ky = pad_feq(iy, n(2), .true.)
      do ix = 1, nx
        kx = pad_feq(ix, n(1), .true.)
        
        g2 = (temp(1)*kx)**2 + (temp(2)*ky)**2 + (temp(3)*kz)**2
        lapl(ix, iy, iz) = - g2 * f(ix, iy, iz)
      end do
    end do
  end do

end subroutine mesh_laplacian_in_FS

! this actually adds to outp
subroutine phase_factor(m, n, vec, inp, outp)
  implicit none
  type(mesh_type), intent(IN) :: m
  integer, intent(in)         :: n(3)
  real(r8), intent(IN)        :: vec(3)
  complex(r8), intent(IN)     :: inp (n(1)/2+1, n(2), n(3))
  complex(r8), intent(inout)  :: outp(n(1)/2+1, n(2), n(3))
  
  complex(r8) :: k(3)
  integer     :: ix, iy, iz, ixx, iyy, izz
  
  k(1:3) = M_zI * ((2.0_r8*M_Pi)/(n(1:3)*m%h(1:3)))
  do iz = 1, n(3)
    izz = pad_feq(iz, n(3), .true.)
    do iy = 1, n(2)
      iyy = pad_feq(iy, n(2), .true.)
      do ix = 1, n(1)/2 + 1
        ixx = pad_feq(ix, n(1), .true.)
        outp(ix, iy, iz) = outp(ix, iy, iz) + &
             exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * inp(ix, iy, iz)
      end do
    end do
  end do
end subroutine phase_factor

