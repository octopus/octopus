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

#include "config_F90.h"
module mesh
use global
use oct_parser
use units
use fft
use atom
  
implicit none

integer, parameter ::     &
     REAL_SPACE = 0,      &
     RECIPROCAL_SPACE = 1

integer, parameter :: &
     SPHERE   = 1,    &
     CYLINDER = 2,    &
     MINIMUM  = 3,    &
     PARALLELEPIPED = 4

type mesh_derivatives_type
  integer :: space ! 0 for FD, 1 for reciprocal-space

  integer :: norder ! order on the discretization of the gradient/laplace operator
  real(r8), pointer :: dgidfj(:) ! the coefficients for the gradient
  real(r8), pointer :: dlidfj(:) ! the coefficients for the laplacian
  
end type mesh_derivatives_type

type mesh_type
  integer  :: box_shape ! 1->sphere, 2->cylinder, 3->sphere around each atom,
                        ! 4->parallelpiped (orthonormal, up to now).
  real(r8) :: h(3)      ! the (constant) spacing between the points
  logical  :: iso
  
  real(r8) :: rsize     ! the radius of the sphere or of the cylinder
  real(r8) :: xsize     ! the length of the cylinder in the x direction
  real(r8) :: lsize(3)  ! half of the length of the parallelepiped in each direction.
  real(r8) :: rlat(3,3)   ! lattice primitive vectors
  real(r8) :: klat(3,3)   ! reciprocal lattice primitive vectors
  real(r8) :: shift(27,3) ! shift to equivalent positions in nearest neighbour primitive cells
  
  integer  :: np        ! number of points in inner mesh
  real(r8) :: vol_pp    ! element of volume for integrations

  integer  :: dnpw, znpw  ! number of plane waves for Fourier space representation.
  real(r8) :: vol_ppw     ! The volume in reciprocal space
  
  ! return x, y and z for each point
  integer, pointer :: Lxyz(:,:)
  ! return points # for each xyz
  integer, pointer :: Lxyz_inv(:,:,:)
  integer, pointer :: ind(:, :, :)

  ! mesh elargement... this is currently only used in solving
  ! the poisson equation with conjugated gradients...
  integer :: nk        ! number of points in outer mesh
  integer, pointer :: Kxyz(:,:)

  ! some other vars
  integer :: nr(3)  ! number of points per radius inside
  integer :: nx(3)  ! the # of points in the whole radius

  type(mesh_derivatives_type) :: d

  integer :: fft_n(3), hfft_n ! we always want to perform FFTs ;)
  integer(POINTER_SIZE) :: dplanf, dplanb, zplanf, zplanb
  real(r8) :: fft_alpha
  integer :: fft_n2(3), hfft_n2
  integer(POINTER_SIZE) :: dplanf2, dplanb2
end type mesh_type

contains

subroutine mesh_init(m, natoms, atom)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)

  integer :: i, k, morder
  logical :: fft_optimize

  sub_name = 'mesh_init'; call push_sub()

  call mesh_init_derivatives_coeff(m)
  
  call mesh_create(m, natoms, atom)

  ! we will probably need ffts in a lot of places

  call oct_parse_logical("FFTOptimize", .true., fft_optimize)
  
  ! only non-periodic directions are optimized
  m%fft_n(:)  = 2*m%nr(:) + 1
  do i = 1, conf%periodic_dim
    m%fft_n(i) = m%fft_n(i) - 1
  end do
  do i = conf%periodic_dim+1, conf%dim ! always ask for an odd number
     if(m%fft_n(i).ne.1 .and. fft_optimize) call oct_fft_optimize(m%fft_n(i), 7, 1)
  end do
  m%hfft_n = m%fft_n(1)/2 + 1
  m%dplanf = int(-1, POINTER_SIZE)
  m%dnpw = m%hfft_n*m%fft_n(2)*m%fft_n(3)
  m%znpw = m%fft_n(1)*m%fft_n(2)*m%fft_n(3)
  m%vol_ppw = m%vol_pp/m%znpw

  call oct_parse_double('DoubleFFTParameter', 2.0_r8, m%fft_alpha)
  if (m%fft_alpha < 1.0_r8 .or. m%fft_alpha > 3.0_r8 ) then
    write(message(1), '(a,f12.5,a)') "Input: '", m%fft_alpha, &
         "' is not a valid DoubleFFTParameter"
    message(2) = '1.0 <= DoubleFFTParameter <= 3.0'
    call write_fatal(2)
  end if

  ! fft box is doubled and optimized only in non-periodic directions
  m%fft_n2 = 1
  do i = conf%periodic_dim+1, conf%dim
    m%fft_n2(i) = 2*nint(m%fft_alpha*maxval(m%nr)) + 1
    if(fft_optimize) call oct_fft_optimize(m%fft_n2(i), 7, 1) ! always ask for an odd number
  end do
  do i = 1, conf%periodic_dim ! for periodic systems
    m%fft_n2(i) = m%fft_n(i)
  end do
  
  m%hfft_n2 = m%fft_n2(1)/2 + 1
  m%dplanf2 = int(-1, POINTER_SIZE)

  call oct_parse_int('DerivativesSpace', REAL_SPACE, m%d%space)
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
    message(1) = 'Info: Derivatives calculated in reciprocal-space'
  end if
  call mesh_alloc_ffts(m, 1)
  call write_info(1)

  call pop_sub()
end subroutine mesh_init

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(4) = (/ &
      'sphere        ', &
      'cylinder      ', &
      'around nuclei ', &
      'parallelepiped'/)

  sub_name = 'mesh_write_info'; call push_sub()
  
  write(message(1), '(a,a,1x)') '  Type = ', bs(m%box_shape)
  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER  .or. m%box_shape == MINIMUM) then
    write(message(2), '(3a,f7.3)') '  Radius  [', trim(units_out%length%abbrev), '] = ', &
                                   m%rsize/units_out%length%factor
  endif
  if(m%box_shape == CYLINDER) then
    write(message(2), '(a,3a,f7.3)') trim(message(2)), ', xlength [', &
         trim(units_out%length%abbrev), '] = ', m%xsize/units_out%length%factor
  end if
  if(m%box_shape == PARALLELEPIPED) then
    write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a)') &
       '  Lengths [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%lsize(1)/units_out%length%factor, ',',                     &
            m%lsize(2)/units_out%length%factor, ',',                     &
            m%lsize(3)/units_out%length%factor, ')'
  endif
  write(message(3),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)') &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%h(1)/units_out%length%factor, ',',                     &
            m%h(2)/units_out%length%factor, ',',                     &
            m%h(3)/units_out%length%factor, ')',                     &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ', &
       m%vol_pp/units_out%length%factor**3
  write(message(4),'(a, i6, a, i6)') '  # inner mesh = ', m%np, &
      '   # outer mesh = ', m%nk
  call write_info(4, unit)

  write(message(1),'(a,f9.3,a)') '  Grid Cutoff ['+trim(units_out%energy%abbrev)+'] = ', &
                                 (M_PI**2/(M_TWO*maxval(m%h)**2))/units_out%energy%factor
      
  if (conf%periodic_dim > 0) then
    write(message(1),'(1x)')
    write(message(2),'(a,3a,a)')'Lattice Primitive Vectors [', trim(units_out%length%abbrev), ']'
    write(message(3),'(a,f8.3)')'    x axis ', &
                    m%rlat(1,1)/units_out%length%factor
    write(message(4),'(a,f8.3)')'    y axis ', &
                    m%rlat(2,2)/units_out%length%factor
    write(message(5),'(a,f8.3)')'    z axis ', &
                    m%rlat(3,3) /units_out%length%factor
    write(message(6),'(a,3a,a)') 'Reciprocal Lattice Primitive Vectors [', &
                    trim(units_out%length%abbrev), '^-1]'
    write(message(7),'(a,f8.3)')'  k_x axis ', &
                    m%klat(1,1)*units_out%length%factor
    write(message(8),'(a,f8.3)')'  k_y axis ', &
                    m%klat(2,2)*units_out%length%factor
    write(message(9),'(a,f8.3)')'  k_z axis ', &
                    m%klat(3,3)*units_out%length%factor
  call write_info(9, unit)
  end if


  call pop_sub()
  return
end subroutine mesh_write_info

! subroutines to get xyzr
subroutine mesh_xyz(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x(conf%dim)

  x(1:conf%dim) = m%Lxyz(1:conf%dim, i)*m%h(1:conf%dim)
end subroutine mesh_xyz

subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x

  x = m%Lxyz(1, i)*m%h(1)
end subroutine mesh_x

subroutine mesh_y(m, i, y)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: y

  y = m%Lxyz(2, i)*m%h(2)
end subroutine mesh_y

subroutine mesh_z(m, i, z)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: z

  z = m%Lxyz(3, i)*m%h(3)
end subroutine mesh_z

subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: r
  real(r8), intent(in), optional :: a(conf%dim)
  real(r8), intent(out), optional :: x(conf%dim)

  real(r8) :: xx(conf%dim)

  call mesh_xyz(m, i, xx)
  if(present(a)) xx = xx - a
  r = sqrt(sum(xx(:)**2))

  if(present(x)) x = xx
end subroutine mesh_r

subroutine mesh_init_derivatives_coeff(m)
  use math, only: weights
  type(mesh_type), intent(inout) :: m

  integer :: k, i, j, morder
  real(r8), allocatable :: cc(:,:,:)

  call oct_parse_int('OrderDerivatives', 4, k)
  m%d%norder = k
  if (k < 1) then
    write(message(1), '(a,i4,a)') "Input: '", k, "' is not a valid OrderDerivatives"
    message(2) = '(1 <= OrderDerivatives)'
    call write_fatal(2)
  end if
  allocate(m%d%dlidfj(-k:k), m%d%dgidfj(-k:k))
  morder = 2*k
  allocate(cc(0:morder, 0:morder, 0:2))
  call weights(2, morder, cc)
  m%d%dgidfj(0) = cc(0, morder, 1)
  m%d%dlidfj(0) = cc(0, morder, 2)
  
  j = 1
  do i = 1, k
    m%d%dgidfj(-i) = cc(j, morder, 1)
    m%d%dlidfj(-i) = cc(j, morder, 2)
    j = j + 1
    m%d%dgidfj( i) = cc(j, morder, 1)
    m%d%dlidfj( i) = cc(j, morder, 2)
    j = j + 1
  end do
  deallocate(cc)

end subroutine mesh_init_derivatives_coeff

subroutine mesh_alloc_ffts(m, i)
  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: i

  sub_name = 'mesh_alloc_ffts'; call push_sub()

  if(i == 1 .and. m%dplanf == int(-1, POINTER_SIZE)) then
    call rfftwnd_f77_create_plan(m%dplanf, conf%dim, m%fft_n, &
         fftw_real_to_complex + fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftwnd_f77_create_plan(m%dplanb, conf%dim, m%fft_n, &
         fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
    call  fftwnd_f77_create_plan(m%zplanf, conf%dim, m%fft_n, &
         fftw_forward,  fftw_measure + fftw_threadsafe)
    call  fftwnd_f77_create_plan(m%zplanb, conf%dim, m%fft_n, &
         fftw_backward, fftw_measure + fftw_threadsafe) 
  else if(i == 2 .and. m%dplanf2 == int(-1, POINTER_SIZE)) then
    if (conf%periodic_dim == conf%dim) then
      message(1)='Info: System is fully periodic: ignoring DoubleFFTParameter'
      call write_info(1)
      write(message(1), '(6x,a,i4,a,i4,a,i4,a)') &
         'box size = (', m%fft_n2(1), ',', m%fft_n2(2), ',',m%fft_n2(3),')'
      call write_info(1)
    else  if (conf%periodic_dim + 1 <= conf%dim) then
      message(1) = "Info: FFTs used in a double box (for poisson | local potential)"
      write(message(2),'(6x,a,i1,a,i1)') &
           'doubled axis: from ',conf%periodic_dim+1,' to ',conf%dim
      write(message(3), '(6x,a,i4,a,i4,a,i4,a)') &
         'box size = (', m%fft_n2(1), ',', m%fft_n2(2), ',',m%fft_n2(3),')'
      write(message(4), '(6x,a,f12.5)') 'alpha = ', m%fft_alpha
      call write_info(4)
    endif
    
    call rfftwnd_f77_create_plan(m%dplanf2, conf%dim, m%fft_n2, &
         fftw_forward, fftw_measure + fftw_threadsafe)
    call rfftwnd_f77_create_plan(m%dplanb2, conf%dim, m%fft_n2, &
         fftw_backward, fftw_measure + fftw_threadsafe)
  end if

  call pop_sub()
end subroutine mesh_alloc_ffts

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
  
  k = M_z0
  k(1:conf%dim) = M_zI * ((2.0_r8*M_Pi)/(n(1:conf%dim)*m%h(1:conf%dim)))
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

! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  sub_name = 'mesh_end'; call push_sub()

  if(associated(m%Lxyz)) then
    deallocate(m%Lxyz, m%Kxyz)
    nullify(m%Lxyz, m%Kxyz)
  end if
  if(associated(m%Lxyz_inv)) then
    deallocate(m%Lxyz_inv, m%ind)
    nullify(m%Lxyz_inv, m%ind)
  end if
  if(m%d%space == 0 .and. associated(m%d%dgidfj)) then
    deallocate(m%d%dgidfj, m%d%dlidfj)
    nullify(m%d%dgidfj, m%d%dlidfj)
  end if

  if(m%dplanf.ne.int(-1, POINTER_SIZE)) then
    call fftw_f77_destroy_plan(m%dplanf)
    call fftw_f77_destroy_plan(m%dplanb)
    call fftw_f77_destroy_plan(m%zplanf)
    call fftw_f77_destroy_plan(m%zplanb)
  end if

  if(m%dplanf2.ne.int(-1, POINTER_SIZE)) then
    call fftw_f77_destroy_plan(m%dplanf2)
    call fftw_f77_destroy_plan(m%dplanb2)
  end if

  call pop_sub()
  return
end subroutine mesh_end

#include "mesh_create.F90"

#include "undef.F90"
#include "real.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_inc.F90"

end module mesh
