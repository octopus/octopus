module mesh
use global
use units
use fft
use atom
  
implicit none

integer, parameter ::     &
     REAL_SPACE = 0,      &
     RECIPROCAL_SPACE = 1

integer, parameter :: &        ! Only SPHERE is permitted in ONE_D.
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

! The type is equal in 1D and in 3D, although in 1D most of the variables
! are not needed.
type mesh_type
  integer  :: box_shape ! 1->sphere
  real(r8) :: h(3)      ! the (constant) spacing between the points
  
  real(r8) :: rsize     ! the radius of the sphere or of the cylinder
  real(r8) :: zsize     ! the length of the cylinder in the z direction
  real(r8) :: lsize(3)  ! the lengths of the parallelepipeds in each direction.
  
  integer  :: np        ! number of points in inner mesh
  real(r8) :: vol_pp    ! element of volume for integrations
  
  ! return x, y and z for each point
  integer, pointer :: Lx(:), Ly(:), Lz(:)
  ! return points # for each xyz
  integer, pointer :: Lxyz_inv(:,:,:)

  ! mesh elargement... this is currently only used in solving
  ! the poisson equation with conjugated gradients...
  integer :: nk        ! number of points in outer mesh
  integer, pointer :: Kx(:), Ky(:), Kz(:)

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
  use math, only: weights

  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: natoms
  type(atom_type), pointer, optional :: atom(:)

  integer :: i, j, k, morder
  real(r8), allocatable :: cc(:,:,:)

  sub_name = 'mesh_init'; call push_sub()

  ! we may still need this for cg hartree
  ! the overhead is small, so we calculate it always :))
  call oct_parse_int(C_string('OrderDerivatives'), 4, k)
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



! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  sub_name = 'mesh_end'; call push_sub()

  if(associated(m%Lx)) then
    deallocate(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
    nullify(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
  end if
  if(associated(m%Lxyz_inv)) then
    deallocate(m%Lxyz_inv)
    nullify(m%Lxyz_inv)
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

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(3) = (/ &
      'sphere       ', &
      'cylinder     ', &
      'around nuclei'/)

  sub_name = 'mesh_write_info'; call push_sub()

  if(m%box_shape == SPHERE .or. m%box_shape == CYLINDER  .or. m%box_shape == MINIMUM) then
    write(message(1), '(a,a,1x,3a,f7.3)') '  Type = ', bs(m%box_shape), &
       ' Radius [', trim(units_out%length%abbrev), '] = ', m%rsize/units_out%length%factor
  endif
  if(m%box_shape == CYLINDER) then
    write(message(1), '(a,3a,f7.3)') trim(message(1)), ', zlength [', &
         trim(units_out%length%abbrev), '] = ', m%zsize/units_out%length%factor
  end if
  if(m%box_shape == PARALLELEPIPED) then
    write(message(1),'(3a, a, f6.3, a, f6.3, a, f6.3, a)') &
       '  Lengths [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%lsize(1)/units_out%length%factor, ',',                     &
            m%lsize(2)/units_out%length%factor, ',',                     &
            m%lsize(3)/units_out%length%factor, ')'
  endif
  write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)') &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',         &
       '(', m%h(1)/units_out%length%factor, ',',                     &
            m%h(2)/units_out%length%factor, ',',                     &
            m%h(3)/units_out%length%factor, ')',                     &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ', &
       (m%h(1)*m%h(2)*m%h(3))/units_out%length%factor**3
  write(message(3),'(a, i6, a, i6)') '  # inner mesh = ', m%np, &
      '   # outer mesh = ', m%nk

  call write_info(3, unit)

  call pop_sub()
  return
end subroutine mesh_write_info

! subroutines to get xyzr
subroutine mesh_xyz(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x(3)

  x(1) = m%Lx(i)*m%H(1)
  x(2) = m%Ly(i)*m%H(2)
  x(3) = m%Lz(i)*m%H(3)
end subroutine mesh_xyz

subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x

  x = m%Lx(i)*m%H(1)
end subroutine mesh_x

subroutine mesh_y(m, i, y)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: y

  y = m%Ly(i)*m%H(2)
end subroutine mesh_y

subroutine mesh_z(m, i, z)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: z

  z = m%Lz(i)*m%H(3)
end subroutine mesh_z

subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: r
  real(r8), intent(in), optional :: a(3)
  real(r8), intent(out), optional :: x(3)

  real(r8) :: xx(3)

  call mesh_xyz(m, i, xx)
  if(present(a)) xx = xx - a
  r = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)

  if(present(x)) x = xx
end subroutine mesh_r

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

#include "undef.F90"
#include "real.F90"
#include "mesh1D_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "mesh1D_inc.F90"
#include "undef.F90"

end module mesh
