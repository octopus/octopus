#include "config.h"

module mesh
use liboct
use global
use units
use fft
use atom
  
implicit none

integer, parameter ::     &
     REAL_SPACE = 0,      &
     RECIPROCAL_SPACE = 1

integer, parameter :: &
     SPHERE   = 1

type mesh_derivatives_type
  integer :: space ! 0 for FD, 1 for reciprocal-space

  integer :: norder ! order on the discretization of the gradient/laplace operator
  real(r8), pointer :: dgidfj(:) ! the coefficients for the gradient
  real(r8), pointer :: dlidfj(:) ! the coefficients for the laplacian
  
end type mesh_derivatives_type

type mesh_type
  integer  :: box_shape ! 1->sphere
  real(r8) :: h         ! the (constant) spacing between the points
  
  real(r8) :: rsize     ! the radius of the sphere
  real(r8) :: zsize     ! the length of the cilinder in the z direction
  
  integer  :: np        ! number of points in inner mesh
  real(r8) :: vol_pp    ! element of volume for integrations
  
  ! return x, y and z for each point
  integer, pointer :: Lx(:)
  ! return points # for each xyz
  integer, pointer :: Lxyz_inv(:)

  ! mesh elargement... this is currently only used in solving
  ! the poisson equation with conjugated gradients...
  integer :: nk        ! number of points in outer mesh
  integer, pointer :: Kx(:)

  ! some other vars
  integer :: nr  ! number of points per radius inside
  integer :: nx  ! the # of points in the whole radius
                 ! (rsize/h) + norder

  type(mesh_derivatives_type) :: d

  integer :: fft_n, hfft_n ! we always want to perform FFTs ;)
  integer(POINTER_SIZE) :: dplanf, dplanb, zplanf, zplanb
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
    message(1) = 'Info: Derivatives calculated in reciprocal-space'
  end if
  call write_info(1)

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
  
  !call mesh3D_create(m, natoms, atom)
  call mesh1D_create(m, natoms, atom)

  ! we will probably need ffts in a lot of places
  ! once again the overhead is small
  m%fft_n  = 2*m%nr + 2
  m%hfft_n = m%fft_n/2 + 1
  call rfftw_f77_create_plan(m%dplanf, m%fft_n, &
       fftw_forward, fftw_measure + fftw_threadsafe)
  call rfftw_f77_create_plan(m%dplanb, m%fft_n, &
       fftw_backward, fftw_measure + fftw_threadsafe)
  call fftw_f77_create_plan(m%zplanf, m%fft_n, &
       fftw_forward, fftw_measure + fftw_threadsafe)
  call fftw_f77_create_plan(m%zplanb, m%fft_n, &
       fftw_backward, fftw_measure + fftw_threadsafe)
!!$  call rfftw3d_f77_create_plan(m%dplanf, m%fft_n, m%fft_n, m%fft_n, &
!!$       fftw_real_to_complex + fftw_forward, fftw_measure + fftw_threadsafe)
!!$  call rfftw3d_f77_create_plan(m%dplanb, m%fft_n, m%fft_n, m%fft_n, &
!!$       fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
!!$  call  fftw3d_f77_create_plan(m%zplanf, m%fft_n, m%fft_n, m%fft_n, &
!!$       fftw_forward,  fftw_measure + fftw_threadsafe)
!!$  call  fftw3d_f77_create_plan(m%zplanb, m%fft_n, m%fft_n, m%fft_n, &
!!$       fftw_backward, fftw_measure + fftw_threadsafe) 

  call pop_sub()
end subroutine mesh_init

! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  sub_name = 'mesh_end'; call push_sub()

  if(associated(m%Lx)) then
    deallocate(m%Lx, m%Kx)
    nullify(m%Lx, m%Kx)
  end if
  if(associated(m%Lxyz_inv)) then
    deallocate(m%Lxyz_inv)
    nullify(m%Lxyz_inv)
  end if
  if(m%d%space == 0 .and. associated(m%d%dgidfj)) then
    deallocate(m%d%dgidfj, m%d%dlidfj)
    nullify(m%d%dgidfj, m%d%dlidfj)
  end if

  call fftw_f77_destroy_plan(m%dplanb)
  call fftw_f77_destroy_plan(m%dplanf)
  call fftw_f77_destroy_plan(m%zplanb)
  call fftw_f77_destroy_plan(m%zplanf)

  call pop_sub()
  return
end subroutine mesh_end

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(3) = (/ &
      'sphere       ', &
      'cilinder     ', &
      'around nuclei'/)

  sub_name = 'mesh_write_info'; call push_sub()

  write(message(1), '(a,a,1x,3a,f7.3)') '  Type = ', bs(m%box_shape), &
       ' Radius [', trim(units_out%length%abbrev), '] = ', m%rsize/units_out%length%factor
  write(message(2),'(3a, f6.3, 1x, 3a, f8.5)') &
       '  Spacing [', trim(units_out%length%abbrev), '] = ', m%h/units_out%length%factor, &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ', &
       (m%h/units_out%length%factor)**3
  write(message(3),'(a, i6, a, i6)') '  # inner mesh = ', m%np, &
      '   # outer mesh = ', m%nk

  call write_info(3, unit)

  call pop_sub()
  return
end subroutine mesh_write_info


subroutine mesh_x(m, i, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: x

  x = m%Lx(i)*m%H
end subroutine mesh_x


subroutine mesh_r(m, i, r, a, x)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: i
  real(r8), intent(out) :: r
  real(r8), intent(in), optional :: a
  real(r8), intent(out), optional :: x

  real(r8) :: xx

  call mesh_x(m, i, xx)
  if(present(a)) xx = xx - a
  r = abs(xx)

  if(present(x)) x = xx
end subroutine mesh_r

#include "mesh1D_create.F90"

#include "undef.F90"
#include "real.F90"
#include "mesh1D_inc.F90"
#include "undef.F90"
#include "complex.F90"
#include "mesh1D_inc.F90"
#include "undef.F90"

end module mesh
