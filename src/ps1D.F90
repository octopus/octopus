module ps
use fdf
use units
use spline
use kb

implicit none

private
public :: ps_type, ps_init, ps_end

type ps_type
  type(spline_type), pointer :: kb(:)   ! Kleynman-Bylander projectors
  type(spline_type), pointer :: dkb(:)  ! derivatives of KB projectors
  type(spline_type), pointer :: Ur(:)   ! atomic wavefunctions
  type(spline_type), pointer :: vps(:)  ! pseudopotential
  type(spline_type) :: vlocal  ! local part
  type(spline_type) :: dvlocal ! derivative of the local part
  type(spline_type) :: core    ! core charge

  integer :: L_max ! maximum value of l to take
  integer :: L_loc ! which component to take as local
  character(len=4) :: icore
  real(r8) :: rc_max
  real(r8) :: vlocal_origin ! local pseudopotential at the orginin

  real(r8), pointer :: dkbcos(:), dknrm(:) ! KB cosinus and norm
end type ps_type

type ps_file
  character(len=2)  :: namatm, icorr
  character(len=3)  :: irel
  character(len=4) :: icore
  character(len=10) :: method(6) 
  character(len=70) :: titleps
  integer :: npotd, npotu, nr
  real(r8) :: b, a, zval
  real(r8), pointer :: rofi(:), vps(:,:), chcore(:), rho_val(:)
  integer :: nrval ! not in file, but very useful :)
end type ps_file

real(r8), parameter :: eps = 1.0e-8_r8
contains

subroutine ps_init(ps, label, z, zval)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: label
  real(r8), intent(inout) :: z
  real(r8), intent(inout) :: zval ! this is inside th ps file

  integer :: i

  sub_name = 'ps_init'; call push_sub()

  allocate(ps%vps(0:0))
  call spline_init(ps%vps(0))
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)
  call spline_init(ps%core)

!!$  ! now we load the necessary information from the ps file
  call ps_load(ps, z, zval)

  call pop_sub()
end subroutine ps_init

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i

  sub_name = 'ps_end'; call push_sub()

  deallocate(ps%vps)
  call spline_end(ps%vps(0))
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  call spline_end(ps%core)  
  
  call pop_sub()
end subroutine ps_end


subroutine ps_load(ps, z, zval)
  type(ps_type), intent(inout) :: ps
  real(r8), intent(in) :: z
  real(r8), intent(out) :: zval

  logical :: found
  integer :: i, iunit
  real(r8), allocatable :: rphi(:,:), eigen(:), rc(:)
  type(ps_file) :: psf


  psf%nr = 1200        ! I am hard-coding this, for the sake of good-old-fortran
  psf%b  = 0.0002
  psf%a   = 0.0125

  psf%nrval = psf%nr + mod((psf%nr + 1), 2) 


  allocate( psf%rofi(psf%nrval), psf%vps(psf%nrval, 0:0) )
  call fetch_data(psf, ps)

  ! now we fit the splines
  call get_splines(psf, ps)

  deallocate(psf%rofi, psf%vps)

! FIX UNITS (Should FIX instead kb.F90)
! Passing from Rydbergs -> Hartree
!  ps%vlocal_origin = ps%vlocal_origin / 2._r8

  return
end subroutine ps_load


subroutine fetch_data(psf, ps)
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(in) :: ps
  
  integer  :: ndown, nup, l, ir
  real(r8) :: r2

  do ir = 1, psf%nrval
     psf%rofi(ir) = psf%b*( exp(psf%a*(ir-1)) - 1)
     psf%vps(ir, 0) = -psf%zval/sqrt(1 + psf%rofi(ir)**2)
  enddo

  return
end subroutine fetch_data


subroutine get_splines(psf, ps)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  
  integer :: l, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  allocate(hato(psf%nrval), derhato(psf%nrval))


  hato = 0.0_r8
  hato(2:psf%nrval) = psf%vps(2:psf%nrval, 0)*psf%rofi(2:psf%nrval) + psf%zval
  hato(1) = 0.0_r8
  
  ! WARNING: Rydbergs -> Hartrees
  ! hato = hato / 2._r8
  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vps(1, 0)

! and the derivative now
  call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

  return
end subroutine get_splines


subroutine derivate_in_log_grid(a, b, nrval, f, dfdr)
  real(r8), intent(in) :: a, b
  integer,  intent(in) :: nrval
  real(r8), intent(IN) :: f(nrval)
  real(r8), intent(out) :: dfdr(nrval)

  real(r8) :: x,y
  integer :: i

  x = 1.0_r8 - exp(-2*a)
  y = 1.0_r8 - exp(-a)

  dfdr(1) = (1/(y*b))*exp(-a)*(f(2)-f(1))  
  do i = 2, nrval-1
    dfdr(i) = (1/(x*b))*exp(-i*a)*(f(i+1)-f(i-1))
  enddo
  dfdr(nrval) = (1/(y*b))*exp(-(nrval-1)*a)*(f(nrval)-f(nrval-1))

  return
end subroutine derivate_in_log_grid



end module ps
