module ps
use io
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

!subroutine ps_init(ps, label, z, lmax, lloc, zval)
subroutine ps_init(ps, label, z, zval)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: label
  !integer, intent(in) :: lmax, lloc
  real(r8), intent(in) :: z
  real(r8), intent(out) :: zval ! this is inside th ps file

  integer :: i

  sub_name = 'ps_init'; call push_sub()

  !ps%L_max = lmax
  !ps%L_loc = lloc
  ps%l_max = 0
  ps%l_loc = 0

  !First we allocate all the stuff
  !allocate(ps%kb(0:ps%L_max), ps%dkb(0:ps%L_max), &
  !     ps%Ur(0:ps%L_max), ps%vps(0:ps%L_max))
  allocate(ps%vps(0:ps%l_max))

  do i = 0, ps%L_max
  !  call spline_init(ps%kb(i))
  !  call spline_init(ps%dkb(i))
  !  call spline_init(ps%Ur(i))
    call spline_init(ps%vps(i))
  end do
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)
  !call spline_init(ps%core)

  !allocate(ps%dkbcos(0:ps%L_max), ps%dknrm(0:ps%L_max))

  ! now we load the necessary information from the ps file
  !call ps_load(ps, trim(label)//'.vps', z, zval)
  call ps_load(ps, z, zval)

  call pop_sub()
end subroutine ps_init

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i

  sub_name = 'ps_end'; call push_sub()

  !if(.not. associated(ps%kb)) return

  do i = 0, ps%L_max
  !  call spline_end(ps%kb(i))
  !  call spline_end(ps%dkb(i))
  !  call spline_end(ps%Ur(i))
    call spline_end(ps%vps(i))
  end do

  !deallocate(ps%kb, ps%dkb, ps%Ur, ps%vps)
  
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  !call spline_end(ps%core)  

  !deallocate(ps%dkbcos, ps%dknrm)
  
  call pop_sub()
end subroutine ps_end

!subroutine ps_load(ps, filename, z, zval)
subroutine ps_load(ps, z, zval)
  type(ps_type), intent(inout) :: ps
  !character(len=*), intent(in) :: filename
  real(r8), intent(in) :: z
  real(r8), intent(in) :: zval

  logical :: found
  integer :: i, iunit
  real(r8), allocatable :: rphi(:,:), eigen(:), rc(:)
  type(ps_file) :: psf

  !read(iunit) psf%namatm, psf%icorr, psf%irel, psf%icore,     &
  !     (psf%method(i),i=1,6), psf%titleps, psf%npotd, psf%npotu, &
  !     psf%nr, psf%b, psf%a, psf%zval
  psf%namatm=''; psf%icorr=''; psf%irel=''; psf%icore='nc';
  psf%method(1:6)=''; psf%titleps=''
  psf%npotd = 1; psf%npotu = 1;
  psf%nr = 1001;
  psf%b  = 0.0002253411_r8
  psf%a  = 0.0125000000_r8
  !zval = psf%zval - Now the zval is read from input file by the user.
  psf%zval = zval
  ps%icore = psf%icore
  !call write_info_about_pseudo_1(stdout, psf, ps, z)

  !ps%L_max = min(psf%npotd - 1, ps%L_max)
  psf%nrval = psf%nr + mod((psf%nr + 1), 2) 

! Reads data file:
!       rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) )
!       vps(1:nrval,0:spec%ps_lmax) : pseudopotential functions

  allocate(psf%rofi(psf%nrval), psf%vps(psf%nrval,0:ps%l_max))
  call read_file_data(iunit, psf, ps)

  call get_splines(psf, ps)

  deallocate(psf%rofi, psf%vps)

  return
end subroutine ps_load

subroutine read_file_data(unit, psf, ps)
  integer, intent(in) :: unit
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(in) :: ps
  
  integer  :: ndown, nup, l, ir
  real(r8) :: r2

  do ir = 1, psf%nrval
     psf%rofi(ir) = psf%b*( exp( psf%a*(ir-1) ) - 1 )
     psf%vps(ir, 0) = -psf%zval/(sqrt(1.0_r8 + psf%rofi(ir)**2))
  enddo

  return
end subroutine read_file_data

!subroutine get_splines(psf, ps, rphi, rc)
subroutine get_splines(psf, ps)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  !real(r8), intent(in) :: rphi(psf%nrval, 0:ps%L_max), rc(0:ps%L_max)
  
  integer :: l, nrc, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  allocate(hato(psf%nrval), derhato(psf%nrval))

!...local part...
  hato = 0.0_r8
  !nrc = nint(log(rc(ps%L_max + 1)/psf%b + 1.0_r8)/psf%a) + 1

  !hato(2:psf%nrval) = psf%vps(2:psf%nrval, ps%L_loc)*psf%rofi(2:psf%nrval) + 2.0_r8*psf%zval
  !hato(1) = 2.0_r8*psf%zval
  hato(1:psf%nrval) = psf%vps(1:psf%nrval,0)*psf%rofi(1:psf%nrval)  + psf%zval

  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vps(1, ps%L_loc)

! and the derivative now
  call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

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
