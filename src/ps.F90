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

subroutine ps_init(ps, label, z, lmax, lloc, zval)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: label
  integer, intent(in) :: lmax, lloc
  real(r8), intent(in) :: z
  real(r8), intent(out) :: zval ! this is inside th ps file

  integer :: i

  sub_name = 'ps_init'; call push_sub()

  ps%L_max = lmax
  ps%L_loc = lloc

  !First we allocate all the stuff
  allocate(ps%kb(0:ps%L_max), ps%dkb(0:ps%L_max), &
       ps%Ur(0:ps%L_max), ps%vps(0:ps%L_max))

  do i = 0, ps%L_max
    call spline_init(ps%kb(i))
    call spline_init(ps%dkb(i))
    call spline_init(ps%Ur(i))
    call spline_init(ps%vps(i))
  end do
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)
  call spline_init(ps%core)

  allocate(ps%dkbcos(0:ps%L_max), ps%dknrm(0:ps%L_max))

  ! now we load the necessary information from the ps file
  call ps_load(ps, trim(label)//'.vps', z, zval)

  call pop_sub()
end subroutine ps_init

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i

  sub_name = 'ps_end'; call push_sub()

  if(.not. associated(ps%kb)) return

  do i = 0, ps%L_max
    call spline_end(ps%kb(i))
    call spline_end(ps%dkb(i))
    call spline_end(ps%Ur(i))
    call spline_end(ps%vps(i))
  end do

  deallocate(ps%kb, ps%dkb, ps%Ur, ps%vps)
  
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  call spline_end(ps%core)  

  deallocate(ps%dkbcos, ps%dknrm)
  
  call pop_sub()
end subroutine ps_end

subroutine ps_load(ps, filename, z, zval)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: filename
  real(r8), intent(in) :: z
  real(r8), intent(out) :: zval

  logical :: found
  integer :: i, iunit
  real(r8), allocatable :: rphi(:,:), eigen(:), rc(:)
  type(ps_file) :: psf

  inquire(file=filename, exist=found)
  if(.not.found) then
     message(1) = "Pseudopotential file '"//trim(filename)//"' not found"
     call write_fatal(1)
  endif
  call io_assign(iunit)
  open(iunit, file=filename, form='unformatted', status='unknown')
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads the header line of the file, with general info about the ps.
!       Writes down info to stdout.

  read(iunit) psf%namatm, psf%icorr, psf%irel, psf%icore,     &
       (psf%method(i),i=1,6), psf%titleps, psf%npotd, psf%npotu, &
       psf%nr, psf%b, psf%a, psf%zval
  zval = psf%zval
  ps%icore = psf%icore
  call write_info_about_pseudo_1(stdout, psf, ps, z)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fix lmax and nrval...

  ps%L_max = min(psf%npotd - 1, ps%L_max)
  psf%nrval = psf%nr + mod((psf%nr + 1), 2) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads data file:
!       rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) )
!       vps(1:nrval,0:spec%ps_lmax) : pseudopotential functions
!       chcore(1:nrval) : core-correction charge distribution
!       rho_val(1:nrval) : pseudo-valence charge distribution.

  allocate(psf%rofi(psf%nrval), psf%vps(psf%nrval, 0:ps%L_max), &
       psf%chcore(1:psf%nrval), psf%rho_val(1:psf%nrval))
  call read_file_data(iunit, psf, ps)
  call io_close(iunit)

  ! get the pseudoatomic eigenfunctions
  allocate(rphi(psf%nrval,0:ps%L_max), eigen(0:ps%L_max))
  call solve_shroedinger(psf, ps, rphi, eigen)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the KB-projector cut-off radii

  allocate(rc(0:ps%L_max + 1))
  call get_cutoff_radii(psf, ps, rphi, rc)
  call write_info_about_pseudo_3(stdout, ps, rc)

  ! now we fit the splines
  call get_splines(psf, ps, rphi, rc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Increase radius a little, just in case, and final info
  ps%rc_max = ps%rc_max * 1.1_r8

  call write_info_about_pseudo_4(stdout, ps, eigen)

  deallocate(rc, rphi, eigen)
  deallocate(psf%rofi, psf%vps, psf%chcore, psf%rho_val)

! FIX UNITS (Should FIX instead kb.F90)
! Passing from Rydbergs -> Hartree
  ps%vlocal_origin = ps%vlocal_origin / 2._r8

  ps%dkbcos = ps%dkbcos / 2._r8
  ps%dknrm  = ps%dknrm  * 2._r8

  return
end subroutine ps_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_file_data(unit, psf, ps)
  integer, intent(in) :: unit
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(in) :: ps
  
  integer  :: ndown, nup, l, ir
  real(r8) :: r2

  read(unit) psf%rofi(2:psf%nrval)
  psf%rofi(1) = 0.0_r8
  do ndown = 1, ps%L_max + 1
    read(unit) l, psf%vps(2:psf%nrval, l)
    if(l /= ndown-1 .and. conf%verbose > 0) then
      message(1) = 'Unexpected angular momentum'
      message(2) = 'Pseudopotential should be ordered by increasing l'
      call write_warning(2)
    end if
    psf%vps(2:psf%nrval,l) = psf%vps(2:psf%nrval,l)/psf%rofi(2:psf%nrval)
    psf%vps(1,l) = psf%vps(2,l)
  end do
  if(ps%L_max + 2 <= psf%npotd)then
    do ndown = ps%L_max + 2, psf%npotd
      read(unit) l
    end do
  end if
  do nup = 1, psf%npotu
    read(unit) l
  end do

  ! read the core correction charge density
  r2 = psf%rofi(2)/(psf%rofi(3) - psf%rofi(2))
  read(unit) (psf%chcore(ir), ir = 2, psf%nrval)
  psf%chcore(1) = psf%chcore(2) - (psf%chcore(3) - psf%chcore(2))*r2

  ! read the pseudo-valence charge density
  read(unit) (psf%rho_val(ir), ir = 2, psf%nrval)
  psf%rho_val(1) = psf%rho_val(2) - (psf%rho_val(3) - psf%rho_val(2))*r2

  ! adjust units from Rydbergs -> Hartree
  ! psf%vps = psf%vps / 2._r8

  return
end subroutine read_file_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solve_shroedinger(psf, ps, rphi, eigen)
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(inout) :: ps
  real(r8), intent(out) :: rphi(psf%nrval,0:ps%L_max), eigen(0:ps%L_max)

  integer :: ir, l, nnode, nprin
  real(r8) :: rpb, ea, vtot, r2, e, z, dr, rmax, f, dsq, a2b4, &
       dnrm, avgv, vphi
  real(r8), allocatable :: s(:), drdi(:), ve(:), hato(:), g(:), y(:), rc(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate parameters for solving Schroedinger equation

  allocate(s(psf%nrval), drdi(psf%nrval), ve(psf%nrval), hato(psf%nrval))
  rpb = psf%b; ea = exp(psf%a)
  do ir = 1, psf%nrval
    drdi(ir) = psf%a*rpb
    s(ir) = sqrt(psf%a*rpb)
    rpb = rpb*ea
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!  ionic pseudopotential if core-correction for hartree
  
  if((ps%icore == 'pche').or.(ps%icore == 'fche')) then
    call vhrtre(psf%chcore, ve, psf%rofi, drdi, s, psf%nrval, psf%a)
!    ve = ve / 2._r8 ! Rydberg -> Hartree conversion
    do l = 0, ps%L_max
      psf%vps(2:psf%nrval, l) = psf%vps(2:psf%nrval, l) + ve(2:psf%nrval)
    end do
    psf%vps(1, l) = psf%vps(2, l)
  end if
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation of the valence screening potential from the density:
!       ve(1:nrval) is the hartree+xc potential created by the pseudo -
!               valence charge distribution.

  call calculate_valence_screening(psf, ps, drdi, s, ve)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check at which radius the asymptotic 2*Zval/r is achieved
!       rc(l) (l=0,spec%ps_lmax, l .ne. spec%ps_lloc) is the radius at which
!               the l-component differ from the one chosen as local.
!       rc(ps%L_loc) is the radius at which this potential differs from its
!               asymptotic behaviour.
!       rc(ps%rc_max) is the maximum value of all the "nonlocal" rc(l)

  allocate(rc(0:ps%L_max))
  call get_asymptotic_radii(psf, ps, rc)
  call write_info_about_pseudo_2(stdout, ps, rc)
  deallocate(rc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation of the pseudo-wave functions.
!       rhpi(1:nrval,0:spec%ps_lmax) : pseudo-wave functions.
!       eigen(0:spec%ps_lmax)        : eigenvalues.

  s(2:psf%nrval) = drdi(2:psf%nrval)*drdi(2:psf%nrval)
  s(1) = s(2)
  a2b4 = 0.25_r8*psf%a**2
  allocate(g(psf%nrval), y(psf%nrval))
  g = 0.0_r8;  y = 0.0_r8

  do l = 0, ps%L_max
    do ir = 2, psf%nrval
      vtot = psf%vps(ir, l) + ve(ir) + dble(l*(l + 1))/(psf%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    
    nnode = 1; nprin = l + 1
    e = -((psf%zval/dble(nprin))**2); z = psf%zval
    dr = -1.0e5_r8; rmax = psf%rofi(psf%nrval)
    call egofv(hato, s, psf%nrval, e, g, y, l, z, psf%a, psf%b, rmax, nprin, nnode, dr)
    eigen(l) = e

    rphi(2:psf%nrval, l) = g(2:psf%nrval) * sqrt(drdi(2:psf%nrval))
    rphi(1, l) = rphi(2, l)
  end do
  deallocate(g, y)

  !  checking normalization of the calculated wave functions
  do l = 0, ps%L_max
    e = sqrt(sum(drdi(2:psf%nrval)*rphi(2:psf%nrval, l)**2))
    e = abs(e - 1.0d0)
           
    if (e > 1.0d-5 .and. conf%verbose > 0) then
      write(message(1), '(a,i2,a)') "Eigenstate for l = ", l , ' is not normalized'
      write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      call write_warning(2)
    end if
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KB-cosines and KB-norms:
!       dkbcos(0:spec%ps_lmax) stores the KB "cosines:"
!               || (v_l - v_local) phi_l ||^2 / < (v_l - v_local)phi_l | phi_l > 
!       dknrm(0:spec%ps_lmax) stores the KB "norms:"
!               1 / || (v_l - v_local) phi_l ||

  do l = 0, ps%L_max
    if(l == ps%L_loc) then
      ps%dkbcos(l) = 0.0_r8; ps%dknrm(l) = 0.0_r8
      cycle
    end if
    dnrm = 0.0_r8
    avgv = 0.0_r8
    do ir = 2, psf%nrval
      vphi = (psf%vps(ir, l) - psf%vps(ir, ps%L_loc))*rphi(ir, l)
      dnrm = dnrm + vphi*vphi*drdi(ir)
      avgv = avgv + vphi*rphi(ir, l)*drdi(ir)
    end do
    ps%dkbcos(l) = dnrm/(avgv + 1.0e-20_r8)
    ps%dknrm(l) = 1.0_r8/(sqrt(dnrm) + 1.0e-20_r8)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ghost analysis.

  call ghost_analysis(psf, ps, ve, s, eigen)

  deallocate(s, drdi, ve, hato)
  return
end subroutine solve_shroedinger

subroutine calculate_valence_screening(psf, ps, drdi, s, ve)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(in) :: ps
  real(r8), intent(IN)  :: drdi(psf%nrval), s(psf%nrval)
  real(r8), intent(out) :: ve(psf%nrval)

  character(len=3) :: xcfunc, xcauth
  integer          :: irelt
  real(r8)         :: e_x, e_c, dx, dc, r2
  real(r8), allocatable :: v_xc(:,:), auxrho(:)

  ve = 0.0_r8
  call vhrtre(psf%rho_val, ve, psf%rofi, drdi, s, psf%nrval, psf%a)
!  ve = ve/2._r8 ! Rydberg -> Hartree

  ! Set the xc functional
  select case(psf%icorr)
  case('ca') 
     xcfunc = 'LDA'
     xcauth = 'PZ'
  case('pw') 
     xcfunc = 'LDA'
     xcauth = 'PW92'
  case('pb') 
     xcfunc = 'GGA'
     xcauth = 'PBE'
  end select
  
  if(psf%irel == 'rel') irelt=1
  if(psf%irel /= 'rel') irelt=0

  allocate(v_xc(psf%nrval, 1), auxrho(psf%nrval))
  auxrho = psf%rho_val
  if(ps%icore /= 'nc  ')  auxrho = auxrho + psf%chcore
  auxrho(2:psf%nrval) = auxrho(2:psf%nrval)/(4.0_r8*M_PI*psf%rofi(2:psf%nrval)**2)

  r2 = psf%rofi(2)/(psf%rofi(3)-psf%rofi(2))
  auxrho(1) = auxrho(2) - (auxrho(3)-auxrho(2))*r2
  call atomxc(xcfunc, xcauth, irelt, psf%nrval, psf%nrval, psf%rofi, &
       1, auxrho, e_x, e_c, dx, dc, v_xc)
  
  ve(1:psf%nrval) = ve(1:psf%nrval) + v_xc(1:psf%nrval, 1)
  deallocate(v_xc, auxrho)

  return
end subroutine calculate_valence_screening

subroutine get_asymptotic_radii(psf, ps, rc)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(IN) :: ps
  real(r8), intent(out) :: rc(0:ps%L_max)

  integer             :: l, ir
  real(r8)            :: dincv

  do l = 0, ps%L_max
    do ir = psf%nrval, 2, -1
      if(l == ps%L_loc) then
        dincv = abs(psf%vps(ir, l)*psf%rofi(ir) + 2.0_r8*psf%zval)
      else
        dincv = abs(psf%vps(ir, l) - psf%vps(ir, ps%L_loc))
      endif
      if(dincv > eps) exit
    end do
    rc(l) = psf%rofi(ir+1)
  end do

  return
end subroutine get_asymptotic_radii

subroutine get_cutoff_radii(psf, ps, rphi, rc)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  real(r8), intent(IN)  :: rphi(psf%nrval, 0:ps%L_max)
  real(r8), intent(out) :: rc(0:ps%L_max)

  integer             :: l, ir
  real(r8)            :: dincv, phi

  ! local part ....
  do ir = psf%nrval, 2, -1
    dincv = abs(psf%vps(ir, ps%L_loc)*psf%rofi(ir) + 2.0_r8*psf%zval)
    if(dincv > eps) exit
  end do
  rc(ps%L_max + 1) = psf%rofi(ir + 1)
  
  ! non-local part....
  ps%rc_max = 0.0d0
  do l = 0, ps%L_max
    if(l == ps%L_loc) then
      rc(l) = 0.0_r8
      cycle
    endif
    do ir = psf%nrval, 2, -1
      phi = (rphi(ir, l)/psf%rofi(ir))*ps%dknrm(l)
      dincv = abs(psf%vps(ir, l) - psf%vps(ir, ps%L_loc))*phi
      if(dincv > eps) exit
    enddo
    rc(l) = psf%rofi(ir + 1)
    ps%rc_max = max(ps%rc_max, rc(l))
  end do
end subroutine get_cutoff_radii

subroutine ghost_analysis(psf, ps, ve, s, eigen)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(in) :: ps
  real(r8), intent(in) :: ve(psf%nrval), s(psf%nrval), eigen(0:ps%L_max)

  integer  :: ir, l, nprin, nnode, ighost
  real(r8) :: vtot, a2b4, z, e, dr, rmax, dnrm, avgv
  real(r8), allocatable :: hato(:), g(:), y(:), elocal(:,:)

  allocate(hato(psf%nrval), g(psf%nrval), y(psf%nrval), elocal(2, 0:ps%L_max))

  ! calculate eigenvalues of the local potential for ghost analysis
  a2b4 = 0.25_r8*psf%a**2
  do l = 0, ps%L_max
    do ir = 2, psf%nrval
      vtot = psf%vps(ir, ps%L_loc) + ve(ir) + dble(l*(l+1))/(psf%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    do nnode = 1, 2
      nprin = l + 1
      e = -(psf%zval/dble(nprin))**2
      z = psf%zval
      dr = -1.0e5_r8
      rmax = psf%rofi(psf%nrval)
      call egofv(hato, s, psf%nrval, e, g, y, l, z, psf%a, psf%b, rmax, nprin, nnode, dr)
      elocal(nnode,l) = e
    end do
  end do

!  calculate KB-cosines and ghost analysis
  do l = 0, ps%L_max
    ighost = -1
    if(ps%dkbcos(l) > 0.0d0) then
      if(eigen(l) > elocal(2, l)) then
        ighost = 1
      end if
    else if(ps%dkbcos(l) < 0d0) then
      if(eigen(l) > elocal(1, l)) then
        ighost = 1
      end if
    end if
    if(ighost >= 0) then
      write(message(1), '(a,i2)') "Ghost state found for l = ", l
      call write_warning(1)
    endif
  end do

  deallocate(hato, g, y, elocal)

  return
end subroutine ghost_analysis

subroutine get_splines(psf, ps, rphi, rc)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  real(r8), intent(in) :: rphi(psf%nrval, 0:ps%L_max), rc(0:ps%L_max)
  
  integer :: l, nrc, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  allocate(hato(psf%nrval), derhato(psf%nrval))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Interpolate the KB-projection functions

  do l = 0, ps%l_max
    if(l == ps%L_loc) cycle

    hato = 0.0d0
    nrc = nint(log(rc(l)/psf%b + 1.0_r8)/psf%a) + 1
    hato(2:nrc) = (psf%vps(2:nrc, l) - psf%vps(2:nrc, ps%L_loc))*rphi(2:nrc, l) &
         *ps%dknrm(l)/psf%rofi(2:nrc)**(l+1)
    hato(1) = hato(2)    
    call spline_fit(psf%nrval, psf%rofi, hato, ps%kb(l))

! and now the derivatives...
    call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
    call spline_fit(psf%nrval, psf%rofi, derhato, ps%dkb(l))
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now the part corresponding to the local and non-local pseudopotential
! where the asymptotic part is substracted 

!...local part...
  hato = 0.0_r8
  nrc = nint(log(rc(ps%L_max + 1)/psf%b + 1.0_r8)/psf%a) + 1

  hato(2:psf%nrval) = psf%vps(2:psf%nrval, ps%L_loc)*psf%rofi(2:psf%nrval) + 2.0_r8*psf%zval
  hato(1) = 2.0_r8*psf%zval
  
  ! WARNING: Rydbergs -> Hartrees
  hato = hato / 2._r8
  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vps(1, ps%L_loc)

! and the derivative now
  call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

! and the non-local parts
  do l = 0 , ps%L_max
    if(l == ps%L_loc) cycle

    hato = 0._r8
    nrc=nint(log(rc(l)/psf%b + 1.0_r8)/psf%a) + 1
    hato(2:psf%nrval) = (psf%vps(2:psf%nrval, l) - psf%vps(2:psf%nrval, ps%L_loc))*&
         psf%rofi(2:psf%nrval)
    hato(1) = 0.0_r8
    
    ! WARNING: Rydbergs -> Hartrees
    hato = hato / 2._r8
    call spline_fit(psf%nrval, psf%rofi, hato, ps%vps(l))
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the table for the pseudo-wavefunction components (using splines)
! with a correct normalization function

  do l = 0 , ps%L_max
    nrc = nint(log(rc(l)/psf%b + 1.0_r8)/psf%a) + 1
    do ir = nrc+ 2, psf%nrval-2
      if ( abs(rphi(ir,l)/psf%rofi(ir)**(l+1)) < eps ) exit
    enddo
    nrc = ir + 1

    hato = 0.0_r8
    hato(2:nrc) = rphi(2:nrc, l)/psf%rofi(2:nrc)**(l + 1)
    hato(1) = hato(2)

    call spline_fit(psf%nrval, psf%rofi, hato, ps%Ur(l))
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  pseudo-core radius and Table with the pseudo-core data

  if(ps%icore /= 'nc  ') then
    nrcore = 0
    do ir = psf%nrval, 2, -1
      chc = psf%chcore(ir)/(4.0_r8*M_PI*(psf%rofi(ir)**2))
      if((chc > eps).and.(nrcore == 0)) then
        nrcore = ir + 1
        exit
      end if
    end do

    hato = 0.0_r8
    hato(2:nrcore) = psf%chcore(2:nrcore)/(4.0d0*M_PI*psf%rofi(2:nrcore)**2)
    hato(1) = hato(2)
    nrc = nint(log(psf%rofi(ir +1)/psf%b + 1.0_r8)/psf%a) + 1
    call spline_fit(psf%nrval, psf%rofi, hato, ps%core)
  end if

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

!---------------------------- Output ------------------------------------
subroutine write_info_about_pseudo_1(unit, psf, ps, z)
  integer, intent(in) :: unit
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(in) :: ps
  real(r8), intent(in) :: z
  
  integer :: i, j
  
  write(message(1),'(a)') ''
  write(message(2),'(a,a,a)')  '**********   Pseudopotential Information for: ', psf%namatm,'   **********'
  write(message(3),'(a,f8.3)') 'Z: ', z
  write(message(4),'(a,f7.3)')   'Valence charge: ', psf%zval
  write(message(5),'(a,i2)')   'Maximum L-component to consider: ', ps%L_max
  write(message(6),'(a,i2)')   'Maximum L-component in file: ', psf%npotd - 1
  write(message(7),'(a,i2)')   'L-component considered as local: ', ps%L_loc
  write(message(8),'(a,a2)')   'Exchange/correlation used in generation: ', psf%icorr
  call write_info(8, unit)
  
  if(conf%verbose >= 999) then
    message(1) = '             Options:'
    message(2) = '               "ca" : Ceperley-Alder (LDA)'
    message(3) = '               "wi  : Wigner (LDA)'
    message(4) = '               "hl" : Hedin-Lundqvist (LDA)'
    message(5) = '               "gl" : Gunnarson-Lundqvist (LDA)'
    message(6) = '               "pb" : Perdew, Burke and Ernzerhof (GGA)'
    message(7) = '               "pw" : PW92 (GGA)'
    call write_info(7, unit)
  endif
  write(message(1),'(a,a3)')   'Relativistic character of calculations: ', psf%irel
  call write_info(1, unit)
  if(conf%verbose >= 999) then
    message(1) = '             Options:'
    message(2) = '               "rel" : Relativistic calculations.'
    message(3) = '               "isp" : Spin polarized calculations'
    message(4) = '               other : Non relativistic, non spin-polarized.'
    call write_info(4, unit)
  endif
  write(message(1),'(a,a4)')   'Type of core corrections: ', ps%icore
  call write_info(1, unit)
  if(conf%verbose >= 999) then
    message(1) = '             Options:'
    message(2) = '               "pcec" : Partial core, xc.'
    message(3) = '               "pche" : Partial core, h+xc.'
    message(4) = '               "fcec" : Full core, xc.'
    message(5) = '               "fche" : Full core, h+xc.'
    message(6) = '               "nc"   : No core correction.'
    call write_info(6, unit)
  endif
  message(1) = 'Valence configuration in calculations:   '
  message(2) = ' (orbital - occupancy - core radius)'
  call write_info(2, unit)
  i = 1
  do while(index(psf%titleps(i:),'/') /= 0)
    j = i
    i = i + index(psf%titleps(i:),'/')
    write(message(1),'(a,a,a)') ' ', psf%titleps(j:i-2)
    call write_info(1, unit)
  enddo
  write(message(1),'(a,i2)')  'Pseudopotential functions for spin-up:   ', psf%npotu
  write(message(2),'(a,i2)')  'Pseudopotential functions for spin-down: ', psf%npotd
  write(message(3),'(a)')     'Radial grid parameters ( R(I) = B*[ EXP(A*(I-1)) -1 ] )'
  write(message(4),'(a,es14.6)') '             A = ', psf%a
  write(message(5),'(a,es14.6)') '             B = ', psf%b
  write(message(6),'(a,i5)')  'Number of radial points: ', psf%nr
  call write_info(6, unit)
  
  return
end subroutine write_info_about_pseudo_1

subroutine write_info_about_pseudo_2(unit, ps, rc) 
  integer, intent(in)  :: unit
  type(ps_type), intent(in) :: ps
  real(r8), intent(in) :: rc(0:ps%L_max)
  
  message(1) = 'Pseudopotential asymptotic radii: ['//trim(units_out%length%abbrev)//']'
  write(message(2),'(10x,10f14.6)') rc/units_out%length%factor
  call write_info(2, unit)
  
  return
end subroutine write_info_about_pseudo_2

subroutine write_info_about_pseudo_3(unit, ps, rc) 
  integer, intent(in) :: unit
  type(ps_type), intent(IN) :: ps
  real(r8), intent(IN) :: rc(0:ps%L_max+1)

  message(1) = 'KB-projectors radii: ['//trim(units_out%length%abbrev)//']'
  write(message(2),'(10x,10f14.6)') rc(0:ps%L_max)/units_out%length%factor
  call write_info(2, unit)

  return
end subroutine write_info_about_pseudo_3

subroutine write_info_about_pseudo_4(unit, ps, eigen)
  integer, intent(in) :: unit
  type(ps_type), intent(in) :: ps
  real(r8), intent(in) :: eigen(0:ps%L_max)

  message(1) = 'KB-cosines: ['//trim(units_out%energy%abbrev)//']'
  ! WARNING: RYDBERGS -> HARTREE (/2._R8)
  write(message(2),'(10x,10f18.6)')  ps%dkbcos(0:ps%L_max)/units_out%energy%factor/2._r8
  message(3) = 'KB-norms:'
  write(message(4),'(10x,10f18.6)')  ps%dknrm(0:ps%L_max)
  message(5) = 'Eigenvalues of pseudo-eigenfunctions ['// &
       trim(units_out%energy%abbrev)//']'
  ! WARNING: RYDBERGS -> HARTREE (/2._R8)
  write(message(6),'(10x,10f18.6)')  eigen(0:ps%L_max)/units_out%energy%factor/2._r8
  write(message(7),'(a,f18.6)')      'Atomic radius: ['//&
       trim(units_out%length%abbrev)//']', ps%rc_max/units_out%length%factor
  message(8) = '*************************************************************'
  message(9) = ''
  call write_info(9, unit)
  
  return
end subroutine write_info_about_pseudo_4

end module ps
