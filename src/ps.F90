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

#include "global.h"

module ps
use global
use lib_oct
use io
use lib_oct_gsl_spline
use logrid
use atomic
use tm
use hgh

implicit none

private
public :: ps_type, ps_init, ps_filter, ps_debug, ps_end

integer, parameter, public :: &
   PS_TM2 = 100, &
   PS_HGH = 101

character(len=3), parameter  :: ps_name(PS_TM2:PS_HGH) = (/"tm2", "hgh"/)

type ps_type
  character(len=10) :: label
  integer           :: flavour

  type(loct_spline_type), pointer :: kb(:, :)   ! Kleynman-Bylander projectors
  type(loct_spline_type), pointer :: so_kb(:, :)
  type(loct_spline_type), pointer :: dkb(:, :)  ! derivatives of KB projectors
  type(loct_spline_type), pointer :: so_dkb(:, :)
  type(loct_spline_type), pointer :: Ur(:, :)   ! atomic wavefunctions
  type(loct_spline_type) :: vlocal    ! local part
  type(loct_spline_type) :: vl        ! local part (once again)
  type(loct_spline_type) :: vlocal_f  ! localized part of local potential 
                                      ! in Fourier space (for periodic)
  type(loct_spline_type) :: dvlocal   ! derivative of the local part
  type(loct_spline_type) :: core      ! core charge
  type(logrid_type) :: g

  integer  :: ispin    ! Consider spin (ispin = 2) or not (ispin = 1)
  integer  :: kbc      ! Number of KB components (1 for TM ps, 3 for HGH)
  FLOAT    :: z, z_val
  integer  :: l_max    ! maximum value of l to take
  integer  :: l_loc    ! which component to take as local
  integer  :: so_l_max ! obvious meaning ;)

  type(valconf) :: conf

  character(len=4) :: icore
  FLOAT :: rc_max
  FLOAT :: vlocal_origin ! local pseudopotential at the orginin

  FLOAT, pointer :: dknrm(:) ! KB norm
  FLOAT, pointer :: so_dknrm(:)
  FLOAT, pointer :: h(:,:,:), k(:, :, :)
end type ps_type

FLOAT, parameter :: eps = CNST(1.0e-8)

contains

subroutine ps_init(ps, label, flavour, z, lmax, lloc, ispin)
  type(ps_type),     intent(out) :: ps
  character(len=10), intent(in)  :: label
  integer         ,  intent(in)  :: flavour
  integer,           intent(in)  :: lmax, lloc, ispin
  FLOAT,             intent(in)  :: z

  type(tm_type)      :: pstm ! In case Troullier-Martins ps are used.
  type(hgh_type)     :: psp  ! In case Hartwigsen-Goedecker-Hutter ps are used.

  integer :: i, j, is

  call push_sub('ps_init')

  ! Makes the directory for debugging.
  if(conf%verbose>999) call loct_mkdir('pseudos')

  ! Sets the flavour, label, and number of spin channels.
  ps%flavour = flavour
  ps%label   = label
  ps%ispin   = ispin

  ! Initialization and processing.
  ASSERT(flavour==PS_TM2.or.flavour==PS_HGH)

  select case(flavour)
  case(PS_TM2)
    call tm_init(pstm, trim(label), ispin)
    ps%l_max = min(pstm%npotd - 1, lmax)   ! Maybe the file has not enough components.
    call valconf_copy(ps%conf, pstm%conf)
    ps%conf%z = z
    ps%conf%p = ps%l_max + 1
    ps%kbc = 1
    ps%l_loc = lloc
    if(ps%l_max == 0) ps%l_loc = 0 ! Vanderbilt is not acceptable if ps%l_max == 0.
    if(pstm%npotu==0) then
      ps%so_l_max = -1
    else
      ps%so_l_max = min(pstm%npotu, ps%l_max)
    endif
    ps%z = z
    call tm_process(pstm, lmax, ps%l_loc)
    ps%g%a = pstm%g%a
    ps%g%b = pstm%g%b
    ps%g%nrval = pstm%g%nrval
    allocate(ps%g%rofi(ps%g%nrval),ps%g%drdi(ps%g%nrval))
    ps%g%rofi = pstm%g%rofi
    ps%g%drdi = pstm%g%drdi
    if(conf%verbose > 999) call tm_debug(pstm)

  case(PS_HGH)
    call hgh_init(psp, trim(label), ispin)
    call valconf_copy(ps%conf, psp%conf)
    ps%kbc = 3
    ps%l_max = psp%l_max
    ps%l_loc = -1
    ps%so_l_max = ps%l_max
    ps%z = z
    call hgh_process(psp)
    ps%g%a = psp%g%a
    ps%g%b = psp%g%b
    ps%g%nrval = psp%g%nrval
    allocate(ps%g%rofi(ps%g%nrval),ps%g%drdi(ps%g%nrval))
    ps%g%rofi = psp%g%rofi
    ps%g%drdi = psp%g%drdi
    if(conf%verbose > 999) call hgh_debug(psp)

  end select

! We allocate all the stuff
  allocate(ps%kb      (0:ps%l_max, ps%kbc),         &
           ps%dkb     (0:ps%l_max, ps%kbc),         &
           ps%ur      (ps%conf%p, ps%ispin),        &
           ps%dknrm   (0:ps%L_max),                 &
           ps%h       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc), &
           ps%k       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))!, &
  
  ps%dknrm   (0:ps%L_max) = M_ZERO
  ps%h       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc) = M_ZERO
  ps%k       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc) = M_ZERO

  if(ps%so_l_max >= 0) then
    allocate(ps%so_kb (0:ps%l_max, ps%kbc), &
             ps%so_dkb(0:ps%l_max, ps%kbc), &
             ps%so_dknrm(0:ps%l_max))
    ps%so_dknrm(0:ps%l_max) = M_ZERO
  endif

  call loct_spline_init(ps%kb)
  call loct_spline_init(ps%dkb)
  if(ps%so_l_max >= 0) call loct_spline_init(ps%so_kb)
  if(ps%so_l_max >= 0) call loct_spline_init(ps%so_dkb)
  call loct_spline_init(ps%vlocal)
  call loct_spline_init(ps%vl)
  call loct_spline_init(ps%dvlocal)
  call loct_spline_init(ps%core)
  if (conf%periodic_dim > 0) then
    call loct_spline_init(ps%vlocal_f)
  end if 

! Now we load the necessary information.
  select case(flavour)
  case(PS_TM2)
    call tm_load(ps, pstm)
    call tm_end(pstm)
  case(PS_HGH)
    call hgh_load(ps, psp)
    call hgh_end(psp)
  end select

  call pop_sub()
end subroutine ps_init

subroutine ps_filter(ps, gmax)
  type(ps_type), intent(inout) :: ps
  FLOAT, intent(in) :: gmax
  FLOAT, parameter :: alpha = CNST(2.0)
  integer :: i, l, k
  type(loct_spline_type) :: vloc
  FLOAT :: r
  FLOAT, allocatable :: y(:)

  call push_sub('ps_filter')

  ! Filter the local part.
  call loct_spline_init(vloc)
  allocate(y(ps%g%nrval))
  y(1) = loct_splint(ps%vl, CNST(0.0)) + ps%z_val*(M_TWO/sqrt(M_PI))*alpha
  do i = 2, ps%g%nrval
     r = ps%g%rofi(i)
     y(i) = loct_splint(ps%vl, r) + ps%z_val*loct_erf(alpha*r)/r
  enddo
  call loct_spline_fit(ps%g%nrval, ps%g%rofi, y, vloc)
  call loct_spline_filter(vloc, fs = (/ CNST(0.75)*gmax, CNST(18.0) /))
  call loct_spline_end(ps%vl)
  y(1) = loct_splint(vloc, CNST(0.0)) - ps%z_val*(M_TWO/sqrt(M_PI))*alpha
  do i = 2, ps%g%nrval
     r = ps%g%rofi(i)
     y(i) = loct_splint(vloc, r) - ps%z_val*loct_erf(alpha*r)/r
  enddo

  call loct_spline_fit(ps%g%nrval, ps%g%rofi, y, ps%vl)

  do l = 0, ps%l_max
     do k = 1, ps%kbc
        call loct_spline_filter(ps%kb(l, k), l, fs = (/ CNST(0.75)*gmax, CNST(18.0) /), &
                                                rs = (/ CNST(2.5), CNST(0.4) /))
     enddo
  enddo

  call loct_spline_end(vloc)
  deallocate(y)
  call pop_sub(); return
end subroutine ps_filter

subroutine ps_debug(ps)
  type(ps_type), intent(IN) :: ps

  ! We will plot also some Fourier transforms.
  type(loct_spline_type), allocatable :: fw(:, :)
  FLOAT, parameter :: gmax = CNST(40.0)

  character(len=30) :: dir
  integer  :: info_unit                            ! A text file with some basic data.
  integer  :: local_unit, dlocal_unit, localw_unit ! The local part, derivative, and FT.
  integer  :: nl_unit, dnl_unit, nlw_unit          ! Nonlocal part
  integer  :: wave_unit                            ! pseudowavefunctions
  integer  :: so_unit, dso_unit, sow_unit          ! The spin-orbit non-local terms.
  integer  :: i, j, k, l, is

  call push_sub('ps_debug')

  ! Opens the files.
  dir = 'pseudos/'//trim(ps%label)
  call loct_mkdir(trim(dir))

  call io_assign(local_unit); call io_assign(dlocal_unit); call io_assign(localw_unit)
  call io_assign(nl_unit)   ; call io_assign(dnl_unit)   ; call io_assign(nlw_unit)
  call io_assign(so_unit)   ; call io_assign(dso_unit)   ; call io_assign(sow_unit)
  call io_assign(info_unit) ;  call io_assign(wave_unit)
  open(    info_unit, file=trim(dir)//'/info')
  open(   local_unit, file=trim(dir)//'/local')
  open(  dlocal_unit, file=trim(dir)//'/local_derivative')
  open(  localw_unit, file=trim(dir)//'/local_ft')
  open(      nl_unit, file=trim(dir)//'/nonlocal')
  open(     dnl_unit, file=trim(dir)//'/nonlocal_derivative')
  open(     nlw_unit, file=trim(dir)//'/nonlocal_ft')
  open(      so_unit, file=trim(dir)//'/so')
  open(     dso_unit, file=trim(dir)//'/so_derivative')
  open(     sow_unit, file=trim(dir)//'/so_ft')
  open(    wave_unit, file=trim(dir)//'/wavefunctions')

  ! Writes down the info.
  write(info_unit,'(a,/)')      ps%label
  write(info_unit,'(a,a,/)')    'Flavour : ', ps_name(ps%flavour)
  write(info_unit,'(a,f6.3)')   'z       : ', ps%z
  write(info_unit,'(a,f6.3,/)') 'zval    : ', ps%z_val
  write(info_unit,'(a,i4)')     'lmax    : ', ps%l_max
  write(info_unit,'(a,i4)')     'lloc    : ', ps%l_loc
  write(info_unit,'(a,i4)')     'so_lmax : ', ps%so_l_max
  write(info_unit,'(a,i4,/)')   'kbc     : ', ps%kbc
  write(info_unit,'(a,f9.5,/)') 'rcmax   : ', ps%rc_max
  write(info_unit,'(/,a,/)')    'h matrix:'
  do l = 0, ps%l_max
     do k = 1, ps%kbc
        write(info_unit,'(3f9.5)') (ps%h(l, k, j), j = 1, ps%kbc)
     enddo
     write(info_unit, '(a)')
  enddo
  write(info_unit,'(/,a,/)')    'k matrix:'
  do l = 0, ps%l_max
     do k = 1, ps%kbc
        write(info_unit,'(3f9.5)') (ps%k(l, k, j), j = 1, ps%kbc)
     enddo
     write(info_unit, '(a)')
  enddo

  ! Local part.
  call loct_spline_print(ps%vlocal, local_unit)
  call loct_spline_print(ps%dvlocal, dlocal_unit)
  allocate(fw(1, 1))
  call loct_spline_init(fw(1, 1))
  call loct_spline_3dft(ps%vlocal, fw(1, 1), gmax = gmax)
  call loct_spline_print(fw(1, 1), localw_unit)
  call loct_spline_end(fw(1, 1))
  deallocate(fw)

  ! Kleinman-Bylander projectors
  call loct_spline_print(ps%kb, nl_unit)
  call loct_spline_print(ps%dkb, dnl_unit)
  allocate(fw(0:ps%l_max, 1:ps%kbc))
  call loct_spline_init(fw)
  do k = 0, ps%l_max
     do j = 1, ps%kbc
        call loct_spline_3dft(ps%kb(k, j), fw(k, j), gmax = gmax)
     enddo
  enddo
  call loct_spline_print(fw, nlw_unit)
  call loct_spline_end(fw)
  deallocate(fw)

  ! Spin-Orbit projectors
  if(ps%so_l_max >= 0) then
    call loct_spline_print(ps%so_kb, so_unit)
    call loct_spline_print(ps%so_dkb, dso_unit)
    allocate(fw(0:ps%so_l_max, 1:ps%kbc))
    call loct_spline_init(fw)
    do k = 0, ps%so_l_max
       do j = 1, ps%kbc
          call loct_spline_3dft(ps%so_kb(k, j), fw(k, j), gmax = gmax)
       enddo
    enddo
    call loct_spline_print(fw, sow_unit)
    call loct_spline_end(fw)
    deallocate(fw)
  endif

  ! Pseudo-wavefunctions
  call loct_spline_print(ps%ur, wave_unit)

  ! Closes files and exits
  call io_close(local_unit); call io_close(dlocal_unit); call io_close(localw_unit)
  call io_close(nl_unit)   ; call io_close(dnl_unit)   ; call io_close(nlw_unit)
  call io_close(so_unit)   ; call io_close(dso_unit)   ; call io_close(sow_unit)
  call io_close(info_unit) ;  call io_close(wave_unit)
  call pop_sub(); return
end subroutine ps_debug

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i, j, is

  call push_sub('ps_end')

  if(.not. associated(ps%kb)) return

  call loct_spline_end(ps%kb)
  call loct_spline_end(ps%dkb)
  call loct_spline_end(ps%ur)
  if(ps%so_l_max>=0) then
     call loct_spline_end(ps%so_kb)
     call loct_spline_end(ps%so_dkb)
  endif

  call loct_spline_end(ps%vlocal)
  call loct_spline_end(ps%dvlocal)
  call loct_spline_end(ps%vl)
  call loct_spline_end(ps%core)  
  if (conf%periodic_dim > 0) then
    call loct_spline_end(ps%vlocal_f)
  end if 

  deallocate(ps%kb, ps%dkb, ps%ur, ps%dknrm, ps%h, ps%k)
  if(ps%so_l_max >=0) deallocate(ps%so_kb, ps%so_dkb)

  call pop_sub()
end subroutine ps_end

subroutine hgh_load(ps, psp)
  type(ps_type),  intent(inout) :: ps
  type(hgh_type), intent(inout) :: psp

  integer :: l, ll
  FLOAT :: x

  call push_sub('hgh_load')

  ! Fixes some components of ps, read in psf
  ps%z_val = psp%z_val
  ps%icore = 'nc'
  if(ps%l_max>=0) then
     ps%rc_max = CNST(1.1) * maxval(psp%kbr(0:ps%l_max)) ! Increase a little.
  else
     ps%rc_max = M_ZERO
  endif
  ps%h = psp%h
  ps%k = psp%k

  ! Fixes the occupations
  if(ps%ispin == 2) then
    do l = 1, ps%conf%p
       ll = ps%conf%l(l)
       x = ps%conf%occ(l, 1)
       ps%conf%occ(l, 1) = min(x, real(2*ll+1, PRECISION))
       ps%conf%occ(l, 2) = x - ps%conf%occ(l, 1)
    enddo
  endif

  ! now we fit the splines
  call get_splines_hgh(psp, ps)
  ps%so_dknrm = ps%dknrm

  call pop_sub()
end subroutine hgh_load

subroutine tm_load(ps, pstm)
  type(ps_type), intent(inout) :: ps
  type(tm_type), intent(inout) :: pstm

  call push_sub('tm_load')

  ! Fixes some components of ps, read in pstm
  ps%z_val = pstm%zval
  ps%icore = pstm%icore
  ps%h(0:ps%l_max, 1, 1) = pstm%dkbcos(0:ps%l_max)
  ps%dknrm (0:ps%l_max) = pstm%dknrm (0:ps%l_max)

  if(ps%so_l_max >= 0) then
    ps%k(1:ps%so_l_max, 1, 1) = pstm%so_dkbcos(1:ps%so_l_max)
    ps%so_dknrm(1:ps%so_l_max) = pstm%so_dknrm(1:ps%so_l_max)
  end if

  ! Increasing radius a little, just in case.
  ! I have hard-coded a larger increase of the cutoff for the filtering.
  ps%rc_max = maxval(pstm%kbr(0:ps%l_max)) * CNST(1.5)

  ! now we fit the splines
  call get_splines_tm(pstm, ps)
  if (conf%periodic_dim > 0) then
    call get_splines_tm_fourier(pstm, ps)
  end if

  ! Passes from Rydbergs to Hartrees.
  ps%h(0:ps%l_max,:,:)    = ps%h(0:ps%l_max,:,:)    / M_TWO
  ps%dknrm(0:ps%L_max)    = ps%dknrm(0:ps%L_max)    * M_TWO
  if(ps%so_l_max >= 0) then
    ps%k(0:ps%l_max,:,:)    = ps%k(0:ps%l_max,:,:)    / M_TWO
    ps%so_dknrm(0:ps%L_max) = ps%so_dknrm(0:ps%L_max) * M_TWO
  end if

  call pop_sub()
end subroutine tm_load

subroutine get_splines_tm(psf, ps)
  type(tm_type), intent(IN)    :: psf
  type(ps_type), intent(inout) :: ps
  
  integer :: is, l, ll, nrc, ir, nrcore
  FLOAT :: chc
  FLOAT, allocatable :: hato(:), derhato(:)

  call push_sub('get_splines_tm')

  allocate(hato(psf%nrval), derhato(psf%nrval))

  ! Interpolate the KB-projection functions
  do l = 0, ps%l_max
    hato = M_ZERO
    nrc = nint(log(psf%kbr(l)/psf%b + M_ONE)/psf%a) + 1
    hato(2:nrc) = (psf%vps(2:nrc, l) - psf%vlocal(2:nrc))*psf%rphi(2:nrc, l, 1) * & 
                   ps%dknrm(l) / psf%rofi(2:nrc)
    hato(1) = hato(2) - ((hato(3)-hato(2))/(psf%rofi(3)-psf%rofi(2)))*psf%rofi(2)    
    call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%kb(l, 1))
    call derivate_in_log_grid(psf%g, hato, derhato)
    call loct_spline_fit(psf%nrval, psf%rofi, derhato, ps%dkb(l, 1))
  end do

  if(ps%so_l_max>=0) then
   do l = 0, ps%l_max
      hato = M_ZERO
      if(l>0 .and. psf%irel=='rel') then
        nrc = psf%g%nrval
        hato(2:nrc) = (psf%vso(2:nrc, l))*psf%rphi(2:nrc, l, 1) * & 
                      ps%so_dknrm(l) / psf%rofi(2:nrc)
        hato(1) = hato(2) - ((hato(3)-hato(2))/(psf%rofi(3)-psf%rofi(2)))*psf%rofi(2)    
      endif
      call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%so_kb(l, 1))
      call derivate_in_log_grid(psf%g, hato, derhato)
      call loct_spline_fit(psf%nrval, psf%rofi, derhato, ps%so_dkb(l, 1))
   end do
  endif

  ! Now the part corresponding to the local pseudopotential
  ! where the asymptotic part is substracted 
  hato = M_ZERO
  nrc = nint(log(psf%kbr(ps%L_max + 1)/psf%b + M_ONE)/psf%a) + 1
  hato(2:psf%nrval) = psf%vlocal(2:psf%nrval)*psf%rofi(2:psf%nrval) + M_TWO*psf%zval
  hato(1) = M_TWO*psf%zval
  ! WARNING: Rydbergs -> Hartrees
  hato = hato / M_TWO
  call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vlocal(1) / M_TWO

  ! and the derivative now
  call derivate_in_log_grid(psf%g, hato, derhato)
  call loct_spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

  hato = psf%vlocal/M_TWO
  call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%vl)

  ! Define the table for the pseudo-wavefunction components (using splines)
  ! with a correct normalization function
  do is = 1, ps%ispin
    do l = 1, ps%conf%p
      ll = ps%conf%l(l)
      nrc = nint(log(psf%kbr(ll)/psf%b + M_ONE)/psf%a) + 1
      do ir = nrc+ 2, psf%nrval-2
        if ( abs(psf%rphi(ir, l-1, 1+is)/psf%rofi(ir)) < eps ) exit
      end do
      nrc = ir + 1
      hato = M_ZERO
      hato(2:nrc) = psf%rphi(2:nrc, l-1, 1 + is)/psf%rofi(2:nrc)
      hato(1) = hato(2)
      call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%ur(l, is))
    end do
  end do

  !  pseudo-core radius and Table with the pseudo-core data
  if(ps%icore /= 'nc  ') then
    nrcore = 0
    do ir = psf%nrval, 2, -1
      chc = psf%chcore(ir)/(M_FOUR*M_PI*(psf%rofi(ir)**2))
      if((chc > eps).and.(nrcore == 0)) then
        nrcore = ir + 1
        exit
      end if
    end do
    hato = M_ZERO
    hato(2:nrcore) = psf%chcore(2:nrcore)/(M_FOUR*M_PI*psf%rofi(2:nrcore)**2)
    hato(1) = hato(2)
    nrc = nint(log(psf%rofi(ir +1)/psf%b + M_ONE)/psf%a) + 1
    call loct_spline_fit(psf%nrval, psf%rofi, hato, ps%core)
  end if

  deallocate(hato, derhato)

  call pop_sub()
end subroutine get_splines_tm

subroutine get_splines_tm_fourier(psf,ps)
  type(tm_type), intent(in) :: psf
  type(ps_type), intent(inout) :: ps

! the G spline is temporary calculated on psf%nrval values
! with a |G| mesh of 0.1 
! This should be thought more accurately.
  real, parameter :: h = CNST(0.1)   ! mesh for the |G| values
  real, parameter :: a_erf = M_TWO   ! see epot
  FLOAT :: modg(psf%nrval), intg(psf%nrval)

  integer :: i, j
  FLOAT   :: d

  call push_sub('get_splines_tm_fourier')

  intg = M_ZERO
  do i = 1, psf%nrval
    modg(i) = (i-1)*h 
    do j = 2, psf%nrval
      d = (loct_splint(ps%vlocal, ps%g%rofi(j)) -    &
           ps%Z_val*(M_ONE - loct_erf(a_erf*ps%g%rofi(j)))) * ps%g%drdi(j)

      if(i == 1) then
        d = d * ps%g%rofi(j)
      else
        d = d * sin(modg(i)*ps%g%rofi(j)) / modg(i)
      end if
      intg(i) = intg(i) + d
    end do
  end do

  call loct_spline_fit(psf%nrval, modg, intg, ps%vlocal_f)

end subroutine get_splines_tm_fourier

subroutine get_splines_hgh(psp, ps)
  type(hgh_type), intent(IN)    :: psp
  type(ps_type),  intent(inout) :: ps

  integer :: l, is, nrc, j
  FLOAT, allocatable :: hato(:), derhato(:)

  call push_sub('get_splines_hgh')

  allocate(hato(psp%g%nrval), derhato(psp%g%nrval))

  ! Interpolate the KB-projection functions
  do l = 0, psp%l_max
    do j = 1, 3
      hato = M_ZERO
      nrc = nint(log(psp%kbr(l)/psp%g%b + M_ONE)/psp%g%a) + 1
      hato(1:nrc) = psp%kb(1:nrc, l, j)
      call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%kb(l, j))
      call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%so_kb(l, j))
      ! and now the derivatives...
      call derivate_in_log_grid(psp%g, hato, derhato)
      call loct_spline_fit(psp%g%nrval, psp%g%rofi, derhato, ps%dkb(l, j))
      call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%so_dkb(l, j))
    end do
  end do

  ! Now the part corresponding to the local pseudopotential
  ! where the asymptotic part is substracted 
  hato(2:psp%g%nrval) = psp%vlocal(2:psp%g%nrval)*psp%g%rofi(2:psp%g%nrval) + psp%z_val
  hato(1) = psp%z_val
  call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psp%vlocal(1)
  ! and the derivative now
  call derivate_in_log_grid(psp%g, hato, derhato)
  call loct_spline_fit(psp%g%nrval, psp%g%rofi, derhato, ps%dvlocal)

  ! Define the table for the pseudo-wavefunction components (using splines)
  ! with a correct normalization function
  do is = 1, ps%ispin
  do l = 1, ps%conf%p
    hato(2:psp%g%nrval) = psp%rphi(2:psp%g%nrval, l)/psp%g%rofi(2:psp%g%nrval)
    hato(1) = hato(2)
    call loct_spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%Ur(l, is))
  end do
  end do

  call pop_sub()
end subroutine get_splines_hgh

end module ps
