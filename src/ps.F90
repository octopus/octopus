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

module ps

use io
use units
use spline
use logrid
use atomic
use tm
use hgh

implicit none

private
public :: ps_type, ps_init, ps_debug, ps_end

type ps_type
  character(len=10) :: label
  character(len=3) :: flavour

  type(spline_type), pointer :: kb(:, :)   ! Kleynman-Bylander projectors
  type(spline_type), pointer :: so_kb(:, :)
  type(spline_type), pointer :: dkb(:, :)  ! derivatives of KB projectors
  type(spline_type), pointer :: so_dkb(:, :)
  type(spline_type), pointer :: Ur(:, :)   ! atomic wavefunctions
  type(spline_type) :: vlocal  ! local part
  type(spline_type) :: dvlocal ! derivative of the local part
  type(spline_type) :: core    ! core charge

  integer  :: ispin ! Consider spin (ispin = 2) or not (ispin = 1)
  integer  :: kbc  ! Number of KB components (1 for TM ps, 3 for HGH)
  real(r8) :: z, z_val
  integer :: l_max ! maximum value of l to take
  integer :: l_loc ! which component to take as local
  integer :: l_max_occ ! maximum l-component which has non-null atomic occupation numbers
  integer :: so_l_max ! obvious meaning ;)
  real(r8), pointer :: occ(:, :)

  character(len=4) :: icore
  real(r8) :: rc_max
  real(r8) :: vlocal_origin ! local pseudopotential at the orginin

  real(r8), pointer :: dknrm(:) ! KB norm
  real(r8), pointer :: so_dknrm(:)
  real(r8), pointer :: h(:,:,:), k(:, :, :)
end type ps_type

real(r8), parameter :: eps = 1.0e-8_r8

contains

subroutine ps_init(ps, label, flavour, z, lmax, lloc, ispin)
  type(ps_type), intent(out) :: ps
  character(len=10), intent(in) :: label
  character(len=3),  intent(in) :: flavour
  integer, intent(in) :: lmax, lloc, ispin
  real(r8), intent(in) :: z

  type(tm_type)      :: pstm ! In case Troullier-Martins ps are used.
  type(hgh_type)     :: psp  ! In case Hartwigsen-Goedecker-Hutter ps are used.

  integer :: i, j, is

  sub_name = 'ps_init'; call push_sub()

  ! Makes the directory for debugging.
  if(conf%verbose>999) call oct_mkdir('pseudos')

  ! Sets the flavour, label, and number of spin channels.
  ps%flavour = flavour
  ps%label   = label
  ps%ispin   = ispin

  ! Initialization and processing.
  select case(flavour(1:2))
  case('tm')
    call tm_init(pstm, trim(label), ispin)
    ps%kbc = 1
    ps%l_max = min(pstm%npotd - 1, lmax)   ! Maybe the file has not enough components.
    ps%l_max_occ = ps%l_max
    ps%l_loc = lloc
    if(ps%l_max == 0) ps%l_loc = 0 ! Vanderbilt is not acceptable if ps%l_max == 0.
    ps%so_l_max = min(pstm%npotu, lmax)
    ps%z = z
    call tm_process(pstm, lmax, ps%l_loc)
    if(conf%verbose > 999) call tm_debug(pstm)
  case('hg')
    call hgh_init(psp, trim(label), ispin)
    ps%kbc = 3
    ps%l_max = psp%l_max
    ps%l_max_occ = max(psp%l_max, 0)
    ps%l_loc = -1
    ps%so_l_max = ps%l_max
    ps%z = z
    call hgh_process(psp)
    if(conf%verbose > 999) call hgh_debug(psp)
  case default
    message(1) = "Unknown pseudopotential type: '"+trim(flavour)+"'"
    call write_fatal(1)
  end select

! We allocate all the stuff
  allocate(ps%kb      (0:ps%l_max, ps%kbc),         &
           ps%dkb     (0:ps%l_max, ps%kbc),         &
           ps%so_kb   (0:ps%l_max, ps%kbc),         &
           ps%so_dkb  (0:ps%l_max, ps%kbc),         &
           ps%ur      (0:ps%L_max_occ, ps%ispin),   &
           ps%dknrm   (0:ps%L_max),                 &
           ps%so_dknrm(0:ps%so_l_max),                 &          
           ps%h       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc), &
           ps%k       (0:ps%l_max, 1:ps%kbc, 1:ps%kbc), &
           ps%occ     (0:ps%l_max_occ, ps%ispin))
           ps%dknrm    = M_ZERO
           ps%so_dknrm = M_ZERO
           ps%h        = M_ZERO
           ps%k        = M_ZERO
           ps%occ      = M_ZERO
  do i = 0, ps%L_max
     do j = 1, ps%kbc
        call spline_init(ps%kb(i, j))
        call spline_init(ps%dkb(i, j))
        call spline_init(ps%so_kb(i, j))
        call spline_init(ps%so_dkb(i, j))
     enddo
     do is = 1, ps%ispin
        call spline_init(ps%ur(i, is))
     enddo
  enddo
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)
  call spline_init(ps%core)

! Now we load the necessary information.
  select case(flavour(1:2))
  case('tm')
    call tm_load(ps, pstm)
    call tm_end(pstm)
  case('hg')
    call hgh_load(ps, psp)
    call hgh_end(psp)
  end select

  call pop_sub()
end subroutine ps_init

subroutine ps_debug(ps)
  type(ps_type), intent(in) :: ps

  ! I think I can hardcode these two numbers.
  integer, parameter  :: npoints = 20001
  real(r8), parameter :: grid = 0.01_r8

  character(len=4)  :: fm
  integer           :: info_unit, local_unit, nonlocal_unit, wave_unit, so_unit, &
                       i, j, k, l, is
  real(r8)          :: r

  sub_name = 'ps_debug'; call push_sub()

  ! Opens the files.
  call oct_mkdir(C_string('pseudos/'+trim(ps%label)))
  call io_assign(info_unit); call io_assign(local_unit)
  call io_assign(nonlocal_unit); call io_assign(wave_unit)
  call io_assign(so_unit)
  open(info_unit, file='pseudos/'+trim(ps%label)+'/info')
  open(local_unit, file='pseudos/'+trim(ps%label)+'/local')
  open(nonlocal_unit, file='pseudos/'+trim(ps%label)+'/nonlocal')
  open(wave_unit, file='pseudos/'+trim(ps%label)+'/wave')
  open(so_unit, file='pseudos/'+trim(ps%label)+'/so')

  ! Writes down the info.
  write(info_unit,'(a,/)')      ps%label
  write(info_unit,'(a,a,/)')    'Flavour : ',ps%flavour
  write(info_unit,'(a,f6.3)')   'z       : ',ps%z
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
  do i=1, npoints
     r = (i-1)*grid
     if(r >= r_small) then
       write(local_unit, *) r, (splint(ps%vlocal,  r) - ps%z_val)/r, &
                              (splint(ps%dvlocal, r)*r - (splint(ps%vlocal, r)-ps%z_val))/r**2
     else
       write(local_unit, *) r, ps%vlocal_origin, 0.0_r8
     end if
  enddo

  ! Kleinman-Bylander projectors
  write(fm,'(i4)') 2*ps%kbc*(ps%l_max+1) + 1; fm = adjustl(fm)
  do i =1, npoints
     r = (i-1)*grid 
     write(nonlocal_unit, '('//trim(fm)//'f16.8)') r, &
           ( (splint(ps%kb(k, j), r), j=1, ps%kbc), k=0, ps%l_max), &
           ( (splint(ps%dkb(k, j), r), j=1, ps%kbc), k=0, ps%l_max)
  enddo

  ! Spin-Orbit projectors
  write(fm,'(i4)') 2*ps%kbc*(ps%l_max+1) + 1; fm = adjustl(fm)
  do i =1, npoints
     r = (i-1)*grid 
     write(so_unit, '('//trim(fm)//'f16.8)') r, &
           ( (splint(ps%so_kb(k, j), r), j=1, ps%kbc), k=0, ps%l_max), &
           ( (splint(ps%so_dkb(k, j), r), j=1, ps%kbc), k=0, ps%l_max)
  enddo

  ! Pseudo-wavefunctions
  write(fm,'(i4)') ps%ispin*(ps%l_max+1)+1; fm = adjustl(fm)
  do i = 1, npoints
     r = (i-1)*grid
     write(wave_unit, '('//trim(fm)//'f16.8)') &
           r, ((splint(ps%ur(l, is), r), l = 0, ps%l_max_occ), is = 1, ps%ispin)
  enddo

  ! Closes files and exits
  call io_close(info_unit); call io_close(local_unit)
  call io_close(nonlocal_unit); call io_close(wave_unit)
  call io_close(so_unit)

  call pop_sub(); return
end subroutine ps_debug

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i, j, is

  sub_name = 'ps_end'; call push_sub()

  if(.not. associated(ps%kb)) return

  do i = 0, ps%L_max
    do j = 1, ps%kbc
      call spline_end(ps%kb(i, j))
      call spline_end(ps%dkb(i, j))
    enddo
    do is = 1, ps%ispin
      call spline_end(ps%Ur(i, is))
    enddo
  end do
  do i = 0, ps%so_L_max
    do j = 1, ps%kbc
      call spline_end(ps%so_kb(i, j))
      call spline_end(ps%so_dkb(i, j))
    end do
  end do

  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  call spline_end(ps%core)  

  deallocate(ps%kb, ps%dkb, ps%so_kb, ps%so_dkb, ps%ur, ps%dknrm, ps%h, ps%k, ps%occ)

  call pop_sub(); return
end subroutine ps_end

subroutine hgh_load(ps, psp)
  type(ps_type), intent(inout)      :: ps
  type(hgh_type), intent(inout)     :: psp

  integer :: l

  sub_name = 'hgh_load'; call push_sub()

  ! Fixes some components of ps, read in psf
  ps%z_val = psp%z_val
  ps%icore = 'nc'
  ps%rc_max = 1.1_r8 * maxval(psp%kbr(0:ps%l_max)) ! Increase a little.
  ps%h = psp%h
  ps%k = psp%k

  ! Fixes the occupations
  ps%occ(0:ps%l_max_occ, 1) = psp%occ(0:ps%l_max_occ)
  if(ps%ispin == 2) then
    do l = 0, psp%l_max_occ
       ps%occ(l, 1) = min(psp%occ(l), real(2*l+1,r8))
       ps%occ(l, 2) = psp%occ(l) - ps%occ(l, 1)
    enddo
  endif

  ! now we fit the splines
  call get_splines_hgh(psp, ps)
  ps%so_dknrm = ps%dknrm

  call pop_sub(); return
end subroutine hgh_load

subroutine tm_load(ps, pstm)
  type(ps_type), intent(inout) :: ps
  type(tm_type), intent(inout) :: pstm

  integer :: l

  sub_name = 'tm_load'; call push_sub()

  ! Fixes some components of ps, read in pstm
  ps%z_val = pstm%zval
  ps%icore = pstm%icore
  ps%h(0:ps%l_max, 1, 1) = pstm%dkbcos(0:ps%l_max)
  ps%dknrm (0:ps%l_max) = pstm%dknrm (0:ps%l_max)

  if(ps%so_l_max >= 0) then
    !ps%k(0:ps%so_l_max, 1, 1) = pstm%so_dkbcos(1:ps%so_l_max+1)
    !ps%so_dknrm(0:ps%so_l_max) = pstm%so_dknrm(1:ps%so_l_max+1)
    ps%k(1:ps%so_l_max, 1, 1) = pstm%so_dkbcos(1:ps%so_l_max)
    ps%so_dknrm(1:ps%so_l_max) = pstm%so_dknrm(1:ps%so_l_max)
  end if

  ! Increasing radius a little, just in case.
  ps%rc_max = maxval(pstm%kbr(0:ps%l_max)) * 1.1_r8

  ! Fixes the occupations
  ps%occ(0:ps%l_max, 1:ps%ispin) = pstm%occ(0:ps%l_max, 1:ps%ispin)

  ! now we fit the splines
  call get_splines_tm(pstm, ps)

  ! Passes from Rydbergs to Hartrees.
  ps%h = ps%h / 2._r8; ps%dknrm = ps%dknrm * 2._r8
  ps%k = ps%k / 2._r8; ps%so_dknrm = ps%so_dknrm * 2._r8

  call pop_sub(); return
end subroutine tm_load

subroutine get_splines_tm(psf, ps)
  type(tm_type), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  
  integer :: is, l, nrc, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  sub_name = 'get_splines_tm'; call push_sub()

  allocate(hato(psf%nrval), derhato(psf%nrval))

  ! Interpolate the KB-projection functions
  do l = 0, ps%l_max
    hato = 0.0_r8
    nrc = nint(log(psf%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
    hato(2:nrc) = (psf%vps(2:nrc, l) - psf%vlocal(2:nrc))*psf%rphi(2:nrc, l, 1) * & 
                   ps%dknrm(l) / psf%rofi(2:nrc)
    hato(1) = hato(2) - ((hato(3)-hato(2))/(psf%rofi(3)-psf%rofi(2)))*psf%rofi(2)    
    call spline_fit(psf%nrval, psf%rofi, hato, ps%kb(l, 1))
    call derivate_in_log_grid(psf%g, hato, derhato)
    call spline_fit(psf%nrval, psf%rofi, derhato, ps%dkb(l, 1))
  end do

  do l = 0, ps%l_max
     hato = 0.0_r8
     if(l>0 .and. psf%irel=='rel') then
       nrc = nint(log(psf%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
       hato(2:nrc) = (psf%vso(2:nrc, l))*psf%rphi(2:nrc, l, 1) * & 
                     ps%so_dknrm(l) / psf%rofi(2:nrc)
       hato(1) = hato(2) - ((hato(3)-hato(2))/(psf%rofi(3)-psf%rofi(2)))*psf%rofi(2)    
     endif
     call spline_fit(psf%nrval, psf%rofi, hato, ps%so_kb(l, 1))
     call derivate_in_log_grid(psf%g, hato, derhato)
     call spline_fit(psf%nrval, psf%rofi, derhato, ps%so_dkb(l, 1))
  end do

  ! Now the part corresponding to the local pseudopotential
  ! where the asymptotic part is substracted 
  hato = 0.0_r8
  nrc = nint(log(psf%kbr(ps%L_max + 1)/psf%b + 1.0_r8)/psf%a) + 1
  hato(2:psf%nrval) = psf%vlocal(2:psf%nrval)*psf%rofi(2:psf%nrval) + 2.0_r8*psf%zval
  hato(1) = 2.0_r8*psf%zval
  ! WARNING: Rydbergs -> Hartrees
  hato = hato / 2._r8
  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vlocal(1) / 2._r8

  ! and the derivative now
  call derivate_in_log_grid(psf%g, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

  ! Define the table for the pseudo-wavefunction components (using splines)
  ! with a correct normalization function
  do is = 1, ps%ispin
    do l = 0 , ps%L_max
      nrc = nint(log(psf%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
      do ir = nrc+ 2, psf%nrval-2
        if ( abs(psf%rphi(ir, l, 1+is)/psf%rofi(ir)) < eps ) exit
      end do
      nrc = ir + 1
      hato = 0.0_r8
      hato(2:nrc) = psf%rphi(2:nrc, l, 1 + is)/psf%rofi(2:nrc)
      hato(1) = hato(2)
      call spline_fit(psf%nrval, psf%rofi, hato, ps%Ur(l, is))
    end do
  end do

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

  deallocate(hato, derhato)

  call pop_sub(); return
end subroutine get_splines_tm

subroutine get_splines_hgh(psp, ps)
  type(hgh_type), intent(in)   :: psp
  type(ps_type), intent(inout) :: ps

  integer :: l, is, nrc, ir, nrcore, j
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  sub_name = 'get_splines_hgh'; call push_sub()

  allocate(hato(psp%g%nrval), derhato(psp%g%nrval))

  ! Interpolate the KB-projection functions
  do l = 0, psp%l_max
  do j = 1, 3
    hato = 0.0_r8
    nrc = nint(log(psp%kbr(l)/psp%g%b + 1.0_r8)/psp%g%a) + 1
    hato(1:nrc) = psp%kb(1:nrc, l, j)
    call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%kb(l, j))
    call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%so_kb(l, j))
    ! and now the derivatives...
    call derivate_in_log_grid(psp%g, hato, derhato)
    call spline_fit(psp%g%nrval, psp%g%rofi, derhato, ps%dkb(l, j))
    call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%so_dkb(l, j))
  end do
  end do

  ! Now the part corresponding to the local pseudopotential
  ! where the asymptotic part is substracted 
  hato = 0.0_r8
  nrc = nint(log(psp%kbr(ps%L_max + 1)/psp%g%b + 1.0_r8)/psp%g%a) + 1
  hato(2:psp%g%nrval) = psp%vlocal(2:psp%g%nrval)*psp%g%rofi(2:psp%g%nrval) + psp%z_val
  hato(1) = psp%z_val
  call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psp%vlocal(1)
  ! and the derivative now
  call derivate_in_log_grid(psp%g, hato, derhato)
  call spline_fit(psp%g%nrval, psp%g%rofi, derhato, ps%dvlocal)

  ! Define the table for the pseudo-wavefunction components (using splines)
  ! with a correct normalization function
  do is = 1, ps%ispin
  do l = 0 , ps%l_max_occ
    nrc = nint(log(psp%kbr(l)/psp%g%b + 1.0_r8)/psp%g%a) + 1
    do ir = nrc+ 2, psp%g%nrval-2
      if ( abs(psp%rphi(ir,l)/psp%g%rofi(ir)) < eps ) exit
    enddo
    nrc = ir + 1
    hato = 0.0_r8
    hato(2:nrc) = psp%rphi(2:nrc, l)/psp%g%rofi(2:nrc)
    hato(1) = hato(2)
    call spline_fit(psp%g%nrval, psp%g%rofi, hato, ps%Ur(l, is))
  end do
  end do

  call pop_sub(); return
end subroutine get_splines_hgh

end module ps
