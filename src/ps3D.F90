
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

use tm
use hgh

implicit none

private
public :: ps_type, ps_init, ps_debug, ps_end

type ps_type
  character(len=10) :: label
  character(len=3) :: flavour

  type(spline_type), pointer :: kb(:, :)   ! Kleynman-Bylander projectors
  type(spline_type), pointer :: dkb(:, :)  ! derivatives of KB projectors
  type(spline_type), pointer :: Ur(:)   ! atomic wavefunctions
  type(spline_type) :: vlocal  ! local part
  type(spline_type) :: dvlocal ! derivative of the local part
  type(spline_type) :: core    ! core charge

  integer  :: kbc  ! Number of KB components (1 for TM ps, 3 for HGH)
  real(r8) :: z, z_val
  integer :: L_max ! maximum value of l to take
  integer :: L_loc ! which component to take as local
  character(len=4) :: icore
  real(r8) :: rc_max
  real(r8) :: vlocal_origin ! local pseudopotential at the orginin

  real(r8), pointer :: dknrm(:) ! KB norm
  real(r8), pointer :: h(:,:,:)
end type ps_type

real(r8), parameter :: eps = 1.0e-8_r8
contains

subroutine ps_init(ps, label, flavour, z, lmax, lloc)
  type(ps_type), intent(out) :: ps
  character(len=10), intent(in) :: label
  character(len=3),  intent(in) :: flavour
  integer, intent(in) :: lmax, lloc
  real(r8), intent(in) :: z

  type(tm_type)      :: pstm ! In case Troullier-Martins ps are used.
  type(ps_st_params) :: psp  ! In case Hartwigsen-Goedecker-Hutter ps are used.

  integer :: i, j

  sub_name = 'ps_init'; call push_sub()

! Sets the flavour and label.
  ps%flavour = flavour
  ps%label   = label

! First of all we read the input files
  select case(flavour(1:2))
  case('tm')
    call tm_init(pstm, trim(label))
    ps%kbc = 1
    ps%L_max = min(pstm%npotd - 1, lmax)   ! Maybe the file has not enough components.
    ps%l_loc = lloc
    ps%z = z
    call tm_process(pstm, lmax, lloc)
  case('hg')
    call ps_ghg_read_file(psp, trim(label))
    ps%kbc = 3
    ps%l_max = psp%l_max
    ps%l_loc = -1
    ps%z = z
  case default
    message(1) = "Unknown pseudopotential type: '"+trim(flavour)+"'"
    call write_fatal(1)
  end select

! We allocate all the stuff
  allocate(ps%kb(0:ps%l_max, ps%kbc), ps%dkb(0:ps%l_max, ps%kbc), ps%Ur(0:ps%L_max))
  do i = 0, ps%L_max
     do j = 1, ps%kbc
        call spline_init(ps%kb(i, j))
        call spline_init(ps%dkb(i, j))
     enddo
     call spline_init(ps%ur(i))
  enddo
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)
  call spline_init(ps%core)
  allocate(ps%dknrm(0:ps%L_max), ps%h(0:ps%l_max, 1:ps%kbc, 1:ps%kbc))

! Now we load the necessary information.
  select case(flavour(1:2))
  case('tm')
    call ps_tm_load(ps, pstm)
    call tm_end(pstm)
  case('hg')
    call ps_hgh_load(ps, psp, trim(label))
  end select

  call pop_sub()
end subroutine ps_init

subroutine ps_debug(ps)
  type(ps_type), intent(in) :: ps

  ! I think I can hardcode these two numbers.
  integer, parameter  :: npoints = 20001
  real(r8), parameter :: grid = 0.01_r8

  character(len=4)  :: fm
  integer           :: info_unit, local_unit, nonlocal_unit, wave_unit, &
                       i, j, k, l
  real(r8)          :: r

  sub_name = 'ps_debug'; call push_sub()

  ! Opens the files.
  call oct_mkdir(C_string(trim(ps%label)))
  call io_assign(info_unit); call io_assign(local_unit)
  call io_assign(nonlocal_unit); call io_assign(wave_unit)
  open(info_unit, file=trim(ps%label)+'/info')
  open(local_unit, file=trim(ps%label)+'/local')
  open(nonlocal_unit, file=trim(ps%label)+'/nonlocal')
  open(wave_unit, file=trim(ps%label)+'/wave')

  ! Writes down the info.
  write(info_unit,'(a,/)')      ps%label
  write(info_unit,'(a,a,/)')    'Flavour : ',ps%flavour
  write(info_unit,'(a,f6.3)')   'z       : ',ps%z
  write(info_unit,'(a,f6.3,/)') 'zval    : ', ps%z_val
  write(info_unit,'(a,i4)')     'lmax    : ', ps%l_max
  write(info_unit,'(a,i4)')     'lloc    : ', ps%l_loc
  write(info_unit,'(a,i4,/)')   'kbc     : ', ps%kbc
  write(info_unit,'(a,f9.5,/)') 'rcmax   : ', ps%rc_max
  write(info_unit,'(/,a,/)')    'h matrix:'
  do l = 0, ps%l_max
     do k = 1, ps%kbc
        write(info_unit,'(3f9.5)') (ps%h(l, k, j), j = 1, ps%kbc)
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
     write(nonlocal_unit, '('//trim(fm)//'e14.6)') r, &
           ( (splint(ps%kb(k, j), r), splint(ps%dkb(k, j), r), k=0, ps%l_max), j=1, ps%kbc)
  enddo

  ! Pseudo-wavefunctions
  write(fm,'(i4)') ps%l_max+2; fm = adjustl(fm)
  do i = 1, npoints
     r = (i-1)*grid
     write(wave_unit, '('//trim(fm)//'e14.6)') &
           r, (splint(ps%ur(l), r), l = 0, ps%l_max)
  enddo

  ! Closes files and exits
  call io_close(info_unit); call io_close(local_unit)
  call io_close(nonlocal_unit); call io_close(wave_unit)

  call pop_sub(); return
end subroutine ps_debug

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i, j

  sub_name = 'ps_end'; call push_sub()

  if(.not. associated(ps%kb)) return

  do i = 0, ps%L_max
     do j = 1, ps%kbc
        call spline_end(ps%kb(i, j))
        call spline_end(ps%dkb(i, j))
     enddo
     call spline_end(ps%Ur(i))
  end do
  
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  call spline_end(ps%core)  

  deallocate(ps%kb, ps%dkb, ps%ur, ps%dknrm, ps%h)

  call pop_sub()
end subroutine ps_end

subroutine ps_hgh_load(ps, psp, filename)
  type(ps_type), intent(inout)      :: ps
  type(ps_st_params), intent(inout) :: psp
  character(len=*), intent(in)      :: filename

  integer  :: iunit, i, ir, l
  logical  :: found
  real(r8) :: ea, rpb

  sub_name = 'ps_hgh_load'; call push_sub()

! Calculate logarithmic grid parameters.
  psp%a = 1.25e-2_r8; psp%b = 4.0e-4_r8
  psp%nrval = 1001
  allocate(psp%s(psp%nrval), psp%drdi(psp%nrval), psp%rofi(psp%nrval))
  rpb = psp%b; ea = exp(psp%a)
  do ir = 1, psp%nrval
    psp%drdi(ir) = psp%a*rpb
    psp%s(ir) = sqrt(psp%a*rpb)
    rpb = rpb*ea
    psp%rofi(ir) = psp%b * ( exp( psp%a * (ir - 1) ) - 1.0_r8 ) 
  end do

! Allocates psp variables
  allocate(psp%vlocal(1:psp%nrval), psp%kb(1:psp%nrval, 0:ps%l_max, 1:3))
  psp%vlocal = 0.0_r8; psp%kb = 0.0_r8

! Fixes some components of ps, read in psf
  ps%z_val = psp%z_val
  ps%icore = 'nc'

! get the pseudoatomic eigenfunctions (WARNING: This is not correctly done yet: "some" wavefunctions
! are obtained, but not the real ones!!!
  allocate(psp%rphi(psp%nrval, 0:ps%l_max), psp%eigen(0:ps%l_max))
  call solve_schroedinger_hgh(psp, ps%l_max)

! Fixes the local potential
  psp%vlocal(1:psp%nrval) = vlocalr(psp%rofi, psp)

! And the projectors
  do l = 0, ps%l_max
     do i = 1, 3
        psp%kb(1:psp%nrval, l, i) = projectorr(psp%rofi, psp, i, l)
     enddo
  enddo

! Define the KB-projector cut-off radii
  call get_cutoff_radii_psp(psp, ps%l_max, ps%l_loc)
       ps%rc_max = maxval(psp%kbr(0:ps%l_max))

! now we fit the splines
  call get_splines_hgh(psp, ps)

! Increase radius a little, just in case.
  ps%rc_max = ps%rc_max * 1.1_r8

! Defines the constant matrix.
  ps%h = psp%h

! Deallocation
  deallocate(psp%rofi, psp%drdi, psp%s, psp%vlocal, psp%kb)

  call pop_sub(); return
end subroutine ps_hgh_load

subroutine ps_tm_load(ps, pstm)
  type(ps_type), intent(inout) :: ps
  type(tm_type), intent(inout) :: pstm

  sub_name = 'ps_tm_load'; call push_sub()

  ! Fixes some components of ps, read in pstm
  ps%z_val = pstm%zval
  ps%icore = pstm%icore
  ps%h(0:ps%l_max, 1, 1) = pstm%dkbcos(0:ps%l_max)
  ps%dknrm (0:ps%l_max) = pstm%dknrm (0:ps%l_max)
  ps%rc_max = maxval(pstm%kbr(0:ps%l_max)) * 1.1_r8
    ! Increasing radius a little, just in case.

  ! now we fit the splines
  call get_splines_tm(pstm, ps)

  ! Passes from Rydbergs to Hartrees.
  ps%h = ps%h / 2._r8; ps%dknrm = ps%dknrm * 2._r8

  call pop_sub(); return
end subroutine ps_tm_load

subroutine get_cutoff_radii_psp(psp, lmax, lloc)
  type(ps_st_params), intent(inout) :: psp
  integer, intent(in)               :: lmax, lloc

  integer  :: ir, l, i
  real(r8) :: dincv, tmp

  sub_name = 'get_cutoff_radii_psp'; call push_sub()

  ! local part ....
  psp%kbr(0:lmax + 1) = 0.0_r8
  do ir = psp%nrval, 2, -1
    dincv = abs(psp%vlocal(ir)*psp%rofi(ir) + psp%z_val)
    if(dincv > eps) exit
  end do
  psp%kbr(lmax + 1) = psp%rofi(ir + 1)

  ! non-local part....
  do l = 0, lmax
    tmp = 0.0_r8
    do i = 1, 3
       do ir = psp%nrval, 2, -1
          dincv = abs(psp%kb(ir, l, i))
          if(dincv > eps) exit
       enddo
       tmp = psp%rofi(ir + 1)
       psp%kbr(l) = max(tmp, psp%kbr(l))
    enddo
  end do

  call pop_sub(); return
end subroutine get_cutoff_radii_psp

subroutine get_splines_tm(psf, ps)
  type(tm_type), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  
  integer :: l, nrc, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  sub_name = 'get_splines_tm'; call push_sub()

  allocate(hato(psf%nrval), derhato(psf%nrval))

! Interpolate the KB-projection functions
  do l = 0, ps%l_max
    !if(l == ps%L_loc) cycle

    hato = 0.0d0
    nrc = nint(log(psf%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
    hato(2:nrc) = (psf%vps(2:nrc, l) - psf%vlocal(2:nrc))*psf%rphi(2:nrc, l) * & 
                   ps%dknrm(l) / psf%rofi(2:nrc)
    !hato(1) = hato(2)
    hato(1) = hato(2) - ((hato(3)-hato(2))/(psf%rofi(3)-psf%rofi(2)))*psf%rofi(2)    
    call spline_fit(psf%nrval, psf%rofi, hato, ps%kb(l, 1))

! and now the derivatives...
    call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
    call spline_fit(psf%nrval, psf%rofi, derhato, ps%dkb(l, 1))
  end do

! Now the part corresponding to the local pseudopotential
! where the asymptotic part is substracted 
!...local part...
  hato = 0.0_r8
  nrc = nint(log(psf%kbr(ps%L_max + 1)/psf%b + 1.0_r8)/psf%a) + 1

  hato(2:psf%nrval) = psf%vlocal(2:psf%nrval)*psf%rofi(2:psf%nrval) + 2.0_r8*psf%zval
  hato(1) = 2.0_r8*psf%zval
  
  ! WARNING: Rydbergs -> Hartrees
  hato = hato / 2._r8
  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vlocal(1) / 2._r8

  ! and the derivative now
  call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

! Define the table for the pseudo-wavefunction components (using splines)
! with a correct normalization function
  do l = 0 , ps%L_max
    nrc = nint(log(psf%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
    do ir = nrc+ 2, psf%nrval-2
      if ( abs(psf%rphi(ir,l)/psf%rofi(ir)**(l+1)) < eps ) exit
    enddo
    nrc = ir + 1

    hato = 0.0_r8
    hato(2:nrc) = psf%rphi(2:nrc, l)/psf%rofi(2:nrc)**(l + 1)
    hato(1) = hato(2)

    call spline_fit(psf%nrval, psf%rofi, hato, ps%Ur(l))
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
  type(ps_st_params), intent(in) :: psp
  type(ps_type), intent(inout) :: ps

  integer :: l, nrc, ir, nrcore, j
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  sub_name = 'get_splines_hgh'; call push_sub()

  allocate(hato(psp%nrval), derhato(psp%nrval))

! Interpolate the KB-projection functions
  do l = 0, ps%l_max
  do j = 1, 3
    hato = 0.0_r8
    nrc = nint(log(psp%kbr(l)/psp%b + 1.0_r8)/psp%a) + 1
    hato(1:nrc) = psp%kb(1:nrc, l, j)
    call spline_fit(psp%nrval, psp%rofi, hato, ps%kb(l, j))
    ! and now the derivatives...
    call derivate_in_log_grid(psp%a, psp%b, psp%nrval, hato, derhato)
    call spline_fit(psp%nrval, psp%rofi, derhato, ps%dkb(l, j))
  end do
  end do

! Now the part corresponding to the local pseudopotential
! where the asymptotic part is substracted 
!...local part...
  hato = 0.0_r8
  nrc = nint(log(psp%kbr(ps%L_max + 1)/psp%b + 1.0_r8)/psp%a) + 1

  hato(2:psp%nrval) = psp%vlocal(2:psp%nrval)*psp%rofi(2:psp%nrval) + psp%z_val
  hato(1) = psp%z_val
  
  call spline_fit(psp%nrval, psp%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psp%vlocal(1)

  ! and the derivative now
  call derivate_in_log_grid(psp%a, psp%b, psp%nrval, hato, derhato)
  call spline_fit(psp%nrval, psp%rofi, derhato, ps%dvlocal)

! Define the table for the pseudo-wavefunction components (using splines)
! with a correct normalization function
  do l = 0 , ps%L_max
    nrc = nint(log(psp%kbr(l)/psp%b + 1.0_r8)/psp%a) + 1
    do ir = nrc+ 2, psp%nrval-2
      if ( abs(psp%rphi(ir,l)/psp%rofi(ir)**(l+1)) < eps ) exit
    enddo
    nrc = ir + 1

    hato = 0.0_r8
    hato(2:nrc) = psp%rphi(2:nrc, l)/psp%rofi(2:nrc)**(l + 1)
    hato(1) = hato(2)

    call spline_fit(psp%nrval, psp%rofi, hato, ps%Ur(l))
  end do

  call pop_sub(); return
end subroutine get_splines_hgh
