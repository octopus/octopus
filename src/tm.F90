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

module tm
use global
use kb
use io
use units

implicit none

private
public :: ps_file,                    &
     read_file_data_bin, read_file_data_ascii, solve_shroedinger, calculate_valence_screening, &
     psf_debug, ps_tm_read_file

type ps_file
  character(len=2)  :: namatm, icorr
  character(len=3)  :: irel
  character(len=4) :: icore
  character(len=10) :: method(6) 
  character(len=70) :: title
  integer :: npotd, npotu, nr
  real(r8) :: b, a, zval
  real(r8), pointer :: s(:), drdi(:)
  real(r8), pointer :: rofi(:), vps(:,:), vlocal(:), rphi(:,:), eigen(:), chcore(:), rho_val(:)
  integer :: nrval ! not in file, but very useful :)
end type ps_file

contains

subroutine ps_tm_read_file(psf, filename)
  type(ps_file), intent(inout) :: psf
  character(len=*), intent(in) :: filename

  character(len=256) :: filename2
  integer :: iunit
  logical :: found

  sub_name = 'ps_tm_read_file'; call push_sub()

  filename2 = trim(filename)//'.vps'
  inquire(file=filename2, exist=found)
  message(1) = "Info: Reading pseudopotential from file:"
  if(found) then
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call write_info(2)
    
    call io_assign(iunit)
    open(iunit, file=filename2, form='unformatted', status='old')
    call read_file_data_bin(iunit, psf)
    call io_close(iunit)
  else
    filename2 = trim(filename)//'.ascii'
    inquire(file=filename2, exist=found)
    if(.not.found) then
      filename2 = SHARE_OCTOPUS//"/PP/TM2/"//trim(filename)//".ascii"
      inquire(file=filename2, exist=found)
      if(.not.found) then
        message(1) = "Pseudopotential file '"//trim(filename)//"{.vps|.ascii}' not found"
        call write_fatal(1)
      end if
    end if

    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call write_info(2)
    
    call io_assign(iunit)
    open(iunit, file=filename2, form='formatted', status='old')
    call read_file_data_ascii(iunit, psf)
    call io_close(iunit)
  end if

  call pop_sub(); return
end subroutine ps_tm_read_file

subroutine calculate_valence_screening(psf, drdi, s, ve)
  type(ps_file), intent(in) :: psf
  real(r8), intent(IN)  :: drdi(psf%nrval), s(psf%nrval)
  real(r8), intent(out) :: ve(psf%nrval)

  character(len=5) :: xcfunc, xcauth
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
  if(psf%icore /= 'nc  ')  auxrho = auxrho + psf%chcore
  auxrho(2:psf%nrval) = auxrho(2:psf%nrval)/(4.0_r8*M_PI*psf%rofi(2:psf%nrval)**2)

  r2 = psf%rofi(2)/(psf%rofi(3)-psf%rofi(2))
  auxrho(1) = auxrho(2) - (auxrho(3)-auxrho(2))*r2
  call atomxc(xcfunc, xcauth, irelt, psf%nrval, psf%nrval, psf%rofi, &
       1, auxrho, e_x, e_c, dx, dc, v_xc)
  
  ve(1:psf%nrval) = ve(1:psf%nrval) + v_xc(1:psf%nrval, 1)
  deallocate(v_xc, auxrho)

  return
end subroutine calculate_valence_screening

subroutine solve_shroedinger(psf)
  type(ps_file), intent(inout) :: psf

  integer :: ir, l, nnode, nprin
  real(r8) :: rpb, ea, vtot, r2, e, z, dr, rmax, f, dsq, a2b4, &
       dnrm, avgv, vphi
  real(r8), allocatable :: s(:), drdi(:), ve(:), hato(:), g(:), y(:)!, rc(:)

  sub_name = 'solve_shroedinger'; call push_sub()

! calculate parameters for solving Schroedinger equation
  allocate(s(psf%nrval), drdi(psf%nrval), ve(psf%nrval), hato(psf%nrval))
  rpb = psf%b; ea = exp(psf%a)
  do ir = 1, psf%nrval
    drdi(ir) = psf%a*rpb
    s(ir) = sqrt(psf%a*rpb)
    rpb = rpb*ea
  end do

!  ionic pseudopotential if core-correction for hartree
  if((psf%icore == 'pche').or.(psf%icore == 'fche')) then
    call vhrtre(psf%chcore, ve, psf%rofi, drdi, s, psf%nrval, psf%a)
    ! ve = ve / 2._r8 ! Rydberg -> Hartree conversion
    do l = 0, psf%npotd - 1
      psf%vps(2:psf%nrval, l) = psf%vps(2:psf%nrval, l) + ve(2:psf%nrval)
    end do
    psf%vps(1, l) = psf%vps(2, l)
  end if

! Calculation of the valence screening potential from the density:
!       ve(1:nrval) is the hartree+xc potential created by the pseudo -
!               valence charge distribution (everything in Rydberts, and bohrs)
  call calculate_valence_screening(psf, drdi, s, ve)

! Calculation of the pseudo-wave functions.
!       rphi(1:nrval,0:spec%ps_lmax) : radial pseudo-wave functions. They are normalized so that
!                                      int dr {rphi^2} = 1. Thus its units are of bohr^(-1/2).
!       eigen(0:spec%ps_lmax)        : eigenvalues, in Rydbergs.
  s(2:psf%nrval) = drdi(2:psf%nrval)*drdi(2:psf%nrval)
  s(1) = s(2)
  a2b4 = 0.25_r8*psf%a**2
  allocate(g(psf%nrval), y(psf%nrval))
  g = 0.0_r8;  y = 0.0_r8

  do l = 0, psf%npotd-1
    do ir = 2, psf%nrval
      vtot = psf%vps(ir, l) + ve(ir) + dble(l*(l + 1))/(psf%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    
    nnode = 1; nprin = l + 1
    e = -((psf%zval/dble(nprin))**2); z = psf%zval
    dr = -1.0e5_r8; rmax = psf%rofi(psf%nrval)
    call egofv(hato, s, psf%nrval, e, g, y, l, z, psf%a, psf%b, rmax, nprin, nnode, dr)
    psf%eigen(l) = e

    psf%rphi(2:psf%nrval, l) = g(2:psf%nrval) * sqrt(drdi(2:psf%nrval))
    psf%rphi(1, l) = psf%rphi(2, l)
  end do
  deallocate(g, y)

  !  checking normalization of the calculated wave functions
  do l = 0, psf%npotd-1! ps%L_max
    e = sqrt(sum(drdi(2:psf%nrval)*psf%rphi(2:psf%nrval, l)**2))
    e = abs(e - 1.0d0)
    if (e > 1.0d-5 .and. conf%verbose > 0) then
      write(message(1), '(a,i2,a)') "Eigenstate for l = ", l , ' is not normalized'
      write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      call write_warning(2)
    end if
  end do

  deallocate(s, drdi, ve, hato)
  call pop_sub; return
end subroutine solve_shroedinger

subroutine read_file_data_bin(unit, psf)
  integer, intent(in) :: unit
  type(ps_file), intent(inout) :: psf
  
  integer  :: ndown, nup, l, ir
  real(r8) :: r2

  sub_name = 'read_file_data_bin'; call push_sub()

  ! Reads the header line of the file, with general info about the ps.
  read(unit) psf%namatm, psf%icorr, psf%irel, psf%icore,     &
       (psf%method(l),l=1,6), psf%title, psf%npotd, psf%npotu, &
       psf%nr, psf%b, psf%a, psf%zval
 
   ! fixes psf%nrval, to obtain an odd number.
  psf%nrval = psf%nr + mod((psf%nr + 1), 2)

  ! Allocates the varibales to psf%nrval:
  allocate(psf%rofi(psf%nrval), psf%vps(psf%nrval, 0:psf%npotd-1), &
       psf%chcore(1:psf%nrval), psf%rho_val(1:psf%nrval))

  ! Reads the radial values, in bohrs
  !   rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) ) [bohr]
  read(unit) psf%rofi(2:psf%nrval)
  psf%rofi(1) = 0.0_r8

  ! Reads the pseudoptential functions, times r, in Rydberg*bohr.
  ! Inmediately afterwards, it is divided by r, so that its final units are Rydbergs
  do ndown = 1, psf%npotd
    read(unit) l, psf%vps(2:psf%nrval, l)
    if(l /= ndown-1 .and. conf%verbose > 0) then
      message(1) = 'Unexpected angular momentum'
      message(2) = 'Pseudopotential should be ordered by increasing l'
      call write_warning(2)
    end if
    psf%vps(2:psf%nrval,l) = psf%vps(2:psf%nrval,l)/psf%rofi(2:psf%nrval)
    psf%vps(1,l) = psf%vps(2,l)
  end do

  ! Skips the spin down wavefunctions.
  do nup = 1, psf%npotu
    read(unit) l
  end do

  ! Reads the core correcction charge density, in bohr^(-3)
  !   chcore(1:nrval) : core-correction charge distribution
  r2 = psf%rofi(2)/(psf%rofi(3) - psf%rofi(2))
  read(unit) (psf%chcore(ir), ir = 2, psf%nrval)
  psf%chcore(1) = psf%chcore(2) - (psf%chcore(3) - psf%chcore(2))*r2

  ! Reads the pseudo-valence charge density, in bohr^(-3)
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  read(unit) (psf%rho_val(ir), ir = 2, psf%nrval)
  psf%rho_val(1) = psf%rho_val(2) - (psf%rho_val(3) - psf%rho_val(2))*r2

  ! adjust units from Rydbergs -> Hartree
  ! psf%vps = psf%vps / 2._r8

  call pop_sub(); return
end subroutine read_file_data_bin

subroutine read_file_data_ascii(unit, psf)
  integer, intent(in) :: unit
  type(ps_file), intent(inout) :: psf
  
  integer  :: ndown, nup, i, l, ir
  real(r8) :: r2
  real(r8), allocatable :: aux(:)
  character(len=70) :: aux_s

  sub_name = 'read_file_data_ascii'; call push_sub()

  ! Reads the header line of the file, with general info about the ps.
  read(unit, 9000) psf%namatm, psf%icorr, psf%irel, psf%icore
  read(unit, 9010) (psf%method(l),l=1,6), psf%title
  read(unit, 9015) psf%npotd, psf%npotu, psf%nr, psf%b, psf%a, psf%zval

9000 format(1x,a2,1x,a2,1x,a3,1x,a4)
9010 format(1x,6a10,/,1x,a70)
9015 format(1x,2i3,i5,3f20.10)

  ! fixes psf%nrval, to obtain an odd number.
  psf%nrval = psf%nr + mod((psf%nr + 1), 2)
  allocate(aux(psf%nr))

  ! Allocates the varibales to psf%nrval:  ! Reads the pseudo-valence charge density, in bohr^(-3)
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  allocate(psf%rofi(psf%nrval), psf%vps(psf%nrval, 0:psf%npotd-1), &
       psf%chcore(1:psf%nrval), psf%rho_val(1:psf%nrval))

  ! Reads the radial values, in bohrs
  !   rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) ) [bohr]
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  psf%rofi(2:psf%nrval) = aux(1:psf%nr)
  psf%rofi(1) = 0.0_r8

8000 format(1x,i2)
9030 format(4(g20.12))
9040 format(1x,a)
  ! Reads the pseudoptential functions, times r, in Rydberg*bohr.
  ! Inmediately afterwards, it is divided by r, so that its final units are Rydbergs
  do ndown = 1, psf%npotd
    read(unit, 9040) aux_s
    read(unit, 8000) l
    read(unit, 9030) (aux(i),i=1, psf%nr)
    if(l /= ndown-1 .and. conf%verbose > 0) then
      message(1) = 'Unexpected angular momentum'
      message(2) = 'Pseudopotential should be ordered by increasing l'
      call write_warning(2)
    end if
    psf%vps(2:psf%nrval, l) = aux(1:psf%nr)/psf%rofi(2:psf%nrval)
    psf%vps(1, l) = psf%vps(2, l)
  end do

  ! Skips the spin down wavefunctions.
  do nup = 1, psf%npotu
    read(unit, 9040) aux_s
    read(unit, 8000) l
    read(unit, 9030) (aux(i),i=1, psf%nr)
  end do

  ! Reads the core correcction charge density, in bohr^(-3)
  !   chcore(1:nrval) : core-correction charge distribution
  r2 = psf%rofi(2)/(psf%rofi(3) - psf%rofi(2))
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  psf%chcore(2:psf%nrval) = aux(1:psf%nr)
  psf%chcore(1) = psf%chcore(2) - (psf%chcore(3) - psf%chcore(2))*r2

  ! Reads the pseudo-valence charge density, in bohr^(-3)  ! Reads the pseudo-valence charge density, in bohr^(-3)
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  psf%rho_val(2:psf%nrval) = aux(1:psf%nr)
  psf%rho_val(1) = psf%rho_val(2) - (psf%rho_val(3) - psf%rho_val(2))*r2

  ! adjust units from Rydbergs -> Hartree
  ! psf%vps = psf%vps / 2._r8

  call pop_sub(); return
end subroutine read_file_data_ascii

subroutine psf_debug(psf, label)
  type(ps_file), intent(in) :: psf

  character(len=*), intent(in) :: label

  integer :: unit, i, l
  real(r8) :: dx

  call io_assign(unit)
  open(unit=unit, file=trim(label)//'.tm.debug')

  do i = 1, psf%nrval
    write(unit, *) &
         psf%rofi(i) / units_out%length%factor,  &
         (psf%vlocal(i)/2.0_r8) /units_out%energy%factor, &
         (psf%vps(i, l), l=0, psf%npotd-1)
  end do
  call io_close(unit)
 
  return
end subroutine psf_debug

end module
