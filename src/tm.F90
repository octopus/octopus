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
use atomic
use io
use units
use logrid

implicit none

private
public :: tm_type, tm_init, tm_end, tm_process, tm_debug

type tm_type
  ! First, the contents of the file.
  character(len=2)  :: namatm, icorr
  character(len=3)  :: irel
  character(len=4) :: icore
  character(len=10) :: method(6) 
  character(len=70) :: title
  integer :: npotd, npotu, nr
  real(r8) :: b, a
  real(r8) :: zval
  real(r8), pointer :: rofi(:),      &
                       vps(:, :),    &
                       vso(:, :),    &
                       vlocal(:),    &
                       rphi(:, :,:), &
                       eigen(:, :),  &
                       chcore(:),    &
                       rho_val(:),   &
                       kbr(:)

  type(valconf) :: conf

  ! Other stuff
  integer :: nrval
  type(logrid_type) :: g
  real(r8), pointer :: dkbcos(:), dknrm(:), &
                       so_dkbcos(:), so_dknrm(:)

  integer :: ispin
end type tm_type

real(r8), parameter :: eps = 1.0e-8_r8

contains

subroutine tm_init(pstm, filename, ispin)
  type(tm_type), intent(inout) :: pstm
  character(len=*), intent(in) :: filename
  integer, intent(in)          :: ispin

  character(len=256) :: filename2
  integer :: iunit, l, n
  logical :: found
  real(r8) :: x, y

  sub_name = 'ps_tm_read_file'; call push_sub()

  ! Sets the spin components
  pstm%ispin = ispin

  ! Find out where the hell the file is.
  filename2 = trim(filename)//'.vps'
  inquire(file=filename2, exist=found)
  message(1) = "Info: Reading pseudopotential from file:"
  if(found) then
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call write_info(2)
    
    call io_assign(iunit)
    open(iunit, file=filename2, form='unformatted', status='old')
    call read_file_data_bin(iunit, pstm)
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
    call read_file_data_ascii(iunit, pstm)
    call io_close(iunit)
  end if

  ! Initializes the logarithmic grid.
  call init_logrid(pstm%g, pstm%a, pstm%b, pstm%nrval)

  ! Allocates some stuff
  allocate(pstm%rphi(pstm%nrval, 0:pstm%npotd-1, 3), & 
           pstm%eigen(0:pstm%npotd-1, 3),            &
           pstm%dkbcos(0:pstm%npotd-1),              &
           pstm%dknrm(0:pstm%npotd-1),               &
           pstm%so_dkbcos(1:pstm%npotu),             &
           pstm%so_dknrm(1:pstm%npotu),              &
           pstm%kbr(0:pstm%npotd))

  ! Fills the valence configuration data.
  call build_valconf(pstm)
  do n = 1, pstm%npotd
     l = pstm%conf%l(n)
     if(ispin==2 .and. pstm%irel.ne.'isp') then
       x = pstm%conf%occ(n, 1)
       pstm%conf%occ(n,1) = min(x, real(2*l+1,r8))
       pstm%conf%occ(n,2) = x - pstm%conf%occ(n,1)
     endif
  enddo

  call pop_sub(); return
end subroutine tm_init

subroutine build_valconf(pstm)
  type(tm_type) :: pstm

  character(len=1) :: char1(6), char2
  integer :: l

  sub_name = 'build_valconf'; call push_sub()

  call valconf_null(pstm%conf)
  pstm%conf%symbol = pstm%namatm
  pstm%conf%p = pstm%npotd
  write(char2,'(i1)') pstm%conf%p

  select case(pstm%irel)
  case('nrl')
    read(pstm%title,'('+char2+'(i1,a1,f5.2,10x))') &
          (pstm%conf%n(l), char1(l), pstm%conf%occ(l,1), l = 1, pstm%npotd)
  case('isp')
    read(pstm%title,'('+char2+'(i1,a1,f4.2,1x,f4,2,6x))') &
          (pstm%conf%n(l), char1(l), pstm%conf%occ(l,1), pstm%conf%occ(l,2), l = 1, pstm%npotd)
  case('rel')
    read(pstm%title,'('+char2+'(i1,a1,f5.2,10x))') &
          (pstm%conf%n(l), char1(l), pstm%conf%occ(l,1), l = 1, pstm%npotd)
  end select

  do l = 1, pstm%npotd
     select case(char1(l))
       case('s'); pstm%conf%l(l) = 0
       case('p'); pstm%conf%l(l) = 1
       case('d'); pstm%conf%l(l) = 2
       case('f'); pstm%conf%l(l) = 3
       case default
         message(1) = 'Error reading pseudopotential file.'
         call write_fatal(1)
     end select
  enddo

  call pop_sub(); return
end subroutine build_valconf

subroutine tm_end(pstm)
  type(tm_type), intent(inout) :: pstm

  sub_name = 'tm_end'; call push_sub()

  deallocate(pstm%rofi, pstm%vps, pstm%vso, pstm%chcore, pstm%rho_val, pstm%vlocal, &
             pstm%rphi, pstm%eigen, pstm%dkbcos, pstm%dknrm, pstm%so_dkbcos, pstm%so_dknrm, pstm%kbr)
  call kill_logrid(pstm%g)

  call pop_sub(); return
end subroutine tm_end

subroutine tm_process(pstm, lmax, lloc)
  type(tm_type), intent(inout) :: pstm
  integer, intent(in) :: lmax, lloc

  sub_name = 'tm_process'; call push_sub()

! get the pseudoatomic eigenfunctions
  call solve_schroedinger(pstm)

! Fixes the local potential. Final argument is the core radius ??!!
  call get_local(pstm, lloc, 3.0_r8)

! calculates kb cosines and norms
  call calculate_kb_cosines(pstm, lloc)

! Ghost analysis.
  call ghost_analysis(pstm, lmax)

! Define the KB-projector cut-off radii
  call get_cutoff_radii(pstm, lloc)

  call pop_sub(); return
end subroutine tm_process

subroutine solve_schroedinger(psf)
  type(tm_type), intent(inout) :: psf

  character(len=3) :: functl
  integer :: iter, ir, is, l, nnode, nprin, irel
  real(r8) :: vtot, diff, r2, e, z, dr, rmax, f, dsq, a2b4, dnrm, avgv, vphi
  real(r8), allocatable :: s(:), ve(:, :), hato(:), g(:), y(:), rho(:, :), prev(:, :)
  real(r8), parameter :: tol = 1.0e-10_r8

  sub_name = 'solve_schroedinger'; call push_sub()

  ! Let us be a bit informative.
  message(1) = '      Calculating atomic pseudo-eigenfunctions for specie '//psf%namatm//'....'
  call write_info(1)

  ! Allocation.
  allocate(s   (psf%g%nrval),            &
           hato(psf%g%nrval),            &
           g   (psf%g%nrval),            &
           y   (psf%g%nrval),            &
           ve  (psf%g%nrval, psf%ispin), &
           rho (psf%g%nrval, psf%ispin), &
           prev(psf%g%nrval, psf%ispin))
           s = 0._r8; hato = 0._r8; g = 0._r8; y = 0._r8; ve = 0._r8; rho = 0._r8; prev = 0._r8

  ! These numerical parameters have to be fixed for egofv to work.  
  s(2:psf%nrval) = psf%g%drdi(2:psf%g%nrval)*psf%g%drdi(2:psf%g%nrval)
  s(1) = s(2)
  a2b4 = 0.25_r8*psf%a**2

  !  ionic pseudopotential if core-correction for hartree
  if((psf%icore == 'pche') .or. (psf%icore == 'fche')) then
    call vhrtre(psf%chcore, ve(:, 1), psf%g%rofi, psf%g%drdi, psf%g%s, psf%g%nrval, psf%g%a)
    do l = 0, psf%npotd - 1
      psf%vps(2:psf%nrval, l) = psf%vps(2:psf%nrval, l) + ve(2:psf%nrval, 1)
    end do
    psf%vps(1, l) = psf%vps(2, l)
  end if

  ! Calculation of the valence screening potential from the density:
  !       ve(1:nrval) is the hartree+xc potential created by the pseudo -
  !               valence charge distribution (everything in Rydberts, and bohrs)
  rho(:, 1) = psf%rho_val(:)
  select case(psf%icorr)
   case('pb');   functl = 'GGA'
   case default; functl = 'LDA'
  end select
  irel = 0; if(psf%irel=='rel') irel = 1
  call atomhxc(functl, irel, psf%g, 1, rho(:, 1:1), ve(:, 1:1), psf%chcore)
 
  ! Calculation of the pseudo-wave functions.
  !       rphi(1:nrval, 1:psf%conf%p, :) : radial pseudo-wave functions. They are normalized so that
  !                                      int dr {rphi^2} = 1. Thus its units are of bohr^(-1/2).
  !       eigen(1:psf%conf%p, :)        : eigenvalues, in Rydbergs.
  do l = 0, psf%npotd-1
    do ir = 2, psf%nrval
      vtot = psf%vps(ir, l) + ve(ir, 1) + dble(l*(l + 1))/(psf%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    
    nnode = 1; nprin = l + 1
    e = -((psf%zval/dble(nprin))**2); z = psf%zval
    dr = -1.0e5_r8; rmax = psf%rofi(psf%nrval)
    call egofv(hato, s, psf%nrval, e, g, y, l, z, psf%a, psf%b, rmax, nprin, nnode, dr)
    psf%eigen(l, 1) = e

    psf%rphi(2:psf%nrval, l, 1) = g(2:psf%nrval) * sqrt(psf%g%drdi(2:psf%nrval))
    psf%rphi(1, l, 1) = psf%rphi(2, l, 1)
  end do

  !  checking normalization of the calculated wave functions
  do l = 0, psf%npotd-1! ps%L_max
    e = sqrt(sum(psf%g%drdi(2:psf%nrval)*psf%rphi(2:psf%nrval, l, 1)**2))
    e = abs(e - 1.0d0)
    if (e > 1.0d-5 .and. conf%verbose > 0) then
      write(message(1), '(a,i2,a)') "Eigenstate for l = ", l , ' is not normalized'
      write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      call write_warning(2)
    end if
  end do

  ! Make a copy of the wavefunctions.
  psf%rphi(:, :, 2) = psf%rphi(:, :, 1)

  ! And now, for the spin polarized case...
  spin_polarized: if(psf%ispin == 2) then

    rho = M_ZERO    ! Here information of previous calculation could be used, but
    prev = M_ZERO   ! to save code lines, let us start from scratch.
    diff = 1.0e5_r8
    iter = 0
    self_consistent: do
      prev = rho
      iter = iter + 1

      spin: do is = 1, psf%ispin
      ang: do l = 0, psf%npotd-1
        do ir = 2, psf%nrval
          vtot = psf%vps(ir, l) + ve(ir, is) + dble(l*(l + 1))/(psf%rofi(ir)**2)
          hato(ir) = vtot*s(ir) + a2b4
        end do
        hato(1) = hato(2)
        nnode = 1; nprin = l + 1
        e = -((psf%zval/dble(nprin))**2); z = psf%zval
        dr = -1.0e5_r8; rmax = psf%rofi(psf%nrval)
        call egofv(hato, s, psf%nrval, e, g, y, l, z, psf%a, psf%b, rmax, nprin, nnode, dr)
        psf%eigen(l, 1 + is) = e

        psf%rphi(2:psf%nrval, l, 1 + is) = g(2:psf%nrval) * sqrt(psf%g%drdi(2:psf%nrval))
        psf%rphi(1, l, 1 + is) = psf%rphi(2, l, 1 + is)
      end do ang
      end do spin

      rho = M_ZERO
      do is = 1, psf%ispin
      do l = 0, psf%npotd-1
         rho(1:psf%g%nrval, is) = rho(1:psf%g%nrval, is) + &
                                  psf%conf%occ(l+1, is)*psf%rphi(1:psf%g%nrval, l, 1 + is)**2
      enddo
      enddo

      diff = M_ZERO
      do is = 1, psf%ispin
         diff = diff + sqrt( sum(psf%g%drdi(2:psf%g%nrval) * &
                               (rho(2:psf%g%nrval, is)-prev(2:psf%g%nrval, is))**2) )
      enddo

      if(diff < tol) exit self_consistent
      if(iter>1) rho = 0.5*rho + 0.5*prev
      !write(message(1),'(a,i4,a,e10.2)') '      Iter =',iter,'; Diff =',diff
      !call write_info(1)
      call atomhxc(functl, irel, psf%g, 2, rho, ve, psf%chcore)
    enddo self_consistent

  endif spin_polarized

  ! Exit this...
  message(1) = '      Done.'; call write_info(1)
  deallocate(s, ve, hato, g, y)

  call pop_sub; return
end subroutine solve_schroedinger

subroutine read_file_data_bin(unit, psf)
  integer, intent(in) :: unit
  type(tm_type), intent(inout) :: psf
  
  integer  :: ndown, nup, l, ir, i
  real(r8) :: r2
  real(r8), allocatable :: aux(:)

  sub_name = 'read_file_data_bin'; call push_sub()

  ! Reads the header line of the file, with general info about the ps.
  read(unit) psf%namatm, psf%icorr, psf%irel, psf%icore,     &
       (psf%method(l),l=1,6), psf%title, psf%npotd, psf%npotu, &
       psf%nr, psf%b, psf%a, psf%zval
 
   ! fixes psf%nrval, to obtain an odd number.
  psf%nrval = psf%nr + mod((psf%nr + 1), 2)
  allocate(aux(psf%nr))

  ! Allocates the varibales to psf%nrval:
  allocate(psf%rofi(psf%nrval), psf%vps(psf%nrval, 0:psf%npotd-1), &
           psf%chcore(1:psf%nrval), psf%rho_val(1:psf%nrval), psf%vso(psf%nrval, psf%npotu))

  ! Reads the radial values, in bohrs
  !   rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) ) [bohr]
  read(unit) (aux(i),i=1, psf%nr)
  psf%rofi(1) = M_ZERO
  do i = 2, psf%nrval
     psf%rofi(i) = aux(i-1)
  enddo

  ! Reads the pseudoptential functions, times r, in Rydberg*bohr.
  ! Inmediately afterwards, it is divided by r, so that its final units are Rydbergs
  do ndown = 1, psf%npotd
    read(unit) l, (aux(i),i=1, psf%nr)
    if(l /= ndown-1 .and. conf%verbose > 0) then
      message(1) = 'Unexpected angular momentum'
      message(2) = 'Pseudopotential should be ordered by increasing l'
      call write_warning(2)
    end if
    do i = 2, psf%nrval
       psf%vps(i, l) = aux(i-1)/psf%rofi(i)
    enddo
    psf%vps(1, l) = psf%vps(2, l)
  end do

  ! Reads --or skips-- the "down" pseudopotentials.
  if(psf%irel=='rel') then
    do nup = 1, psf%npotu
       read(unit) l, (aux(i),i=1, psf%nr)
       if(l /= nup) then
         message(1) = 'Unexpected angular momentum'
         message(2) = 'Pseudopotential should be ordered by increasing l'
         call write_warning(2)
       end if
       do i = 2, psf%nrval
          psf%vso(2:psf%nrval,l) = aux(i-1)/psf%rofi(i)
       enddo
       psf%vso(1,l) = psf%vso(2,l)
    end do       
  else
    psf%vso = M_ZERO
    do nup = 1, psf%npotu
       read(unit) l
    end do
  endif

  ! Reads the core correcction charge density, in bohr^(-3)
  !   chcore(1:nrval) : core-correction charge distribution
  r2 = psf%rofi(2)/(psf%rofi(3) - psf%rofi(2))
  read(unit) (aux(i),i=1, psf%nr)
  do i = 2, psf%nrval
     psf%chcore(i) = aux(i-1)
  enddo
  psf%chcore(1) = psf%chcore(2) - (psf%chcore(3) - psf%chcore(2))*r2

  ! Reads the pseudo-valence charge density, in bohr^(-3)
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  read(unit) (aux(i),i=1, psf%nr)
  do i = 2, psf%nrval
     psf%rho_val(i) = aux(i-1)
  enddo
  psf%rho_val(1) = psf%rho_val(2) - (psf%rho_val(3) - psf%rho_val(2))*r2

  ! adjust units from Rydbergs -> Hartree
  ! psf%vps = psf%vps / 2._r8

  deallocate(aux)
  call pop_sub(); return
end subroutine read_file_data_bin

subroutine read_file_data_ascii(unit, psf)
  integer, intent(in) :: unit
  type(tm_type), intent(inout) :: psf
  
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
           psf%chcore(1:psf%nrval), psf%rho_val(1:psf%nrval), psf%vso(psf%nrval, psf%npotu) )

  ! Reads the radial values, in bohrs
  !   rofi(1:nrval) : radial values ( rofi(i) = b*( exp(a*(i-1)) - 1 ) ) [bohr]
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  psf%rofi(1) = M_ZERO
  do i = 2, psf%nrval
     psf%rofi(i) = aux(i-1)
  enddo

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
    do i = 2, psf%nrval
       psf%vps(i, l) = aux(i-1)/psf%rofi(i)
    enddo
    psf%vps(1, l) = psf%vps(2, l)
  end do

  ! Reads --or skips-- the "down" pseudopotentials.
  if(psf%irel=='rel') then
    do nup = 1, psf%npotu
       read(unit, 9040) aux_s
       read(unit, 8000) l
       read(unit, 9030) (aux(i),i=1, psf%nr)
       if(l /= nup) then
         message(1) = 'Unexpected angular momentum'
         message(2) = 'Pseudopotential should be ordered by increasing l'
         call write_warning(2)
       end if
       do i = 2, psf%nrval
          psf%vso(i, l) = aux(i-1)/psf%rofi(i)
       enddo
       psf%vso(1, l) = psf%vso(2, l)
    end do       
  else
    psf%vso = M_ZERO
    do nup = 1, psf%npotu
       read(unit, 9040) aux_s
       read(unit, 8000) l
       read(unit, 9030) (aux(i),i=1, psf%nr)
    end do
  endif

  ! Reads the core correcction charge density, in bohr^(-3)
  !   chcore(1:nrval) : core-correction charge distribution
  r2 = psf%rofi(2)/(psf%rofi(3) - psf%rofi(2))
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  do i = 2, psf%nrval
     psf%chcore(i) = aux(i-1)
  enddo
  psf%chcore(1) = psf%chcore(2) - (psf%chcore(3) - psf%chcore(2))*r2

  ! Reads the pseudo-valence charge density, in bohr^(-3)
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  !   rho_val(1:nrval) : pseudo-valence charge distribution
  read(unit, 9040) aux_s
  read(unit, 9030) (aux(i),i=1, psf%nr)
  do i = 2, psf%nrval
     psf%rho_val(i) = aux(i-1)
  enddo
  psf%rho_val(1) = psf%rho_val(2) - (psf%rho_val(3) - psf%rho_val(2))*r2

  ! adjust units from Rydbergs -> Hartree
  ! psf%vps = psf%vps / 2._r8

  call pop_sub(); return
end subroutine read_file_data_ascii

subroutine calculate_kb_cosines(pstm, lloc)
  type(tm_type), intent(inout) :: pstm
  integer, intent(in)          :: lloc

  integer :: ir, l
  real(r8) :: dnrm, avgv, vphi

  sub_name = 'calculate_kb_cosines'; call push_sub()

! KB-cosines and KB-norms:
!       dkbcos(0:spec%ps_lmax) stores the KB "cosines:"
!               || (v_l - v_local) phi_l ||^2 / < (v_l - v_local)phi_l | phi_l >  [Rydberg]
!       dknrm(0:spec%ps_lmax) stores the KB "norms:"
!               1 / || (v_l - v_local) phi_l || [1/Rydberg]
  do l = 0, pstm%npotd-1
    if(l == lloc) then
      pstm%dkbcos(l) = 0.0_r8; pstm%dknrm(l) = 0.0_r8
      cycle
    end if
    dnrm = 0.0_r8
    avgv = 0.0_r8
    do ir = 2, pstm%g%nrval
      vphi = (pstm%vps(ir, l) - pstm%vlocal(ir))*pstm%rphi(ir, l, 1)
      dnrm = dnrm + vphi*vphi*pstm%g%drdi(ir)
      avgv = avgv + vphi*pstm%rphi(ir, l, 1)*pstm%g%drdi(ir)
    end do
    pstm%dkbcos(l) = dnrm/(avgv + 1.0e-20_r8)
    pstm%dknrm(l) = 1.0_r8/(sqrt(dnrm) + 1.0e-20_r8)
  end do

  do l = 1, pstm%npotu
    dnrm = 0.0_r8
    avgv = 0.0_r8
    do ir = 2, pstm%g%nrval
      vphi = pstm%vso(ir, l)*pstm%rphi(ir, l, 1)
      dnrm = dnrm + vphi*vphi*pstm%g%drdi(ir)
      avgv = avgv + vphi*pstm%rphi(ir, l, 1)*pstm%g%drdi(ir)
    end do
    pstm%so_dkbcos(l) = dnrm/(avgv + 1.0e-20_r8)
    pstm%so_dknrm(l) =  1.0_r8/(sqrt(dnrm) + 1.0e-20_r8)
  end do

  call pop_sub; return
end subroutine calculate_kb_cosines

subroutine ghost_analysis(pstm, lmax)
  type(tm_type), intent(in) :: pstm
  integer, intent(in)       :: lmax

  character(len=3) :: functl
  integer :: ir, l, nnode, nprin, ighost, irel
  real(r8) :: vtot, a2b4, z, e, dr, rmax, dnrm, avgv
  real(r8), allocatable :: ve(:), s(:), hato(:), g(:), y(:), elocal(:,:)

  sub_name = 'ghost_analysis'; call push_sub()

  allocate(ve(pstm%nrval), s(pstm%nrval), hato(pstm%nrval), g(pstm%nrval), y(pstm%nrval), elocal(2, 0:lmax))

  select case(pstm%icorr)
   case('pb');   functl = 'GGA'
   case default; functl = 'LDA'
  end select
  irel = 0; if(pstm%irel=='rel') irel = 1
  call atomhxc(functl, irel, pstm%g, 1, pstm%rho_val, ve(:), pstm%chcore)

  s(2:pstm%nrval) = pstm%g%drdi(2:pstm%nrval)**2
  s(1) = s(2)
  ! calculate eigenvalues of the local potential for ghost analysis
  a2b4 = 0.25_r8*pstm%a**2
  do l = 0, lmax
    do ir = 2, pstm%g%nrval
      !vtot = psf%vps(ir, ps%L_loc) + ve(ir) + dble(l*(l+1))/(psf%rofi(ir)**2)
      vtot = pstm%vlocal(ir) + ve(ir) + dble(l*(l+1))/(pstm%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    do nnode = 1, 2
      nprin = l + 1
      e = -(pstm%zval/dble(nprin))**2
      z = pstm%zval
      dr = -1.0e5_r8
      rmax = pstm%rofi(pstm%nrval)
      call egofv(hato, s, pstm%nrval, e, g, y, l, z, pstm%a, pstm%b, rmax, nprin, nnode, dr)
      elocal(nnode,l) = e
    end do
  end do

! Ghost analysis
  do l = 0, lmax
    ighost = -1
    if(pstm%dkbcos(l) > 0.0d0) then
      if(pstm%eigen(l, 1) > elocal(2, l)) then
        ighost = 1
      end if
    else if(pstm%dkbcos(l) < 0d0) then
      if(pstm%eigen(l, 1) > elocal(1, l)) then
        ighost = 1
      end if
    end if
    if(ighost >= 0) then
      write(message(1), '(a,i2)') "Ghost state found for l = ", l
      call write_warning(1)
    endif
  end do

  deallocate(hato, g, y, elocal, s, ve)

  call pop_sub; return
end subroutine ghost_analysis

subroutine get_cutoff_radii(pstm, lloc)
  type(tm_type), intent(inout) :: pstm
  integer, intent(in)          :: lloc

  integer             :: l, ir
  real(r8)            :: dincv, phi
  real(r8), parameter :: threshold = 1.0e-6_r8

  sub_name = 'get_cutoff_radii'; call push_sub()

  ! local part ....
  do ir = pstm%g%nrval, 2, -1
    dincv = abs(pstm%vlocal(ir)*pstm%rofi(ir) + 2.0_r8*pstm%zval)
    if(dincv > threshold) exit
  enddo
  pstm%kbr(pstm%npotd) = pstm%rofi(ir + 1)
  
  ! non-local part....
  do l = 0, pstm%npotd-1
    if(l == lloc) then
      pstm%kbr(l) = 0.0_r8
      cycle
    endif
    do ir = pstm%g%nrval, 2, -1
      phi = (pstm%rphi(ir, l, 1)/pstm%rofi(ir))*pstm%dknrm(l)
      dincv = abs((pstm%vps(ir, l) - pstm%vlocal(ir))*phi)
      if(dincv > threshold) exit
      phi = (pstm%rphi(ir, l, 1)/pstm%rofi(ir))*pstm%dknrm(l)
      if(pstm%irel=='rel' .and. l>0 .and. l<=pstm%npotu) then
        dincv = abs((pstm%vso(ir, l))*phi)
        if(dincv > threshold) exit
      endif
    enddo
    pstm%kbr(l) = pstm%rofi(ir + 1)
  end do

  call pop_sub(); return
end subroutine get_cutoff_radii

subroutine get_local(psf, l_loc, rcore)
  type(tm_type), intent(inout) :: psf
  integer, intent(in)          :: l_loc
  real(r8), intent(in)         :: rcore

  integer :: ir
  real(r8) :: a, b, qtot
  real(r8), allocatable :: rho(:)

  sub_name = 'get_local'; call push_sub()

  allocate(psf%vlocal(psf%nrval))
  if(l_loc >= 0) then
    write(message(1), '(a,i2,a)') "Info: l = ", l_loc, " component used as local potential"
    call write_info(1)

    psf%vlocal(1:psf%nrval) = psf%vps(1:psf%nrval, l_loc)
  else if(l_loc == -1) then
    message(1) = "Info: Vanderbilt function local potential"
    call write_info(1)

    a = 1.82_r8 / rcore
    b = 1.0_r8
    allocate(rho(psf%nrval))

    do ir = 1, psf%nrval
      rho(ir) = exp( -( sinh(a*b*psf%rofi(ir)) / sinh(b) )**2 )
      rho(ir) = 4.0_r8 * M_Pi * rho(ir) * psf%rofi(ir)**2
    end do
    qtot = sum(rho(2:psf%nrval)*psf%g%drdi(2:psf%nrval))
    rho(:) = rho(:)*(psf%zval/qtot)

    call vhrtre(-rho, psf%vlocal, psf%rofi, psf%g%drdi, psf%g%s, psf%g%nrval, psf%g%a)
    psf%vlocal(1) = psf%vlocal(2)

    deallocate(rho)
  endif

  call pop_sub()
end subroutine get_local

subroutine tm_debug(pstm)
  type(tm_type), intent(in) :: pstm

  integer :: loc_unit, kbp_unit, dat_unit, wav_unit, so_unit, i, l, is

  sub_name = 'tm_debug'; call push_sub()

  ! Opens files.
  call oct_mkdir(C_string('pseudos/'+'tm2.'+trim(pstm%namatm)))
  call io_assign(loc_unit); call io_assign(wav_unit)
  call io_assign(dat_unit); call io_assign(kbp_unit)
  call io_assign(so_unit) ; call io_assign(so_unit)
  open(loc_unit, file = 'pseudos/'+'tm2.'+trim(pstm%namatm)+'/'+'local')
  open(dat_unit, file = 'pseudos/'+'tm2.'+trim(pstm%namatm)+'/'+'info')
  open(kbp_unit, file = 'pseudos/'+'tm2.'+trim(pstm%namatm)+'/'+'nonlocal')
  open(wav_unit, file = 'pseudos/'+'tm2.'+trim(pstm%namatm)+'/'+'wave')
  open(so_unit,  file = 'pseudos/'+'tm2.'+trim(pstm%namatm)+'/'+'so')

  ! First of all, writes down the info.
  write(dat_unit,'(a,/)') pstm%namatm
  write(dat_unit,'(a)')   'Valence configuration and core radii:'
  write(dat_unit,'(a,/)')   pstm%title
  write(dat_unit,'(a)')   'Descriptive string of the generation process:'
  write(dat_unit,'(6a,/)')   (pstm%method(i), i = 1, 6)
  write(dat_unit,'(a)')   'Relativistic character: '+pstm%irel
  write(dat_unit,'(a)')   'XC:                     '+pstm%icorr
  write(dat_unit,'(a)')   'Core correction:        '+pstm%icore
  write(dat_unit,'(/,a,i4)') 'Maximum L: ', pstm%npotd-1
  write(dat_unit,'(/,a)')   'Occupations:'
  do is = 1, pstm%ispin
     write(dat_unit,'(4x,4f14.6)') pstm%conf%occ(1:pstm%conf%p, is)
  enddo
  write(dat_unit,'(/,a)') 'Eigenvalues: '
  write(dat_unit,'(4x,4f14.6)') pstm%eigen(0:pstm%npotd-1, 1)/2
  if(pstm%ispin==2) then
    write(dat_unit,'(a)') 'Eigenvalues in spin-polarized calculations:'
    write(dat_unit,'(4x,4f14.6)') pstm%eigen(0:pstm%npotd-1, 2)/2
    write(dat_unit,'(4x,4f14.6)') pstm%eigen(0:pstm%npotd-1, 3)/2
  endif
  write(dat_unit,'(/,a)') 'Cutoff radii: '
  write(dat_unit,'(4x,5f14.6)') pstm%kbr(0:pstm%npotd)
  write(dat_unit,'(/,a)') 'KB norms: '
  write(dat_unit,'(4x,4e14.4)') pstm%dknrm(0:pstm%npotd-1)*2
  write(dat_unit,'(/,a)') 'KB cosines: '
  write(dat_unit,'(4x,4e14.6)') pstm%dkbcos(0:pstm%npotd-1)/2
  if(pstm%irel=='rel') then
    write(dat_unit,'(/,a)') 'SO-KB norms: '
    write(dat_unit,'(4x,4e14.4)') pstm%so_dknrm(1:pstm%npotu)*2
    write(dat_unit,'(/,a)') 'SO-KB cosines: '
    write(dat_unit,'(4x,4e14.6)') pstm%so_dkbcos(1:pstm%npotu)/2
  endif

  ! Now the local part
  do i = 1, pstm%g%nrval
     write(loc_unit,*) pstm%rofi(i), pstm%vlocal(i)/2
  enddo

  ! Nonlocal components.
  do i = 1, pstm%g%nrval
     write(kbp_unit,'(5es14.5)') pstm%rofi(i), (pstm%vps(i, l)/2, l = 0, pstm%npotd-1)
  enddo

  if(pstm%irel=='rel') then
    do i = 1, pstm%g%nrval
       write(so_unit,'(5es14.5)') pstm%rofi(i), (pstm%vso(i, l)/2, l = 1, pstm%npotu)
    enddo
  endif
  
  ! Pseudo wave-functions
  do i = 1, pstm%nrval
     write(wav_unit,'(5es14.5)') pstm%rofi(i), ((pstm%rphi(i, l, is), l = 0, pstm%npotd-1), is = 1, pstm%ispin)
  enddo

  ! Closes files and exits
  call io_close(loc_unit); call io_close(wav_unit)
  call io_close(dat_unit); call io_close(kbp_unit)
  call io_close(so_unit)

  call pop_sub(); return
end subroutine tm_debug

end module
