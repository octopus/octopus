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

use io
use units
use spline
use kb
use tm
use hgh

implicit none

private
public :: ps_type, ps_init, ps_end

type ps_type
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
  real(r8) :: rc(0:3)
  real(r8), pointer :: kbr(:), eigen(:)
  real(r8) :: rc_max
  real(r8) :: vlocal_origin ! local pseudopotential at the orginin

  real(r8), pointer :: dkbcos(:), dknrm(:) ! KB cosinus and norm
  real(r8), pointer :: h(:,:,:)
end type ps_type

real(r8), parameter :: eps = 1.0e-8_r8
contains

subroutine ps_init(ps, label, flavour, z, lmax, lloc, debug)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: label, flavour
  integer, intent(in) :: lmax, lloc
  real(r8), intent(in) :: z
  logical, intent(in) :: debug

  type(ps_file)      :: psf ! In case Troullier-Martins ps are used.
  type(ps_st_params) :: psp ! In case Hartwigsen-Goedecker-Hutter ps are used.

  integer :: i, j

  sub_name = 'ps_init'; call push_sub()

! Sets the flavour
  ps%flavour = trim(flavour)

! First of all we read the input files
  select case(flavour(1:2))
  case('tm')
    call ps_tm_read_file(psf, trim(label))
    ps%kbc = 1
    ps%L_max = min(psf%npotd - 1, lmax)   ! Maybe the file has not enough components.
    ps%l_loc = lloc
  case('hg')
    call ps_ghg_read_file(psp, trim(label))
    ps%kbc = 3
    ps%l_max = psp%l_max
    ps%l_loc = -1
  case default
    message(1) = "Unknown pseudopotential type: '"+trim(flavour)+"'"
    call write_fatal(1)
  end select

! Fixes the nuclear charge (which is probably useless)
  ps%z = z

! We allocate all the stuff
  allocate(ps%kb(0:ps%l_max, ps%kbc), ps%dkb(0:ps%l_max, ps%kbc))
  allocate(ps%Ur(0:ps%L_max))
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

  allocate(ps%dkbcos(0:ps%L_max), ps%dknrm(0:ps%L_max), ps%kbr(0:ps%l_max+1))
  allocate(ps%h(0:ps%l_max, 1:ps%kbc, 1:ps%kbc))

! Now we load the necessary information.
  select case(flavour(1:2))
  case('tm')
    call ps_tm_load(ps, psf, trim(label), debug)
  case('hg')
    call ps_hgh_load(ps, psp, trim(label), debug)
  end select

  call pop_sub()
end subroutine ps_init

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

  deallocate(ps%kb, ps%dkb, ps%Ur)
  
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  call spline_end(ps%core)  

  deallocate(ps%dkbcos, ps%dknrm, ps%kbr)
  deallocate(ps%h)

  call pop_sub()
end subroutine ps_end

subroutine ps_hgh_load(ps, psp, filename, debug)
  type(ps_type), intent(inout)      :: ps
  type(ps_st_params), intent(inout) :: psp
  character(len=*), intent(in)      :: filename
  logical, intent(in)               :: debug

  integer  :: iunit, i, ir, l
  logical  :: found
  real(r8) :: ea, rpb

  sub_name = 'ps_hgh_load'; call push_sub()

! Reads the cut-off radii used to generate the ps, in bohrs.
  ps%rc = 0.0_r8
  do i = 0, ps%l_max
     ps%rc(i) = psp%rc(i)
  enddo

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
  allocate(ps%eigen(0:ps%l_max))
  ps%eigen(0:ps%l_max) = psp%eigen(0:ps%l_max)

! Fixes the local potential
  psp%vlocal(1:psp%nrval) = vlocalr(psp%rofi, psp)

! And the projectors
  do l = 0, ps%l_max
     do i = 1, 3
        psp%kb(1:psp%nrval, l, i) = projectorr(psp%rofi, psp, i, l)
     enddo
  enddo

! Define the KB-projector cut-off radii
  call get_cutoff_radii_psp(psp, ps)

! now we fit the splines
  call get_splines_hgh(psp, ps)

! Increase radius a little, just in case, and final info
  ps%rc_max = ps%rc_max * 1.1_r8

! Defines the constant matrix.
  ps%h = psp%h

! And outputs some nice info about it...
  call write_pseudo_info_hgh(stdout, psp, ps)

! Debugging?
  if(debug) call hgh_debug(psp, filename)

! Deallocation
  deallocate(psp%rofi, psp%drdi, psp%s, psp%vlocal, psp%kb)

  call pop_sub(); return
end subroutine ps_hgh_load

subroutine ps_tm_load(ps, psf, filename, debug)
  type(ps_type), intent(inout) :: ps
  type(ps_file), intent(inout) :: psf
  character(len=*), intent(in) :: filename
  logical, intent(in)          :: debug

  integer :: i, ir, j, l, iunit
  real(r8) :: ea, rpb

  sub_name = 'ps_tm_load'; call push_sub()

! Reads the cut-off radii used to generate the ps, in bohrs.
  ps%rc = 0.0_r8; i = 1; l = 0
  do while(index(psf%title(i:),'/') /= 0)
    j = i
    i = i + index(psf%title(i:),'/')
    read(psf%title(j+12:i-2),*) ps%rc(l)
    l = l + 1
  enddo

! calculate logarithmic grid parameters.
  allocate(psf%s(psf%nrval), psf%drdi(psf%nrval))
  rpb = psf%b; ea = exp(psf%a)
  do ir = 1, psf%nrval
    psf%drdi(ir) = psf%a*rpb
    psf%s(ir) = sqrt(psf%a*rpb)
    rpb = rpb*ea
  end do

! Fixes some components of ps, read in psf
  ps%z_val = psf%zval
  ps%icore = psf%icore

! get the pseudoatomic eigenfunctions
  allocate(psf%rphi(psf%nrval, 0:psf%npotd-1), psf%eigen(0:psf%npotd-1))
  call solve_shroedinger(psf)
  allocate(ps%eigen(0:ps%l_max))
  ps%eigen(0:ps%l_max) = psf%eigen(0:ps%l_max)

! Fixes the local potential
  call get_local(psf, ps%l_loc, maxval(ps%rc) )

! calculates kb cosines and norms
  call calculate_kb_cosines(psf, ps)

! Ghost analysis.
  call ghost_analysis(psf, ps)

! Define the KB-projector cut-off radii
  call get_cutoff_radii_tm(psf, ps)

! now we fit the splines
  call get_splines_tm(psf, ps)

! Increase radius a little, just in case, and final info
  ps%rc_max = ps%rc_max * 1.1_r8

! FIX UNITS (Should FIX instead kb.F90)
! Passing from Rydbergs -> Hartree
  ps%vlocal_origin = ps%vlocal_origin / 2._r8
  ps%eigen = ps%eigen / 2._r8
  ps%dkbcos = ps%dkbcos / 2._r8
  ps%dknrm  = ps%dknrm  * 2._r8

! Defines the constant matrix.
  ps%h(0:ps%l_max, 1, 1) = ps%dkbcos(0:ps%l_max)

! And outputs some nice info about it...
  call write_pseudo_info_tm(stdout, psf, ps)

! Debugging?
  if(debug) call psf_debug(psf, filename)

! Deallocate stuff
  deallocate(ps%eigen)
  deallocate(psf%rofi, psf%vps, psf%chcore, psf%rho_val)

  call pop_sub(); return
end subroutine ps_tm_load

subroutine calculate_kb_cosines(psf, ps)
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(inout) :: ps

  integer :: ir, l
  real(r8) :: dnrm, avgv, vphi

  sub_name = 'calculate_kb_cosines'; call push_sub()

! KB-cosines and KB-norms:
!       dkbcos(0:spec%ps_lmax) stores the KB "cosines:"
!               || (v_l - v_local) phi_l ||^2 / < (v_l - v_local)phi_l | phi_l >  [Rydberg]
!       dknrm(0:spec%ps_lmax) stores the KB "norms:"
!               1 / || (v_l - v_local) phi_l || [1/Rydberg]
  do l = 0, ps%L_max
    if(l == ps%L_loc) then
      ps%dkbcos(l) = 0.0_r8; ps%dknrm(l) = 0.0_r8
      cycle
    end if
    dnrm = 0.0_r8
    avgv = 0.0_r8
    do ir = 2, psf%nrval
      vphi = (psf%vps(ir, l) - psf%vlocal(ir))*psf%rphi(ir, l)
      dnrm = dnrm + vphi*vphi*psf%drdi(ir)
      avgv = avgv + vphi*psf%rphi(ir, l)*psf%drdi(ir)
    end do
    ps%dkbcos(l) = dnrm/(avgv + 1.0e-20_r8)
    ps%dknrm(l) = 1.0_r8/(sqrt(dnrm) + 1.0e-20_r8)
  end do

  call pop_sub; return
end subroutine calculate_kb_cosines

subroutine ghost_analysis(psf, ps)
  type(ps_file), intent(inout) :: psf
  type(ps_type), intent(inout) :: ps

  integer :: ir, l, nnode, nprin, ighost
  real(r8) :: vtot, a2b4, z, e, dr, rmax, dnrm, avgv
  real(r8), allocatable :: ve(:), s(:), hato(:), g(:), y(:), elocal(:,:)

  sub_name = 'ghost_analysis'; call push_sub()

  allocate(ve(psf%nrval), s(psf%nrval), hato(psf%nrval), g(psf%nrval), y(psf%nrval), elocal(2, 0:ps%L_max))
! Calculation of the valence screening potential from the density:
!       ve(1:nrval) is the hartree+xc potential created by the pseudo -
!               valence charge distribution (everything in Rydberts, and bohrs)
  call calculate_valence_screening(psf, psf%drdi, psf%s, ve)

  s(2:psf%nrval) = psf%drdi(2:psf%nrval)**2
  s(1) = s(2)
  ! calculate eigenvalues of the local potential for ghost analysis
  a2b4 = 0.25_r8*psf%a**2
  do l = 0, ps%L_max
    do ir = 2, psf%nrval
      !vtot = psf%vps(ir, ps%L_loc) + ve(ir) + dble(l*(l+1))/(psf%rofi(ir)**2)
      vtot = psf%vlocal(ir) + ve(ir) + dble(l*(l+1))/(psf%rofi(ir)**2)
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

! Ghost analysis
  do l = 0, ps%L_max
    ighost = -1
    if(ps%dkbcos(l) > 0.0d0) then
      if(psf%eigen(l) > elocal(2, l)) then
        ighost = 1
      end if
    else if(ps%dkbcos(l) < 0d0) then
      if(psf%eigen(l) > elocal(1, l)) then
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

subroutine get_local(psf, l_loc, rcore)
  type(ps_file), intent(inout) :: psf
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
    qtot = sum(rho(2:psf%nrval)*psf%drdi(2:psf%nrval))
    rho(:) = rho(:)*(psf%zval/qtot)

    call vhrtre(-rho, psf%vlocal, psf%rofi, psf%drdi, psf%s, psf%nrval, psf%a)
    psf%vlocal(1) = psf%vlocal(2)

    deallocate(rho)
  endif

  call pop_sub()
end subroutine get_local

subroutine get_cutoff_radii_psp(psp, ps)
  type(ps_st_params), intent(in) :: psp
  type(ps_type), intent(inout)   :: ps

  integer  :: ir, l, i
  real(r8) :: dincv, tmp

  sub_name = 'get_cutoff_radii_psp'; call push_sub()

  ! local part ....
  ps%kbr(0:ps%l_max + 1) = 0.0_r8
  do ir = psp%nrval, 2, -1
    dincv = abs(psp%vlocal(ir)*psp%rofi(ir) + psp%z_val)
    if(dincv > eps) exit
  end do
  ps%kbr(ps%L_max + 1) = psp%rofi(ir + 1)

  ! non-local part....
  ps%rc_max = 0.0d0
  do l = 0, ps%L_max
    tmp = 0.0_r8
    do i = 1, 3
       do ir = psp%nrval, 2, -1
          dincv = abs(psp%kb(ir, l, i))
          if(dincv > eps) exit
       enddo
       tmp = psp%rofi(ir + 1)
       ps%kbr(l) = max(tmp, ps%kbr(l))
    enddo
    ps%rc_max = max(ps%rc_max, ps%kbr(l))
  end do

  call pop_sub(); return
end subroutine get_cutoff_radii_psp

subroutine get_cutoff_radii_tm(psf, ps)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps

  integer             :: l, ir
  real(r8)            :: dincv, phi

  sub_name = 'get_cutoff_radii_tm'; call push_sub()

  ! local part ....
  do ir = psf%nrval, 2, -1
    dincv = abs(psf%vlocal(ir)*psf%rofi(ir) + 2.0_r8*psf%zval)
    if(dincv > eps) exit
  end do
  ps%kbr(ps%L_max + 1) = psf%rofi(ir + 1)
  
  ! non-local part....
  ps%rc_max = 0.0d0
  do l = 0, ps%L_max
    if(l == ps%L_loc) then
      ps%kbr(l) = 0.0_r8
      cycle
    endif
    do ir = psf%nrval, 2, -1
      phi = (psf%rphi(ir, l)/psf%rofi(ir))*ps%dknrm(l)
      dincv = abs((psf%vps(ir, l) - psf%vlocal(ir))*phi)
      if(dincv > eps) exit
    enddo
    ps%kbr(l) = psf%rofi(ir + 1)
    ps%rc_max = max(ps%rc_max, ps%kbr(l))
  end do

  call pop_sub(); return
end subroutine get_cutoff_radii_tm

subroutine get_splines_tm(psf, ps)
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(inout) :: ps
  
  integer :: l, nrc, ir, nrcore
  real(r8) :: chc
  real(r8), allocatable :: hato(:), derhato(:)

  sub_name = 'get_splines_tm'; call push_sub()

  allocate(hato(psf%nrval), derhato(psf%nrval))

! Interpolate the KB-projection functions
  do l = 0, ps%l_max
    if(l == ps%L_loc) cycle

    hato = 0.0d0
    nrc = nint(log(ps%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
    hato(2:nrc) = (psf%vps(2:nrc, l) - psf%vlocal(2:nrc))*psf%rphi(2:nrc, l) * ps%dknrm(l) / psf%rofi(2:nrc)
    hato(1) = hato(2)    
    call spline_fit(psf%nrval, psf%rofi, hato, ps%kb(l, 1))

! and now the derivatives...
    call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
    call spline_fit(psf%nrval, psf%rofi, derhato, ps%dkb(l, 1))
  end do

! Now the part corresponding to the local pseudopotential
! where the asymptotic part is substracted 
!...local part...
  hato = 0.0_r8
  nrc = nint(log(ps%kbr(ps%L_max + 1)/psf%b + 1.0_r8)/psf%a) + 1

  hato(2:psf%nrval) = psf%vlocal(2:psf%nrval)*psf%rofi(2:psf%nrval) + 2.0_r8*psf%zval
  hato(1) = 2.0_r8*psf%zval
  
  ! WARNING: Rydbergs -> Hartrees
  hato = hato / 2._r8
  call spline_fit(psf%nrval, psf%rofi, hato, ps%vlocal)
  ps%vlocal_origin = psf%vlocal(1)

  ! and the derivative now
  call derivate_in_log_grid(psf%a, psf%b, psf%nrval, hato, derhato)
  call spline_fit(psf%nrval, psf%rofi, derhato, ps%dvlocal)

! Define the table for the pseudo-wavefunction components (using splines)
! with a correct normalization function
  do l = 0 , ps%L_max
    nrc = nint(log(ps%kbr(l)/psf%b + 1.0_r8)/psf%a) + 1
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
    nrc = nint(log(ps%kbr(l)/psp%b + 1.0_r8)/psp%a) + 1
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
  nrc = nint(log(ps%kbr(ps%L_max + 1)/psp%b + 1.0_r8)/psp%a) + 1

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
    nrc = nint(log(ps%kbr(l)/psp%b + 1.0_r8)/psp%a) + 1
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

! Output
subroutine write_pseudo_info_tm(unit, psf, ps)
  integer, intent(in)       :: unit
  type(ps_file), intent(in) :: psf
  type(ps_type), intent(in) :: ps

  integer :: l, i, j, k

  write(message(1),'(a)') ''
  write(message(2),'(a,a,a)')  '**********   Pseudopotential Information for: ', psf%namatm,'   **********'
  write(message(3),'(a)')      '**********   FLAVOUR: TROULLIER-MARTINS.'
  write(message(4),'(a,f10.4)')'Z    : ', ps%z
  write(message(5),'(a,f10.4)')'Z_val: ', ps%z_val
  write(message(6),'(a,a2)')   'Exchange/correlation used in generation: ', psf%icorr
  write(message(7),'(a,a3)')   'Relativistic character of calculations: ', psf%irel
  write(message(8),'(a,a4)')   'Type of core corrections: ', ps%icore
  write(message(9),'(a)')      'Signature of pseudopotential: '
  write(message(10), '(6x,6a10)') (psf%method(l), l=1,6)
  message(11)  = 'Valence configuration in calculations:   '
  message(12) = ' (orbital - occupancy - core radius)'
  call write_info(12, unit)
  i = 1; k = 1
  do while(index(psf%title(i:),'/') /= 0)
    j = i
    i = i + index(psf%title(i:),'/')
    write(message(k),'(a,a)') '  ', psf%title(j:i-2)
    k = k + 1
  enddo
  call write_info(k-1, unit)
  write(message(1),'(a,i2)')       'Maximum L-component to consider: ', ps%L_max
  write(message(2),'(a,i2)')     'Maximum L-component in file: ', psf%npotd - 1
  write(message(3),'(a)')        'Radial grid parameters ( R(I) = B*[ EXP(A*(I-1)) -1 ] )'
  write(message(4),'(a,es14.6)') '             A = ', psf%a
  write(message(5),'(a,es14.6)') '             B = ', psf%b
  write(message(6),'(a,i5)')     'Number of radial points: ', psf%nr
  write(message(7),'(a,i5)')     'nrval: ', psf%nrval
  write(message(8),'(a,a1,a,a1,4f9.5)') 'PS-generation cut-off radii: ', '[', &
                                   trim(units_out%length%abbrev), ']', ps%rc(0:3)/units_out%length%factor
  write(message(9),'(a,a1,a,a1,4f9.5)')'KB-spheres radii:            ', '[', &
                                   trim(units_out%length%abbrev), ']', ps%kbr(0:ps%l_max)/units_out%length%factor
  call write_info(9, unit)
  message(1) = 'KB-cosines: ['//trim(units_out%energy%abbrev)//']'
  write(message(2),'(10x,10f18.6)')  ps%dkbcos(0:ps%L_max)/units_out%energy%factor
  message(3) = 'KB-norms: [1/'//trim(units_out%energy%abbrev)//']'
  write(message(4),'(10x,10f18.6)')  ps%dknrm(0:ps%L_max)
  message(5) = 'Eigenvalues of pseudo-eigenfunctions ['// &
       trim(units_out%energy%abbrev)//']'
  write(message(6),'(10x,10f18.6)')  ps%eigen(0:ps%L_max)/units_out%energy%factor
  write(message(7),'(a,f18.6)')      'Atomic radius: ['//&
       trim(units_out%length%abbrev)//']', ps%rc_max/units_out%length%factor
  message(8) = '*************************************************************'
  message(9) = ''
  call write_info(9, unit)

  return
end subroutine write_pseudo_info_tm

subroutine write_pseudo_info_hgh(unit, psp, ps)
  integer, intent(in)       :: unit
  type(ps_st_params), intent(in) :: psp
  type(ps_type), intent(in) :: ps

  integer :: l, i, j, k

  write(message(1),'(a)') ''
  write(message(2),'(a,a,a)')  '**********   Pseudopotential Information for: ', psp%atom_name,'   **********'
  write(message(3),'(a)')      '**********   FLAVOUR: HARTWIGSEN-GOEDECKER-HUTTER.'
  write(message(4),'(a,f10.4)')'Z    : ', ps%z
  write(message(5),'(a,f10.4)')'Z_val: ', ps%z_val
  write(message(6),'(a,i2)')       'Maximum L-component to consider: ', ps%L_max
  call write_info(6, unit)
  write(message(1),'(a,a1,a,a1,4f9.5)') 'PS-generation cut-off radii: ', '[', &
                                   trim(units_out%length%abbrev), ']', ps%rc(0:3)/units_out%length%factor
  write(message(2),'(a,a1,a,a1,4f9.5)')'KB-spheres radii:            ', '[', &
                                   trim(units_out%length%abbrev), ']', ps%kbr(0:ps%l_max)/units_out%length%factor
  call write_info(2, unit)
  write(message(1),'(a,f18.6)')      'Atomic radius: ['//&
       trim(units_out%length%abbrev)//']', ps%rc_max/units_out%length%factor
  message(2) = '*************************************************************'
  message(3) = ''
  call write_info(3, unit)

  return
end subroutine write_pseudo_info_hgh







