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

module hgh
! For information about the Hartwinger-Goedecker-Hutter pseudopotentials, take a look at:
!  (1) S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996).
!  (2) C. Hartwinger, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998).
use global
use liboct
use io
use units
use kb
use logrid

implicit none

private
public  :: hgh_type, hgh_init, hgh_process, hgh_debug, hgh_end

! Next data type contains:
!   (a) the pseudopotential parameters, as read from a *.hgh file,
!   (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
type hgh_type
  ! HGH parameters.
  character(len=5) :: atom_name
  integer          :: z_val
  real(r8)         :: rlocal
  real(r8)         :: rc(0:3)
  real(r8)         :: c(1:4)
  real(r8)         :: h(0:3, 1:3, 1:3)
  real(r8)         :: k(0:3, 1:3, 1:3)

  integer          :: l_max     ! Maximum l for the Kleinmann-Bylander component.
  integer          :: l_max_occ ! Maximum l-value of the occupied wavefunctions
  real(r8), pointer :: occ(:)   ! Occupations of the pseudo-wavefunctions

  ! Logarithmic grid parameters
  type(logrid_type) :: g

  real(r8), pointer :: vlocal(:) ! Local potential
  real(r8), pointer :: kb(:,:,:) ! KB projectors
  real(r8), pointer :: kbr(:)    ! KB radii
  real(r8), pointer :: rphi(:,:), eigen(:)
end type hgh_type

real(r8), parameter :: eps = 1.0e-8_r8

interface vlocalr
  module procedure vlocalr_scalar, vlocalr_vector
end interface
interface projectorr
  module procedure projectorr_scalar, projectorr_vector
end interface

contains

subroutine hgh_init(psp, filename, ispin)
  type(hgh_type), intent(inout)     :: psp
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: ispin

  integer :: iunit, i
  logical :: found
  character(len=256) :: filename2

  sub_name = 'hgh_init'; call push_sub()

  filename2 = filename//'.hgh'
  inquire(file=filename2, exist=found)
  if(.not.found) then
    filename2 = SHARE_OCTOPUS//"/PP/HGH/"//filename//".hgh"
    inquire(file=filename2, exist=found)
    if(.not.found) then
      message(1) = "Pseudopotential file '"//trim(filename)//".hgh' not found!"
      call write_fatal(1)
    end if
  end if

  message(1) = "Info: Reading pseudopotential from file:"
  write(message(2), '(6x,3a)') "'", trim(filename2), "'"
  call write_info(2)

  call io_assign(iunit)
  open(iunit, file=filename2, form='formatted', status='old')
  i = load_params(iunit, psp)
  if(i .ne. 0) then
    message(1) = 'Error reading hgh file!'
    call write_fatal(1)
  endif
  call io_close(iunit)

  ! Finds out psp%l_max
  psp%l_max = 0
  do while(psp%rc(psp%l_max)>0.01_r8)
     psp%l_max = psp%l_max + 1
  end do
  psp%l_max = psp%l_max - 1

  ! Fixes psp%l_max_occ. It is equal to l_max, except for H and He, where
  ! l_max = -1 (no nonlocal components).
  psp%l_max_occ = max(0, psp%l_max)

  ! Initializes the logarithmic grid. Parameters are hard-coded.
       !psp%g%a = 1.25e-2_r8; psp%g%b = 4.0e-4_r8
       !psp%g%nrval = 1001
       psp%g%a = 3e-2_r8; psp%g%b = 4.0e-4_r8
       psp%g%nrval = 431
  call init_logrid(psp%g, psp%g%a, psp%g%b, psp%g%nrval)

  ! Allocation of stuff.
  allocate(psp%vlocal(1:psp%g%nrval), &
           psp%kbr(0:psp%l_max+1),    &
           psp%occ(0:psp%l_max_occ))
  if(psp%l_max >= 0) then
    allocate(psp%kb(1:psp%g%nrval, 0:psp%l_max, 1:3))
  endif
  allocate(psp%rphi(psp%g%nrval, 0:psp%l_max_occ), psp%eigen(0:psp%l_max_occ))
  psp%vlocal = 0.0_r8; psp%kb = 0.0_r8; psp%rphi = 0.0_r8


  ! Fixes the occupations.
  psp%occ = hgh_occs(psp%atom_name, psp%l_max_occ)

  call pop_sub(); return
end subroutine hgh_init

subroutine hgh_end(psp)
  type(hgh_type), intent(inout) :: psp

  sub_name = 'hgh_end'; call push_sub()

  deallocate(psp%kbr, psp%occ, psp%vlocal, psp%rphi, psp%eigen)
  if(associated(psp%kb)) deallocate(psp%kb)
  call kill_logrid(psp%g)

  call pop_sub(); return
end subroutine hgh_end

subroutine hgh_process(psp)
  type(hgh_type), intent(inout) :: psp

  integer :: l, i

  sub_name = 'hgh_process'; call push_sub()

! Fixes the local potential
  psp%vlocal(1:psp%g%nrval) = vlocalr(psp%g%rofi, psp)

! And the projectors
  do l = 0, psp%l_max
     do i = 1, 3
       psp%kb(1:psp%g%nrval, l, i) = projectorr(psp%g%rofi, psp, i, l)
     enddo
  enddo

! get the pseudoatomic eigenfunctions (WARNING: This is not correctly done yet: "some" wavefunctions
! are obtained, but not the real ones!!!
  call solve_schroedinger(psp)

! Define the KB-projector cut-off radii
  call get_cutoff_radii(psp)

  call pop_sub(); return
end subroutine hgh_process

function load_params(unit, params)
  integer, intent(in)             :: unit        ! where to read from
  type(hgh_type), intent(out)     :: params      ! obvious
  integer                         :: load_params ! 0 if success, 
                                                 ! 1 otherwise.

  integer :: i, iostat, j, k
  character(len=200) :: line

  sub_name = 'load_params'; call push_sub()

! Set initially everything to zero.
  params%c(1:4) = 0.0_r8; params%rlocal = 0.0_r8;
  params%rc = 0.0_r8; params%h = 0.0_r8; params%k = 0.0_r8

! Reads the file in a hopefully smart way
  iostat = 1; j = 5
  read(unit,'(a)') line
  do while((iostat .ne. 0) .and. (j > 0))
    j = j - 1
    read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal, params%c(1:j)
  end do
  if(j<1) read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal
  if( iostat.ne.0 ) then
    load_params = 1
    call pop_sub(); return
  endif

  read(unit,'(a)', iostat = iostat) line
  if(iostat .ne. 0) then
    load_params = 0
    call pop_sub(); return
  endif
  iostat = 1; j = 4
  do while((iostat .ne. 0) .and. (j > 0))
     j = j - 1
     read(line, *, iostat=iostat) params%rc(0), (params%h(0, i, i), i = 1, j)
  enddo
  if(j < 0) then
    load_params = 2
    call pop_sub(); return
  endif

  kloop: do k = 1, 3
     read(unit, '(a)', iostat = iostat) line
     if(iostat .ne. 0) exit kloop
     iostat = 1; j = 4
     do while((iostat .ne. 0) .and. (j > 0))
        j = j - 1
        read(line, *, iostat = iostat) params%rc(k), (params%h(k, i, i), i = 1, j)
     enddo
     if(params%rc(k) == 0.0_r8) exit kloop
     read(unit, '(a)') line
     iostat = 1; j = 4
     do while((iostat .ne. 0) .and. (j>0))
        j = j - 1
        read(line, *, iostat = iostat) (params%k(k, i, i), i = 1, 3)
     enddo
  enddo kloop

! Fill in the rest of the parameters matrices...
  params%h(0, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(3.0_r8/5.0_r8)       * params%h(0, 2, 2)
  params%h(0, 1, 3) =  (1.0_r8/2.0_r8) * sqrt(5.0_r8/21.0_r8)      * params%h(0, 3, 3)
  params%h(0, 2, 3) = -(1.0_r8/2.0_r8) * sqrt(100.0_r8/63.0_r8)    * params%h(0, 3, 3)
  params%h(1, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(5.0_r8/7.0_r8)       * params%h(1, 2, 2)
  params%h(1, 1, 3) =  (1.0_r8/6.0_r8) * sqrt(35.0_r8/11.0_r8)     * params%h(1, 3, 3)
  params%h(1, 2, 3) = -(1.0_r8/6.0_r8) * (14.0_r8 / sqrt(11.0_r8)) * params%h(1, 3, 3)
  params%h(2, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(7.0_r8/9.0_r8)       * params%h(2, 2, 2)
  params%h(2, 1, 3) =  (1.0_r8/2.0_r8) * sqrt(63.0_r8/143.0_r8)    * params%h(2, 3, 3)
  params%h(2, 2, 3) = -(1.0_r8/2.0_r8) * (18.0_r8/sqrt(143.0_r8))  * params%h(2, 3, 3)

  params%k(0, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(3.0_r8/5.0_r8)       * params%k(0, 2, 2)
  params%k(0, 1, 3) =  (1.0_r8/2.0_r8) * sqrt(5.0_r8/21.0_r8)      * params%k(0, 3, 3)
  params%k(0, 2, 3) = -(1.0_r8/2.0_r8) * sqrt(100.0_r8/63.0_r8)    * params%k(0, 3, 3)
  params%k(1, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(5.0_r8/7.0_r8)       * params%k(1, 2, 2)
  params%k(1, 1, 3) =  (1.0_r8/6.0_r8) * sqrt(35.0_r8/11.0_r8)     * params%k(1, 3, 3)
  params%k(1, 2, 3) = -(1.0_r8/6.0_r8) * (14.0_r8 / sqrt(11.0_r8)) * params%k(1, 3, 3)
  params%k(2, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(7.0_r8/9.0_r8)       * params%k(2, 2, 2)
  params%k(2, 1, 3) =  (1.0_r8/2.0_r8) * sqrt(63.0_r8/143.0_r8)    * params%k(2, 3, 3)
  params%k(2, 2, 3) = -(1.0_r8/2.0_r8) * (18.0_r8/sqrt(143.0_r8))  * params%k(2, 3, 3)

! Parameters are symmetric.
  do k = 0, 3
     do i = 1, 3
        do j = i + 1, 3
           params%h(k, j, i) = params%h(k, i, j)
           params%k(k, j, i) = params%k(k, i, j)
        enddo
     enddo
  enddo

  load_params = 0
  call pop_sub(); return
end function load_params

subroutine get_cutoff_radii(psp)
  type(hgh_type), intent(inout)     :: psp

  integer  :: ir, l, i
  real(r8) :: dincv, tmp

  sub_name = 'get_cutoff_radii_psp'; call push_sub()

  ! local part ....
  psp%kbr(0:psp%l_max + 1) = 0.0_r8
  do ir = psp%g%nrval, 2, -1
    dincv = abs(psp%vlocal(ir)*psp%g%rofi(ir) + psp%z_val)
    if(dincv > eps) exit
  end do
  psp%kbr(psp%l_max + 1) = psp%g%rofi(ir + 1)

  ! non-local part....
  do l = 0, psp%l_max
    tmp = 0.0_r8
    do i = 1, 3
       do ir = psp%g%nrval, 2, -1
          dincv = abs(psp%kb(ir, l, i))
          if(dincv > eps) exit
       enddo
       tmp = psp%g%rofi(ir + 1)
       psp%kbr(l) = max(tmp, psp%kbr(l))
    enddo
  end do

  call pop_sub(); return
end subroutine get_cutoff_radii

! Local pseudopotential, both in real and reciprocal space.
function vlocalr_scalar(r, p)
  type(hgh_type), intent(in)     :: p
  real(r8), intent(in)           :: r
  real(r8)                       :: vlocalr_scalar

  real(r8) :: r1, r2, r4, r6

  r1 = r/p%rlocal; r2 = r1**2; r4 = r2**2; r6 = r4*r2

  if(r < 1e-7_r8) then
     vlocalr_scalar = - (2.0_r8 * p%z_val)/(sqrt(2.0_r8*M_Pi)*p%rlocal) + p%c(1)
     return
  endif

  vlocalr_scalar = - (p%z_val/r)*oct_erf(r1/sqrt(2.0_r8))   &
                   + exp( -(1.0_r8/2.0_r8)*r2 ) *    &
                   ( p%c(1) + p%c(2)*r2 + p%c(3)*r4 + p%c(4)*r6 )

end function vlocalr_scalar

function vlocalr_vector(r, p)
  type(hgh_type), intent(in)      :: p
  real(r8), intent(in)            :: r(:)
  real(r8), pointer               :: vlocalr_vector(:)

  integer :: i

  allocate(vlocalr_vector(size(r)))
  do i = 1, size(r)
     vlocalr_vector(i) = vlocalr_scalar(r(i), p)
  enddo

end function vlocalr_vector

function vlocalg(g, p)
  type(hgh_type), intent(in)     :: p
  real(r8), intent(in)           :: g
  real(r8)                       :: vlocalg

  real(r8) :: g1, g2, g4, g6

  g1 = g*p%rlocal; g2 = g1*g1; g4 = g2*g2; g6 = g4*g2

  vlocalg = -(4.0_r8*M_Pi*p%z_val/g**2) * exp( -g2/2.0_r8) +            &
            sqrt(8.0_r8*M_Pi**3) * p%rlocal**3 * exp( -g2/2.0_r8) *    &
            ( p%c(1) + p%c(2)*(3.0_r8 - g2) + p%c(3)*(15.0_r8 - 10.0_r8*g2 + g4) + &
              p%c(4)*(105.0_r8 -105.0_r8*g2 + 21.0_r8*g4 - g6) )

end function vlocalg

function projectorr_scalar(r, p, i, l)
  type(hgh_type), intent(in)     :: p
  real(r8), intent(in)           :: r
  integer, intent(in)            :: i, l
  real(r8)                       :: projectorr_scalar

  real(r8) :: x, y, rr

  x = l + real(4*i-1, r8)/2
  y = oct_gamma(x); x = sqrt(y)
  if(l==0 .and. i==1) then
    rr = 1._r8
  else
    rr = r ** (l + 2*(i-1))
  endif
 
  projectorr_scalar = sqrt(2.0_r8) * rr * exp(-r**2/(2.0_r8*p%rc(l)**2)) / &
              (  p%rc(l)**(l + real(4*i-1,r8)/2) * x )  

end function projectorr_scalar

function projectorr_vector(r, p, i, l)
  type(hgh_type), intent(in)     :: p
  real(r8), intent(in)           :: r(:)
  integer, intent(in)            :: i, l
  real(r8), pointer              :: projectorr_vector(:)

  integer :: j

  allocate(projectorr_vector(size(r)))
  do j=1, size(r)
     projectorr_vector(j) = projectorr_scalar(r(j), p, i, l)
  enddo
 
end function projectorr_vector

function projectorg(g, p, i, l)
  type(hgh_type), intent(in)     :: p
  real(r8), intent(in)           :: g
  integer, intent(in)            :: i, l
  real(r8)                       :: projectorg

  !real(r8), external :: gamma
  real(r8) :: pif, ex

  pif = M_Pi**(5.0_r8/4.0_r8)

  ex = exp( (1.0_r8/2.0_r8)*(g*p%rc(l))**2 )

  select case(l)
  case(0)
    select case(i)
    case(1)
      projectorg = ( 4.0_r8*sqrt(2.0_r8*p%rc(0)**3)*pif ) / ex
    case(2)
      projectorg = ( sqrt(8.0_r8*2*p%rc(0)**3/15.0_r8)*pif * &
           (3.0_r8 - (g*p%rc(0))**2) ) / ex
    case(3)
      projectorg = ( 16.0_r8*sqrt(2.0_r8*p%rc(0)**3/105.0_r8) * pif * &
           (15.0_r8 - 10.0_r8*g**2*p%rc(0)**2 + g**4*p%rc(0)**2) ) / (3.0_r8*ex)
    end select

  case(1)
    select case(i)
    case(1)
      projectorg = ( 8.0_r8*sqrt(p%rc(1)**5/3.0_r8)*pif*g ) / ex
    case(2)
      projectorg = ( 16.0_r8*sqrt(p%rc(1)**5/105.0_r8)* pif * g * & 
           ( 5.0_r8 - (g*p%rc(1))**2 ) ) / ex
    case(3)
      projectorg = ( 32.0_r8*sqrt(p%rc(1)**5/1155.0_r8)* pif * g * &
           ( 35.0_r8 - 14.0_r8*g**2*p%rc(1)**2 + (g*p%rc(1))**4 ) ) / &
           (3.0_r8*ex)
    end select

  case(2)
    select case(i)
    case(1)
      projectorg = ( 8.0_r8 * sqrt(2.0_r8*p%rc(2)**7/15.0_r8) * pif * g**2 ) / ex
    case(2)
      projectorg = ( 16.0_r8 * sqrt(2.0_r8*p%rc(2)**7/105.0_r8) * pif * g**2 * &
           (7.0_r8 - g**2*p%rc(2)**2) ) / (3.0_r8*ex)
    case(3)
      projectorg = 0.0_r8 ! ??
    end select
    
  case(3)    
    ! This should be checked. Probably will not be needed in an near future...
  end select

end function projectorg

subroutine solve_schroedinger(psp)
  type(hgh_type), intent(inout)     :: psp
  
  integer :: iter, ir, l, nnode, nprin, i, j, irr
  real(r8) :: vtot, r2, e, z, dr, rmax, f, dsq, a2b4, diff, nonl
  real(r8), allocatable :: s(:), hato(:), g(:), y(:), prev(:), rho(:), ve(:)
  real(r8), parameter :: tol = 1.0e-4_r8

  sub_name = 'solve_schroedinger'; call push_sub()

  ! Allocations...
  allocate(s(psp%g%nrval), hato(psp%g%nrval), g(psp%g%nrval), y(psp%g%nrval), prev(psp%g%nrval), &
           rho(psp%g%nrval), ve(psp%g%nrval))
  hato = 0.0_r8; g = 0.0_r8;  y = 0.0_r8; rho = 0.0_r8; ve = 0.0_r8

  ! These numerical parameters have to be fixed for egofv to work.
  s(2:psp%g%nrval) = psp%g%drdi(2:psp%g%nrval)*psp%g%drdi(2:psp%g%nrval)
  s(1) = s(2)
  a2b4 = 0.25_r8*psp%g%a**2

  ! Let us be a bit informative.
  if(conf%verbose > 20 .and. mpiv%node == 0) then
    message(1) = '      Calculating atomic pseudo-eigenfunctions for specie '//psp%atom_name//'....'
    call write_info(1)
  endif

  ! "Double" self consistent loop: nonlocal and xc parts have to be calculated
  ! self-consistently.
  diff = 1e5_r8
  iter = 0
  self_consistent: do while( diff > tol )
    prev = rho
    iter = iter + 1
    do l = 0, psp%l_max_occ
      do ir = 2, psp%g%nrval
        vtot = 2*psp%vlocal(ir) + ve(ir) + dble(l*(l + 1))/(psp%g%rofi(ir)**2)
        nonl = 0.0_r8
        if(iter>1 .and. psp%l_max >=0 .and. psp%rphi(ir, l) > 1.0e-7) then
           do i = 1, 3
           do j = 1, 3
              do irr = 2, psp%g%nrval
                 nonl = nonl + psp%h(l, i, j)*psp%kb(ir, l, i)* & 
                        psp%g%drdi(irr)*psp%g%rofi(irr)*psp%rphi(irr, l)*psp%kb(irr,l,j)
              enddo
           enddo
           enddo
           nonl = 2*nonl/psp%rphi(ir, l)*psp%g%rofi(ir)
        endif
        vtot = vtot + nonl
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = hato(2)
      nnode = 1; nprin = l + 1
      if(iter == 1) then
         e = -((psp%z_val/dble(nprin))**2); z = psp%z_val
      else
         e = psp%eigen(l); z = psp%z_val
      endif
      dr = -1.0e5_r8; rmax = psp%g%rofi(psp%g%nrval)
      call egofv(hato, s, psp%g%nrval, e, g, y, l, z, psp%g%a, psp%g%b, rmax, nprin, nnode, dr)
      psp%eigen(l) = e

      psp%rphi(2:psp%g%nrval, l) = g(2:psp%g%nrval) * sqrt(psp%g%drdi(2:psp%g%nrval))
      psp%rphi(1, l) = psp%rphi(2, l)
    end do
    rho = 0.0_r8
    do l = 0, psp%l_max_occ
       rho = rho + psp%occ(l)*psp%rphi(1:psp%g%nrval, l)**2
    enddo
    if(iter>1) rho = 0.5*rho + 0.5*prev
    diff = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*(rho(2:psp%g%nrval)-prev(2:psp%g%nrval))**2))
    !if(conf%verbose>20 .and. mpiv%node == 0) then
    !  write(message(1),'(a,i4,a,e10.2)') '      Iter =',iter,'; Diff =',diff
    !  call write_info(1)
    !endif
    call calculate_valence_screening(psp, rho, ve)

  enddo self_consistent

  if(conf%verbose > 20 .and. mpiv%node == 0) then
    message(1) = '      Done.'
    call write_info(1)
  endif

  !  checking normalization of the calculated wave functions
  do l = 0, psp%l_max_occ
    e = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*psp%rphi(2:psp%g%nrval, l)**2))
    e = abs(e - 1.0d0)
    if (e > 1.0d-5 .and. conf%verbose > 0) then
      write(message(1), '(a,i2,a)') "Eigenstate for l = ", l , ' is not normalized'
      write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      call write_warning(2)
    end if
  end do

  ! Output in Ha and not in stupid Rydbergs.
  psp%eigen = psp%eigen / 2.0_r8

  ! Deallocations.
  deallocate(s, hato, g, y, rho, prev)

  call pop_sub; return
end subroutine solve_schroedinger

subroutine calculate_valence_screening(psp, rhoval, ve)
  type(hgh_type), intent(in) :: psp
  real(r8), intent(in)       :: rhoval(1:psp%g%nrval)
  real(r8), intent(out)      :: ve(1:psp%g%nrval)

  character(len=5) :: xcfunc, xcauth
  integer          :: irelt
  real(r8)         :: e_x, e_c, dx, dc, r2
  real(r8), allocatable :: v_xc(:,:), auxrho(:)

  sub_name = 'valence_screening'; call push_sub()

  ! Set this to zero...
  ve = 0.0_r8

  ! Allocates the auxiliary stuff
  allocate(v_xc(psp%g%nrval, 1), auxrho(psp%g%nrval))
  auxrho = rhoval

  ! Gets the hartree potential first.
  call vhrtre(rhoval, ve, psp%g%rofi, psp%g%drdi, psp%g%s, psp%g%nrval, psp%g%a)

  ! We will assume this for the moment (GGA pseudos could also be got from Goedecker...)
  xcfunc = 'LDA'
  xcauth = 'PZ'
  irelt = 1

  ! Fix aux_rho (just rho/(4 \pi r^2))
  auxrho(2:psp%g%nrval) = auxrho(2:psp%g%nrval)/(4.0_r8*M_PI*psp%g%rofi(2:psp%g%nrval)**2)
  r2 = psp%g%rofi(2)/(psp%g%rofi(3)-psp%g%rofi(2))
  auxrho(1) = auxrho(2) - (auxrho(3)-auxrho(2))*r2

  ! And gets v_xc.
  call atomxc(xcfunc, xcauth, irelt, psp%g%nrval, psp%g%nrval, psp%g%rofi, &
              1, auxrho, e_x, e_c, dx, dc, v_xc)

  ! Sums the vxc to the hartree.
  ve(1:psp%g%nrval) = ve(1:psp%g%nrval) + v_xc(1:psp%g%nrval, 1)

  ! And deallocates.
  deallocate(v_xc, auxrho)

  call pop_sub(); return
end subroutine calculate_valence_screening

function hgh_occs(label, lmax)
  character(len=5), intent(in) :: label
  integer, intent(in)          :: lmax
  real(r8)                     :: hgh_occs(0:lmax)

  select case(label)
    case('H');  hgh_occs(0)   = 1
    case('He'); hgh_occs(0)   = 2
    case('Li'); hgh_occs(0:1) = (/ 1, 0 /)
    case('Be'); hgh_occs(0:1) = (/ 2, 0 /)
    case('B');  hgh_occs(0:1) = (/ 2, 1 /)
    case('C');  hgh_occs(0:1) = (/ 2, 2 /)
    case('N');  hgh_occs(0:1) = (/ 2, 3 /)
    case('O');  hgh_occs(0:1) = (/ 2, 4 /)
    case('F');  hgh_occs(0:1) = (/ 2, 5 /)
    case('Ne'); hgh_occs(0:1) = (/ 2, 6 /)
    case('Na'); hgh_occs(0:1) = (/ 1, 0 /)
    case('Mg'); hgh_occs(0:1) = (/ 2, 0 /)
    case('Al'); hgh_occs(0:1) = (/ 2, 1 /)
    case('Si'); hgh_occs(0:1) = (/ 2, 2 /)
    case('P');  hgh_occs(0:1) = (/ 2, 3 /)
    case('S');  hgh_occs(0:1) = (/ 2, 4 /)
    case('Cl'); hgh_occs(0:1) = (/ 2, 5 /)
    case('Ar'); hgh_occs(0:1) = (/ 2, 6 /)
    case('K');  hgh_occs(0:2) = (/ 1, 0, 0 /)
    case('Ca'); hgh_occs(0:2) = (/ 2, 0, 0 /)
    case('Sc'); hgh_occs(0:2) = (/ 2, 0, 1 /)
    case('Ti'); hgh_occs(0:2) = (/ 2, 0, 2 /)
    case('Au'); hgh_occs(0:2) = (/ 1, 0, 0 /)
    case('Au_sc'); hgh_occs(0:2) = (/1, 0, 10 /)
    case('Hg'); hgh_occs(0:2) = (/ 2, 0, 0 /)
    case default
    if(mpiv%node == 0) then
      message(1) = 'Specie '//trim(label)//' not included in occupation-numbers data base.'
      message(2) = 'Add it or complain to the lazy developer of this piece of code:'
      message(3) = 'alberto.castro@tddft.org'
      call write_fatal(3)
    endif
  end select

end function hgh_occs

subroutine hgh_debug(psp)
  type(hgh_type), intent(in) :: psp

  integer :: hgh_unit, loc_unit, dat_unit, kbp_unit, wav_unit, i, l, k

  sub_name = 'hgh_debug'; call push_sub()

  ! Opens files.
  call oct_mkdir(C_string('pseudos/'+'hgh.'+trim(psp%atom_name)))
  call io_assign(hgh_unit); call io_assign(loc_unit); call io_assign(wav_unit)
  call io_assign(dat_unit); call io_assign(kbp_unit)
  open(hgh_unit, file = 'pseudos/'+'hgh.'+trim(psp%atom_name)+'/'+'hgh')
  open(loc_unit, file = 'pseudos/'+'hgh.'+trim(psp%atom_name)+'/'+'local')
  open(dat_unit, file = 'pseudos/'+'hgh.'+trim(psp%atom_name)+'/'+'info')
  open(kbp_unit, file = 'pseudos/'+'hgh.'+trim(psp%atom_name)+'/'+'nonlocal')
  open(wav_unit, file = 'pseudos/'+'hgh.'+trim(psp%atom_name)+'/'+'wave')

  ! Writes down the input file, to be checked agains SHARE_OCTOPUS/PP/HGH/ATOM_NAME.hgh
  write(hgh_unit,'(a5,i6,5f12.6)') psp%atom_name, psp%z_val, psp%rlocal, psp%c(1:4)
  write(hgh_unit,'(  11x,4f12.6)') psp%rc(0), (psp%h(0,i,i), i = 1, 3)
  do k = 1, 3
     write(hgh_unit,'(  11x,4f12.6)') psp%rc(k), (psp%h(k, i, i), i = 1, 3)
     write(hgh_unit,'(  23x,4f12.6)')            (psp%k(k, i, i), i = 1, 3)
  enddo

  ! Writes down some info.
  write(dat_unit,'(a,i3)')        'lmax  = ', psp%l_max
  if(psp%l_max >= 0) then
     write(dat_unit,'(a,4f14.6)') 'kbr   = ', psp%kbr
  endif
  write(dat_unit,'(a,4f14.6)')    'eigen = ', psp%eigen

  ! Writes down local part.
  do i = 1, psp%g%nrval
     write(loc_unit, *) psp%g%rofi(i), psp%vlocal(i)
  enddo

  ! Writes down nonlocal part.
  if(psp%l_max >=0) then
    do i = 1, psp%g%nrval
       write(kbp_unit, *) psp%g%rofi(i), ( (psp%kb(i, l, k) ,k = 1, 3),l = 0, psp%l_max)
    enddo
  endif

  ! And the pseudo-wavefunctions.
  do i = 1, psp%g%nrval
     write(wav_unit, *) psp%g%rofi(i), (psp%rphi(i, l), l = 0, psp%l_max_occ)
  enddo

  ! Closes files and exits
  call io_close(hgh_unit); call io_close(loc_unit); call io_close(wav_unit)
  call io_close(dat_unit); call io_close(kbp_unit)

  call pop_sub(); return
end subroutine hgh_debug

end module hgh
