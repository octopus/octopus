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

implicit none

private
public  :: &! DATA TYPES:
           ps_st_params,          &
           &! PROCEDURES:
           load_params,           &
           vlocalr,               &
           vlocalg,               &
           projectorr,            &
           projectorg,            &
           ps_ghg_read_file,      &
           hgh_debug,             &
           solve_schroedinger_hgh

! Next data type contains (a) the pseudopotential parameters, as read from a *.hgh file,
! (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
type ps_st_params
! HGH parameters.
  character(len=2) :: atom_name
  integer          :: z_val
  real(r8)         :: rlocal
  real(r8)         :: rc(0:3)
  real(r8)         :: c(1:4)
  real(r8)         :: h(0:3, 1:3, 1:3)
  real(r8)         :: k(1:3, 1:3, 1:3)
! Maximum l
  integer          :: l_max
! Logarithmic grid parameters
  real(r8)         :: a, b
  integer          :: nrval
  real(r8), pointer :: s(:), drdi(:), rofi(:)
! Local potential
  real(r8), pointer :: vlocal(:)
! KB projectors
  real(r8), pointer :: kb(:,:,:)
! Pseudo wave functions, pseudo eigenvalues...
  real(r8), pointer :: rphi(:,:), eigen(:)
end type ps_st_params

interface vlocalr
  module procedure vlocalr_scalar, vlocalr_vector
end interface
interface projectorr
  module procedure projectorr_scalar, projectorr_vector
end interface

contains

subroutine ps_ghg_read_file(psp, filename)
  type(ps_st_params), intent(inout) :: psp
  character(len=*), intent(in)      :: filename

  integer :: iunit, i
  logical :: found

  sub_name = 'ps_ghg_read_file'; call push_sub()

  inquire(file=filename//'.hgh', exist=found)
  if(found) then
    call io_assign(iunit)
    open(iunit, file=filename//'.hgh', form='formatted')
    i = load_params(iunit, psp)
    if(i .ne. 0) then
      message(1) = 'Error reading hgh file'
      call write_fatal(1)
    endif
    call io_close(iunit)
  else
    message(1) = "Pseudopotential file '"//trim(filename)//"hgh' not found"
    call write_fatal(1)
  end if
!!$      ! This should go someday...
!!$      write(*,'(a)') '****'
!!$      write(*,'(a2,i5,5f14.8)') psf%atom_name, psf%z_val, psf%rlocal, psf%c(1:4)
!!$      write(*,'(7x,4f14.8)')    psf%rc(0), (psf%h(0, i, i), i=1, 3)
!!$      write(*,'(7x,4f14.8)')    psf%rc(1), (psf%h(1, i, i), i=1, 3)
!!$      write(*,'(21x, 3f14.8)')             (psf%k(1, i, i), i=1, 3)
!!$      stop

  psp%l_max = 0
  do while(psp%rc(psp%l_max)>0.01_r8)
     psp%l_max = psp%l_max + 1
  end do
  psp%l_max = psp%l_max - 1

  call pop_sub(); return
end subroutine ps_ghg_read_file

! Loads all the necessary parameters.
function load_params(unit, params)
  integer, intent(in)             :: unit        ! where to read from
  type(ps_st_params), intent(out) :: params      ! obvious
  integer                         :: load_params ! 0 if success, 
                                                 ! 1 otherwise.

  integer :: i, iostat, j, k
  character(len=200) :: line

  sub_name = 'load_params'; call push_sub()

! Set initially everything to zero.
  params%c(1:4 ) = 0.0_r8; params%rlocal = 0.0_r8;
  params%rc = 0.0_r8; params%h = 0.0_r8; params%k = 0.0_r8

! Reads the file in a hopefully smart way
  iostat = 1; j = 5
  read(unit,'(a)') line
  do while((iostat .ne. 0) .and. (j > 0))
     j = j - 1
     read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal, params%c(1:j)
  enddo
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
  params%h(0, 1, 2) = - (1.0_r8/2.0_r8) * sqrt(3.0_r8/5.0_r8) * params%h(0, 2, 2)
  params%h(0, 1, 3) = (1.0_r8/2.0_r8) * sqrt(5.0_r8/21.0_r8) * params%h(0, 3, 3)
  params%h(0, 2, 3) = - (1.0_r8/2.0_r8) * sqrt(100.0_r8/63.0_r8) * params%h(0, 3, 3)
  params%h(1, 1, 2) = - (1.0_r8/2.0_r8) * sqrt(5.0_r8/7.0_r8) * params%h(1, 2, 2)
  params%h(1, 1, 3) = (1.0_r8/6.0_r8) * sqrt(35.0_r8/11.0_r8) * params%h(1, 3, 3)
  params%h(1, 2, 3) = -(1.0_r8/6.0_r8) * (14.0_r8 / sqrt(11.0_r8)) * params%h(1, 3, 3)
  params%h(2, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(7.0_r8/9.0_r8) * params%h(2, 2, 2)
  params%h(2, 1, 3) = (1.0_r8/2.0_r8) * sqrt(63.0_r8/143.0_r8) * params%h(2, 3, 3)
  params%h(2, 2, 3) = -(1.0_r8/2.0_r8) * (18.0_r8/sqrt(143.0_r8)) * params%h(2, 3, 3)

  params%k(0, 1, 2) = - (1.0_r8/2.0_r8) * sqrt(3.0_r8/5.0_r8) * params%k(0, 2, 2)
  params%k(0, 1, 3) = (1.0_r8/2.0_r8) * sqrt(5.0_r8/21.0_r8) * params%k(0, 3, 3)
  params%k(0, 2, 3) = - (1.0_r8/2.0_r8) * sqrt(100.0_r8/63.0_r8) * params%k(0, 3, 3)
  params%k(1, 1, 2) = - (1.0_r8/2.0_r8) * sqrt(5.0_r8/7.0_r8) * params%k(1, 2, 2)
  params%k(1, 1, 3) = (1.0_r8/6.0_r8) * sqrt(35.0_r8/11.0_r8) * params%k(1, 3, 3)
  params%k(1, 2, 3) = -(1.0_r8/6.0_r8) * (14.0_r8 / sqrt(11.0_r8)) * params%k(1, 3, 3)
  params%k(2, 1, 2) = -(1.0_r8/2.0_r8) * sqrt(7.0_r8/9.0_r8) * params%k(2, 2, 2)
  params%k(2, 1, 3) = (1.0_r8/2.0_r8) * sqrt(63.0_r8/143.0_r8) * params%k(2, 3, 3)
  params%k(2, 2, 3) = -(1.0_r8/2.0_r8) * (18.0_r8/sqrt(143.0_r8)) * params%k(2, 3, 3)

! Parameters are symmetric.
  do k = 0, 3
     do i = 1, 3
        do j = i + 1, 3
           params%h(k, j, i) = params%h(k, i, j)
           params%k(k, j, i) = params%k(k, i, j)
        enddo
     enddo
  enddo

!!$  do k = 0, 3
!!$     write(*, '(a)')
!!$     do i = 1, 3
!!$        write(*, *) (params%h(k, i, j), j = 1, 3)
!!$     enddo
!!$  enddo
!!$  stop

  load_params = 0
  call pop_sub(); return
end function load_params

! Local pseudopotential, both in real and reciprocal space.
function vlocalr_scalar(r, p)
  type(ps_st_params), intent(in) :: p
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
  type(ps_st_params), intent(in)  :: p
  real(r8), intent(in)            :: r(:)
  real(r8), pointer               :: vlocalr_vector(:)

  integer :: i

  allocate(vlocalr_vector(size(r)))
  do i=1, size(r)
     vlocalr_vector(i) = vlocalr_scalar(r(i), p)
  enddo

end function vlocalr_vector

function vlocalg(g, p)
  type(ps_st_params), intent(in) :: p
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
  type(ps_st_params), intent(in) :: p
  real(r8), intent(in)           :: r
  integer, intent(in)            :: i, l
  real(r8)                       :: projectorr_scalar

  real(r8) :: x, y

  x = l + real(4*i-1, r8)/2
  y = oct_gamma(x); x = sqrt(y)
 
  projectorr_scalar = sqrt(2.0_r8) * r ** (l + 2*(i-1)) * exp(-r**2/(2.0_r8*p%rc(l)**2)) / &
              (  p%rc(l)**(l + real(4*i-1,r8)/2) * x )  

end function projectorr_scalar

function projectorr_vector(r, p, i, l)
  type(ps_st_params), intent(in) :: p
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
  type(ps_st_params), intent(in) :: p
  real(r8), intent(in)           :: g
  integer, intent(in)            :: i, l
  real(r8)                       :: projectorg

  !real(r8), external :: gamma
  real(r8) :: pif, ex

  pif = M_Pi**(5.0_r8/4.0_r8)

  select case(l)
  case(0)

    ex = exp( (1.0_r8/2.0_r8)*(g*p%rc(0))**2 )
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

    ex = exp( (1.0_r8/2.0_r8)*(g*p%rc(1))**2 )
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

    ex = exp( (1.0_r8/2.0_r8)*(g*p%rc(2))**2 )
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

    !!! This should be checked. Probably will not be needed in an near future...
  end select

end function projectorg

subroutine hgh_debug(psp, filename)
  type(ps_st_params), intent(in) :: psp
  character(len=*)               :: filename

  integer :: iunit, ir, i, l

  call io_assign(iunit)
  open(unit=iunit, file=trim(filename)//'.hgh.debug')

   do ir = 1, psp%nrval
      write(iunit, *) psp%rofi(ir) / units_out%length%factor, &
                      psp%vlocal(ir) / units_out%energy%factor, &
                      ( psp%kb(ir, l, 1), l = 0, psp%l_max), &
                      ( psp%kb(ir, l, 2), l = 0, psp%l_max), &
                      ( psp%kb(ir, l, 2), l = 0, psp%l_max)
   end do

  call io_close(iunit)

end subroutine hgh_debug

subroutine solve_schroedinger_hgh(psp, lmax)
  type(ps_st_params), intent(inout) :: psp
  integer, intent(in)               :: lmax

  integer :: ir, l, nnode, nprin
  real(r8) :: rpb, ea, vtot, r2, e, z, dr, rmax, f, dsq, a2b4
  real(r8), allocatable :: s(:), hato(:), g(:), y(:)

  sub_name = 'solve_schroedinger_hgh'; call push_sub()

! Allocations...
  allocate(s(psp%nrval), hato(psp%nrval))

! Calculation of the pseudo-wave functions.
  s(2:psp%nrval) = psp%drdi(2:psp%nrval)*psp%drdi(2:psp%nrval)
  s(1) = s(2)
  a2b4 = 0.25_r8*psp%a**2
  allocate(g(psp%nrval), y(psp%nrval))
  g = 0.0_r8;  y = 0.0_r8

  do l = 0, lmax
    do ir = 2, psp%nrval
      vtot = psp%vlocal(ir) + dble(l*(l + 1))/(psp%rofi(ir)**2)
      hato(ir) = vtot*s(ir) + a2b4
    end do
    hato(1) = hato(2)
    
    nnode = 1; nprin = l + 1
    e = -((psp%z_val/dble(nprin))**2); z = psp%z_val
    dr = -1.0e5_r8; rmax = psp%rofi(psp%nrval)
    call egofv(hato, s, psp%nrval, e, g, y, l, z, psp%a, psp%b, rmax, nprin, nnode, dr)
    psp%eigen(l) = e

    psp%rphi(2:psp%nrval, l) = g(2:psp%nrval) * sqrt(psp%drdi(2:psp%nrval))
    psp%rphi(1, l) = psp%rphi(2, l)
  end do
  deallocate(g, y)

  !  checking normalization of the calculated wave functions
  do l = 0, lmax
    e = sqrt(sum(psp%drdi(2:psp%nrval)*psp%rphi(2:psp%nrval, l)**2))
    e = abs(e - 1.0d0)
    if (e > 1.0d-5 .and. conf%verbose > 0) then
      write(message(1), '(a,i2,a)') "Eigenstate for l = ", l , ' is not normalized'
      write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      call write_warning(2)
    end if
  end do

  deallocate(s, hato)
  call pop_sub; return
end subroutine solve_schroedinger_hgh

end module hgh
