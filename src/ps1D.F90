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

implicit none

private
public :: ps_type, ps_init, ps_end

! This public type should be kept equal to the 3D case, although most
! of the components are unnecessary.
! Some of the things should be simplified eventually...
type ps_type
  character(len=3) :: flavour

  type(spline_type), pointer :: kb(:,:)   ! Kleynman-Bylander projectors
  type(spline_type), pointer :: dkb(:,:)  ! derivatives of KB projectors
  type(spline_type), pointer :: Ur(:)   ! atomic wavefunctions
  type(spline_type), pointer :: vps(:)  ! pseudopotential
  type(spline_type) :: vlocal  ! local part
  type(spline_type) :: dvlocal ! derivative of the local part
  type(spline_type) :: core    ! core charge

  integer :: kbc
  integer :: L_max ! maximum value of l to take
  integer :: L_loc ! which component to take as local
  character(len=4) :: icore
  real(r8) :: rc_max
  real(r8) :: vlocal_origin ! local pseudopotential at the orginin

  real(r8), pointer :: dkbcos(:), dknrm(:) ! KB cosinus and norm
  real(r8), pointer :: h(:,:,:)
end type ps_type

real(r8), parameter :: eps = 1.0e-8_r8

contains

subroutine ps_init(ps, label, z, zval, pot)
  type(ps_type), intent(inout) :: ps
  character(len=*), intent(in) :: label
  real(r8), intent(in) :: z
  real(r8), intent(in) :: zval
  character(len=*), intent(in) :: pot

  integer :: i

  sub_name = 'ps_init'; call push_sub()

  ps%l_max = 0 ! In 1D, both thigs are zero (fully local
  ps%l_loc = 0 ! potential)

  !First we allocate all the stuff
  allocate(ps%vps(0:0))
  call spline_init(ps%vps(0))
  call spline_init(ps%vlocal)
  call spline_init(ps%dvlocal)

  ! now we load the necessary information from the ps file
  call ps_load(ps, z, zval, pot)

  call pop_sub()
end subroutine ps_init

subroutine ps_end(ps)
  type(ps_type), intent(inout) :: ps

  integer :: i

  sub_name = 'ps_end'; call push_sub()

  call spline_end(ps%vps(0))
  call spline_end(ps%vlocal)
  call spline_end(ps%dvlocal)
  
  call pop_sub()
end subroutine ps_end

subroutine ps_load(ps, z, zval, pot)
  type(ps_type), intent(inout) :: ps
  real(r8), intent(in) :: z
  real(r8), intent(in) :: zval
  character(len=*), intent(in) :: pot

  integer :: i, ir
  real(r8), allocatable :: hato(:), derhato(:)

  integer :: npotd, npotu, nr
  real(r8) :: b, a
  real(r8), pointer :: rofi(:), vps(:)
  integer :: nrval

  nr = 1001;
  b  = 0.0002253411_r8
  a  = 0.0125000000_r8
  nrval = nr + mod((nr + 1), 2) 

  ps%icore = 'nc'

  allocate(rofi(nrval), vps(nrval))

  do ir = 1, nrval
     rofi(ir) = b*( exp( a*(ir-1) ) - 1 )
     vps(ir) = oct_parse_potential(0.0_r8, 0.0_r8, 0.0_r8, rofi(ir), pot)
  enddo

  allocate(hato(nrval), derhato(nrval))

  hato = 0.0_r8
  hato(1:nrval) = vps(1:nrval)*rofi(1:nrval)  + zval
  call spline_fit(nrval, rofi, hato, ps%vlocal)
  ps%vlocal_origin = vps(1)

  call derivate_in_log_grid(a, b, nrval, hato, derhato)
  call spline_fit(nrval, rofi, derhato, ps%dvlocal)

  deallocate(rofi, vps, hato, derhato)

  return
end subroutine ps_load
