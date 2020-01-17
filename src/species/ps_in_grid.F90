!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module ps_in_grid_oct_m
  use atomic_oct_m
  use global_oct_m
  use logrid_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                   &
    ps_in_grid_t,             &
    ps_in_grid_init,          &
    ps_in_grid_end,           &
    ps_in_grid_vlocal,        &
    ps_in_grid_kb_projectors, &
    ps_in_grid_kb_cosines,    &
    ps_in_grid_cutoff_radii,  &
    ps_in_grid_check_rphi,    &
    first_point_extrapolate

  type ps_in_grid_t
    ! Components are public by default
    type(logrid_t) :: g                !< log grid where the pseudos are defined

    FLOAT          :: zval             !< valence charge

    integer        :: no_l_channels    !< number of l channels to consider
    FLOAT, pointer :: vps(:, :)        !< the pseudopotential (l=0 .. no_l_channels-1)
    FLOAT, pointer :: KB(:,:)          !< Kleinman-Bylander projectors
    FLOAT, pointer :: dkbcos(:)        !< Kleinman-Bylander cosine
    FLOAT, pointer, private :: dknorm(:)        !< Kleinman-Bylander norm
    FLOAT, pointer :: kb_radius(:)     !< radius of KB projectors

    integer        :: so_no_l_channels
    FLOAT, pointer :: so_vps(:,:)      !< spin-orbit components (l=1 .. so_no_l_channels)
    FLOAT, pointer, private :: so_KB(:,:)       !< Kleinman-Bylander projectors
    FLOAT, pointer, private :: so_dkbcos(:)     !< Kleinman-Bylander cosine
    FLOAT, pointer, private :: so_dknorm(:)     !< Kleinman-Bylander norm
    FLOAT, pointer, private :: so_kb_radius(:)  !< radius of KB projectors

    FLOAT, pointer :: vlocal(:)        !< local part of the pseudopotential
    FLOAT, pointer :: rphi(:, :,:)     !< pseudo wavefunctions

    logical        :: core_corrections
    FLOAT, pointer :: chcore(:)        !< core charge density

  end type ps_in_grid_t

contains

  subroutine ps_in_grid_init(ps, flavor, a, b, nrval, no_l, so_no_l)
    type(ps_in_grid_t), intent(out) :: ps
    integer,            intent(in)  :: flavor, nrval
    FLOAT,              intent(in)  :: a, b
    integer,            intent(in)  :: no_l, so_no_l

    PUSH_SUB(ps_in_grid_init)

    ! initialize logarithmic grid
    call logrid_init(ps%g, flavor, a, b, nrval)

    ! copy data
    ps%   no_l_channels =    no_l
    ps%so_no_l_channels = so_no_l

    ! Allocate some stuff
    SAFE_ALLOCATE(ps%vps    (1:ps%g%nrval, 1:no_l))
    SAFE_ALLOCATE(ps%chcore (1:ps%g%nrval))
    SAFE_ALLOCATE(ps%vlocal (1:ps%g%nrval))

    SAFE_ALLOCATE(ps%rphi  (1:ps%g%nrval, 1:no_l, 1:3))
    SAFE_ALLOCATE(ps%KB    (1:ps%g%nrval, 1:no_l))
    SAFE_ALLOCATE(ps%dkbcos(1:no_l))
    SAFE_ALLOCATE(ps%dknorm(1:no_l))
    SAFE_ALLOCATE(ps%kb_radius(1:no_l+1))

    if(so_no_l > 0) then
      SAFE_ALLOCATE(ps%so_vps(1:ps%g%nrval, 1:so_no_l))
      SAFE_ALLOCATE(ps%so_KB (1:ps%g%nrval, 1:so_no_l))
      SAFE_ALLOCATE(ps%so_dkbcos(1:so_no_l))
      SAFE_ALLOCATE(ps%so_dknorm(1:so_no_l))
      SAFE_ALLOCATE(ps%so_kb_radius(1:so_no_l))
    end if

    POP_SUB(ps_in_grid_init)
  end subroutine ps_in_grid_init


  ! ---------------------------------------------------------
  subroutine ps_in_grid_end(ps)
    type(ps_in_grid_t), intent(inout) :: ps

    PUSH_SUB(ps_in_grid_end)

    SAFE_DEALLOCATE_P(ps%vps)
    SAFE_DEALLOCATE_P(ps%chcore)
    SAFE_DEALLOCATE_P(ps%vlocal)

    SAFE_DEALLOCATE_P(ps%rphi)
    SAFE_DEALLOCATE_P(ps%KB)
    SAFE_DEALLOCATE_P(ps%dkbcos)
    SAFE_DEALLOCATE_P(ps%dknorm)
    SAFE_DEALLOCATE_P(ps%kb_radius)

    if(ps%so_no_l_channels > 0) then
      SAFE_DEALLOCATE_P(ps%so_vps)
      SAFE_DEALLOCATE_P(ps%so_KB)
      SAFE_DEALLOCATE_P(ps%so_dkbcos)
      SAFE_DEALLOCATE_P(ps%so_dknorm)
      SAFE_DEALLOCATE_P(ps%so_kb_radius)
    end if

    call logrid_end(ps%g)

    POP_SUB(ps_in_grid_end)
  end subroutine ps_in_grid_end


  ! ---------------------------------------------------------
  subroutine ps_in_grid_vlocal(ps, l_loc, rcore, namespace)
    type(ps_in_grid_t), intent(inout) :: ps
    integer,            intent(in)    :: l_loc
    FLOAT,              intent(in)    :: rcore
    type(namespace_t),  intent(in)    :: namespace

    integer :: ir
    FLOAT :: a, b, qtot
    FLOAT, allocatable :: rho(:)

    PUSH_SUB(ps_in_grid_vlocal)

    if(l_loc >= 0) then
      ps%vlocal(:) = ps%vps(:, l_loc+1)

    else if(l_loc == -1) then
      if(ps%g%flavor /= LOGRID_PSF) then
        message(1) = "For the moment, Vanderbilt local potentials are only possible with tm grids."
        call messages_fatal(1, namespace=namespace)
      end if

      a = CNST(1.82) / rcore
      b = M_ONE
      SAFE_ALLOCATE(rho(1:ps%g%nrval))

      do ir = 1, ps%g%nrval
        rho(ir) = exp( -( sinh(a*b*ps%g%rofi(ir)) / sinh(b) )**2 )
        rho(ir) = M_FOUR * M_PI * rho(ir) * ps%g%rofi(ir)**2
      end do
!      do ir =1,2
      qtot = sum(rho(:)*ps%g%drdi(:))
      rho(:) = - rho(:)*(ps%zval/qtot)

      call vhrtre(rho, ps%vlocal, ps%g%rofi, ps%g%drdi, ps%g%s, ps%g%nrval, ps%g%a)
      ps%vlocal(1) = ps%vlocal(2)

      SAFE_DEALLOCATE_A(rho)
    end if

    POP_SUB(ps_in_grid_vlocal)
  end subroutine ps_in_grid_vlocal


  ! ---------------------------------------------------------
  !> KB-projectors
  !!   kb = (vps - vlocal) |phi> * dknorm
  subroutine ps_in_grid_kb_projectors(ps)
    type(ps_in_grid_t), intent(inout) :: ps

    integer :: l

    PUSH_SUB(ps_in_grid_kb_projectors)

    do l = 1, ps%no_l_channels
      ps%KB(2:, l) = (ps%vps(2:, l) - ps%vlocal(2:))*(ps%rphi(2:, l, 1)/ps%g%rofi(2:))*ps%dknorm(l)

      ps%KB(1, l) = first_point_extrapolate(ps%g%rofi, ps%KB(:, l))
    end do

    do l = 1, ps%so_no_l_channels
      ps%so_KB(2:, l) = ps%so_vps(2:, l)*(ps%rphi(2:, l, 1)/ps%g%rofi(2:))*ps%dknorm(l)

      ps%so_KB(1, l) = first_point_extrapolate(ps%g%rofi, ps%so_KB(:, l))
    end do

    POP_SUB(ps_in_grid_kb_projectors)
  end subroutine ps_in_grid_kb_projectors


  ! ---------------------------------------------------------
  !> KB-cosines and KB-norms:
  !!       dkbcos stores the KB "cosines:"
  !!               || (v_l - v_local) phi_l ||^2 / < (v_l - v_local)phi_l | phi_l >  [Rydberg]
  !!       dknorm stores the KB "norms:"
  !!               1 / || (v_l - v_local) phi_l || [1/Rydberg]
  subroutine ps_in_grid_kb_cosines(ps, lloc)
    type(ps_in_grid_t), intent(inout) :: ps
    integer,            intent(in)    :: lloc

    integer :: ir, l
    FLOAT :: dnrm, avgv, vphi

    PUSH_SUB(ps_in_grid_kb_cosines)

    do l = 1, ps%no_l_channels
      if(l-1 == lloc) then
        ps%dkbcos(l) = M_ZERO
        ps%dknorm(l) = M_ZERO
        cycle
      end if

      dnrm = M_ZERO
      avgv = M_ZERO
      do ir = 1, ps%g%nrval
        vphi = (ps%vps(ir, l) - ps%vlocal(ir))*ps%rphi(ir, l, 1)
        dnrm = dnrm + vphi*vphi*ps%g%drdi(ir)
        avgv = avgv + vphi*ps%rphi(ir, l, 1)*ps%g%drdi(ir)
      end do
      ps%dkbcos(l) = dnrm/(avgv + CNST(1.0e-20))
      ps%dknorm(l) = M_ONE/(sqrt(dnrm) + CNST(1.0e-20))
    end do

    do l = 1, ps%so_no_l_channels
      dnrm = M_ZERO
      avgv = M_ZERO
      do ir = 1, ps%g%nrval
        vphi = ps%so_vps(ir, l)*ps%rphi(ir, l, 1)
        dnrm = dnrm + vphi*vphi*ps%g%drdi(ir)
        avgv = avgv + vphi*ps%rphi(ir, l, 1)*ps%g%drdi(ir)
      end do
      ps%so_dkbcos(l) = dnrm/(avgv + CNST(1.0e-20))
      ps%so_dknorm(l) = M_ONE/(sqrt(dnrm) + CNST(1.0e-20))
    end do

    POP_SUB(ps_in_grid_kb_cosines)
  end subroutine ps_in_grid_kb_cosines


  ! ---------------------------------------------------------
  subroutine ps_in_grid_cutoff_radii(ps, lloc)
    type(ps_in_grid_t), intent(inout) :: ps
    integer,            intent(in)    :: lloc

    integer          :: l, ir
    FLOAT            :: dincv, phi
    FLOAT, parameter :: threshold = CNST(1.0e-6)

    PUSH_SUB(ps_in_grid_cutoff_radii)

    ! local part ....
    do ir = ps%g%nrval-1, 2, -1
      dincv = abs(ps%vlocal(ir)*ps%g%rofi(ir) + M_TWO*ps%zval)
      if(dincv > threshold) exit
    end do
    ps%kb_radius(ps%no_l_channels+1) = ps%g%rofi(ir + 1)

    ! non-local part....
    do l = 1, ps%no_l_channels
      if(l-1 == lloc) then
        ps%kb_radius(l) = M_ZERO
        cycle
      end if

      do ir = ps%g%nrval-1, 2, -1
        phi = (ps%rphi(ir, l, 1)/ps%g%rofi(ir))*ps%dknorm(l)
        dincv = abs((ps%vps(ir, l) - ps%vlocal(ir))*phi)

        if(dincv > threshold) exit

        phi = (ps%rphi(ir, l, 1)/ps%g%rofi(ir))*ps%dknorm(l)
      end do

      ps%kb_radius(l) = ps%g%rofi(ir + 1)
    end do

    ! now the SO part
    do l = 1, ps%so_no_l_channels
      do ir = ps%g%nrval-1, 2, -1
        phi = (ps%rphi(ir, l, 1)/ps%g%rofi(ir))*ps%so_dknorm(l)
        dincv = abs((ps%so_vps(ir, l))*phi)

        if(dincv > threshold) exit
      end do

      ps%so_kb_radius(l) = ps%g%rofi(ir + 1)
    end do

    POP_SUB(ps_in_grid_cutoff_radii)
  end subroutine ps_in_grid_cutoff_radii


  ! ---------------------------------------------------------
  !> checks normalization of the pseudo wavefunctions
  subroutine ps_in_grid_check_rphi(ps, namespace)
    type(ps_in_grid_t), intent(in) :: ps
    type(namespace_t),  intent(in) :: namespace

    integer :: l
    FLOAT   :: nrm

    PUSH_SUB(ps_in_grid_check_rphi)

    !  checking normalization of the wavefunctions
    do l = 1, ps%no_l_channels
      nrm = sqrt(sum(ps%g%drdi(:)*ps%rphi(:, l, 1)**2))
      nrm = abs(nrm - M_ONE)
      if (nrm > CNST(1.0e-5)) then
        write(message(1), '(a,i2,a)') "Eigenstate for l = ", l-1, ' is not normalized.'
        write(message(2), '(a, f12.6,a)') '(abs(1 - norm) = ', nrm, ')'
        call messages_warning(2, namespace=namespace)
      end if
    end do

    POP_SUB(ps_in_grid_check_rphi)
  end subroutine ps_in_grid_check_rphi

  ! ---------------------------------------------------------

  FLOAT function first_point_extrapolate(x, y, high_order) result(y0)
    FLOAT,             intent(in)    :: x(:)
    FLOAT,             intent(in)    :: y(:)
    logical, optional, intent(in)    :: high_order

    FLOAT :: x1, x2, x3
    FLOAT :: y1, y2, y3

    x1 = x(2) - x(1)
    x2 = x(3) - x(1)
    x3 = x(4) - x(1)
    y1 = y(2)
    y2 = y(3)
    y3 = y(4)

    if(optional_default(high_order, .false.)) then

      y0 = y1*x2*x3*(x2 - x3) + y2*x1*x3*(x3 - x1) + y3*x1*x2*(x1 - x2);
      y0 = y0/((x1 - x2)*(x1 - x3)*(x2 - x3));

    else

      y0 = y1 - (y2 - y1)*x1/(x2 - x1)

    end if

  end function first_point_extrapolate

end module ps_in_grid_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
