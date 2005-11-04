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
!!
!! $Id$

#include "global.h"

module xc
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use lib_basic_alg
  use lib_adv_alg
  use mesh
  use grid
  use functions
  use poisson
  use states
  use lib_xc
  use xc_functl
  use mpi_mod

  implicit none

  private
  public ::             &
    xc_type,            &
    xc_init,            &
    xc_end,             &
    xc_write_info,      &
    xc_get_vxc,         &
    xc_get_vxc_and_axc, &
    xc_get_fxc,         &
    xc_get_kxc


  type xc_type
    logical :: cdft

    integer :: family                   ! the families present
    type(xc_functl_type) :: functl(2,2) ! (1,:) => exchange,    (2,:) => correlation
                                        ! (:,1) => unpolarized, (:,2) => polarized
    type(xc_functl_type) :: j_functl    ! current-depent part of the functional

    ! the meta-GGA can be implemented in two ways
    integer :: mGGA_implementation      ! 1 => as a GGA like functional
    ! 2 => using the OEP method
  end type xc_type


  FLOAT, parameter :: tiny      = CNST(1.0e-12)
  FLOAT, parameter :: denom_eps = CNST(1.0e-20) ! added to denominators to avoid overflows...

contains

  ! ---------------------------------------------------------
  subroutine xc_write_info(xcs, iunit)
    type(xc_type), intent(in) :: xcs
    integer,       intent(in) :: iunit

    integer :: i

    if (xcs%cdft .and. iand(xcs%family, XC_FAMILY_LCA) /= 0) then
      write(iunit,'(a)') " Current-dependent exchange and correlation:"
      call xc_functl_write_info(xcs%j_functl, iunit)

      write(iunit,'(a)') " Auxiliary exchange and correlation functionals:"
    else
      write(iunit,'(a)') " Exchange and correlation:"
    end if

    do i = 1, 2
      call xc_functl_write_info(xcs%functl(i, 1), iunit)
    end do

  end subroutine xc_write_info


  ! ---------------------------------------------------------
  subroutine xc_init(xcs, ndim, spin_channels, cdft)
    type(xc_type), intent(out) :: xcs
    integer,       intent(in)  :: ndim
    integer,       intent(in)  :: spin_channels
    logical,       intent(in)  :: cdft

    integer :: i

    call push_sub('xc.xc_init')

    xcs%cdft   = cdft  ! make a copy of flag indicating the use of current-dft

    ! get current-dependent functional
    call xc_j_functl_init (xcs%j_functl, cdft, spin_channels)
    xcs%family = xcs%j_functl%family

    if (xcs%family == XC_FAMILY_LCA .or. xcs%family == 0) then
      !we also need xc functionals that do not depend on the current
      !get both spin polarized and unpolarized
      do i = 1, 2
        call xc_functl_init_exchange   (xcs%functl(1,i), ndim, i)
        call xc_functl_init_correlation(xcs%functl(2,i), ndim, i)
      end do
      xcs%family = ior(xcs%family, xcs%functl(1,1)%family)
      xcs%family = ior(xcs%family, xcs%functl(2,1)%family)

      if (iand(xcs%family, XC_FAMILY_LCA).ne.0 .and. &
        iand(xcs%family, XC_FAMILY_MGGA + XC_FAMILY_OEP).ne.0) then
        message(1) = "LCA functional can only be used along with LDA or GGA functionals"
        call write_fatal(1)
      end if

      if(iand(xcs%family, XC_FAMILY_MGGA).ne.0) then
        call loct_parse_int(check_inp('MGGAimplementation'), 1, xcs%mGGA_implementation)
        if(xcs%mGGA_implementation.ne.1.and.xcs%mGGA_implementation.ne.2) then
          message(1) = 'MGGAimplementation can only assume the values:'
          message(2) = '  1 : GEA implementation'
          message(3) = '  2 : OEP implementation'
          call write_fatal(3)
        end if
      end if

    end if

    call pop_sub()

  end subroutine xc_init


  ! ---------------------------------------------------------
  subroutine xc_end(xcs)
    type(xc_type), intent(inout) :: xcs

    integer :: i

    if (xcs%cdft) then
      call xc_functl_end(xcs%j_functl)
    end if
    do i = 1, 2
      call xc_functl_end(xcs%functl(1,i))
      call xc_functl_end(xcs%functl(2,i))
    end do
    xcs%family = 0

  end subroutine xc_end


#include "xc_vxc.F90"
#include "xc_axc.F90"
#include "xc_fxc.F90"
#include "xc_kxc.F90"

end module xc
