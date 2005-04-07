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

#include "global.h"

module v_ks
  use global
  use states
  use lib_basic_alg
  use functions
  use mesh
  use mesh_function
  use geometry
  use poisson
  use lib_xc
  use xc
  use xc_OEP
  use hamiltonian

  implicit none

  private
  public ::              &
     v_ks_type,          &
     v_ks_init,          &
     v_ks_end,           &
     v_ks_write_info,    &
     dh_calc_vhxc,       &
     zh_calc_vhxc

  type v_ks_type
    logical :: ip_app

    integer           :: xc_family  ! the xc stuff
    type(xc_type)     :: xc
    type(xc_OEP_type) :: oep
  end type v_ks_type

contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, m, geo, d)
    type(v_ks_type), intent(out) :: ks
    type(mesh_type),        intent(inout) :: m
    type(geometry_type),    intent(in)    :: geo
    type(states_dim_type),  pointer       :: d

    call push_sub('v_ks_init');

    ! Should we treat the particles as independent?
    call loct_parse_logical("NonInteractingElectrons", .false., ks%ip_app)
    if(ks%ip_app) then
      message(1) = 'Info: Treating the electrons as non-interacting'
      call write_info(1)

    else
      ! initilize hartree and xc modules
      message(1) = "Info: Init Hartree"
      call write_info(1)
      call poisson_init(m)

      call xc_init(ks%xc, geo%nlcc, d%spin_channels, d%cdft)
      call xc_oep_init(ks%oep, ks%xc%family, d)

      if(conf%verbose >= VERBOSE_NORMAL) call xc_write_info(ks%xc, stdout)
    end if

    call pop_sub()
  end subroutine v_ks_init


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks)
    type(v_ks_type), intent(inout) :: ks

    call push_sub('v_ks_end');

    if(.not.ks%ip_app) then
      call xc_end(ks%xc)
      call poisson_end()
    end if

    call pop_sub();
  end subroutine v_ks_end


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_type), intent(in) :: ks
    integer,         intent(in) :: iunit 

    call push_sub('v_ks_write_info');

    call xc_write_info(ks%xc, iunit)
    if(iand(ks%xc_family, XC_FAMILY_OEP).ne.0) then
      call xc_oep_write_info(ks%oep, iunit)
    end if

    call pop_sub()
  end subroutine v_ks_write_info

#include "undef.F90"
#include "real.F90"
#include "v_ks_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "v_ks_inc.F90"

end module v_ks
