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

module v_ks
  use global
  use messages
  use syslabels
  use states
  use lib_oct_parser
  use io
  use lib_basic_alg
  use functions
  use mesh
  use grid
  use mesh_function
  use poisson
  use lib_xc
  use xc
  use xc_OEP
  use hamiltonian

  implicit none

  private
  public ::             &
    v_ks_type,          &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    dv_ks_calc,         &
    zv_ks_calc

  integer, parameter :: &
    sic_none   = 1,     &  ! no self interaction correction
    sic_pz     = 2,     &  ! SIC a la Perdew Zunger (OEP way)
    sic_amaldi = 3         ! Amaldi correction term

  type v_ks_type
    logical :: ip_app

    integer           :: xc_family  ! the xc stuff
    integer           :: sic_type   ! what kind of Self Interaction Correction to apply
    type(xc_type)     :: xc
    type(xc_OEP_type) :: oep
  end type v_ks_type


contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(gr, ks, d)
    type(v_ks_type),        intent(out)   :: ks
    type(grid_type),        intent(inout) :: gr
    type(states_dim_type),  pointer       :: d

    call push_sub('v_ks.v_ks_init');

    !%Variable NonInteractingElectrons
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Sometimes it may be helpful to treat the electrons as non-interacting particles,
    !% i.e., not to take into account Hartree and exchange-correlation effects between
    !% the electrons. This variable may be used to toogle this behavior on and off
    !%Option no 0
    !% Electrons are treated as *interacting* particles
    !%Option yes 1
    !% Electrons are handled as *non-interacting* paticles
    !%End
    call loct_parse_logical(check_inp('NonInteractingElectrons'), .false., ks%ip_app)

    if(ks%ip_app) then
      message(1) = 'Info: Treating the electrons as non-interacting'
      call write_info(1)

    else
      ! initilize hartree and xc modules
      message(1) = "Info: Init Hartree"
      call write_info(1)
      call poisson_init(gr)

      message(1) = "Info: Init Exchange-Correlation"
      call write_info(1)
      call xc_init(ks%xc, NDIM, d%spin_channels, d%cdft)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA).ne.0) then

        !%Variable SICorrection
        !%Type integer
        !%default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which Self Interaction Correction to use. Note that
        !% this correction will be applyed to the functional chosen by 'XFunctional' and
        !% 'CFunctional'
        !%Option sic_none 1
        !% No Self Interaction Correction
        !%Option sic_pz 2
        !% SIC a Perdew Zunger, hadled by the OEP technique
        !%Option sic_amaldi 3
        !% Amaldi correction term (NOT WORKING)
        !%End
        call loct_parse_int(check_inp('SICorrection'), sic_none, ks%sic_type)
        if((ks%sic_type.ne.sic_none).and.(ks%sic_type.ne.sic_pz).and.(ks%sic_type.ne.sic_amaldi)) &
          call input_error('SICorrection')

        ! Perdew Zunger corrections
        if(ks%sic_type == sic_pz) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
      end if

      call xc_oep_init(ks%oep, ks%xc_family, gr%m, d)

      call v_ks_write_info(ks, stdout)
    end if

    call pop_sub()
  end subroutine v_ks_init


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks)
    type(v_ks_type), intent(inout) :: ks

    call push_sub('v_ks.v_ks_end');

    if(.not.ks%ip_app) then
      call xc_oep_end(ks%oep)
      call xc_end(ks%xc)

      call poisson_end()
    end if

    call pop_sub();
  end subroutine v_ks_end


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_type), intent(in) :: ks
    integer,         intent(in) :: iunit

    call push_sub('v_ks.v_ks_write_info');

#ifdef HAVE_MPI
    if(mpiv%node == 0) then
#endif
      write(iunit,'(/,a)') stars
      call xc_write_info(ks%xc, iunit)

      select case(ks%sic_type)
      case(sic_pz)
        write(iunit, '(/,2x,a)') 'using Perdew-Zunger SI corrections'
      case(sic_amaldi)
        write(iunit, '(/,2x,a)') 'using Amaldi correction term'
      end select

      if(iand(ks%xc_family, XC_FAMILY_OEP).ne.0) then
        write(iunit, '(1x)')
        call xc_oep_write_info(ks%oep, iunit)
      end if
      write(iunit,'(a,/)') stars

#ifdef HAVE_MPI
    end if
#endif

    call pop_sub()
  end subroutine v_ks_write_info

#include "undef.F90"
#include "real.F90"
#include "v_ks_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "v_ks_inc.F90"

end module v_ks
