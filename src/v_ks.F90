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

module v_ks_m
  use global_m
  use mpi_m
  use messages_m
  use datasets_m
  use profiling_m
  use states_m
  use lib_oct_parser_m
  use io_m
  use lib_basic_alg_m
  use functions_m
  use mesh_m
  use grid_m
  use mesh_function_m
  use poisson_m
  use lib_xc_m
  use xc_m
  use xc_OEP_m
  use hamiltonian_m
  use varinfo_m
  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    dv_ks_calc,         &
    zv_ks_calc

  integer, parameter :: &
    sic_none   = 1,     &  ! no self interaction correction
    sic_pz     = 2,     &  ! SIC a la Perdew Zunger (OEP way)
    sic_amaldi = 3         ! Amaldi correction term

  type v_ks_t
    logical :: ip_app

    integer           :: xc_family  ! the xc stuff
    integer           :: sic_type   ! what kind of Self Interaction Correction to apply
    type(xc_t)     :: xc
    type(xc_OEP_t) :: oep
  end type v_ks_t


contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(gr, ks, d)
    type(v_ks_t),        intent(out)   :: ks
    type(grid_t),        intent(inout) :: gr
    type(states_dim_t),  pointer       :: d

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
      call messages_print_stress(stdout, "Hartree")
      call poisson_init(gr)
      call messages_print_stress(stdout)

      call xc_init(ks%xc, NDIM, d%spin_channels, d%cdft)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA).ne.0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
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
        call loct_parse_int(check_inp('SICCorrection'), sic_none, ks%sic_type)
        if(.not.varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')

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
    type(v_ks_t), intent(inout) :: ks

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
    type(v_ks_t), intent(in) :: ks
    integer,         intent(in) :: iunit

    call push_sub('v_ks.v_ks_write_info');

    if(mpi_grp_is_root(mpi_world)) then
      call messages_print_stress(iunit, "Exchange-Correlation")
      call xc_write_info(ks%xc, iunit)

      write(iunit, '(1x)')
      call messages_print_var_option(iunit, 'SICCorrection', ks%sic_type)

      if(iand(ks%xc_family, XC_FAMILY_OEP).ne.0) then
        write(iunit, '(1x)')
        call xc_oep_write_info(ks%oep, iunit)
      end if
      call messages_print_stress(iunit)
    end if

    call pop_sub()
  end subroutine v_ks_write_info

#include "undef.F90"
#include "real.F90"
#include "v_ks_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "v_ks_inc.F90"

end module v_ks_m
