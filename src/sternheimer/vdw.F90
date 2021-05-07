!! Copyright (C) 2006 Hyllios
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

module vdw_oct_m
  use em_resp_calc_oct_m
  use gauss_legendre_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use linear_response_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multisystem_basic_oct_m
  use parser_oct_m
  use pert_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use sternheimer_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m

  implicit none

  private
  public :: &
       vdw_run

contains

  ! ---------------------------------------------------------
  subroutine vdw_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(in)    :: from_scratch

    PUSH_SUB(vdw_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = vdw not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call vdw_run_legacy(system, from_scratch)
    end select

    POP_SUB(vdw_run)
  end subroutine vdw_run

  ! ---------------------------------------------------------
  subroutine vdw_run_legacy(sys, fromScratch)
    type(electrons_t),    intent(inout) :: sys
    logical,              intent(in)    :: fromScratch

    type(lr_t) :: lr(MAX_DIM, 1)
    type(sternheimer_t)     :: sh

    integer :: dir, i, iunit, gauss_start, ndir
    CMPLX :: omega
    FLOAT :: domega, pol, c3, c6, cat

    integer :: gaus_leg_n
    FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
    FLOAT, parameter :: omega0 = CNST(0.3)

    type(restart_t) :: restart_dump

    PUSH_SUB(vdw_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if(sys%space%is_periodic()) then
      call messages_not_implemented('Van der Waals calculation for periodic system')
    end if

    if (sys%kpoints%use_symmetries) call messages_experimental("KPoints symmetries with CalculationMode = vdw")

    call input()
    call init_()
    call sternheimer_init(sh, sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ks%xc, sys%mc, wfs_are_cplx = .true.)

    if(gauss_start == 1 .and. mpi_grp_is_root(mpi_world)) then
      iunit = io_open(VDW_DIR//'vdw_c6', sys%namespace, action='write')
      write(iunit, '(a,i3)') '# npoints = ', gaus_leg_n
      write(iunit, '(a1,a12,2a20)') '#', 'omega', 'domega', 'pol'
      call io_close(iunit)
    end if

    do i = gauss_start, gaus_leg_n
      omega  = M_zI*omega0*(M_ONE - gaus_leg_points(i))/(M_ONE + gaus_leg_points(i))
      domega = gaus_leg_weights(i) * omega0 * (M_TWO)/(M_ONE + gaus_leg_points(i))**2

      pol = get_pol(omega)
      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(VDW_DIR//'vdw_c6', sys%namespace, action='write', position='append')
        write(iunit, '(3es20.12)') aimag(omega), domega, pol
        call io_close(iunit)
      end if

      c3  = c3  + M_THREE/M_PI * domega * pol 
      c6  = c6  + M_THREE/M_PI * domega * pol**2
      cat = cat + M_THREE/M_PI * domega * pol**3
    end do

    if((gauss_start  <=  gaus_leg_n).and.mpi_grp_is_root(mpi_world)) then
      iunit = io_open(VDW_DIR//'vdw_c6', sys%namespace, action='write', position='append')
      write(iunit, '(1x)')
      write(iunit, '(a,es20.12)') "C_3  [a.u. ] = ", c3
      write(iunit, '(a,es20.12)') "C_6  [a.u. ] = ", c6
      write(iunit, '(a,es20.12)') "C_AT [a.u. ] = ", cat
      write(iunit, '(1x)')

      write(iunit, '(3a,es20.12)') "C_3  [", &
        trim(units_abbrev(units_out%energy * units_out%length**sys%space%dim)), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**sys%space%dim, c3)
      write(iunit, '(3a,es20.12)') "C_6  [", &
        trim(units_abbrev(units_out%energy * units_out%length**(2*sys%space%dim))), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**(2*sys%space%dim), c6)
      write(iunit, '(3a,es20.12)') "C_AT [", &
        trim(units_abbrev(units_out%energy * units_out%length**(3*sys%space%dim))), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**(3*sys%space%dim), cat)

      call io_close(iunit)
    end if

    call sternheimer_end(sh)
    call end_()

    POP_SUB(vdw_run_legacy)
  contains

    ! --------------------------------------------------------------------
    subroutine input()
      integer :: equiv_axes

      PUSH_SUB(vdw_run_legacy.input)

      !%Variable vdWNPoints
      !%Type integer
      !%Default 6
      !%Section Linear Response::Polarizabilities
      !%Description
      !% How many points to use in the Gauss-Legendre integration to obtain the
      !% van der Waals coefficients.
      !%End
      call messages_obsolete_variable(sys%namespace, 'vdW_npoints', 'vdWNPoints')
      call parse_variable(sys%namespace, 'vdWNPoints', 6, gaus_leg_n)

      ! \todo symmetry stuff should be general
      call parse_variable(sys%namespace, 'TDPolarizationEquivAxes', 0, equiv_axes)

      select case(equiv_axes)
      case(3);      ndir = 1
      case(2);      ndir = min(2, sys%space%dim)
      case default; ndir = min(3, sys%space%dim)
      end select

      POP_SUB(vdw_run_legacy.input)
    end subroutine input


    ! --------------------------------------------------------------------
    subroutine init_()
      integer :: ierr, iunit, ii
      logical :: file_exists
      character(len=80) :: dirname
      FLOAT :: iomega, domega, pol

      type(restart_t) :: restart_load, gs_restart

      PUSH_SUB(vdw_run_legacy.init_)

      ! make some space for static polarizability
      gaus_leg_n = gaus_leg_n + 1

      ! get Gauss-Legendre points
      SAFE_ALLOCATE(gaus_leg_points (1:gaus_leg_n))
      SAFE_ALLOCATE(gaus_leg_weights(1:gaus_leg_n))

      call gauss_legendre_points(gaus_leg_n-1, gaus_leg_points, gaus_leg_weights)
      c3  = M_ZERO
      c6  = M_ZERO
      cat = M_ZERO
      gauss_start = 1
      gaus_leg_points (gaus_leg_n) = CNST(0.99999)
      gaus_leg_weights(gaus_leg_n) = M_ZERO

      ! FIXME: this should be part of the restart framework
      ! check if we can restart
      inquire(file=VDW_DIR//'vdw_c6', exist=file_exists)
      if(.not.fromScratch .and. file_exists) then
        iunit = io_open(VDW_DIR//'vdw_c6', sys%namespace, action='read')
        read(iunit, '(a12,i3)', iostat=ierr) dirname, ii
        if(ii /= gaus_leg_n) then
          message(1) = "Invalid restart of van der Waals calculation."
          message(2) = "The number of points in the Gauss-Legendre integration changed."
          write(message(3), '(i3,a,i3,a)') gaus_leg_n, " (input) != ", ii, "(restart)"
          call messages_fatal(3)
        end if
        read(iunit,*) ! skip comment line
        do
          read(iunit, *, iostat=ierr) iomega, domega, pol
          if(ierr /= 0) exit
          gauss_start = gauss_start + 1
          c3  = c3  + M_THREE/M_PI * domega * pol 
          c6  = c6  + M_THREE/M_PI * domega * pol**2
          cat = cat + M_THREE/M_PI * domega * pol**3
        end do
        call io_close(iunit)
      end if

      ! we always need complex response
      call restart_init(gs_restart, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr == 0) then
        call states_elec_look_and_load(gs_restart, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
          is_complex = .true.)
        call restart_end(gs_restart)
      else
        message(1) = "Previous gs calculation required."
        call messages_fatal(1)
      end if

      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian for linear response.'
      call messages_info(1)
      call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)

      do dir = 1, ndir
        call lr_init(lr(dir,1))
        call lr_allocate(lr(dir,1), sys%st, sys%gr%mesh)
      end do

      ! load wavefunctions
      if (.not. fromScratch) then
        call restart_init(restart_load, sys%namespace, RESTART_VDW, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)

        do dir = 1, ndir
          write(dirname,'(a,i1,a)') "wfs_", dir, "_1_1"
          call restart_open_dir(restart_load, dirname, ierr)
          if (ierr == 0) then
            call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr, lr=lr(dir,1))
          end if
          if(ierr /= 0) then
            message(1) = "Unable to read response wavefunctions from '"//trim(dirname)//"'."
            call messages_warning(1)
          end if
          call restart_close_dir(restart_load)
        end do

        call restart_end(restart_load)
      end if

      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir(VDW_DIR, sys%namespace)               ! output data
      end if

      call restart_init(restart_dump, sys%namespace, RESTART_VDW, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)

      POP_SUB(vdw_run_legacy.init_)
    end subroutine init_

    ! --------------------------------------------------------------------
    subroutine end_()
      integer :: dir

      PUSH_SUB(vdw_run_legacy.end_)

      SAFE_DEALLOCATE_A(gaus_leg_points)
      SAFE_DEALLOCATE_A(gaus_leg_weights)

      do dir = 1, ndir
        call lr_dealloc(lr(dir, 1))
      end do

      call restart_end(restart_dump)

      POP_SUB(vdw_run_legacy.end_)
    end subroutine end_


    ! --------------------------------------------------------------------
    FLOAT function get_pol(omega)
      CMPLX, intent(in) :: omega

      CMPLX        :: alpha(1:MAX_DIM, 1:MAX_DIM)
      type(pert_t) :: perturbation

      PUSH_SUB(vdw_run_legacy.get_pol)

      call pert_init(perturbation, sys%namespace, PERTURBATION_ELECTRIC, sys%gr, sys%ions)
      do dir = 1, ndir
        write(message(1), '(3a,f7.3)') 'Info: Calculating response for the ', index2axis(dir), &
          '-direction and imaginary frequency ', units_from_atomic(units_out%energy, aimag(omega))
        call messages_info(1)   

        call pert_setup_dir(perturbation, dir)
        call zsternheimer_solve(sh, sys%namespace, sys%space, sys%gr, sys%kpoints, sys%st, sys%hm, sys%ks%xc, sys%mc, sys%ions, &
          lr(dir, :), 1, omega, perturbation, restart_dump, em_rho_tag(TOFLOAT(omega),dir), em_wfs_tag(dir,1))
      end do

      call zcalc_polarizability_finite(sys%namespace, sys%space, sys%gr, sys%st, sys%hm, sys%ions, lr(:,:), 1, perturbation, &
        alpha(:,:), ndir = ndir)

      get_pol = M_ZERO
      do dir = 1, ndir
        get_pol = get_pol + TOFLOAT(alpha(dir, dir))
      end do
      do dir = ndir+1, sys%space%dim
        get_pol = get_pol + TOFLOAT(alpha(ndir, ndir))
      end do

      get_pol = get_pol / TOFLOAT(sys%space%dim)

      call pert_end(perturbation)
      POP_SUB(vdw_run_legacy.get_pol)
    end function get_pol

  end subroutine vdw_run_legacy

end module vdw_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
