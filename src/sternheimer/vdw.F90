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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: em_resp.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module vdw_m
  use datasets_m
  use em_resp_m
  use em_resp_calc_m
  use gauss_legendre_m
  use global_m
  use grid_m
  use output_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use linear_response_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use mpi_m
  use parser_m
  use pert_m
  use poisson_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use unit_m
  use unit_system_m
  use utils_m

  implicit none

  private
  public :: &
       vdw_run

contains

  ! ---------------------------------------------------------
  subroutine vdw_run(sys, hm, fromScratch)
    type(system_t),         intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(lr_t) :: lr(MAX_DIM, 1)
    type(sternheimer_t)     :: sh

    integer :: dir, i, iunit, gauss_start, ndir
    CMPLX :: omega
    FLOAT :: domega, pol, c3, c6, cat

    integer :: gaus_leg_n
    FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
    FLOAT, parameter :: omega0 = CNST(0.3)

    PUSH_SUB(vdw_run)

    if(simul_box_is_periodic(sys%gr%sb)) then
      message(1) = "Van der Waals calculation for periodic system not implemented."
      call messages_fatal(1)
    endif

    call input()
    call init_()
    call sternheimer_init(sh, sys, hm, "Pol", .true.)

    if(gauss_start == 1 .and. mpi_grp_is_root(mpi_world)) then
      iunit = io_open(VDW_DIR//'vdw_c6', action='write')
      write(iunit, '(a,i3)') '# npoints = ', gaus_leg_n
      write(iunit, '(a1,a12,2a20)') '#', 'omega', 'domega', 'pol'
      call io_close(iunit)
    end if

    do i = gauss_start, gaus_leg_n
      omega  = M_zI*omega0*(M_ONE - gaus_leg_points(i))/(M_ONE + gaus_leg_points(i))
      domega = gaus_leg_weights(i) * omega0 * (M_TWO)/(M_ONE + gaus_leg_points(i))**2

      pol = get_pol(omega)
      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(VDW_DIR//'vdw_c6', action='write', position='append')
        write(iunit, '(3es20.12)') aimag(omega), domega, pol
        call io_close(iunit)
      end if

      c3  = c3  + M_THREE/M_PI * domega * pol 
      c6  = c6  + M_THREE/M_PI * domega * pol**2
      cat = cat + M_THREE/M_PI * domega * pol**3
    end do

    if((gauss_start .le. gaus_leg_n).and.mpi_grp_is_root(mpi_world)) then
      iunit = io_open(VDW_DIR//'vdw_c6', action='write', position='append')
      write(iunit, '(1x)')
      write(iunit, '(a,es20.12)') "C_3  [a.u. ] = ", c3
      write(iunit, '(a,es20.12)') "C_6  [a.u. ] = ", c6
      write(iunit, '(a,es20.12)') "C_AT [a.u. ] = ", cat
      write(iunit, '(1x)')

      write(iunit, '(3a,es20.12)') "C_3  [", &
        trim(units_abbrev(units_out%energy * units_out%length**sys%gr%mesh%sb%dim)), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**sys%gr%mesh%sb%dim, c3)
      write(iunit, '(3a,es20.12)') "C_6  [", &
        trim(units_abbrev(units_out%energy * units_out%length**(2*sys%gr%mesh%sb%dim))), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**(2*sys%gr%mesh%sb%dim), c6)
      write(iunit, '(3a,es20.12)') "C_AT [", &
        trim(units_abbrev(units_out%energy * units_out%length**(3*sys%gr%mesh%sb%dim))), "] = ", &
        units_from_atomic(units_out%energy * units_out%length**(3*sys%gr%mesh%sb%dim), cat)

      call io_close(iunit)
    end if

    call sternheimer_end(sh)
    call end_()

    POP_SUB(vdw_run)
  contains

    ! --------------------------------------------------------------------
    subroutine input()
      integer :: equiv_axes

      PUSH_SUB(vdw_run.input)

      !%Variable vdW_npoints
      !%Type integer
      !%Section Linear Response::Polarizabilities
      !%Description
      !% How many points to use in the Gauss-Legendre integration to obtain the
      !% van der Waals coefficients.
      !%End
      call  parse_integer(datasets_check('vdW_npoints'), 6, gaus_leg_n)

      ! \todo symmetry stuff should be general
      call parse_integer(datasets_check('TDPolarizationEquivAxes'), 0, equiv_axes)

      select case(equiv_axes)
      case(3);      ndir = 1
      case(2);      ndir = min(2, sys%gr%mesh%sb%dim)
      case default; ndir = min(3, sys%gr%mesh%sb%dim)
      end select

      POP_SUB(vdw_run.input)
    end subroutine input


    ! --------------------------------------------------------------------
    subroutine init_()
      integer :: ierr, iunit, ii
      logical :: file_exists
      character(len=80) :: dirname
      FLOAT :: iomega, domega, pol

      PUSH_SUB(vdw_run.init_)

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

      ! check if we can restart
      inquire(file=VDW_DIR//'vdw_c6', exist=file_exists)
      if(.not.fromScratch .and. file_exists) then
        iunit = io_open(VDW_DIR//'vdw_c6', action='read')
        read(iunit, '(a12,i3)', iostat=ierr) dirname, ii
        if(ii .ne. gaus_leg_n) then
          message(1) = "Invalid restart of van der Waals calculation."
          message(2) = "The number of points in the Gauss-Legendre integration changed."
          write(message(3), '(i3,a,i3,a)') gaus_leg_n, " (input) != ", ii, "(restart)"
          call messages_fatal(3)
        end if
        read(iunit,*) ! skip comment line
        do
          read(iunit, *, iostat=ierr) iomega, domega, pol
          if(ierr.ne.0) exit
          gauss_start = gauss_start + 1
          c3  = c3  + M_THREE/M_PI * domega * pol 
          c6  = c6  + M_THREE/M_PI * domega * pol**2
          cat = cat + M_THREE/M_PI * domega * pol**3
        end do
        call io_close(iunit)
      end if

      ! we always need complex response
      call restart_look_and_read(sys%st, sys%gr, sys%geo, is_complex = .true.)

      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian for linear response.'
      call messages_info(1)
      call system_h_setup(sys, hm)

      do dir = 1, ndir
        call lr_init(lr(dir,1))
        call lr_allocate(lr(dir,1), sys%st, sys%gr%mesh)

        ! load wavefunctions
        if(.not.fromScratch) then
          write(dirname,'(a,i1,a)') VDW_DIR//"wfs_", dir, "_1_1"
          call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
            ierr, lr=lr(dir,1))
          
          if(ierr.ne.0) then
            message(1) = "Could not load response wavefunctions from '"//trim(tmpdir)//dirname
            call messages_warning(1)
          end if
        end if
      end do

      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir(trim(tmpdir)//VDW_DIR, is_tmp=.true.) ! restart
        call io_mkdir(VDW_DIR)               ! output data
      endif

      POP_SUB(vdw_run.init_)
    end subroutine init_

    ! --------------------------------------------------------------------
    subroutine end_()
      integer :: dir

      PUSH_SUB(vdw_run.end_)

      SAFE_DEALLOCATE_A(gaus_leg_points)
      SAFE_DEALLOCATE_A(gaus_leg_weights)

      do dir = 1, ndir
        call lr_dealloc(lr(dir, 1))
      end do

      POP_SUB(vdw_run.end_)
    end subroutine end_


    ! --------------------------------------------------------------------
    FLOAT function get_pol(omega)
      CMPLX, intent(in) :: omega

      CMPLX        :: alpha(1:MAX_DIM, 1:MAX_DIM)
      type(pert_t) :: perturbation

      PUSH_SUB(vdw_run.get_pol)

      call pert_init(perturbation, PERTURBATION_ELECTRIC, sys%gr, sys%geo)
      do dir = 1, ndir
        write(message(1), '(3a,f7.3)') 'Info: Calculating response for the ', index2axis(dir), &
          '-direction and imaginary frequency ', units_from_atomic(units_out%energy, aimag(omega))
        call messages_info(1)   

        call pert_setup_dir(perturbation, dir)
        call zsternheimer_solve(sh, sys, hm, lr(dir, :), 1,  omega, perturbation, &
             VDW_DIR, em_rho_tag(real(omega),dir), em_wfs_tag(dir,1))
      end do

      call zcalc_polarizability_finite(sys, hm, lr(:,:), 1, perturbation, alpha(:,:), ndir = ndir)

      get_pol = M_ZERO
      do dir = 1, ndir
        get_pol = get_pol + alpha(dir, dir)
      end do
      do dir = ndir+1, sys%gr%mesh%sb%dim
        get_pol = get_pol + alpha(ndir, ndir)
      end do

      get_pol = get_pol / real(sys%gr%mesh%sb%dim)

      call pert_end(perturbation)
      POP_SUB(vdw_run.get_pol)
    end function get_pol

  end subroutine vdw_run

end module vdw_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
