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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: em_resp.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"
#define RESTART_DIR "restart_pol_lr/"

module vdw_m
  use datasets_m
  use em_resp_calc_m
  use gauss_legendre_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lib_basic_alg_m
  use lib_oct_m
  use lib_oct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use mpi_m
  use output_m
  use poisson_m
  use restart_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m

  use pol_lr_m

  implicit none

  private
  public :: &
       vdw_run

contains

  ! ---------------------------------------------------------
  subroutine vdw_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(lr_t) :: lr(MAX_DIM, 1)
    type(sternheimer_t)     :: sh

    integer :: dir, i, iunit, gauss_start
    CMPLX :: omega
    FLOAT :: domega, pol, c3, c6, cat

    integer :: gaus_leg_n
    FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
    FLOAT, parameter :: omega0 = CNST(0.3)

    call push_sub('vdw.vdw_run')

    call sternheimer_init(sh, sys, h, "Pol", hermitian=.false.)
    call input()
    call init()

    if(gauss_start == 1 .and. mpi_grp_is_root(mpi_world)) then
      iunit = io_open('linear/vdw_c6', action='write')
      write(iunit, '(a,i3)') '# npoints = ', gaus_leg_n
      write(iunit, '(a1,a12,2a20)') '#', 'omega', 'domega', 'pol'
      call io_close(iunit)
    end if

    do i = gauss_start, gaus_leg_n
      omega  = M_zI*omega0*(M_ONE - gaus_leg_points(i))/(M_ONE + gaus_leg_points(i))
      domega = gaus_leg_weights(i) * omega0 * (M_TWO)/(M_ONE + gaus_leg_points(i))**2

      pol = get_pol(omega)
      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open('linear/vdw_c6', action='write', position='append')
        write(iunit, '(3es20.12)') aimag(omega), domega, pol
        call io_close(iunit)
      end if

      c3  = c3  + M_THREE/M_PI * domega * pol 
      c6  = c6  + M_THREE/M_PI * domega * pol**2
      cat = cat + M_THREE/M_PI * domega * pol**3
    end do

    if((gauss_start .le. gaus_leg_n).and.mpi_grp_is_root(mpi_world)) then
      iunit = io_open('linear/vdw_c6', action='write', position='append')
      write(iunit, '(1x)')
      write(iunit, '(a,f18.8)') "C_3  [a.u.  ] = ", c3
      write(iunit, '(a,f18.8)') "C_6  [a.u.  ] = ", c6
      write(iunit, '(a,f18.8)') "C_AT [a.u.  ] = ", cat
      write(iunit, '(1x)')

      write(iunit, '(5a,i1,a,f18.8)') "C_3  [", trim(units_out%energy%abbrev), " ",  &
        trim(units_out%length%abbrev), "^", sys%NDIM, "] = ", &
        c3/(units_out%energy%factor * units_out%length%factor**sys%NDIM)
      write(iunit, '(5a,i1,a,f18.8)') "C_6  [", trim(units_out%energy%abbrev), " ",  &
        trim(units_out%length%abbrev), "^", 2*sys%NDIM, "] = ", &
        c6/(units_out%energy%factor * units_out%length%factor**(2*sys%NDIM))
      write(iunit, '(5a,i1,a,f18.8)') "C_AT [", trim(units_out%energy%abbrev), " ",  &
        trim(units_out%length%abbrev), "^", 3*sys%NDIM, "] = ", &
        cat/(units_out%energy%factor * units_out%length%factor**(3*sys%NDIM))

      call io_close(iunit)
    end if

    deallocate(gaus_leg_points)
    deallocate(gaus_leg_weights)
    call sternheimer_end(sh)

    call pop_sub()
  contains

    ! --------------------------------------------------------------------
    subroutine input()
      !%Variable vdW_npoints
      !%Type integers
      !%Section Linear Response::Polarizabilities
      !%Description
      !% How many points to use in the Gauss-Legendre integration to obtain the
      !% van der Waals coefficients
      !%
      !%End
      call  loct_parse_int(check_inp('vdW_npoints'), 6, gaus_leg_n)

    end subroutine input


    ! --------------------------------------------------------------------
    subroutine init()
      integer :: ierr, iunit, ii
      logical :: file_exists
      character(len=80) :: dirname
      FLOAT :: iomega, domega, pol

      ! get gauss legendre points
      ALLOCATE(gaus_leg_points(gaus_leg_n), gaus_leg_n)
      ALLOCATE(gaus_leg_weights(gaus_leg_n), gaus_leg_n)

      call gauss_legendre_points(gaus_leg_n, gaus_leg_points, gaus_leg_weights)
      c3  = M_ZERO
      c6  = M_ZERO
      cat = M_ZERO
      gauss_start = 1

      ! check if we can restart
      inquire(file='linear/vdw_c6', exist=file_exists)
      if(.not.fromScratch .and. file_exists) then
        iunit = io_open('linear/vdw_c6', action='read')
        read(iunit, '(a12,i3)', iostat=ierr) dirname, ii
        if(ii .ne. gaus_leg_n) then
          message(1) = "Invalid restart of van der Waals calculation"
          message(2) = "The number of points in the Gauss-Legendre integration changed"
          write(message(3), '(i3,a,i3,a)') gaus_leg_n, " (input) != ", ii, "(restart)"
          call write_fatal(3)
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
      call read_wfs(sys%st, sys%gr, sys%geo, .true.)

      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian for linear response'
      call write_info(1)
      call system_h_setup(sys, h)

      do dir = 1, sys%NDIM
        call init_lr_wfs(sys%st, sys%gr, lr(dir,1))

        ! load wave-functions
        if(.not.fromScratch) then
          write(dirname,'(a,i1,a)') RESTART_DIR//"wfs", dir, "_1_1"
          call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
            ierr, lr=lr(dir,1))
          
          if(ierr.ne.0) then
            message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//dirname
            call write_warning(1)
          end if
        end if
      end do

      call io_mkdir(trim(tmpdir)//RESTART_DIR)
      call io_mkdir('linear/')
    end subroutine init


    ! --------------------------------------------------------------------
    FLOAT function get_pol(omega)
      CMPLX, intent(in) :: omega

      CMPLX :: alpha(1:MAX_DIM, 1:MAX_DIM)

      do dir = 1, sys%NDIM
        write(message(1), '(a,i1,a,f7.3)') 'Info: Calculating response for direction ', dir, &
          ' and imaginary frequency ', aimag(omega)/units_out%energy%factor
        call write_info(1)   

        call zsternheimer_solve(sh, sys, h, lr(dir, :), dir,  omega, sys%gr%m%x(:,dir), &
             RESTART_DIR, em_rho_tag(real(omega),dir), em_wfs_tag(dir,1))
      end do

      call zlr_calc_polarizability(sys, lr(:,:), alpha(:,:))

      get_pol = M_ZERO
      do dir = 1, sys%NDIM
        get_pol = get_pol + alpha(dir, dir)
      end do
      get_pol = get_pol / real(sys%NDIM)

    end function get_pol

  end subroutine vdw_run

end module vdw_m
