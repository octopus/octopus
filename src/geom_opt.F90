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

module geom_opt_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use units_m
  use mesh_m
  use external_pot_m
  use v_ks_m
  use hamiltonian_m
  use geometry_m
  use states_m
  use system_m
  use scf_m
  use restart_m
  use varinfo_m

  implicit none

  private
  public :: geom_opt_run

  type geom_opt_t
    integer  :: method
    FLOAT :: step
    FLOAT :: tol
    integer  :: max_iter

    FLOAT :: f
    FLOAT, pointer :: x(:), df(:)

  end type geom_opt_t

contains

  ! ---------------------------------------------------------
  subroutine geom_opt_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(scf_t)               :: scfv
    type(mesh_t),     pointer :: m    ! shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo

    type(geom_opt_t) :: g_opt

    integer :: i, ierr
    FLOAT, allocatable :: x(:)

    call init_()

    ! load wave-functions
    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    call scf_init(sys%gr, scfv, sys%st, h)

    ALLOCATE(x(3*geo%natoms), 3*geo%natoms)
    do i = 0, geo%natoms - 1
      x(3*i + 1) = geo%atom(i + 1)%x(1)
      x(3*i + 2) = geo%atom(i + 1)%x(2)
      x(3*i + 3) = geo%atom(i + 1)%x(3)
    end do

    select case(g_opt%method)
    case(1)
      i = steepest_descents(x)
    end select

    if(i == 0) then
      message(1) = "Info: Minimum found"
      call write_info(1)
    else
      message(1) = "Did not reach the minimum!"
      message(2) = " (the geometry can make some sense though - do not dispair!)"
      call write_warning(2)
    end if

    ! print out geometry
    do i = 0, geo%natoms - 1
      geo%atom(i+1)%x(1) = x(3*i + 1)
      geo%atom(i+1)%x(2) = x(3*i + 2)
      geo%atom(i+1)%x(3) = x(3*i + 3)
    end do
    call atom_write_xyz(".", "min", geo)

    deallocate(x)

    call scf_end(scfv)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('geom_opt.geom_opt_run')

      ! allocate wfs
      ALLOCATE(sys%st%X(psi)(sys%gr%m%np_part, sys%st%d%dim, sys%st%nst, sys%st%d%nik), sys%gr%m%np_part*sys%st%d%dim*sys%st%nst*sys%st%d%nik)

      ! shortcuts
      m   => sys%gr%m
      geo => sys%gr%geo
      st  => sys%st

      !%Variable GOMethod
      !%Type integer
      !%Default steep
      !%Section Geometry Optimization
      !%Description
      !% Method by which the minimization is performed.
      !%Option steep 1
      !% simple steepest descent.
      !%End
      call loct_parse_int(check_inp('GOMethod'), 1, g_opt%method)
      if(.not.varinfo_valid_option('GOMethod', g_opt%method)) call input_error('GOMethod')
      call messages_print_var_option(stdout, "GOMethod", g_opt%method)

      !%Variable GOTolerance
      !%Type float
      !%Default 0.0001 a.u.
      !%Section Geometry Optimization
      !%Description
      !% Convergence criterium to stop the minimization. In units of force; minimization
      !% is stopped when all forces on ions are smaller.
      !%End
      call loct_parse_float(check_inp('GOTolerance'), CNST(0.0001)/units_inp%force%factor, g_opt%tol)
      g_opt%tol = g_opt%tol*units_inp%force%factor
      
      !%Variable GOStep
      !%Type float
      !%Default 0.5
      !%Section Geometry Optimization
      !%Description
      !% Initial step for the geometry optimizer.
      !%End
      call loct_parse_float(check_inp('GOStep'), M_HALF, g_opt%step) ! WARNING: in some weird units

      !%Variable GOMaxIter
      !%Type integer
      !%Default 200
      !%Section Geometry Optimization
      !%Description
      !% Even if previous convergence criterium is not satisfied, minimization will stop
      !% after this number of iterations.
      !%End
      call loct_parse_int(check_inp('GOMaxIter'), 200, g_opt%max_iter)
      if(g_opt%max_iter <= 0) then
        message(1) = "GoMaxIter has to be larger than 0"
        call write_fatal(1)
      end if

      call pop_sub()
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      deallocate(sys%st%X(psi))
    end subroutine end_


    ! ---------------------------------------------------------
    subroutine geom_calc_point(x, f, df)
      FLOAT, intent(in)  :: x(3*geo%natoms)
      FLOAT, intent(out) :: f, df(3*geo%natoms)

      integer :: i

      do i = 0, geo%natoms - 1
        geo%atom(i+1)%x(1) = x(3*i + 1)
        geo%atom(i+1)%x(2) = x(3*i + 2)
        geo%atom(i+1)%x(3) = x(3*i + 3)
      end do
      call atom_write_xyz(".", "work-min", geo)

      call epot_generate(h%ep, m, sys%gr%sb, geo, st, h%reltype)
      call X(states_calc_dens) (st, m%np, st%rho)
      call X(v_ks_calc) (sys%gr, sys%ks, h, st, calc_eigenval=.true.)
      call hamiltonian_energy(h, st, geo%eii, -1)

      ! do scf calculation
      call scf_run(scfv, sys%gr, st, sys%ks, h, sys%outp)

      ! store results
      f = h%etot

      do i = 0, geo%natoms - 1
        df(3*i + 1) = - geo%atom(i+1)%f(1)
        df(3*i + 2) = - geo%atom(i+1)%f(2)
        df(3*i + 3) = - geo%atom(i+1)%f(3)
      end do

    end subroutine geom_calc_point


    ! ---------------------------------------------------------
    integer function steepest_descents(x)
      FLOAT, intent(inout) :: x(:)

      FLOAT, allocatable :: x1(:), df(:), df1(:)
      FLOAT :: f, f1
      integer :: iter, count

      ALLOCATE( x1(3*geo%natoms), 3*geo%natoms)
      ALLOCATE( df(3*geo%natoms), 3*geo%natoms)
      ALLOCATE(df1(3*geo%natoms), 3*geo%natoms)

      count = 0
      steepest_descents = 1

      ! get initial point
      call geom_calc_point(x, f, df)

      do iter = 1, g_opt%max_iter
        x1 = x - g_opt%step * df

        call geom_calc_point(x1, f1, df1)

        if(f1 < f) then
          f = f1; x = x1; df = df1
          g_opt%step = 2*g_opt%step
          count = count + 1

          if(maxval(abs(df)) < g_opt%tol) then
            steepest_descents = 0
            exit
          end if
        else
          ! try with a smaller step
          g_opt%step = g_opt%step/2
          count = count - 1
        end if

        if(count < -5) then
          ! too many subdivisions
          steepest_descents = 2
          exit
        end if

        write(message(1), '(a,i5,a)') "Info: geom_opt (iter = ", iter, ")"
        write(message(2), '(6x,2(a,f16.10))') "energy = ", f/units_out%energy%factor, &
          " max force = ", maxval(abs(df))/units_out%force%factor
        write(message(3), '(6x,2(a,f16.10))') "step   = ", g_opt%step, &
          "       tol = ", g_opt%tol/units_out%force%factor
        call write_info(3)
      end do

      deallocate(x1, df, df1)
    end function steepest_descents

  end subroutine geom_opt_run
end module geom_opt_m
