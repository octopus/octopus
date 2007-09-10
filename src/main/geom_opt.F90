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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module geom_opt_m
  use lib_oct_m
  use datasets_m
  use external_pot_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m
  use v_ks_m
  use varinfo_m
  use specie_pot_m
  use lcao_m

  implicit none

  private
  public :: geom_opt_run

  type geom_opt_t
    integer  :: method
    FLOAT :: step
    FLOAT :: tolgrad
    FLOAT :: toldr
    integer  :: max_iter

!!$    FLOAT :: f
!!$    FLOAT, pointer :: x(:), df(:)
  end type geom_opt_t

  integer, parameter ::             &
    MINMETHOD_STEEPEST_DESCENT = 1, &
    MINMETHOD_FR_CG            = 2, &
    MINMETHOD_PR_CG            = 3, &
    MINMETHOD_BFGS             = 4


  type(geometry_t),    pointer :: geo
  type(hamiltonian_t), pointer :: hamilt
  type(system_t),      pointer :: syst
  type(scf_t)                  :: scfv
  type(mesh_t),        pointer :: m    ! shortcuts
  type(states_t),      pointer :: st

  integer                      :: geom_iter

contains

  ! ---------------------------------------------------------
  subroutine geom_opt_run(sys, h)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: h

    type(geom_opt_t) :: g_opt
    integer :: i, ierr, lcao_start, lcao_start_default
    real(8), allocatable :: x(:)
    type(lcao_t) :: lcao_data
    FLOAT :: energy

    call init_()

    ! load wave-functions

    call restart_read(trim(tmpdir)//'gs', sys%st, sys%gr, sys%geo, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)

      ! Randomly generate the initial wave-functions
      call states_generate_random(sys%st, sys%gr%m)

      ! We do not compute the density from the random wave-functions. 
      ! Instead, we try to get a better guess for the density
      call guess_density(sys%gr%m, sys%gr%sb, sys%geo, sys%st%qtot, sys%st%d%nspin, &
           sys%st%d%spin_channels, sys%st%rho)

      ! setup Hamiltonian (we do not call system_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call v_ks_calc(sys%gr, sys%ks, h, sys%st, calc_eigenval=.true.) ! get potentials
      call states_fermi(sys%st, sys%gr%m)                                ! occupations
      call hamiltonian_energy(h, sys%gr, sys%geo, sys%st, -1)             ! total energy

      lcao_start_default = LCAO_START_FULL
      if(sys%geo%only_user_def) lcao_start_default = LCAO_START_NONE

      call loct_parse_int(check_inp('LCAOStart'), lcao_start_default, lcao_start)
      if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
      call messages_print_var_option(stdout, 'LCAOStart', lcao_start)

      if (lcao_start > LCAO_START_NONE) then

        lcao_data%state = 0 ! Uninitialized here.
        call lcao_init(sys%gr, sys%geo, lcao_data, sys%st, h)
        if(lcao_data%state == 1) then
          write(message(1),'(a,i4,a)') 'Info: Performing initial LCAO calculation with ', &
               lcao_data%st%nst,' orbitals.'
          call write_info(1)
          
          call lcao_wf(lcao_data, sys%st, sys%gr, h)
          call lcao_end(lcao_data, sys%st%nst)

          !Just populate again the states, so that the eigenvalues are properly written
          call states_fermi(sys%st, sys%gr%m)
          call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb)

          if (lcao_start == LCAO_START_FULL) then
            ! Update the density and the Hamiltonian
            call system_h_setup(sys, h)
          end if

        end if

      end if

    else

      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call system_h_setup(sys, h)

    end if

    call scf_init(sys%gr, sys%geo, scfv, sys%st, h)

    geom_iter = 0 

    ALLOCATE(x(3*geo%natoms), 3*geo%natoms)
    do i = 0, geo%natoms - 1
      x(3*i + 1) = geo%atom(i + 1)%x(1)
      x(3*i + 2) = geo%atom(i + 1)%x(2)
      x(3*i + 3) = geo%atom(i + 1)%x(3)
    end do

    energy = loct_minimize(g_opt%method, 3*geo%natoms, x(1), &
         real(g_opt%step, 8), real(g_opt%tolgrad, 8), real(g_opt%toldr, 8), g_opt%max_iter, calc_point)

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

      call states_allocate_wfns(sys%st, sys%gr%m)

      ! shortcuts
      m   => sys%gr%m
      geo => sys%geo
      st  => sys%st
      hamilt => h
      syst => sys

      !%Variable GOMethod
      !%Type integer
      !%Default steep
      !%Section Geometry Optimization
      !%Description
      !% Method by which the minimization is performed.
      !%Option steep 1
      !% simple steepest descent.
      !%Option cg_fr 2
      !% Fletcher-Reeves conjugate gradient algorithm. The
      !% conjugate gradient algorithm proceeds as a succession of line
      !% minimizations. The sequence of search directions is used to build
      !% up an approximation to the curvature of the function in the
      !% neighborhood of the minimum. 
      !%Option cg_pr 3
      !% Polak-Ribiere conjugate gradient algorithm.
      !%Option cg_bfgs 4
      !% Vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) conjugate gradient algorithm.
      !% It is a quasi-Newton method which builds up an approximation to the second 
      !% derivatives of the function f using the difference between successive gradient
      !% vectors.  By combining the first and second derivatives the algorithm is able 
      !% to take Newton-type steps towards the function minimum, assuming quadratic 
      !% behavior in that region.
      !%Option cg_bfgs2 5
      !% The bfgs2 version of this minimizer is the most efficient version available, 
      !% and is a faithful implementation of the line minimization scheme described in 
      !% Fletcher's _Practical Methods of Optimization_, Algorithms 2.6.2 and 2.6.4.

      !%End
      call loct_parse_int(check_inp('GOMethod'), MINMETHOD_STEEPEST_DESCENT, g_opt%method)
      if(.not.varinfo_valid_option('GOMethod', g_opt%method)) call input_error('GOMethod')
      call messages_print_var_option(stdout, "GOMethod", g_opt%method)

      !%Variable GOTolerance
      !%Type float
      !%Default 0.0001 a.u.
      !%Section Geometry Optimization
      !%Description
      !% Convergence criterium to stop the minimization. In units of force; minimization
      !% is stopped when all forces on ions are smaller.
      !% Used in conjunction with GOMinimumMove. If GOTolerance = 0, this criterium is ignored.
      !%End
      call loct_parse_float(check_inp('GOTolerance'), CNST(0.0001)/units_inp%force%factor, g_opt%tolgrad)
      g_opt%tolgrad = g_opt%tolgrad*units_inp%force%factor
      
      !%Variable GOMinimumMove
      !%Type float
      !%Default 0.0 a.u.
      !%Section Geometry Optimization
      !%Description
      !% Convergence criterium to stop the minimization. In units of length; minimization
      !% is stopped when all species coordinates change less than GOMinimumMove.
      !% Used in conjunction with GOTolerance. If GOMinimumMove = 0, this criterium is ignored.
      !%End
      call loct_parse_float(check_inp('GOMinimumMove'), CNST(0.0)/units_inp%length%factor, g_opt%toldr)
      g_opt%toldr = g_opt%toldr*units_inp%length%factor
      
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
      call states_deallocate_wfns(sys%st)
      nullify(m)
      nullify(geo)
      nullify(st)
      nullify(hamilt)
      nullify(syst)
    end subroutine end_

  end subroutine geom_opt_run


  subroutine calc_point(n, x, f, getgrad, df)
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    real(8), intent(out) :: f
    integer, intent(in)  :: getgrad
    real(8), intent(out) :: df(n)

    integer :: i
    character(len=256) :: c_geom_iter, title

    do i = 0, geo%natoms - 1
      geo%atom(i+1)%x(1) = x(3*i + 1)
      geo%atom(i+1)%x(2) = x(3*i + 2)
      geo%atom(i+1)%x(3) = x(3*i + 3)
    end do

    call atom_write_xyz(".", "work-min", geo)

    call epot_generate(hamilt%ep, syst%gr, syst%geo, syst%mc, syst%st, hamilt%reltype)
    call states_calc_dens(st, m%np, st%rho)
    call v_ks_calc(syst%gr, syst%ks, hamilt, st, calc_eigenval=.true.)
    call hamiltonian_energy(hamilt, syst%gr, geo, st, -1)

    ! do scf calculation
    call scf_run(scfv, syst%gr, geo, st, syst%ks, hamilt, syst%outp)

    ! store results
    f = hamilt%etot

    if(getgrad .eq. 1) then
      do i = 0, geo%natoms - 1
        df(3*i + 1) = -geo%atom(i+1)%f(1)
        df(3*i + 2) = -geo%atom(i+1)%f(2)
        df(3*i + 3) = -geo%atom(i+1)%f(3)
      end do
    end if

    write(c_geom_iter, '(a,i4.4)') "go.", geom_iter
    write(title, '(f16.10)') f/units_out%energy%factor
    call atom_write_xyz("geom", trim(c_geom_iter), geo, trim(title))
    geom_iter = geom_iter + 1

        write(message(1),'(a)')
        write(message(2), '(6x,2(a,f16.10))')  "Energy = ", f/units_out%energy%factor
        call write_info(2)
        if(getgrad .eq. 1) then
          write(message(1),'(6x,2(a,f16.10))') "Max force = ", maxval(abs(df))/units_out%force%factor
          call write_info(1)
        end if

  end subroutine calc_point

end module geom_opt_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
