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

module geom_opt
  use global
  use mesh
  use hamiltonian
  use geometry
  use states
  use scf
  
  implicit none

  private
  public :: geom_opt_run
  
  type geom_opt_type
    integer  :: method
    FLOAT :: step
    FLOAT :: tol
    integer  :: max_iter
    
    FLOAT :: f
    FLOAT, pointer :: x(:), df(:)
    
  end type geom_opt_type

contains

  subroutine geom_opt_run(scf, m, f_der, st, geo, h, outp)
    type(scf_type),         intent(inout) :: scf
    type(mesh_type),        intent(IN)    :: m
    type(f_der_type),       intent(inout) :: f_der
    type(states_type),      intent(inout) :: st
    type(geometry_type),    intent(inout) :: geo
    type(hamiltonian_type), intent(inout) :: h
    type(output_type),      intent(IN)    :: outp
    
    type(geom_opt_type) :: g_opt

    integer :: i
    FLOAT, allocatable :: x(:)
    !integer, external :: loct_geom_opt
    
    call push_sub('geom_opt_run')
 
    call geo_init()

    allocate(x(3*geo%natoms))
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
    
    call pop_sub()
    return
    
  contains
    subroutine geo_init()
      call loct_parse_int("GOMethod", 1, g_opt%method)
      if(g_opt%method < 1 .or. g_opt%method >1) then
        message(1) = "'GOMethod' can only take the values:"
        message(2) = "   1 = Steepest descent"
        call write_fatal(2)
      end if

      call loct_parse_float("GOTolerance", CNST(0.0001)/units_inp%force%factor, g_opt%tol)
      g_opt%tol = g_opt%tol*units_inp%force%factor

      ! WARNING: in some weird units
      call loct_parse_float("GOStep", M_HALF, g_opt%step)

      call loct_parse_int("GOMaxIter", 200, g_opt%max_iter)
      if(g_opt%max_iter <= 0) then
        message(1) = "GoMaxIter has to be larger than 0"
        call write_fatal(1)
      end if

    end subroutine geo_init

    subroutine geom_calc_point(x, f, df)
      FLOAT, intent(IN)  :: x(3*geo%natoms)
      FLOAT, intent(out) :: f, df(3*geo%natoms)
      
      integer :: i
      
      do i = 0, geo%natoms - 1
        geo%atom(i+1)%x(1) = x(3*i + 1)
        geo%atom(i+1)%x(2) = x(3*i + 2)
        geo%atom(i+1)%x(3) = x(3*i + 3)
      end do
      call atom_write_xyz(".", "work-min", geo)

      call epot_generate(h%ep, m, st, geo, h%Vpsl, h%reltype)
      call X(calcdens) (st, m%np, st%rho)
      call X(h_calc_vhxc) (h, m, f_der, st, calc_eigenval=.true.)
      call hamiltonian_energy(h, st, geo%eii, -1)
  
      ! do scf calculation
      call scf_run(scf, m, f_der, st, geo, h, outp)

      ! store results
      f = h%etot

      do i = 0, geo%natoms - 1
        df(3*i + 1) = - geo%atom(i+1)%f(1)
        df(3*i + 2) = - geo%atom(i+1)%f(2)
        df(3*i + 3) = - geo%atom(i+1)%f(3)
      end do

    end subroutine geom_calc_point
    
    integer function steepest_descents(x)
      FLOAT, intent(inout) :: x(:)

      FLOAT, allocatable :: x1(:), df(:), df1(:)
      FLOAT :: f, f1
      integer :: iter, count

      allocate(x1(3*geo%natoms), df(3*geo%natoms), df1(3*geo%natoms))

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
end module geom_opt
