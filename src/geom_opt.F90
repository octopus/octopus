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

#include "config_F90.h"

module geom_opt
  use liboct
  use system
  use hamiltonian
  use scf
  
  implicit none

  private
  public :: geom_opt_run
  
  type geom_opt_type
    integer  :: method
    real(r8) :: step
    real(r8) :: tol
    integer  :: max_iter
    
    real(r8) :: f
    real(r8), pointer :: x(:), df(:)
    
  end type geom_opt_type

contains
  
  subroutine geom_opt_run(scf, sys, h)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(scf_type), intent(inout) :: scf
    
    type(geom_opt_type) :: geo

    integer :: i
    real(r8), allocatable :: x(:)
    integer, external :: oct_geom_opt
    
    sub_name = 'geom_opt_run'; call push_sub()
    
    call geo_init()

    allocate(x(3*sys%natoms))
    do i = 0, sys%natoms - 1
      x(3*i + 1) = sys%atom(i + 1)%x(1)
      x(3*i + 2) = sys%atom(i + 1)%x(2)
      x(3*i + 3) = sys%atom(i + 1)%x(3)
    end do

    select case(geo%method)
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
    do i = 0, sys%natoms - 1
      sys%atom(i+1)%x(1) = x(3*i + 1)
      sys%atom(i+1)%x(2) = x(3*i + 2)
      sys%atom(i+1)%x(3) = x(3*i + 3)
    end do
    call system_write_xyz(".", "min", sys)
    
    deallocate(x)
    
    call pop_sub()
    return
    
  contains
    subroutine geo_init()
      call oct_parse_int(C_string("GOMethod"), 1, geo%method)
      if(geo%method < 1 .or. geo%method >1) then
        message(1) = "'GOMethod' can only take the values:"
        message(2) = "   1 = Steepest descent"
        call write_fatal(2)
      end if

      call oct_parse_double(C_string("GOTolerance"), 0.0001_r8/units_inp%force%factor, geo%tol)
      geo%tol = geo%tol*units_inp%force%factor

      ! WARNING: in some weird units
      call oct_parse_double(C_string("GOStep"), 0.5_r8, geo%step)

      call oct_parse_int(C_string("GOMaxIter"), 200, geo%max_iter)
      if(geo%max_iter <= 0) then
        message(1) = "GoMaxIter has to be larger than 0"
        call write_fatal(1)
      end if

    end subroutine geo_init

    subroutine geom_calc_point(x, f, df)
      real(r8), intent(IN)  :: x(3*sys%natoms)
      real(r8), intent(out) :: f, df(3*sys%natoms)
      
      integer :: i
      
      do i = 0, sys%natoms - 1
        sys%atom(i+1)%x(1) = x(3*i + 1)
        sys%atom(i+1)%x(2) = x(3*i + 2)
        sys%atom(i+1)%x(3) = x(3*i + 3)
      end do
      call system_write_xyz(".", "work-min", sys)

      call generate_external_pot(h, sys)
      call R_FUNC(calcdens) (sys%st, sys%m%np, sys%st%rho)
      call R_FUNC(hamiltonian_setup) (h, sys%m, sys%st, sys)
      call hamiltonian_energy(h, sys%st, sys%eii, -1)
  
      ! do scf calculation
      call scf_run(scf, sys, h)

      ! store results
      f = h%etot

      do i = 0, sys%natoms - 1
        df(3*i + 1) = - sys%atom(i+1)%f(1)
        df(3*i + 2) = - sys%atom(i+1)%f(2)
        df(3*i + 3) = - sys%atom(i+1)%f(3)
      end do

    end subroutine geom_calc_point
    
    integer function steepest_descents(x)
      real(r8), intent(inout) :: x(:)

      real(r8), allocatable :: x1(:), df(:), df1(:)
      real(r8) :: f, f1
      integer :: iter, count

      allocate(x1(3*sys%natoms), df(3*sys%natoms), df1(3*sys%natoms))

      count = 0
      steepest_descents = 1

      ! get initial point
      call geom_calc_point(x, f, df)

      do iter = 1, geo%max_iter
        x1 = x - geo%step * df

        call geom_calc_point(x1, f1, df1)

        if(f1 < f) then
          f = f1; x = x1; df = df1
          geo%step = 2*geo%step
          count = count + 1

          if(maxval(abs(df)) < geo%tol) then
            steepest_descents = 0
            exit
          end if
        else
          ! try with a smaller step
          geo%step = geo%step/2
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
        write(message(3), '(6x,2(a,f16.10))') "step   = ", geo%step, &
             "       tol = ", geo%tol/units_out%force%factor
        call write_info(3)
      end do
       
      deallocate(x1, df, df1)
    end function steepest_descents

  end subroutine geom_opt_run
end module geom_opt
