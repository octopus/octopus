#include "config_F90.h"

module geom_opt
use liboct
use io
use units
use system
use hamiltonian
use scf

implicit none

private
public :: geom_opt_run

type geom_opt_type
  integer  :: method
  real(r8) :: tol
  integer  :: max_iter

  real(r8) :: f
  real(r8), pointer :: x(:), df(:)

  type(scf_type),         pointer :: scf
  type(system_type),      pointer :: sys
  type(hamiltonian_type), pointer :: h

end type geom_opt_type

type(geom_opt_type) :: geo

contains

subroutine geom_opt_run(scf, sys, h)
  type(system_type), intent(inout), target :: sys
  type(hamiltonian_type), intent(inout), target :: h
  type(scf_type), intent(inout), target :: scf

  integer :: i
  real(r8), allocatable :: x(:)
  integer, external :: oct_geom_opt

  sub_name = 'geom_opt_run'; call push_sub()
  
  call geo_init()

  allocate(x(3*geo%sys%natoms))
  do i = 0, geo%sys%natoms - 1
    x(3*i + 1) = geo%sys%atom(i + 1)%x(1)
    x(3*i + 2) = geo%sys%atom(i + 1)%x(2)
    x(3*i + 3) = geo%sys%atom(i + 1)%x(3)
  end do  

  i = oct_geom_opt(x, 3*geo%sys%natoms, geo%method, geo%tol, geo%max_iter, &
       geom_calc_point)

  if(i.ne.0) then
    message(1) = "Info: Minimum found"
    call write_info(1)
  else
    message(1) = "Did not reach the minimum!"
    message(2) = " (the geometry can make some sense though - do not dispair!)"
    call write_warning(2)
  end if

  ! print out geometry
  do i = 0, geo%sys%natoms - 1
    sys%atom(i+1)%x(1) = x(3*i + 1)
    sys%atom(i+1)%x(2) = x(3*i + 2)
    sys%atom(i+1)%x(3) = x(3*i + 3)
  end do
  call system_write_xyz("min", sys)
  
  deallocate(x)

  call pop_sub()
  return

  contains
    subroutine geo_init()
      call oct_parse_int(C_string("GOMethod"), 1, geo%method)
      if(geo%method < 1 .or. geo%method >4) then
        message(1) = "'GOMethod' can only take the values:"
        message(5) = "   1 = Steepest descent"
        message(3) = "   2 = Polak-Ribiere conjugate gradient"
        message(2) = "   3 = Fletcher-Reeves conjugate gradient"
        message(4) = "   4 = Broyden-Fletcher-Goldfarb-Shanno conjugate gradient"
        call write_fatal(5)
      end if

      call oct_parse_double(C_string("GOTolerance"), 0.0001_r8/units_inp%force%factor, geo%tol)
      geo%tol = geo%tol*units_inp%force%factor

      call oct_parse_int(C_string("GOMaxIter"), 200, geo%max_iter)
      if(geo%max_iter <= 0) then
        message(1) = "GoMaxIter has to be larger than 0"
        call write_fatal(1)
      end if

      ! store variables in global area
      geo%scf => scf
      geo%sys => sys
      geo%h   => h

    end subroutine geo_init

end subroutine geom_opt_run

subroutine geom_calc_point(x, f, df)
  real(r8), intent(IN)  :: x(3*geo%sys%natoms)
  real(r8), intent(out) :: f, df(3*geo%sys%natoms)

  integer :: i
  integer, save :: iter = 0

  do i = 0, geo%sys%natoms - 1
    geo%sys%atom(i+1)%x(1) = x(3*i + 1)
    geo%sys%atom(i+1)%x(2) = x(3*i + 2)
    geo%sys%atom(i+1)%x(3) = x(3*i + 3)
  end do
  call system_write_xyz("work-min", geo%sys)

  ! generate external potential
  call generate_external_pot(geo%h, geo%sys)
  
  ! setup hamiltonian
  call R_FUNC(calcdens) (geo%sys%st, geo%sys%m%np, geo%sys%st%rho)
  call R_FUNC(hamiltonian_setup)    (geo%h, geo%sys)
  
  ! update hamiltonian and eigenvalues (fermi is *not* called)
  call R_FUNC(hamiltonian_eigenval) (geo%h, geo%sys,  1, geo%sys%st%nst)
  call hamiltonian_energy           (geo%h, geo%sys, -1)
  
  ! do scf calculation
  call scf_run(geo%scf, geo%sys, geo%h)

  ! store results
  f = geo%h%etot

  do i = 0, geo%sys%natoms - 1
    df(3*i + 1) = - geo%sys%atom(i+1)%f(1)
    df(3*i + 2) = - geo%sys%atom(i+1)%f(2)
    df(3*i + 3) = - geo%sys%atom(i+1)%f(3)
  end do

  iter = iter + 1
  write(message(1), '(a,i5,a,f16.10)') "Info: geom_opt iter = ", iter, &
       " energy = ", f/units_out%energy%factor
  call write_info(1)

end subroutine geom_calc_point

end module geom_opt
