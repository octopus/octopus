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

module output
  use global
  use messages
  use syslabels
  use lib_oct
  use lib_oct_parser
  use io
  use units
  use mesh
  use simul_box
  use mesh_function
  use cube_function
  use functions
#if defined(HAVE_NETCDF)
  use netcdf
#endif
implicit none

  private
  public :: output_type, &
            output_init, &
            output_fill_how, &
            dinput_function, zinput_function, &
            doutput_function, zoutput_function
            

type output_type
  logical :: what(8)
  integer :: how    ! how to output

  integer :: iter   ! output every iter
  logical :: duringscf
  
  integer :: wfs(32) ! which wfs to output
end type output_type

integer, public, parameter :: &
     output_potential  =  1, &
     output_density    =  2, &
     output_wfs        =  3, &
     output_ELF        =  4, &
     output_ELF_FS     =  5, &
     output_geometry   =  6, &
     output_wfs_sqmod  =  7, &   
     output_something  =  8   ! this one should be the last

integer, parameter, private :: &
     output_axis_x     =    1, &
     output_axis_y     =    2, &
     output_axis_z     =    4, &
     output_plane_x    =    8, &
     output_plane_y    =   16, &
     output_plane_z    =   32, &
     output_dx         =   64, &
     output_dx_cdf     =  128, &
     output_plain      =  256, &
     output_mesh_index =  512, &
     output_gnuplot    = 1024


! doutput_kind => real variables; zoutput_kind => complex variables.
integer, parameter, private :: &
     doutput_kind =  1, &
     zoutput_kind = -1

contains

subroutine output_init(outp)
  type(output_type), intent(out) :: outp

  integer :: i
  logical :: l
  character(len=80) :: owf

  call loct_parse_logical(check_inp('OutputKSPotential'), .false., outp%what(output_potential))
  call loct_parse_logical(check_inp('OutputDensity'),     .false., outp%what(output_density))
  call loct_parse_logical(check_inp('OutputWfs'),         .false., outp%what(output_wfs))
  call loct_parse_logical(check_inp('OutputELF'),         .false., outp%what(output_elf))
  call loct_parse_logical(check_inp('OutputELF_FS'),      .false., outp%what(output_elf_FS))
  call loct_parse_logical(check_inp('OutputGeometry'),    .false., outp%what(output_geometry))
  call loct_parse_logical(check_inp('OutputWfsSqMod'),    .false., outp%what(output_wfs_sqmod))
 
  outp%what(output_something) = .false.
  do i = 1, output_something - 1
    outp%what(output_something) = outp%what(output_something).or.outp%what(i)
  end do

  if(outp%what(output_wfs)) then
    call loct_parse_string(check_inp('OutputWfsNumber'), "1-1024", owf)
    call loct_wfs_list(owf, outp%wfs)
  end if

  if(outp%what(output_something)) then
    outp%how = 0
    call loct_parse_logical(check_inp('OutputMeshIndex'), .false., l)
    if(l) outp%how = ior(outp%how, output_mesh_index)
    call loct_parse_logical(check_inp('OutputAxisX'), .false., l)
    if(l) outp%how = ior(outp%how, output_axis_x)
    if(conf%dim > 1) then
      call loct_parse_logical(check_inp('OutputAxisY'), .false., l)
      if(l) outp%how = ior(outp%how, output_axis_y)
      call loct_parse_logical(check_inp('OutputPlaneZ'), .false., l)
      if(l) outp%how = ior(outp%how, output_plane_z)
      call loct_parse_logical(check_inp('OutputGnuplotMode'), .false., l)
      if(l) outp%how = ior(outp%how, output_gnuplot)
      if(conf%dim > 2) then
        call loct_parse_logical(check_inp('OutputAxisZ'), .false., l)
        if(l) outp%how = ior(outp%how, output_axis_z)
        call loct_parse_logical(check_inp('OutputPlaneX'), .false., l)
        if(l) outp%how = ior(outp%how, output_plane_x)
        call loct_parse_logical(check_inp('OutputPlaneY'), .false., l)
        if(l) outp%how = ior(outp%how, output_plane_y)
        call loct_parse_logical(check_inp('OutputDX'), .false., l)
        if(l) outp%how = ior(outp%how, output_dx)
#if defined(HAVE_NETCDF)
        call loct_parse_logical(check_inp('OutputNETCDF'), .false., l)
        if(l) then
          outp%how = ior(outp%how, output_dx_cdf)
        end if
#endif
      end if
    end if

  end if
  
  ! this is always needed in a time-dependent calculation
  call loct_parse_int(check_inp('OutputEvery'), 1000, outp%iter)

  call loct_parse_logical(check_inp('OutputDuringSCF'), .false., outp%duringscf)
end subroutine output_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use this function to quickly plot functions for debugging purposes:
! call doutput_function(output_fill_how("AxisX_and_PlaneX_and_DX", &
!                       ".", "func", m, func, M_ONE, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function output_fill_how(where) result(how)
  character(len=*), intent(in) :: where
  how = 0
  if(index(where, "AxisX").ne.0)     how = ior(how, output_axis_x)
  if(index(where, "AxisY").ne.0)     how = ior(how, output_axis_y)
  if(index(where, "AxisZ").ne.0)     how = ior(how, output_axis_z)
  if(index(where, "PlaneX").ne.0)    how = ior(how, output_plane_x)
  if(index(where, "PlaneY").ne.0)    how = ior(how, output_plane_y)
  if(index(where, "PlaneZ").ne.0)    how = ior(how, output_plane_z)
  if(index(where, "DX").ne.0)        how = ior(how, output_dx)
  if(index(where, "Plain").ne.0)     how = ior(how, output_plain)
  if(index(where, "MeshIndex").ne.0) how = ior(how, output_mesh_index)
  if(index(where, "Gnuplot").ne.0)   how = ior(how, output_gnuplot)
#if defined(HAVE_NETCDF)
  if(index(where, "NETCDF").ne.0) how = ior(how, output_dx_cdf)
#endif
end function output_fill_how

#if defined(HAVE_NETCDF)
subroutine ncdf_error(func, status, filename, ierr)
    character(len=*), intent(in) :: func
    integer, intent(in) :: status
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr

    if(status .eq. NF90_NOERR) return
    write(message(1),'(3a)') "NETCDF error in function '" , trim(func) , "'"
    write(message(2),'(3a)') "(reading/writing ", trim(filename) , ")"
    write(message(3), '(6x,a,a)')'Error code = ', trim(nf90_strerror(status))
    call write_warning(3)
    ierr = 5
end subroutine ncdf_error
#endif


#include "undef.F90"
#include "real.F90"
#include "out_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "out_inc.F90"
#include "undef.F90"

end module output
