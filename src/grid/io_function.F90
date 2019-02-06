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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module io_function_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use distributed_oct_m
  use geometry_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use vtk_oct_m

  implicit none

  private
  public ::                       &
    io_function_read_how,         &
    io_function_fill_how,         &
    write_canonicalized_xyz_file, &
    dio_function_input,           &
    zio_function_input,           &
    dio_function_output,          &
    zio_function_output,          &
    io_function_output_vector,    &
    dio_function_output_global,   &
    zio_function_output_global

  !> doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private ::  &
    DOUTPUT_KIND      =    1,     &
    ZOUTPUT_KIND      =   -1

  !> index to label mapping
  character(len=3), parameter ::  &
    index2label(3) = (/ 're ', 'im ', 'abs' /)

  type(profile_t), save :: read_prof, write_prof

  interface io_function_output_vector
    module procedure dio_function_output_vector, zio_function_output_vector
  end interface io_function_output_vector

contains

  ! ---------------------------------------------------------
  subroutine io_function_read_how(sb, how, ignore_error)
    type(simul_box_t), intent(in)  :: sb
    integer(8),        intent(out) :: how
    logical, optional, intent(in)  :: ignore_error !> Ignore error check. Used when called from some external utility.

    PUSH_SUB(io_function_read_how)

    how = 0_8
    
    call messages_obsolete_variable('OutputHow', 'OutputFormat')
    
    !%Variable OutputFormat
    !%Type flag
    !%Default 0
    !%Section Output
    !%Description
    !% Describes the format of the output files (see <tt>Output</tt>).
    !% Example: <tt>axis_x + plane_x + dx</tt>
    !%Option axis_x bit(0) 
    !% The values of the function on the <i>x</i> axis are printed. The string <tt>.y=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_y bit(1)
    !% The values of the function on the <i>y</i> axis are printed. The string <tt>.x=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_z bit(2)
    !% The values of the function on the <i>z</i> axis are printed. The string <tt>.x=0,y=0</tt> is appended
    !% to previous file names.
    !%Option plane_x bit(3)
    !% A plane slice at <i>x</i> = 0 is printed. The string <tt>.x=0</tt> is appended
    !% to previous file names.
    !%Option plane_y bit(4)
    !% A plane slice at <i>y</i> = 0 is printed. The string <tt>.y=0</tt> is appended
    !% to previous file names.
    !%Option plane_z bit(5)
    !% A plane slice at <i>z</i> = 0 is printed. The string <tt>.z=0</tt> is appended to
    !% previous file names.
    !%Option mesh_index bit(8)
    !% Generates output files of a given quantity (density, wavefunctions, ...) which include
    !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
    !% useful for small meshes and debugging purposes.
    !% The output can also be used to display the mesh directly. A Gnuplot script for mesh visualization
    !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>.
    !%Option boundary_points bit(12)
    !% This option includes the output of the mesh enlargement. Default is without.
    !% Supported only by <tt>binary</tt>, <tt>axis</tt>, <tt>plane</tt>, <tt>mesh_index</tt>,
    !% and <tt>matlab</tt> formats.
    !% Not all types of <tt>Output</tt> will have this information available. Not supported when parallel in domains.
    !%Option binary bit(13)
    !% Plain binary, new format.
    !%Option xyz bit(15)
    !% Geometry will be output in XYZ format. Does not affect other outputs.
    !%Option cube bit(16)
    !% Generates output in the <a href=http://paulbourke.net/dataformats/cube>cube file format</a>.
    !% Available only in 3D. Only writes the real part of complex functions.
    !%Option vtk bit(20)
    !% Generates output in <a href=http://www.vtk.org/VTK/img/file-formats.pdf>VTK legacy format</a>.
    !%End
    call parse_variable('OutputFormat', 0, how)
    if(.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
      call messages_input_error('OutputFormat')
    end if

    if(how  ==  0 .and. .not. optional_default(ignore_error, .false.)) then
      write(message(1), '(a)') 'Must specify output method with variable OutputFormat.'
      call messages_fatal(1, only_root_writes = .true.)
     end if

    ! some modes are not available in some circumstances
    if(sb%dim == 1) then
      if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Y) /= 0) then
        message(1) = "OutputFormat = axis_y not available with Dimensions = 1."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Z) /= 0) then
        message(1) = "OutputFormat = plane_z not available with Dimensions = 1."
        call messages_fatal(1)
      end if
    end if

    if(sb%dim <= 2) then
      if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Z) /= 0) then
        message(1) = "OutputFormat = axis_z not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_X) /= 0) then
        message(1) = "OutputFormat = plane_x not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Y) /= 0) then
        message(1) = "OutputFormat = plane_y not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__CUBE) /= 0) then
        message(1) = "OutputFormat = cube not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
    end if

    POP_SUB(io_function_read_how)
  end subroutine io_function_read_how

  ! -------------------------------------------------------------------
  !> Use this function to quickly plot functions for debugging purposes:
  !! call dio_function_output(io_function_fill_how("AxisX_and_PlaneX_and_DX"), &
  !                       ".", "func", mesh, sb, func, M_ONE, ierr)
  ! -------------------------------------------------------------------
  integer(8) function io_function_fill_how(where) result(how)
    character(len=*), intent(in) :: where

    PUSH_SUB(io_function_fill_how)

    how = 0
    if(index(where, "AxisX")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_X)
    if(index(where, "AxisY")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_Y)
    if(index(where, "AxisZ")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_Z)
    IF(index(where, "PlaneX")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_X)
    if(index(where, "PlaneY")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Y)
    if(index(where, "PlaneZ")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Z)
    if(index(where, "PlaneZ")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Z)
    if(index(where, "Binary")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__BINARY)
    if(index(where, "MeshIndex") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__MESH_INDEX)
    if(index(where, "XYZ")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__XYZ)
    if(index(where, "Cube")      /= 0) how = ior(how, OPTION__OUTPUTFORMAT__CUBE)
    if(index(where, "VTK")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__VTK)

    POP_SUB(io_function_fill_how)
  end function io_function_fill_how

  ! ---------------------------------------------------------
  !> Write canonicalized xyz file with atom labels and positions in Angstroms.
  !> Includes information about simulation box and periodicity when applicable.
  !> This differs from a normal xyz file by including information about box
  !> shape and always using Angstroms.
  subroutine write_canonicalized_xyz_file(dir, fname, geo, mesh)
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: fname
    type(geometry_t), intent(in) :: geo
    type(mesh_t),     intent(in) :: mesh

    integer :: iunit
    integer :: idir
    FLOAT :: position
    integer :: iatom

    PUSH_SUB(write_canonicalized_xyz_file)

    call io_mkdir(dir)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xyz', action='write', position='asis')

    write(iunit, '(i6)') geo%natoms
    call simul_box_write_short_info(mesh%sb, iunit)

    ! xyz-style labels and positions:
    do iatom=1, geo%natoms
      write(iunit, '(10a)', advance='no') geo%atom(iatom)%label
      do idir=1, 3
        if(idir <= mesh%sb%dim) then
          position = geo%atom(iatom)%x(idir)
        else
          position = M_ZERO
        end if
        write(iunit, '(xf11.6)', advance='no') units_from_atomic(unit_angstrom, position)
      end do
      write(iunit, '()')
    end do

    call io_close(iunit)

    POP_SUB(write_canonicalized_xyz_file)
  end subroutine write_canonicalized_xyz_file

#include "undef.F90"
#include "real.F90"
#include "io_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "io_function_inc.F90"
#include "undef.F90"

end module io_function_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
