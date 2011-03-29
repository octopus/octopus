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

module io_function_m
  use cube_function_m
  use datasets_m
  use geometry_m
  use global_m
  use index_m
  use io_m
  use io_binary_m
  use io_csv_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_debug_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif
  use par_vec_m
  use parser_m
  use profiling_m
  use simul_box_m
  use species_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private
  public ::                       &
    io_function_read_how,         &
    io_function_fill_how,         &
    write_xsf_geometry,           &
    write_xsf_geometry_file,      &
    dio_function_input,           &
    zio_function_input,           &
    dio_function_output,          &
    zio_function_output,          &
    dio_function_out_text,        &
    zio_function_out_text

  integer, parameter, public ::   &
    C_OUTPUT_HOW_AXIS_X          =      1,    &
    C_OUTPUT_HOW_AXIS_Y          =      2,    &
    C_OUTPUT_HOW_AXIS_Z          =      4,    &
    C_OUTPUT_HOW_PLANE_X         =      8,    &
    C_OUTPUT_HOW_PLANE_Y         =     16,    &
    C_OUTPUT_HOW_PLANE_Z         =     32,    &
    C_OUTPUT_HOW_DX              =     64,    &
    C_OUTPUT_HOW_NETCDF          =    128,    &
    C_OUTPUT_HOW_MESH_INDEX      =    512,    &
    C_OUTPUT_HOW_XCRYSDEN        =   1024,    &
    C_OUTPUT_HOW_MATLAB          =   2048,    &
    C_OUTPUT_HOW_MESHGRID        =   4096,    &
    C_OUTPUT_HOW_BOUNDARY_POINTS =   8192,    &
    C_OUTPUT_HOW_BINARY          =  16384,    &
    C_OUTPUT_HOW_ETSF            =  32768,    &
    C_OUTPUT_HOW_XYZ             =  65536,    &
    C_OUTPUT_HOW_CUBE            = 131072

  ! doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private ::  &
    DOUTPUT_KIND      =    1,     &
    ZOUTPUT_KIND      =   -1

  ! index to label mapping
  character(len=3), parameter ::  &
    index2label(3) = (/ 're ', 'im ', 'abs' /)

  type(profile_t), save :: read_prof, write_prof

contains

  ! ---------------------------------------------------------
  subroutine io_function_read_how(sb, how)
    type(simul_box_t), intent(in)  :: sb
    integer,           intent(out) :: how

    PUSH_SUB(io_function_read_how)

    how = 0

    !%Variable OutputHow
    !%Type flag
    !%Default 0
    !%Section Output
    !%Description
    !% Describes the format of the output files (see <tt>Output</tt>).
    !% Example: <tt>axis_x + plane_x + dx</tt>
    !%Option axis_x 1
    !% The values of the function on the <i>x</i> axis are printed. The string <tt>.y=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_y 2
    !% The values of the function on the <i>y</i> axis are printed. The string <tt>.x=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_z 4
    !% The values of the function on the <i>z</i> axis are printed. The string <tt>.x=0,y=0</tt> is appended
    !% to previous file names.
    !%Option plane_x 8
    !% A plane slice at <i>x</i> = 0 is printed. The string <tt>.x=0</tt> is appended
    !% to previous file names.
    !%Option plane_y 16
    !% A plane slice at <i>y</i> = 0 is printed. The string <tt>.y=0</tt> is appended
    !% to previous file names.
    !%Option plane_z 32
    !% A plane slice at <i>z</i> = 0 is printed. The string <tt>.z=0</tt> is appended to
    !% previous file names.
    !%Option dx 64
    !% For printing three-dimensional information, the open-source program
    !% visualization tool OpenDX (<tt>http://www.opendx.org/</tt>) can be used. The string
    !% <tt>.dx</tt> is appended to previous file names.
    !%Option netcdf 128
    !% Outputs in NetCDF (<tt>http://www.unidata.ucar.edu/packages/netcdf/</tt>) format. This file
    !% can then be read, for example, by OpenDX. The string <tt>.ncdf</tt> is appended to previous file names.
    !% Requires the NetCDF library.
    !%Option mesh_index 512
    !% Generates output files of a given quantity (density, wavefunctions, ...) which include
    !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
    !% useful for small meshes and debugging purposes.
    !% The output can also be used to display the mesh directly. A Gnuplot script for mesh vizualization
    !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>.
    !%Option xcrysden 1024
    !% A format for printing structures and three-dimensional information, which can be visualized by
    !% the free open-source program XCrySDen (<tt>http://www.xcrysden.org/</tt>). The string
    !% <tt>.xsf</tt> is appended to previous file names. Note that lattice vectors and coordinates are as
    !% specified by <tt>UnitsOutput</tt>.
    !%Option matlab 2048
    !% In combination with <tt>plane_x</tt>, <tt>plane_y</tt> and
    !% <tt>plane_z</tt>, this option produces output files which are
    !% suitable for 2D Matlab functions like <tt>mesh()</tt>,
    !% <tt>surf()</tt>, or <tt>waterfall()</tt>. To load these files
    !% into Matlab you can use, <i>e.g.</i>
    !%<tt>
    !%   >> density = load('static/density-1.x=0.matlab.abs');
    !%   >> mesh(density);
    !%</tt>
    !%Option meshgrid 4096
    !% Outputs in Matlab mode the internal mesh in a format similar to
    !%<tt>
    !%   >> [x,y] = meshgrid(-2:.2:2,-1:.15:1)
    !%</tt>
    !% The <i>x</i> meshgrid is contained in a file <tt>*.meshgrid.x</tt> and the <i>y</i>-grid can be found in
    !% <tt>*.meshgrid.y</tt>.
    !%Option boundary_points 8192
    !% This option includes the output of the mesh enlargement. Default is without.
    !%Option binary 16384
    !% Plain binary, new format.
    !%Option etsf 32768
    !% ETSF file format (<tt>http://www.etsf.eu/resources/software/standardization_project</tt>).
    !% Requires the ETSF_IO library.
    !%Option xyz 65536
    !% Geometry will be output in XYZ format. Does not affect other outputs.
    !%Option cube 131072
    !% Generates output in the cube file format (<tt>http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/</tt>)
    !%End
    call parse_integer(datasets_check('OutputHow'), 0, how)
    if(.not.varinfo_valid_option('OutputHow', how, is_flag=.true.)) then
      call input_error('OutputHow')
    end if

    if(how .eq. 0) then
      write(message(1), '(a)') 'Must specify output method with variable OutputHow.'
      call messages_fatal(1, only_root_writes = .true.)
     endif

    ! some modes are not available in some circumstances, so we reset how
    if(sb%dim == 1) how = iand(how, not(C_OUTPUT_HOW_AXIS_Y + C_OUTPUT_HOW_PLANE_Z))
    if(sb%dim <= 2) how = iand(how, not(C_OUTPUT_HOW_AXIS_Z + &
      C_OUTPUT_HOW_PLANE_X + C_OUTPUT_HOW_PLANE_Y + C_OUTPUT_HOW_DX + C_OUTPUT_HOW_CUBE))

#if !defined(HAVE_NETCDF)
    if (iand(how, C_OUTPUT_HOW_NETCDF) .ne. 0) then
      message(1) = 'Octopus was compiled without NetCDF support.'
      message(2) = 'It is not possible to write output in NetCDF format.'
      call messages_fatal(2)
    end if
#endif
#if !defined(HAVE_ETSF_IO)
    if (iand(how, C_OUTPUT_HOW_ETSF) .ne. 0) then
      message(1) = 'Octopus was compiled without ETSF_IO support.'
      message(2) = 'It is not possible to write output in ETSF format.'
      call messages_fatal(2)
    end if
#endif

    POP_SUB(io_function_read_how)
  end subroutine io_function_read_how

  ! -------------------------------------------------------------------
  ! Use this function to quickly plot functions for debugging purposes:
  ! call dio_function_output(io_function_fill_how("AxisX_and_PlaneX_and_DX"), &
  !                       ".", "func", mesh, sb, func, M_ONE, ierr)
  ! -------------------------------------------------------------------
  integer function io_function_fill_how(where) result(how)
    character(len=*), intent(in) :: where

    PUSH_SUB(io_function_fill_how)

    how = 0
    if(index(where, "AxisX")     .ne. 0) how = ior(how, C_OUTPUT_HOW_AXIS_X)
    if(index(where, "AxisY")     .ne. 0) how = ior(how, C_OUTPUT_HOW_AXIS_Y)
    if(index(where, "AxisZ")     .ne. 0) how = ior(how, C_OUTPUT_HOW_AXIS_Z)
    IF(INDEX(WHERE, "PlaneX")    .ne. 0) how = ior(how, C_OUTPUT_HOW_PLANE_X)
    if(index(where, "PlaneY")    .ne. 0) how = ior(how, C_OUTPUT_HOW_PLANE_Y)
    if(index(where, "PlaneZ")    .ne. 0) how = ior(how, C_OUTPUT_HOW_PLANE_Z)
    if(index(where, "DX")        .ne. 0) how = ior(how, C_OUTPUT_HOW_DX)
    if(index(where, "XCrySDen")  .ne. 0) how = ior(how, C_OUTPUT_HOW_XCRYSDEN)
    if(index(where, "Binary")    .ne. 0) how = ior(how, C_OUTPUT_HOW_BINARY)
    if(index(where, "MeshIndex") .ne. 0) how = ior(how, C_OUTPUT_HOW_MESH_INDEX)
    if(index(where, "XYZ")       .ne. 0) how = ior(how, C_OUTPUT_HOW_XYZ)
#if defined(HAVE_NETCDF)
    if(index(where, "NETCDF")    .ne. 0) how = ior(how, C_OUTPUT_HOW_NETCDF)
#endif
    if(index(where, "Cube")      .ne. 0) how = ior(how, C_OUTPUT_HOW_CUBE)

    POP_SUB(io_function_fill_how)
  end function io_function_fill_how

  ! ---------------------------------------------------------
  subroutine write_xsf_geometry_file(dir, fname, geo, sb, offset, write_forces)
    character(len=*),   intent(in) :: dir, fname
    type(geometry_t),   intent(in) :: geo
    type(simul_box_t),  intent(in) :: sb
    FLOAT,    optional, intent(in) :: offset(:)
    logical,  optional, intent(in) :: write_forces

    integer iunit
    character(len=6) position
    FLOAT, allocatable:: offset_(:)
    logical write_forces_

    if( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_xsf_geometry_file)

    call io_mkdir(dir)
    position = 'asis'
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xsf', action='write', position=position)

    if(.not. present(write_forces)) then
      write_forces_ = .false.
    else
      write_forces_ = write_forces
    endif

    SAFE_ALLOCATE(offset_(1:sb%dim))
    if(.not. present(offset)) then
      offset_(1:sb%dim) = 0
    else
      ASSERT(ubound(offset, 1) >= sb%dim)
      offset_(1:sb%dim) = offset(1:sb%dim)
    endif

    call write_xsf_geometry(iunit, geo, sb, offset_(:), write_forces_)

    call io_close(iunit)

    POP_SUB(write_xsf_geometry_file)
  end subroutine write_xsf_geometry_file

  ! ---------------------------------------------------------
! for format specification see:
! http://www.xcrysden.org/doc/XSF.html#__toc__11
  subroutine write_xsf_geometry(iunit, geo, sb, offset, write_forces)
    integer,            intent(in) :: iunit
    type(geometry_t),   intent(in) :: geo
    type(simul_box_t),  intent(in) :: sb
    FLOAT,    optional, intent(in) :: offset(:)
    logical,  optional, intent(in) :: write_forces

    integer idir, idir2, iatom
    logical write_forces_

    PUSH_SUB(write_xsf_geometry)

    if(present(write_forces)) then
      write_forces_ = write_forces
    else
      write_forces_ = .false.
    endif

    if(simul_box_is_periodic(sb)) then
      select case(sb%periodic_dim)
        case(3)
          write(iunit, '(a)') 'CRYSTAL'
        case(2)
          write(iunit, '(a)') 'SLAB'
        case(1)
          write(iunit, '(a)') 'POLYMER'
      end select

      write(iunit, '(a)') 'PRIMVEC'

      do idir = 1, sb%dim
        write(iunit, '(3f12.6)') (units_from_atomic(units_out%length, &
          sb%rlattice(idir2, idir)), idir2 = 1, sb%dim)
      enddo

      write(iunit, '(a)') 'PRIMCOORD'
      write(iunit, '(i10, a)') geo%natoms, ' 1'
    else
      write(iunit, '(a)') 'ATOMS'
    endif

    ! BoxOffset should be considered here
    do iatom = 1, geo%natoms
      write(iunit, '(a10, 3f12.6)', advance='no') trim(geo%atom(iatom)%label), &
        (units_from_atomic(units_out%length, geo%atom(iatom)%x(idir) - offset(idir)), idir = 1, sb%dim)
      if(write_forces_) then
        write(iunit, '(5x, 3f12.6)', advance='no') (units_from_atomic(units_out%force, geo%atom(iatom)%f(idir)), idir = 1, sb%dim)
      endif
      write(iunit, '()')
    enddo

    POP_SUB(write_xsf_geometry)
  end subroutine write_xsf_geometry


#if defined(HAVE_NETCDF)
  ! ---------------------------------------------------------
  subroutine ncdf_error(func, status, filename, ierr)
    character(len=*), intent(in)    :: func
    integer,          intent(in)    :: status
    character(len=*), intent(in)    :: filename
    integer,          intent(inout) :: ierr

    PUSH_SUB(ncdf_error)

    if(status .eq. NF90_NOERR) then
    POP_SUB(ncdf_error)
      return
    endif

    write(message(1),'(3a)') "NETCDF error in function '" , trim(func) , "'"
    write(message(2),'(3a)') "(reading/writing ", trim(filename) , ")"
    write(message(3), '(6x,a,a)')'Error code = ', trim(nf90_strerror(status))
    call messages_warning(3)
    ierr = 5

    POP_SUB(ncdf_error)
  end subroutine ncdf_error
#endif

  ! ---------------------------------------------------------
  subroutine transpose3(in, out)
    FLOAT, intent(in)  :: in(:, :, :)
    FLOAT, intent(out) :: out(:, :, :)
    integer :: ix, iy, iz

    PUSH_SUB(transpose3)

    do ix = lbound(in, 1), ubound(in, 1)
      do iy = lbound(in, 2), ubound(in, 2)
        do iz = lbound(in, 3), ubound(in, 3)
          out(iz, iy, ix) = in(ix, iy, iz)
        end do
      end do
    end do

    POP_SUB(transpose3)
  end subroutine transpose3

#include "undef.F90"
#include "real.F90"
#include "io_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "io_function_inc.F90"
#include "undef.F90"

end module io_function_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
