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
  use par_vec
  use mpi_mod
  use mpi_debug_mod

  implicit none

  private
  public ::                      &
    output_type,                 &
    output_init,                 &
    output_fill_how,             &
    dinput_function,             &
    zinput_function,             &
    doutput_function,            &
    zoutput_function,            &
    iopar_open,                  &
    iopar_close,                 &
    iopar_read,                  &
    iopar_backspace

  integer, parameter, public  :: &
    output_potential  =    1,    &
    output_density    =    2,    &
    output_wfs        =    3,    &
    output_ELF        =    4,    &
    output_ELF_FS     =    5,    &
    output_geometry   =    6,    &
    output_wfs_sqmod  =    7,    &
    output_something  =    8     ! this one should be the last

  integer, parameter, private :: &
    output_axis_x     =    1,    &
    output_axis_y     =    2,    &
    output_axis_z     =    4,    &
    output_plane_x    =    8,    &
    output_plane_y    =   16,    &
    output_plane_z    =   32,    &
    output_dx         =   64,    &
    output_dx_cdf     =  128,    &
    output_plain      =  256,    &
    output_mesh_index =  512,    &
    output_gnuplot    = 1024

  ! doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private :: &
    doutput_kind      =    1,    &
    zoutput_kind      =   -1

  type output_type
    logical :: what(8)
    integer :: how               ! how to output

    integer :: iter              ! output every iter
    logical :: duringscf

    character(len=80) :: wfs_list
  end type output_type


contains

  ! ---------------------------------------------------------
  subroutine output_init(sb, outp)
    type(simul_box_type), intent(in)  :: sb
    type(output_type),    intent(out) :: outp

    integer :: i
    logical :: l

    !%Variable OutputKSPotential
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% Prints out Kohn-Sham potential, separated by parts. File names would be "v0" for 
    !% the local part, "vc" for the classical potential (if it exists), "vh" for the
    !% Hartree potential, and "vxc-x" for each of the exchange and correlation potentials
    !% of a give spin channel, where "x" stands for the spin channel.
    !%End
    call loct_parse_logical(check_inp('OutputKSPotential'), .false., outp%what(output_potential))

    !%Variable OutputDensity
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% Prints out the density. The output file is called "density-i", where "i" stands for 
    !% the spin channel.
    !%End
    call loct_parse_logical(check_inp('OutputDensity'),     .false., outp%what(output_density))

    !%Variable OutputWfs
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% Prints out wave-functions. Which wavefunctions are to be printed is specified
    !% by the variable <tt>OutputWfsNumber</tt> -- see below. The output file is called
    !% "wf-k-p-i", where k stands for the <i>k</i> number, p for the state, and
    !% i for the spin channel.
    !%End
    call loct_parse_logical(check_inp('OutputWfs'),         .false., outp%what(output_wfs))

    !%Variable OutputELF
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% Prints out the electron localization function, ELF. The output file is called
    !% "elf-i", where i stands for the spin channel.
    !%End
    call loct_parse_logical(check_inp('OutputELF'),         .false., outp%what(output_elf))

    call loct_parse_logical(check_inp('OutputELF_FS'),      .false., outp%what(output_elf_FS))

    !%Variable OutputGeometry
    !%Type logical
    !%Default no
    !%Section Output
    !%Description    
    !% If true @code{octopus} outputs a XYZ file called 
    !% "geometry.xyz" containing the coordinates of the atoms
    !% treated within Quantum Mechanics. If point charges were defined
    !% in the PDB file (see <tt>PDBCoordinates</tt>), they will be output
    !% in the file "geometry_classical.xyz".
    !%End
    call loct_parse_logical(check_inp('OutputGeometry'),    .false., outp%what(output_geometry))

    !%Variable OutputWfSqMod
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% Prints out squared module of wave-functions. 
    !% The output file is called "sqm-wf-k-p-i",
    !% where k stands for the <i>k</i> number, p for the state,
    !% and i for the spin channel.
    !%End

    call loct_parse_logical(check_inp('OutputWfsSqMod'),    .false., outp%what(output_wfs_sqmod))

    outp%what(output_something) = .false.
    do i = 1, output_something - 1
      outp%what(output_something) = outp%what(output_something).or.outp%what(i)
    end do

    if(outp%what(output_wfs)) then
      !%Variable OutputWfsNumber
      !%Type string
      !%Default "1-1024"
      !%Section Output
      !%Description
      !% Which wavefunctions to print, in list form, i.e., "1-5" to print the first
      !% five states, "2,3" to print the second and the third state, etc.
      !%End
      call loct_parse_string(check_inp('OutputWfsNumber'), "1-1024", outp%wfs_list)
    end if

    if(outp%what(output_something)) then
      outp%how = 0
      !%Variable OutputMeshIndex
      !%Type logical
      !%Default no
      !%Section Output
      !%Description
      !% Generates output files of a given quantity (Density, Wfs, ...) which include
      !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
      !% useful for small meshes and debugging purposes.
      !% The output can also be used to display the mesh directly. A gnuplot script for mesh vizualization
      !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>
      !%End
      call loct_parse_logical(check_inp('OutputMeshIndex'), .false., l)
      if(l) outp%how = ior(outp%how, output_mesh_index)

      !%Variable OutputAxisX
      !%Type logical
      !%Default no
      !%Section Output
      !%Description
      !% The values of the function on the <math>x</math> axis are printed. The string ".y=0,z=0" is appended
      !% to previous file names.
      !%End

      call loct_parse_logical(check_inp('OutputAxisX'), .false., l)
      if(l) outp%how = ior(outp%how, output_axis_x)
      if(sb%dim > 1) then
        !%Variable OutputAxisY
        !%Type logical
        !%Default no
        !%Section Output
        !%Description
        !% The values of the function on the <math>y</math> axis are printed. The string ".x=0,z=0" is appended
        !% to previous file names.
        !%End
        call loct_parse_logical(check_inp('OutputAxisY'), .false., l)
        if(l) outp%how = ior(outp%how, output_axis_y)

        !%Variable OutputPlaneZ
        !%Type logical
        !%Default no
        !%Section Output
        !%Description
        !% A plane slice at <math>y=0</math> is printed. The string ".z=0" is appended to
        !% previous file names.
        !%End
        call loct_parse_logical(check_inp('OutputPlaneZ'), .false., l)
        if(l) outp%how = ior(outp%how, output_plane_z)
        call loct_parse_logical(check_inp('OutputGnuplotMode'), .false., l)
        if(l) outp%how = ior(outp%how, output_gnuplot)
        if(sb%dim > 2) then

          !%Variable OutputAxisZ
          !%Type logical
          !%Default no
          !%Section Output
          !%Description
          !% The values of the function on the <math>z</math> axis are printed. The string ".x=0,y=0" is appended
          !% to previous file names.
          !%End
          call loct_parse_logical(check_inp('OutputAxisZ'), .false., l)
          if(l) outp%how = ior(outp%how, output_axis_z)

          !%Variable OutputPlaneX
          !%Type logical
          !%Default no
          !%Section Output
          !%Description
          !% A plane slice at @math{x=0} is printed. The string ``.x=0'' is appended
          !% to previous file names.
          !%End
          call loct_parse_logical(check_inp('OutputPlaneX'), .false., l)
          if(l) outp%how = ior(outp%how, output_plane_x)

          !%Variable OutputPlaneY
          !%Type logical
          !%Default no
          !%Section Output
          !%Description
          !% A plane slice at <math>y=0</math> is printed. The string ".y=0" is appended
          !% to previous file names.
          !%End
          call loct_parse_logical(check_inp('OutputPlaneY'), .false., l)
          if(l) outp%how = ior(outp%how, output_plane_y)

          !%Variable OutputDX
          !%Type logical
          !%Default no
          !%Section Output
          !%Description
          !% For printing all the three dimensional information, the open source program
          !% visualization tool OpenDX (http://www.opendx.org/) is used. The string
          !% ".dx" is appended to previous file names.
          !%End
          call loct_parse_logical(check_inp('OutputDX'), .false., l)
          if(l) outp%how = ior(outp%how, output_dx)

#if defined(HAVE_NETCDF)
          !%Variable OutputNETCDF
          !%Type logical
          !%Default no
          !%Section Output
          !%Description
          !% Outputs in NetCDF (http://www.unidata.ucar.edu/packages/netcdf/) format. This file
          !% can then be read, for example, by OpenDX. The string ".ncdf" is appended to previous file names. 
          !%End
          call loct_parse_logical(check_inp('OutputNETCDF'), .false., l)
          if(l) then
            outp%how = ior(outp%how, output_dx_cdf)
          end if
#endif
        end if
      end if

    end if

    !%Variable OutputEvery
    !%Type integer
    !%Default 1000
    !%Section Time Dependent::TD Output
    !%Description
    !% The output is saved when the iteration number is a multiple of the
    !% <tt>OutputEvery</tt> variable.
    !%End
    call loct_parse_int(check_inp('OutputEvery'), 1000, outp%iter)

    call loct_parse_logical(check_inp('OutputDuringSCF'), .false., outp%duringscf)
  end subroutine output_init

  ! -------------------------------------------------------------------
  ! Use this function to quickly plot functions for debugging purposes:
  ! call doutput_function(output_fill_how("AxisX_and_PlaneX_and_DX", &
  !                       ".", "func", m, func, M_ONE, ierr)
  ! -------------------------------------------------------------------
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
  ! ---------------------------------------------------------
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


  ! ---------------------------------------------------------
  integer function iopar_open(m, file, action, status, form, position, die) &
    result(iunit)
    type(mesh_type),  intent(in)           :: m
    character(len=*), intent(in)           :: file, action
    character(len=*), intent(in), optional :: status, form, position
    logical,          intent(in), optional :: die

#if defined(HAVE_MPI)
    integer :: mpi_err
#endif

    call push_sub('out.iopar_open')

    if(mpi_grp_is_root(m%mpi_grp)) then
      iunit = io_open(file, action, status, form, position, die)
    end if

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      call MPI_Bcast(iunit, 1, MPI_INTEGER, m%vp%root, m%vp%comm, mpi_err)
    end if
#endif

    call pop_sub()
  end function iopar_open


  ! ---------------------------------------------------------
  subroutine iopar_close(m, iunit)
    type(mesh_type), intent(in)    :: m
    integer,         intent(inout) :: iunit

#if defined(HAVE_MPI)
    integer :: mpi_err
#endif

    call push_sub('out.iopar_close')

    if(mpi_grp_is_root(m%mpi_grp)) then
      call io_close(iunit)
    end if

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      call MPI_Bcast(iunit, 1, MPI_INTEGER, m%vp%root, m%vp%comm, mpi_err)
    end if
#endif

    call pop_sub()
  end subroutine iopar_close


  ! ---------------------------------------------------------
  subroutine iopar_read(m, iunit, line, ierr)
    type(mesh_type),  intent(in)  :: m
    integer,          intent(in)  :: iunit
    character(len=*), intent(out) :: line
    integer,          intent(out) :: ierr

#if defined(HAVE_MPI)
    integer :: mpi_err
#endif

    call push_sub('out.iopar_read')

    if(mpi_grp_is_root(m%mpi_grp)) then
      read(iunit, '(a)', iostat=ierr) line
    end if

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      call MPI_Bcast(ierr, 1, MPI_INTEGER, m%vp%root, m%vp%comm, mpi_err)
      call MPI_Bcast(line, len(line), MPI_CHARACTER, m%vp%root, m%vp%comm, mpi_err)
    end if
#endif

    call pop_sub()
  end subroutine iopar_read


  ! ---------------------------------------------------------
  subroutine iopar_backspace(m, iunit)
    type(mesh_type),  intent(in)  :: m
    integer,          intent(in)  :: iunit

    call push_sub('out.iopar_read')

    if(mpi_grp_is_root(m%mpi_grp)) then
      backspace(iunit)
    end if

    call pop_sub()

  end subroutine iopar_backspace


#include "undef.F90"
#include "real.F90"
#include "out_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "out_inc.F90"
#include "undef.F90"

end module output
