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
  use box_hypercube_oct_m
  use box_parallelepiped_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use distributed_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use io_csv_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use namespace_oct_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use vtk_oct_m

  implicit none

  private
  public ::                         &
    io_function_read_what_how_when, &
    io_function_fill_how,           &
    write_bild_forces_file,         &
    write_canonicalized_xyz_file,   &
    write_xsf_geometry,             &
    write_xsf_geometry_file,        &
    dio_function_input,             &
    zio_function_input,             &
    dio_function_output,            &
    zio_function_output,            &
    io_function_output_vector,      &
    dio_function_output_global,     &
    zio_function_output_global,     &
    io_function_output_vector_BZ,   &
    io_function_output_global_BZ

#if defined(HAVE_NETCDF)
 public ::                          &
    dout_cf_netcdf,                 &
    zout_cf_netcdf
#endif

  !> doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private ::    &
    DOUTPUT_KIND      =    1,       &
    ZOUTPUT_KIND      =   -1

  !> index to label mapping
  character(len=3), parameter ::    &
    index2label(3) = (/ 're ', 'im ', 'abs' /)

  type(profile_t), save :: read_prof, write_prof

  interface io_function_output_vector
    module procedure dio_function_output_vector, zio_function_output_vector
  end interface io_function_output_vector

  interface io_function_output_vector_BZ
    module procedure dio_function_output_vector_BZ, zio_function_output_vector_BZ
  end interface io_function_output_vector_BZ

  interface io_function_output_global_BZ
    module procedure dio_function_output_global_BZ, zio_function_output_global_BZ
  end interface io_function_output_global_BZ


contains

  ! ---------------------------------------------------------
  subroutine io_function_read_what_how_when(namespace, space, what, how, output_interval, &
    what_tag_in, how_tag_in, output_interval_tag_in, ignore_error)
    type(namespace_t), intent(in)           :: namespace
    type(space_t),     intent(in)           :: space
    logical,           intent(inout)        :: what(MAX_OUTPUT_TYPES)
    integer(8),        intent(out)          :: how(0:MAX_OUTPUT_TYPES)
    integer,           intent(out)          :: output_interval(0:MAX_OUTPUT_TYPES) 
    character(len=*),  optional, intent(in) :: what_tag_in
    character(len=*),  optional, intent(in) :: how_tag_in
    character(len=*),  optional, intent(in) :: output_interval_tag_in
    logical, optional, intent(in)           :: ignore_error !> Ignore error check. Used when called from some external utility.

    type(block_t) :: blk
    integer :: ncols, nrows, iout, column_index, what_i
    integer(8) :: what_no_how(10)       !> these kinds of Output do not have a how
    character(len=80) :: what_tag
    character(len=80) :: how_tag
    character(len=80) :: output_interval_tag
    character(len=80) :: output_column_marker
    PUSH_SUB(io_function_read_what_how_when)

    how = 0_8
    output_interval = 0_8
    ncols = 0
    what_tag = optional_default(what_tag_in, 'Output')
    how_tag = optional_default(how_tag_in, 'OutputFormat')
    output_interval_tag = optional_default(output_interval_tag_in, 'OutputInterval')
    
    call messages_obsolete_variable(namespace, 'OutputHow', 'OutputFormat')
    call messages_obsolete_variable(namespace, 'Output_KPT', 'Output')
    call messages_obsolete_variable(namespace, 'OutputLDA_U' , 'Output')
    call messages_obsolete_variable(namespace, 'OutputEvery', 'OutputInterval/RestartWriteInterval')
    
    !%Variable Output
    !%Type block
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. 
    !% Each output must be in a separate row. Optionally individual output formats and output intervals can be defined
    !% for each row or they can be read separately from <tt>OutputFormat</tt> and <tt>OutputInterval</tt> variables
    !% in the input file.
    !% The output files are written at the end of the run into the output directory for the
    !% relevant kind of run (<i>e.g.</i> <tt>static</tt> for <tt>CalculationMode = gs</tt>).
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% (available for <tt>CalculationMode</tt> = <tt>gs</tt>, <tt>unocc</tt>,  <tt>td</tt>, and <tt>opt_control</tt>)
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir</tt>.
    !% For linear-response run modes, the derivatives of many quantities can be printed, as listed in
    !% the options below. Indices in the filename are labelled as follows:
    !% <tt>sp</tt> = spin (or spinor component), <tt>k</tt> = <i>k</i>-point, <tt>st</tt> = state/band.
    !% There is no tag for directions, given as a letter. The perturbation direction is always
    !% the last direction for linear-response quantities, and a following +/- indicates the sign of the frequency.
    !%
    !% Example (minimal):
    !% <br><br><tt>%Output
    !% <br>&nbsp;&nbsp;density
    !% <br>&nbsp;&nbsp;potential
    !% <br>%<br></tt>
    !%
    !% Example (with OutputFormat):
    !% <br><br><tt>%Output
    !% <br>&nbsp;&nbsp;density   | cube + axis_z
    !% <br>&nbsp;&nbsp;potential | cube
    !% <br>%<br></tt>
    !%
    !% Example (with OutputFormat, incomplete):
    !% <br><br><tt>%Output
    !% <br>&nbsp;&nbsp;density   | cube + axis_z
    !% <br>&nbsp;&nbsp;potential
    !% <br>%<br></tt>
    !%
    !% Example (tagged):
    !% <br><br><tt>%Output
    !% <br>&nbsp;&nbsp;density   | "output_format" | cube + axis_z | "output_interval" | 50
    !% <br>&nbsp;&nbsp;potential | "output_format" | cube          | "output_interval" | 20
    !% <br>%<br></tt>
    !%
    !% Example (tagged, incomplete):
    !% <br><br><tt>%Output
    !% <br>&nbsp;&nbsp;density   | "output_format"   | cube + axis_z 
    !% <br>&nbsp;&nbsp;potential | "output_interval" | 20
    !% <br>%<br></tt>
    !% Missing information for the incomplete blocks will be parsed form the out-of-block
    !% definitions. It is also possible to mix the order of columns in the tagged format. 
    !% See <tt>OutputFormat</tt>, and <tt>OutputInterval</tt>.
    !%Option potential 1
    !% Outputs Kohn-Sham potential, separated by parts. File names are <tt>v0</tt> for 
    !% the local part of the ionic potential, <tt>vc</tt> for the classical potential (if it exists),
    !% <tt>vh</tt> for the Hartree potential, <tt>vks</tt> for the local part of the Kohn-Sham potential, and
    !% <tt>vxc-</tt> for the exchange-correlation potentials. For <tt>vks</tt> and <tt>vxc</tt>,
    !% a suffix for spin is added in the spin-polarized case.
    !%Option density 2
    !% Outputs density. The output file is called <tt>density-</tt>, or <tt>lr_density-</tt> in linear response.
    !%Option wfs 3
    !% Outputs wavefunctions. Which wavefunctions are to be printed is specified
    !% by the variable <tt>OutputWfsNumber</tt> -- see below. The output file is called
    !% <tt>wf-</tt>, or <tt>lr_wf-</tt> in linear response.
    !%Option wfs_sqmod 4
    !% Outputs modulus squared of the wavefunctions. 
    !% The output file is called <tt>sqm-wf-</tt>. For linear response, the filename is <tt>sqm_lr_wf-</tt>.
    !%Option geometry 5
    !% Outputs file containing the coordinates of the atoms treated within quantum mechanics.
    !% If <tt>OutputFormat = xyz</tt>, the file is called <tt>geometry.xyz</tt>; a
    !% file <tt>crystal.xyz</tt> is written with a supercell geometry if the system is periodic;
    !% if point charges were defined in the PDB file (see <tt>PDBCoordinates</tt>), they will be output
    !% in the file <tt>geometry_classical.xyz</tt>.
    !% If <tt>OutputFormat = xcrysden</tt>, a file called <tt>geometry.xsf</tt> is written.
    !%Option current 6
    !% Outputs the total current density. The output file is called <tt>current-</tt>.
    !% For linear response, the filename is <tt>lr_current-</tt>.
    !%Option ELF 7
    !% Outputs electron localization function (ELF). The output file is called <tt>elf-</tt>,
    !% or <tt>lr_elf-</tt> in linear response, in which case the associated function D is also written,
    !% as <tt>lr_elf_D-</tt>. Only in 2D and 3D.
    !%Option ELF_basins 8
    !% Outputs basins of attraction of the ELF. The output file is called
    !% <tt>elf_rs_basins.info</tt>. Only in 2D and 3D.
    !%Option Bader 9
    !% Outputs Laplacian of the density which shows lone pairs, bonded charge concentrations
    !% and regions subject to electrophilic or nucleophilic attack.
    !% See RF Bader, <i>Atoms in Molecules: A Quantum Theory</i> (Oxford Univ. Press, Oxford, 1990).
    !%Option el_pressure 10
    !% Outputs electronic pressure. See Tao, Vignale, and Tokatly, <i>Phys Rev Lett</i> <b>100</b>, 206405 (2008).
    !%Option matrix_elements 11
    !% Outputs a series of matrix elements of the Kohn-Sham states. What is output can
    !% be controlled by the <tt>OutputMatrixElements</tt> variable.
    !%Option pol_density 12
    !% Outputs dipole-moment density <tt>dipole_density-</tt>, or polarizability density <tt>alpha_density-</tt>
    !% in linear response. If <tt>ResponseMethod = finite_differences</tt>, the hyperpolarizability density
    !% <tt>beta_density-</tt> is also printed.
    !%Option mesh_r 13
    !% Outputs values of the coordinates over the grid. Files
    !% will be called <tt>mesh_r-</tt> followed by the direction.
    !%Option kinetic_energy_density 14
    !% Outputs kinetic-energy density, defined as:
    !%
    !% <math>\tau_\sigma(\vec{r}) = \sum_{i=1}^{N_\sigma} 
    !%  \left| \vec{\nabla} \phi_{i\sigma}(\vec{r}) \right|^2\,. </math>
    !%
    !% The index <math>\sigma</math> is the spin index for the spin-polarized case,
    !% or if you are using spinors. For spin-unpolarized calculations, you
    !% get the total kinetic-energy density. The previous expression assumes full 
    !% or null occupations. If fractional occupation numbers, each term in the sum
    !% is weighted by the occupation. Also, if we are working with an infinite 
    !% system, all <i>k</i>-points are summed up, with their corresponding weights. The
    !% files will be called <tt>tau-sp1</tt> and <tt>tau-sp2</tt>, if the spin-resolved kinetic
    !% energy density is produced (runs in spin-polarized and spinors mode), or
    !% only <tt>tau</tt> if the run is in spin-unpolarized mode.
    !%Option dos 15
    !% Outputs density of states. See <tt>DOSEnergyMax</tt>, <tt>DOSEnergyMin</tt>, <tt>DOSEnergyPoints</tt>,
    !% and <tt>DOSGamma</tt>.
    !%Option tpa 16
    !% Outputs transition-potential approximation (TPA) matrix elements, using <math>\vec{q}</math>-vector specified
    !% by <tt>MomentumTransfer</tt>.
    !%Option forces 17
    !% Outputs file <tt>forces.xsf</tt> containing structure and forces on the atoms as 
    !% a vector associated with each atom, which can be visualized with XCrySDen.
    !%Option wfs_fourier 18
    !% (Experimental) Outputs wavefunctions in Fourier space. This is
    !% only implemented for the ETSF file format output. The file will
    !% be called <tt>wfs-pw-etsf.nc</tt>.  
    !%Option xc_density 19
    !% Outputs the XC density, which is the charge density that
    !% generates the XC potential. (This is <math>-1/4\pi</math> times
    !% the Laplacian of the XC potential). The files are called <tt>nxc</tt>.
    !%Option PES_wfs 20
    !% Outputs the photoelectron wavefunctions. The file name is <tt>pes_wfs-</tt>  
    !% plus the orbital number.
    !%Option PES_density 21
    !% Outputs the photolectron density. Output file is <tt>pes_dens-</tt> plus spin species if
    !% spin-polarized calculation is performed. 
    !%Option PES 22
    !% Outputs the time-dependent photoelectron spectrum.
    !%Option BerkeleyGW 23
    !% Output for a run with <a href=http://www.berkeleygw.org>BerkeleyGW</a>.
    !% See <tt>Output::BerkeleyGW</tt> for further specification.
    !%Option delta_perturbation 24
    !% Outputs the "kick", or time-delta perturbation applied to compute optical response in real time.
    !%Option external_td_potential 25
    !% Outputs the (scalar) time-dependent potential.
    !%Option mmb_wfs 26
    !% Triggers the ModelMB wavefunctions to be output for each state.
    !%Option mmb_den 27
    !% Triggers the ModelMB density matrix to be output for each state, and the particles
    !% specified by the <tt>DensitytoCalc</tt> block. Calculates, and outputs, the reduced density
    !% matrix. For the moment the trace is made over the second dimension, and
    !% the code is limited to 2D. The idea is to model <i>N</i> particles in 1D as an
    !% <i>N</i>-dimensional non-interacting problem, then to trace out <i>N</i>-1 coordinates.
    !%Option potential_gradient 28
    !% Prints the gradient of the potential.
    !%Option energy_density 29
    !% Outputs the total energy density to a file called
    !% <tt>energy_density</tt>.
    !%Option heat_current 30
    !% Outputs the total heat current density. The output file is
    !% called <tt>heat_current-</tt>.
    !%Option photon_correlator 31
    !% Outputs the electron-photon correlation function. The output file is
    !% called <tt>photon_correlator</tt>.
    !%Option J_flow 32
    !% todo: document J_flow option!
    !%Option current_kpt 33
    !% Outputs the current density resolved in momentum space. The output file is called <tt>current_kpt-</tt>.
    !%Option density_kpt 34
    !% Outputs the electronic density resolved in momentum space. 
    !%Option occ_matrices 35
    !% Outputs the occupation matrices of LDA+U
    !%Option effectiveU 36
    !% Outputs the value of the effectiveU for each atoms 
    !%Option magnetization 37
    !% Outputs file containing structure and magnetization of the localized subspace 
    !% on the atoms as a vector associated with each atom, which can be visualized.
    !% For the moment, it only works if a +U is added on one type of orbital per atom. 
    !%Option local_orbitals 38
    !% Outputs the localized orbitals that form the correlated subspace
    !%Option kanamoriU 39
    !% Outputs the Kanamori interaction parameters U, U`, and J.
    !% These parameters are not determined self-consistently, but are taken from the 
    !% occupation matrices and Coulomb integrals comming from a standard +U calculation.
    !%Option xc_torque 40
    !% Outputs the exchange-correlation torque. Only for the spinor case and in the 3D case.
    !%End

    !%Variable OutputFormat
    !%Type flag
    !%Default 0
    !%Section Output
    !%Description
    !% Describes the format of the output files.
    !% This variable can also be defined inside the <tt>Output</tt> block.
    !% See <tt>Output</tt>.
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
    !%Option dx bit(6)
    !% For printing three-dimensional information, the open-source program
    !% visualization tool <a href=http://www.opendx.org>OpenDX</a> can be used. The string
    !% <tt>.dx</tt> is appended to previous file names. Available only in 3D.
    !%Option netcdf bit(7)
    !% Outputs in <a href=http://www.unidata.ucar.edu/packages/netcdf>NetCDF</a> format. This file
    !% can then be read, for example, by OpenDX. The string <tt>.ncdf</tt> is appended to previous file names.
    !% Requires the NetCDF library. Only writes the real part of complex functions.
    !%Option mesh_index bit(8)
    !% Generates output files of a given quantity (density, wavefunctions, ...) which include
    !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
    !% useful for small meshes and debugging purposes.
    !% The output can also be used to display the mesh directly. A Gnuplot script for mesh visualization
    !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>.
    !%Option xcrysden bit(9)
    !% A format for printing structures and three-dimensional information, which can be visualized by
    !% the free open-source program <a href=http://www.xcrysden.org>XCrySDen</a> and others. The string
    !% <tt>.xsf</tt> is appended to previous file names. Note that lattice vectors and coordinates are as
    !% specified by <tt>UnitsOutput</tt>. Available in 2D and 3D.
    !%Option matlab bit(10)
    !% In combination with <tt>plane_x</tt>, <tt>plane_y</tt> and
    !% <tt>plane_z</tt>, this option produces output files which are
    !% suitable for 2D Matlab functions like <tt>mesh()</tt>,
    !% <tt>surf()</tt>, or <tt>waterfall()</tt>. To load these files
    !% into Matlab you can use, <i>e.g.</i>
    !%<tt>
    !%   >> density = load('static/density-1.x=0.matlab.abs');
    !%   >> mesh(density);
    !%</tt>
    !%Option meshgrid bit(11)
    !% Outputs in Matlab mode the internal mesh in a format similar to
    !%<tt>
    !%   >> [x,y] = meshgrid(-2:.2:2,-1:.15:1)
    !%</tt>
    !% The <i>x</i> meshgrid is contained in a file <tt>*.meshgrid.x</tt> and the <i>y</i>-grid can be found in
    !% <tt>*.meshgrid.y</tt>.
    !%Option boundary_points bit(12)
    !% This option includes the output of the mesh enlargement. Default is without.
    !% Supported only by <tt>binary</tt>, <tt>axis</tt>, <tt>plane</tt>, <tt>mesh_index</tt>,
    !% and <tt>matlab</tt> formats.
    !% Not all types of <tt>Output</tt> will have this information available. Not supported when parallel in domains.
    !%Option binary bit(13)
    !% Plain binary, new format.
    !%Option etsf bit(14)
    !% <a href=http://www.etsf.eu/resources/software/standardization_project>ETSF file format</a>.
    !% Requires the ETSF_IO library. Applies only to <tt>Output = density</tt>, <tt>geometry</tt>,
    !% <tt>wfs</tt>, and/or <tt>wfs_fourier</tt>.
    !%Option xyz bit(15)
    !% Geometry will be output in XYZ format. Does not affect other outputs.
    !%Option cube bit(16)
    !% Generates output in the <a href=http://paulbourke.net/dataformats/cube>cube file format</a>.
    !% Available only in 3D. Only writes the real part of complex functions.
    !%Option bild bit(19)
    !% Generates output in <a href=http://plato.cgl.ucsf.edu/chimera/docs/UsersGuide/bild.html>BILD format</a>.
    !%Option vtk bit(20)
    !% Generates output in <a href=http://www.vtk.org/VTK/img/file-formats.pdf>VTK legacy format</a>.
    !%Option integrate_xy bit(21)
    !% Integrates the function in the x-y plane and the result on the <i>z</i> axis is printed.
    !%Option integrate_xz bit(22)
    !% Integrates the function in the x-z plane and the result on the <i>y</i> axis is printed
    !%Option integrate_yz bit(23)
    !% Integrates the function in the y-z plane and the result on the <i>x</i> axis is printed
    !%Option ascii bit(24)
    !% Plain text format regardless of dimensionality. For the moment only employed by the oct-phototoelectron_spectrum
    !% post-processing utility.
    !%End

    !%Variable OutputInterval
    !%Type integer
    !%Default 50
    !%Section Output
    !%Description
    !% The output requested by variable <tt>Output</tt> is written
    !% to the directory <tt>OutputIterDir</tt>
    !% when the iteration number is a multiple of the <tt>OutputInterval</tt> variable.
    !% Subdirectories are named Y.X, where Y is <tt>td</tt>, <tt>scf</tt>, or <tt>unocc</tt>, and
    !% X is the iteration number. To use the working directory, specify <tt>"."</tt>
    !% (Output of restart files is instead controlled by <tt>RestartWriteInterval</tt>.)
    !% Must be >= 0. If it is 0, then no output is written. For <tt>gs</tt> and <tt>unocc</tt>
    !% calculations, <tt>OutputDuringSCF</tt> must be set too for this output to be produced.
    !% This variable can also be defined inside the <tt>Output</tt> block.
    !% See <tt>Output</tt>.
    !%End


    what_no_how = (/ OPTION__OUTPUT__MATRIX_ELEMENTS, OPTION__OUTPUT__BERKELEYGW, OPTION__OUTPUT__DOS, &
      OPTION__OUTPUT__TPA, OPTION__OUTPUT__MMB_DEN, OPTION__OUTPUT__J_FLOW, OPTION__OUTPUT__OCC_MATRICES, &
      OPTION__OUTPUT__EFFECTIVEU, OPTION__OUTPUT__MAGNETIZATION, OPTION__OUTPUT__KANAMORIU /)

    if (parse_block(namespace, what_tag, blk) == 0) then
      nrows = parse_block_n(blk)
      do iout = 0, nrows - 1
        ncols = max(ncols , parse_block_cols(blk, iout))
      end do

      if (ncols == 1) then
        !new format, Type 0
        !%Output
        !  density
        !  wfs
        !%
        do iout = 1, nrows
          call parse_block_integer(blk, iout - 1, 0, what_i)
          if (.not. varinfo_valid_option(what_tag, what_i)) then
            call messages_input_error(namespace, what_tag)
          end if
          if (what_i > 0) then
            what(what_i) = .true.
            call parse_variable(namespace, output_interval_tag, 50, output_interval(what_i))
            if (((what_tag == 'Output') .and. (.not. any(what_no_how == what_i)))&
              .or. (what_tag /= 'Output')) then
              call parse_variable(namespace, how_tag, 0, how(what_i))
              if (.not. varinfo_valid_option(how_tag, how(what_i), is_flag=.true.)) then
                call messages_input_error(namespace, how_tag)
              end if
            end if
          end if
        end do
      else if (ncols == 2) then
        !new format, Type 1
        !%Output
        !  density | cube + axis_z
        !  wfs     | cube
        !%

        do iout = 1, nrows
          call parse_block_integer(blk, iout - 1, 0, what_i)
          if (.not. varinfo_valid_option(what_tag, what_i)) then
            call messages_input_error(namespace, what_tag)
          end if
          if (what_i > 0) then 
            what(what_i) = .true.
            call parse_variable(namespace, output_interval_tag, 50, output_interval(what_i))
            if (((what_tag == 'Output') .and. (.not. any(what_no_how == what_i)))&
              .or. (what_tag /= 'Output')) then
              call parse_block_integer(blk, iout - 1, 1, how(what_i))
              if (how(what_i) == 0) call parse_variable(namespace, how_tag, 0, how(what_i))
              if (.not. varinfo_valid_option(how_tag, how(what_i), is_flag=.true.)) then
                call messages_input_error(namespace, how_tag)
              end if
            end if
          end if
        end do

      else 
        !new format, Type 2 (tagged)
        !%Output
        !  density | "output_interval" | 10   | "output_format"   | cube + axis_z
        !  wfs     | "output_format"   | cube | "output_interval" | 50
        !%
        !
        ! OR
        !
        !%Output
        !  density | "output_interval" | 10   | "output_format" | cube + axis_z
        !  wfs     | "output_format"   | cube
        !%
        do iout = 1, nrows
          call parse_block_integer(blk, iout - 1, 0, what_i)
          if (.not. varinfo_valid_option(what_tag, what_i)) then
            call messages_input_error(namespace, what_tag)
          end if
          if (what_i > 0) then
            what(what_i) = .true.
            do column_index = 0, 1  
              call parse_block_string(blk, iout - 1, 1 + column_index * 2, output_column_marker)
              if (output_column_marker == 'output_interval') then
                call parse_block_integer(blk, iout - 1, 2 + column_index * 2, output_interval(what_i))
              else if (output_column_marker == 'output_format') then
                if (((what_tag == 'Output') .and. (.not. any(what_no_how == what_i)))&
                  .or. (what_tag /= 'Output')) then
                  call parse_block_integer(blk, iout - 1, 2 + column_index * 2, how(what_i))
                  if (.not. varinfo_valid_option(how_tag, how(what_i), is_flag=.true.)) then
                    call messages_input_error(namespace, how_tag)
                  end if
                end if
              else if (len_trim(output_column_marker) /= 0) then
                ! Unknown output_column_marker
                call messages_input_error(namespace, what_tag)
              else
                ! no output_column_marker -> full output info is not in this block
                if (output_interval(what_i) == 0) then
                  call parse_variable(namespace, output_interval_tag, 50, output_interval(what_i))
                end if
                if (how(what_i) == 0) then
                  if (((what_tag == 'Output') .and. (.not. any(what_no_how == what_i)))&
                    .or. (what_tag /= 'Output')) then
                    call parse_variable(namespace, how_tag, 0, how(what_i))
                    if (.not. varinfo_valid_option(how_tag, how(what_i), is_flag=.true.)) then
                      call messages_input_error(namespace, how_tag)
                    end if
                  end if
                end if
              end if
            end do
          end if
        end do
      end if
      call parse_block_end(blk)
    else

      call messages_variable_is_block(namespace, what_tag)

      ! Output block does not exist but we may have OutputHow/OutputInterval 
      call parse_variable(namespace, how_tag, 0, how(0))
      call parse_variable(namespace, output_interval_tag, 50, output_interval(0))
    endif

    
    do what_i = lbound(what, 1), ubound(what, 1)
      if (what_tag == 'Output') then
        if (what(what_i) .and. (.not. any(what_no_how == what_i))) then
          if (.not. varinfo_valid_option(how_tag, how(what_i), is_flag=.true.)) then
            call messages_input_error(namespace, how_tag)
          end if

          if (how(what_i) == 0 .and. .not. optional_default(ignore_error, .false.)) then
            write(message(1), '(a)') 'Must specify output method with variable OutputFormat.'
            call messages_fatal(1, only_root_writes = .true.)
          end if

          ! some modes are not available in some circumstances
          if (space%dim == 1) then
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__AXIS_Y) /= 0) then
              message(1) = "OutputFormat = axis_y not available with Dimensions = 1."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__PLANE_Z) /= 0) then
              message(1) = "OutputFormat = plane_z not available with Dimensions = 1."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__XCRYSDEN) /= 0) then
              message(1) = "OutputFormat = xcrysden not available with Dimensions = 1."
              call messages_fatal(1)
            end if
          end if

          if (space%dim <= 2) then
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__AXIS_Z) /= 0) then
              message(1) = "OutputFormat = axis_z not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__PLANE_X) /= 0) then
              message(1) = "OutputFormat = plane_x not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__PLANE_Y) /= 0) then
              message(1) = "OutputFormat = plane_y not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__INTEGRATE_XY) /= 0) then
              message(1) = "OutputFormat = integrate_xy not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__INTEGRATE_XZ) /= 0) then
              message(1) = "OutputFormat = integrate_xz not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__INTEGRATE_YZ) /= 0) then
              message(1) = "OutputFormat = integrate_yz not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__DX) /= 0) then
              message(1) = "OutputFormat = dx not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
            if (bitand(how(what_i), OPTION__OUTPUTFORMAT__CUBE) /= 0) then
              message(1) = "OutputFormat = cube not available with Dimensions <= 2."
              call messages_fatal(1)
            end if
          end if

  #if !defined(HAVE_NETCDF)
          if (bitand(how(what_i), OPTION__OUTPUTFORMAT__NETCDF) /= 0) then
            message(1) = 'Octopus was compiled without NetCDF support.'
            message(2) = 'It is not possible to write output in NetCDF format.'
            call messages_fatal(2)
          end if
  #endif
  #if !defined(HAVE_ETSF_IO)
          if (bitand(how(what_i), OPTION__OUTPUTFORMAT__ETSF) /= 0) then
            message(1) = 'Octopus was compiled without ETSF_IO support.'
            message(2) = 'It is not possible to write output in ETSF format.'
            call messages_fatal(2)
          end if
  #endif


        end if
      end if

      if (output_interval(what_i) < 0) then
        message(1) = "OutputInterval must be >= 0."
        call messages_fatal(1, namespace=namespace)
      end if
    end do
  
    POP_SUB(io_function_read_what_how_when)
  end subroutine io_function_read_what_how_when

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
    if(index(where, "IntegrateXY") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_XY)
    if(index(where, "IntegrateXZ") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_XZ)
    if(index(where, "IntegrateYZ") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_YZ)
    if(index(where, "PlaneZ")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Z)
    if(index(where, "DX")        /= 0) how = ior(how, OPTION__OUTPUTFORMAT__DX)
    if(index(where, "XCrySDen")  /= 0) how = ior(how, OPTION__OUTPUTFORMAT__XCRYSDEN)
    if(index(where, "Binary")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__BINARY)
    if(index(where, "MeshIndex") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__MESH_INDEX)
    if(index(where, "XYZ")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__XYZ)
#if defined(HAVE_NETCDF)
    if(index(where, "NETCDF")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__NETCDF)
#endif
    if(index(where, "Cube")      /= 0) how = ior(how, OPTION__OUTPUTFORMAT__CUBE)
    if(index(where, "VTK")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__VTK)

    POP_SUB(io_function_fill_how)
  end function io_function_fill_how

  ! ---------------------------------------------------------
  subroutine write_bild_forces_file(dir, fname, ions, namespace)
    character(len=*),   intent(in) :: dir, fname
    type(ions_t),       intent(in) :: ions
    type(namespace_t),  intent(in) :: namespace

    integer :: iunit, iatom, idir
    FLOAT, allocatable :: forces(:,:), center(:,:)
    character(len=20) frmt

    if (.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_bild_forces_file)

    call io_mkdir(dir, namespace)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.bild', namespace, action='write', &
      position='asis')

    write(frmt,'(a,i0,a)')'(a,2(', ions%space%dim,'f16.6,1x))'

    SAFE_ALLOCATE(forces(1:ions%space%dim, 1:ions%natoms))
    SAFE_ALLOCATE(center(1:ions%space%dim, 1:ions%natoms))
    center = units_from_atomic(units_out%length, ions%pos)
    forces = units_from_atomic(units_out%force, ions%tot_force)
    write(iunit, '(a)')'.comment : force vectors in ['//trim(units_abbrev(units_out%force))//']'
    write(iunit, *)
    write(iunit, '(a)')'.color red'
    write(iunit, *)
    do iatom = 1, ions%natoms
      write(iunit, '(a,1x,i4,1x,a2,1x,a6,1x,f10.6,a)') '.comment :', iatom, trim(ions%atom(iatom)%label), & 
                         'force:', norm2(forces(:, iatom)),'['//trim(units_abbrev(units_out%force))//']'
      write(iunit,fmt=trim(frmt)) '.arrow', (center(idir, iatom), idir = 1, ions%space%dim), &
                                 (center(idir, iatom) + forces(idir, iatom), idir = 1, ions%space%dim)
      write(iunit,*)
    end do

    SAFE_DEALLOCATE_A(forces)
    SAFE_DEALLOCATE_A(center)

    call io_close(iunit)

    POP_SUB(write_bild_forces_file)
  end subroutine write_bild_forces_file

  ! ---------------------------------------------------------
  !> Write canonicalized xyz file with atom labels and positions in Angstroms.
  !> Includes information about simulation box and periodicity when applicable.
  !> This differs from a normal xyz file by including information about box
  !> shape and always using Angstroms.
  subroutine write_canonicalized_xyz_file(dir, fname, ions, sb, namespace)
    character(len=*),  intent(in) :: dir
    character(len=*),  intent(in) :: fname
    type(ions_t),      intent(in) :: ions
    type(simul_box_t), intent(in) :: sb
    type(namespace_t), intent(in) :: namespace

    integer :: iunit
    integer :: idir
    FLOAT :: position
    integer :: iatom

    PUSH_SUB(write_canonicalized_xyz_file)

    call io_mkdir(dir, namespace)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xyz', namespace, action='write', position='asis')

    write(iunit, '(i6)') ions%natoms
    write(iunit, '(a,a,a)', advance='no') trim(ions%space%short_info()), '; ', trim(sb%short_info(unit_angstrom))
    if (ions%space%is_periodic()) then
      write(iunit, '(a,a)') '; ', trim(ions%latt%short_info(unit_angstrom))
    else
      write(iunit, '()')
    end if

    ! xyz-style labels and positions:
    do iatom = 1, ions%natoms
      write(iunit, '(10a)', advance='no') ions%atom(iatom)%label
      do idir = 1, 3
        if(idir <= ions%space%dim) then
          position = ions%pos(idir, iatom)
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

  ! ---------------------------------------------------------
  subroutine write_xsf_geometry_file(dir, fname, ions, mesh, namespace, write_forces)
    character(len=*),   intent(in) :: dir, fname
    type(ions_t),       intent(in) :: ions
    type(mesh_t),       intent(in) :: mesh
    type(namespace_t),  intent(in) :: namespace
    logical,  optional, intent(in) :: write_forces

    integer :: iunit
    FLOAT, allocatable :: forces(:,:)
    logical :: write_forces_

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_xsf_geometry_file)

    call io_mkdir(dir, namespace)
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xsf', namespace, action='write', position='asis')

    if(.not. present(write_forces)) then
      write_forces_ = .false.
    else
      write_forces_ = write_forces
    end if

    if(write_forces_) then
      SAFE_ALLOCATE(forces(1:ions%space%dim, 1:ions%natoms))
      forces = units_from_atomic(units_out%force, ions%tot_force)
      call write_xsf_geometry(iunit, ions, mesh, forces = forces)
      SAFE_DEALLOCATE_A(forces)
    else
      call write_xsf_geometry(iunit, ions, mesh)
    end if

    call io_close(iunit)

    POP_SUB(write_xsf_geometry_file)
  end subroutine write_xsf_geometry_file

  ! ---------------------------------------------------------
  !> for format specification see:
  !! http://www.xcrysden.org/doc/XSF.html#__toc__11
  subroutine write_xsf_geometry(iunit, ions, mesh, forces, index)
    integer,           intent(in) :: iunit
    type(ions_t),      intent(in) :: ions
    type(mesh_t),      intent(in) :: mesh
    FLOAT,   optional, intent(in) :: forces(1:ions%space%dim, 1:ions%natoms)
    integer, optional, intent(in) :: index !< for use in writing animated files

    integer :: idir, idir2, iatom, index_
    character(len=7) :: index_str
    FLOAT :: offset(3)

    PUSH_SUB(write_xsf_geometry)

    if(present(index)) then
      write(index_str, '(a,i6)') ' ', index
      index_ = index
    else
      write(index_str, '(a)') ''
      index_ = 1
    end if

    ! The corner of the cell is always (0,0,0) to XCrySDen
    ! so the offset is applied to the atomic coordinates.
    ! Along periodic dimensions the offset is -1/2 in reduced coordinates, as
    ! our origin is at the center of the cell instead of being at the edge.
    offset(1:ions%space%dim) = mesh%sb%latt%red_to_cart(spread(-M_HALF, 1, ions%space%dim))
    ! Offset in aperiodic directions:
    do idir = ions%space%periodic_dim + 1, 3
      offset(idir) = -(mesh%idx%ll(idir) - 1)/2 * mesh%spacing(idir)
    end do

    if(ions%space%is_periodic()) then
      if(index_ == 1) then
        select case(ions%space%periodic_dim)
        case(3)
          write(iunit, '(a)') 'CRYSTAL'
        case(2)
          write(iunit, '(a)') 'SLAB'
        case(1)
          write(iunit, '(a)') 'POLYMER'
        end select
      end if

      write(iunit, '(a)') 'PRIMVEC'//trim(index_str)

      do idir = 1, ions%space%dim
        write(iunit, '(3f12.6)') (units_from_atomic(units_out%length, ions%latt%rlattice(idir2, idir)), idir2 = 1, ions%space%dim)
      end do

      write(iunit, '(a)') 'PRIMCOORD'//trim(index_str)
      write(iunit, '(i10, a)') ions%natoms, ' 1'
    else
      write(iunit, '(a)') 'ATOMS'//trim(index_str)
    end if

    ! BoxOffset should be considered here
    do iatom = 1, ions%natoms
      write(iunit, '(a10, 3f12.6)', advance='no') trim(ions%atom(iatom)%label), &
        (units_from_atomic(units_out%length, ions%pos(idir, iatom) - offset(idir)), idir = 1, ions%space%dim)
      if(present(forces)) then
        write(iunit, '(5x, 3f12.6)', advance='no') forces(:, iatom)
      end if
      write(iunit, '()')
    end do

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

    if(status  ==  NF90_NOERR) then
    POP_SUB(ncdf_error)
      return
    end if

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

end module io_function_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
