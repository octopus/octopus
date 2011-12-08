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
!! $Id: output.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

module output_m
  use basins_m
  use cube_function_m
  use cube_m
  use datasets_m
  use density_m
  use derivatives_m
  use elf_m
#if defined(HAVE_ETSF_IO)
  use etsf_io
  use etsf_io_tools
#endif
  use fft_m
  use fourier_shell_m
  use fourier_space_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kpoints_m
  use linear_response_m
  use loct_m
  use loct_math_m
  use magnetic_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use modelmb_density_m
  use modelmb_density_matrix_m
  use modelmb_exchange_syms_m
  use mpi_m
  use output_me_m
  use parser_m
  use profiling_m
  use simul_box_m
  use solids_m
  use species_m
  use states_m
  use states_dim_m
  use states_io_m
  use symm_op_m
  use symmetries_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m
#if defined(HAVE_BERKELEYGW)
  use wfn_rho_vxc_io_m
#endif
  use young_m
  use xc_m

  implicit none

  private
  public ::              &
    output_t,            &
    output_bgw_t,        &
    output_init,         &
    output_states,       &
    output_modelmb,      &
    output_hamiltonian,  &
    output_all,          &
    output_current_flow, &
    doutput_lr,          &
    zoutput_lr

  type output_bgw_t
    integer           :: nbands
    integer           :: vxc_diag_nmin
    integer           :: vxc_diag_nmax
    integer           :: vxc_offdiag_nmin
    integer           :: vxc_offdiag_nmax
    logical           :: complex
    character(len=80) :: wfn_filename
  end type output_bgw_t

  type output_t
    ! General output variables:
    integer :: what                ! what to output
    integer :: how                 ! how to output

    type(output_me_t) :: me        ! this handles the output of matrix elements

    ! These variables fine-tune the output for some of the possible output options:
    integer :: iter                ! output every iter
    logical :: duringscf

    character(len=80) :: wfs_list  ! If output_wfs, this list decides which wavefunctions to print.

    type(mesh_plane_t) :: plane    ! This is to calculate the current flow across a plane
    type(mesh_line_t)  :: line     ! or though a line (in 2D)

    type(output_bgw_t) :: bgw      ! parameters for BerkeleyGW output

  end type output_t

  integer, parameter, public ::             &
    C_OUTPUT_POTENTIAL       =        1,    &
    C_OUTPUT_DENSITY         =        2,    &
    C_OUTPUT_WFS             =        4,    &
    C_OUTPUT_WFS_SQMOD       =        8,    &
    C_OUTPUT_GEOMETRY        =       16,    &
    C_OUTPUT_CURRENT         =       32,    &
    C_OUTPUT_ELF             =       64,    &
    C_OUTPUT_ELF_BASINS      =      128,    &
    C_OUTPUT_ELF_FS          =      256,    &
    C_OUTPUT_BADER           =      512,    &
    C_OUTPUT_EL_PRESSURE     =     1024,    &
    C_OUTPUT_MATRIX_ELEMENTS =     2048,    &
    C_OUTPUT_POL_DENSITY     =     4096,    &
    C_OUTPUT_R               =     8192,    &
    C_OUTPUT_KED             =    16384,    &
    C_OUTPUT_J_FLOW          =    32768,    &
    C_OUTPUT_DOS             =    65536,    &
    C_OUTPUT_TPA             =   131072,    &
    C_OUTPUT_DENSITY_MATRIX  =   262144,    &
    C_OUTPUT_MODELMB         =   524288,    &
    C_OUTPUT_FORCES          =  1048576,    &
    C_OUTPUT_WFS_FOURIER     =  2097152,    &
    C_OUTPUT_XC_DENSITY      =  4194304,    &
    C_OUTPUT_PES_WFS         =  8388608,    &
    C_OUTPUT_PES_DENSITY     = 16777216,    &
    C_OUTPUT_PES             = 33554432,    &
    C_OUTPUT_BERKELEYGW      = 67108864

contains

  subroutine output_init(sb, nst, outp)
    type(simul_box_t),    intent(in)  :: sb
    integer,              intent(in)  :: nst
    type(output_t),       intent(out) :: outp

    type(block_t) :: blk
    FLOAT :: norm
    character(len=80) :: nst_string, default

    PUSH_SUB(output_init)

    !%Variable Output
    !%Type flag
    !%Default no
    !%Section Output
    !%Description
    !% Specifies what to print. The output files go into the <tt>static</tt> directory, except when
    !% running a time-dependent simulation, when the directory <tt>td.XXXXXXX</tt> is used. For
    !% linear-response run modes, the derivatives of many quantities can be printed, as listed in
    !% the options below; the files will be printed in the directory
    !% for the run mode. Indices in the filename are labelled as follows:
    !% <tt>sp</tt> = spin, <tt>k</tt> = <i>k</i>-point, <tt>st</tt> = state/band,
    !% There is no tag for directions, given as a letter. The perturbation direction is always
    !% the last direction for linear-response quantities, and a following +/- indicates the sign of the frequency.
    !% Example: <tt>density + potential</tt>
    !%Option potential 1
    !% Outputs Kohn-Sham potential, separated by parts. File names are <tt>v0</tt> for 
    !% the local part, <tt>vc</tt> for the classical potential (if it exists), <tt>vh</tt> for the
    !% Hartree potential, and <tt>vxc-</tt> for the exchange-correlation potentials.
    !%Option density 2
    !% Outputs density. The output file is called <tt>density-</tt>, or <tt>lr_density-</tt> in linear response.
    !%Option wfs 4
    !% Outputs wavefunctions. Which wavefunctions are to be printed is specified
    !% by the variable <tt>OutputWfsNumber</tt> -- see below. The output file is called
    !% <tt>wf-</tt>, or <tt>lr_wf-</tt> in linear response.
    !%Option wfs_fourier 2097152
    !% (Experimental) Outputs wavefunctions in Fourier space. This is
    !% only implemented for the ETSF file format output. The file will
    !% be called <tt>static/wfs-pw-etsf.nc</tt>.  
    !%Option wfs_sqmod 8
    !% Outputs modulus squared of the wavefunctions. 
    !% The output file is called <tt>sqm-wf-</tt>. For linear response, the filename is <tt>sqm_lr_wf-</tt>.
    !%Option geometry 16
    !% Outputs file containing the coordinates of the atoms treated within quantum mechanics.
    !% If <tt>OutputHow = xyz</tt>, the file is called <tt>geometry.xyz</tt>; a
    !% file <tt>crystal.xyz</tt> is written with a supercell geometry if the system is periodic;
    !% if point charges were defined in the PDB file (see <tt>PDBCoordinates</tt>), they will be output
    !% in the file <tt>geometry_classical.xyz</tt>.
    !% If <tt>OutputHow = xcrysden</tt>, a file called <tt>geometry.xsf</tt> is written.
    !%Option current 32
    !% Outputs paramagnetic current density. The output file is called <tt>current-</tt>.
    !% For linear response, the filename is <tt>lr_current-</tt>.
    !%Option ELF 64
    !% Outputs electron localization function (ELF). The output file is called <tt>elf-</tt>,
    !% or <tt>lr_elf-</tt> in linear response, in which case the associated function D is also written,
    !% as <tt>lr_elf_D-</tt>. Only in 2D and 3D.
    !%Option ELF_basins 128
    !% Outputs basins of attraction of the ELF. The output file is called
    !% <tt>elf_rs_basins.info</tt>. Only in 2D and 3D.
    !%Option ELF_FS 256
    !% Outputs electron localization function in Fourier space (experimental). The output file is called
    !% <tt>elf_FS-</tt>. Only in 2D and 3D.
    !%Option Bader 512
    !% Outputs Laplacian of the density which shows lone pairs, bonded charge concentrations
    !% and regions subject to electrophilic or nucleophilic attack.
    !% See RF Bader, <i>Atoms in Molecules: A Quantum Theory</i> (Oxford Univ. Press, Oxford, 1990).
    !%Option el_pressure 1024
    !% Outputs electronic pressure. See Tao, Vignale, and Tokatly, <i>Phys Rev Lett</i> <b>100</b>, 206405 (2008).
    !%Option matrix_elements 2048
    !% Outputs a series of matrix elements of the Kohn-Sham states. What is output can
    !% be controlled by the <tt>OutputMatrixElements</tt> variable.
    !%Option pol_density 4096
    !% Outputs dipole-moment density <tt>dipole_density-</tt>, or polarizability density <tt>alpha_density-</tt>
    !% in linear response. If <tt>ResponseMethod = finite_differences</tt>, the hyperpolarizability density
    !% <tt>beta_density-</tt> is also printed.
    !%Option mesh_r 8192
    !% Outputs values of the coordinates over the grid. Files
    !% will be in the <tt>exec/</tt> directory.
    !%Option kinetic_energy_density 16384
    !% Outputs kinetic-energy density, defined as:
    !%
    !% <math>\tau_\sigma(\vec{r}) = \sum_{i=1}^{N_\sigma} 
    !%  \vert \nabla \phi_{i\sigma}(\vec{r}) \vert^2\,. </math>
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
    !%Option dos 65536
    !% Outputs density of states.
    !%Option tpa 131072
    !% Outputs transition-potential approximation (TPA) matrix elements.
    !%Option density_matrix 262144
    !% Calculates, and outputs, the reduced density
    !% matrix. For the moment the trace is made over the second dimension, and
    !% the code is limited to 2D. The idea is to model <i>N</i> particles in 1D as an
    !% <i>N</i>-dimensional non-interacting problem, then to trace out <i>N</i>-1 coordinates.
    !%Option modelmb 524288
    !% This flag turns on the output for model many-body calculations, and
    !% triggers the density, wavefunctions, or density matrix to be output for the
    !% particles described in the <tt>DescribeParticlesModelMB</tt> block. Which
    !% quantities will be output depends on the simultaneous presence of <tt>wfs</tt>,
    !% <tt>density</tt>, etc.
    !%Option forces 1048576
    !% Outputs file <tt>forces.xsf</tt> containing structure and forces on the atoms as 
    !% a vector associated with each atom, which can be visualized with XCrySDen.
    !%Option xc_density 4194304
    !% Outputs the XC density, which is the charge density that
    !% generates the XC potential. (This is <math>-1/4\pi</math> times
    !% the Laplacian of the XC potential). The files are called <tt>nxc</tt>.
    !%Option PES_wfs 8388608
    !% Outputs the photoelectron wavefunctions. The file name is <tt>pes_wfs-</tt>  
    !% plus the orbital number.
    !%Option PES_density 16777216
    !% Outputs the photolectron density. Output file is <tt>pes_dens-</tt> plus spin species if
    !% spin-polarized calculation is performed. 
    !%Option PES 33554432   
    !% Outputs the time-dependent photoelectron spectrum.
    !%Option BerkeleyGW 67108864
    !% Output for a run with BerkeleyGW (<tt>www.berkeleygw.org</tt>). See <tt>Output::BerkeleyGW</tt> for further specification.
    !%End
    call parse_integer(datasets_check('Output'), 0, outp%what)

    if(iand(outp%what, C_OUTPUT_ELF_FS) .ne. 0) then
      call messages_experimental("ELF in Fourier space")
    endif

    if(iand(outp%what, C_OUTPUT_WFS_FOURIER) .ne. 0) then
      call messages_experimental("Wave-functions in Fourier space")
    endif

    ! cannot calculate the ELF in 1D
    if(iand(outp%what, C_OUTPUT_ELF) .ne. 0 .or. iand(outp%what, C_OUTPUT_elf_basins) .ne. 0 &
       .or. iand(outp%what, C_OUTPUT_ELF_FS) .ne. 0) then
       if(sb%dim .ne. 2 .and. sb%dim .ne. 3) then
         outp%what = iand(outp%what, not(C_OUTPUT_ELF + C_OUTPUT_ELF_BASINS + C_OUTPUT_ELF_FS))
         write(message(1), '(a)') 'Cannot calculate ELF except in 2D and 3D.'
         call messages_warning(1)
       endif
    endif

    if(.not.varinfo_valid_option('Output', outp%what, is_flag=.true.)) then
      call input_error('Output')
    end if

    if(iand(outp%what, C_OUTPUT_MODELMB) .ne. 0 .or. iand(outp%what, C_OUTPUT_DENSITY_MATRIX) .ne. 0) then
      call messages_experimental("Model many-body and density matrix")
    endif

    if(iand(outp%what, C_OUTPUT_MODELMB) .ne. 0) then
      write(message(1),'(a)') 'Model many-body quantities will be output, according to the presence of'
      write(message(2),'(a)') '  wfs, density, or density_matrix in Output.'
      call messages_info(2)
    end if

    if(iand(outp%what, C_OUTPUT_DENSITY_MATRIX) .ne. 0) then
      write(message(1),'(a)') 'Info: The density matrix will be calculated, traced'
      write(message(2),'(a)') 'over the second dimension, diagonalized, and output.'
      call messages_info(2)
      if(iand(outp%what, C_OUTPUT_MODELMB) .eq. 0) then
        write(message(1),'(a)') 'Note that density matrix only works for model MB calculations for the moment.'
        call messages_info(1)
      end if
      ! NOTES:
      !   could be made into block to be able to specify which dimensions to trace
      !   in principle all combinations are interesting, but this means we need to
      !   be able to output density matrices for multiple particles or multiple
      !   dimensions. The current 1D 1-particle case is simple.
    end if

    if(iand(outp%what, C_OUTPUT_WFS) .ne. 0  .or.  iand(outp%what, C_OUTPUT_WFS_SQMOD) .ne. 0 ) then

      !%Variable OutputWfsNumber
      !%Type string
      !%Default all states
      !%Section Output
      !%Description
      !% Which wavefunctions to print, in list form: <i>i.e.</i>, "1-5" to print the first
      !% five states, "2,3" to print the second and the third state, etc.
      !% If more states are specified than available, extra ones will be ignored.
      !%End

      write(nst_string,'(i6)') nst
      write(default,'(a,a)') "1-", trim(adjustl(nst_string))
      call parse_string(datasets_check('OutputWfsNumber'), default, outp%wfs_list)
    end if

    if(parse_block(datasets_check('CurrentThroughPlane'), blk) == 0) then
      outp%what = ior(outp%what, C_OUTPUT_J_FLOW)

      !%Variable CurrentThroughPlane
      !%Type block
      !%Section States
      !%Description
      !% At the end of the ground-state calculation, the code can calculate
      !% the steady-state current in the ground state 
      !% traversing a user-defined portion of a plane, as specified by this block.
      !% In the format below, <tt>origin</tt> is a point in the plane.
      !% <tt>u</tt> and <tt>v</tt> are the (dimensionless) lattice vectors defining the plane;
      !% they will be normalized by the code. <tt>spacing</tt> is the fineness of the mesh
      !% on the plane. Integers <tt>nu</tt> and <tt>mu</tt> are the length and
      !% width of the portion of the plane, in units of <tt>spacing</tt>.
      !% Thus, the grid points included in the plane are
      !% <tt>x_ij = origin + i*spacing*u + j*spacing*v</tt>,
      !% for <tt>nu <= i <= mu </tt> and <tt>nv <= j <= mv</tt>.
      !% Analogously, in the 2D case, the current flow is calculated through a line;
      !% in the 1D case, the current flow is calculated through a point.
      !%
      !% Example (3D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 0.0  # origin
      !% <br>&nbsp;&nbsp; 0.0 | 1.0 | 0.0  # u
      !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 1.0  # v
      !% <br>&nbsp;&nbsp; 0.2              # spacing
      !% <br>&nbsp;&nbsp; 0 | 50           # nu | mu
      !% <br>&nbsp;&nbsp; -50 | 50         # nv | mv
      !% <br>%</tt>
      !%
      !% Example (2D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0 | 0.0        # origin
      !% <br>&nbsp;&nbsp; 1.0 | 0.0        # u
      !% <br>&nbsp;&nbsp; 0.2              # spacing
      !% <br>&nbsp;&nbsp; 0 | 50           # nu | mu
      !% <br>%</tt>
      !%
      !% Example (1D):
      !%
      !% <tt>%CurrentThroughPlane
      !% <br>&nbsp;&nbsp; 0.0              # origin
      !% <br>%</tt>
      !%
      !%End
        
      select case(sb%dim)
      case(3)

        call parse_block_float(blk, 0, 0, outp%plane%origin(1), units_inp%length)
        call parse_block_float(blk, 0, 1, outp%plane%origin(2), units_inp%length)
        call parse_block_float(blk, 0, 2, outp%plane%origin(3), units_inp%length)
        call parse_block_float(blk, 1, 0, outp%plane%u(1))
        call parse_block_float(blk, 1, 1, outp%plane%u(2))
        call parse_block_float(blk, 1, 2, outp%plane%u(3))
        call parse_block_float(blk, 2, 0, outp%plane%v(1))
        call parse_block_float(blk, 2, 1, outp%plane%v(2))
        call parse_block_float(blk, 2, 2, outp%plane%v(3))
        call parse_block_float(blk, 3, 0, outp%plane%spacing, units_inp%length)
        call parse_block_integer(blk, 4, 0, outp%plane%nu)
        call parse_block_integer(blk, 4, 1, outp%plane%mu)
        call parse_block_integer(blk, 5, 0, outp%plane%nv)
        call parse_block_integer(blk, 5, 1, outp%plane%mv)

        norm = sqrt(sum(outp%plane%u(1:3)**2))
        if(norm < M_EPSILON) then
          write(1, '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        endif
        outp%plane%u(1:3) = outp%plane%u(1:3) / norm

        norm = sqrt(sum(outp%plane%v(1:3)**2))
        if(norm < M_EPSILON) then
          write(1, '(a)') 'v-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        endif
        outp%plane%v(1:3) = outp%plane%v(1:3) / norm

        outp%plane%n(1) = outp%plane%u(2)*outp%plane%v(3) - outp%plane%u(3)*outp%plane%v(2)
        outp%plane%n(2) = outp%plane%u(3)*outp%plane%v(1) - outp%plane%u(1)*outp%plane%v(3)
        outp%plane%n(3) = outp%plane%u(1)*outp%plane%v(2) - outp%plane%u(2)*outp%plane%v(1)

      case(2)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)
        call parse_block_float(blk, 0, 1, outp%line%origin(2), units_inp%length)
        call parse_block_float(blk, 1, 0, outp%line%u(1))
        call parse_block_float(blk, 1, 1, outp%line%u(2))
        call parse_block_float(blk, 2, 0, outp%line%spacing, units_inp%length)
        call parse_block_integer(blk, 3, 0, outp%line%nu)
        call parse_block_integer(blk, 3, 1, outp%line%mu)

        norm = sqrt(sum(outp%line%u(1:2)**2))
        if(norm < M_EPSILON) then
          write(1, '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        endif
        outp%line%u(1:2) = outp%line%u(1:2) / norm

        outp%line%n(1) = -outp%line%u(2)
        outp%line%n(2) =  outp%line%u(1)

      case(1)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)

      end select
    end if

    if(iand(outp%what, C_OUTPUT_MATRIX_ELEMENTS) .ne. 0) then
      call output_me_init(outp%me, sb)
    end if

    if(iand(outp%what, C_OUTPUT_BERKELEYGW) .ne. 0) then
      call output_berkeleygw_init(nst, outp%bgw)
    end if

    !%Variable OutputEvery
    !%Type integer
    !%Default 50
    !%Section Output
    !%Description
    !% The output is saved when the iteration number is a multiple of the
    !% <tt>OutputEvery</tt> variable. This works for the ground-state and
    !% time-dependent runs.
    !%End
    call parse_integer(datasets_check('OutputEvery'), 50, outp%iter)

    !%Variable OutputDuringSCF
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% If this variable is set to yes, during a ground-state run,
    !% <tt>Octopus</tt> output will be written after every self-consistent
    !% iteration to a directory called <tt>scf.nnnn/</tt> (with
    !% <tt>nnnn</tt> the iteration number).
    !%End

    call parse_logical(datasets_check('OutputDuringSCF'), .false., outp%duringscf)

    if(outp%what .ne. 0 .and. outp%what .ne. C_OUTPUT_MATRIX_ELEMENTS) call io_function_read_how(sb, outp%how)

    POP_SUB(output_init)
  end subroutine output_init


  ! ---------------------------------------------------------
  subroutine output_all(outp, gr, geo, st, hm, xc, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    type(xc_t),           intent(in)    :: xc
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir

    integer :: idir
    FLOAT   :: offset(1:MAX_DIM)
    
    PUSH_SUB(output_all)
    
    call output_states(st, gr, geo, dir, outp)
    call output_hamiltonian(hm, gr%der, dir, outp, geo)
    call output_localization_funct(st, hm, gr, dir, outp, geo)
    call output_current_flow(gr, st, dir, outp, geo)

    if(iand(outp%what, C_OUTPUT_GEOMETRY) .ne. 0) then
      if(iand(outp%how, C_OUTPUT_HOW_XCRYSDEN) .ne. 0) then        
        call write_xsf_geometry_file(dir, "geometry", geo, gr%mesh)
      endif
      if(iand(outp%how, C_OUTPUT_HOW_XYZ) .ne. 0) then
        call atom_write_xyz(dir, "geometry", geo, gr%sb%dim)
        if(simul_box_is_periodic(gr%sb)) &
          call periodic_write_crystal(gr%sb, geo, dir)
      endif
    end if

    if(iand(outp%what, C_OUTPUT_FORCES) .ne. 0) then
      call write_xsf_geometry_file(dir, "forces", geo, gr%mesh, write_forces = .true.)
    endif

    if(iand(outp%what, C_OUTPUT_MATRIX_ELEMENTS) .ne. 0) then
      call output_me(outp%me, dir, st, gr, geo, hm)
    end if

    if (iand(outp%how, C_OUTPUT_HOW_ETSF) .ne. 0) then
      call output_etsf(st, gr, geo, dir, outp)
    end if

    if (iand(outp%what, C_OUTPUT_BERKELEYGW) .ne. 0) then
      call output_berkeleygw(outp%bgw, dir, st, gr, xc, geo)
    end if
    
    POP_SUB(output_all)
  end subroutine output_all


  ! ---------------------------------------------------------
  subroutine output_localization_funct(st, hm, gr, dir, outp, geo)
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: dir
    type(output_t),         intent(in)    :: outp
    type(geometry_t),       intent(in)    :: geo

    FLOAT, allocatable :: f_loc(:,:)
    character(len=256) :: fname
    integer :: is, ierr, imax
    type(mpi_grp_t) :: mpi_grp

    PUSH_SUB(output_localization_funct)
    
    mpi_grp = gr%mesh%mpi_grp
    if(st%parallel_in_states) mpi_grp = st%mpi_grp
    if(st%d%kpt%parallel) mpi_grp = st%d%kpt%mpi_grp

    ! if SPIN_POLARIZED, the ELF contains one extra channel: the total ELF
    imax = st%d%nspin
    if(st%d%ispin == SPIN_POLARIZED) imax = 3

    SAFE_ALLOCATE(f_loc(1:gr%mesh%np, 1:imax))

    ! First the ELF in real space
    if(iand(outp%what, C_OUTPUT_ELF) .ne. 0 .or. iand(outp%what, C_OUTPUT_ELF_BASINS) .ne. 0) then
      ASSERT(gr%mesh%sb%dim .ne. 1)

      call elf_calc(st, gr, f_loc)
      
      ! output ELF in real space
      if(iand(outp%what, C_OUTPUT_ELF) .ne. 0) then
        write(fname, '(a)') 'elf_rs'
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,imax), unit_one, ierr, is_tmp = .false., geo = geo)
        ! this quantity is dimensionless

        if(st%d%ispin .ne. UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,i1)') 'elf_rs-sp', is
            call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
              f_loc(:, is), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
            ! this quantity is dimensionless
          end do
        end if
      end if

      if(iand(outp%what, C_OUTPUT_ELF_BASINS) .ne. 0) &
        call out_basins(f_loc(:,1), "elf_rs_basins")
    end if

    ! Second, ELF in Fourier space.
    if(iand(outp%what, C_OUTPUT_ELF_FS) .ne. 0) then
      call elf_calc_fs(st, gr, f_loc)
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'elf_fs-sp', is
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,is), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
        ! this quantity is dimensionless
      end do
    end if

    ! Now Bader analysis
    if(iand(outp%what, C_OUTPUT_BADER) .ne. 0) then
      do is = 1, st%d%nspin
        call dderivatives_lapl(gr%der, st%rho(:,is), f_loc(:,is))
        write(fname, '(a,i1)') 'bader-sp', is
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,is), units_out%length**(-2 - gr%sb%dim), ierr, is_tmp = .false., &
          geo = geo, grp = mpi_grp)

        write(fname, '(a,i1)') 'bader_basins-sp', is
        call out_basins(f_loc(:,1), fname)
      end do
    end if

    ! Now the pressure
    if(iand(outp%what, C_OUTPUT_EL_PRESSURE) .ne. 0) then
      call calc_electronic_pressure(st, hm, gr, f_loc(:,1))
      call dio_function_output(outp%how, dir, "el_pressure", gr%mesh, &
        f_loc(:,1), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
      ! this quantity is dimensionless
    end if

    SAFE_DEALLOCATE_A(f_loc)

    POP_SUB(output_localization_funct)

  contains
    ! ---------------------------------------------------------
    subroutine out_basins(ff, filename)
      FLOAT,            intent(in)    :: ff(:)
      character(len=*), intent(in)    :: filename

      character(len=256) :: fname
      type(basins_t)     :: basins
      integer            :: iunit

      PUSH_SUB(output_localization_funct.out_basins)

      call basins_init(basins, gr%mesh)
      call basins_analyze(basins, gr%mesh, st%d%nspin, ff(:), st%rho, CNST(0.01))

      call dio_function_output(outp%how, dir, trim(filename), gr%mesh, &
        real(basins%map, REAL_PRECISION), unit_one, ierr, is_tmp = .false., geo = geo)
      ! this quantity is dimensionless

      write(fname,'(4a)') trim(dir), '/', trim(filename), '.info'
      iunit = io_open(file=trim(fname), action = 'write')
      call basins_write(basins, gr%mesh, iunit)
      call io_close(iunit)

      call basins_end(basins)

      POP_SUB(output_localization_funct.out_basins)
    end subroutine out_basins

  end subroutine output_localization_funct

  
  ! ---------------------------------------------------------
  subroutine calc_electronic_pressure(st, hm, gr, pressure)
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(grid_t),           intent(inout) :: gr
    FLOAT,                  intent(out)   :: pressure(:)

    FLOAT, allocatable :: rho(:,:), lrho(:), tau(:,:)
    FLOAT   :: p_tf, dens
    integer :: is, ii

    PUSH_SUB(calc_electronic_pressure)

    SAFE_ALLOCATE( rho(1:gr%mesh%np_part, 1:st%d%nspin))
    SAFE_ALLOCATE(lrho(1:gr%mesh%np))
    SAFE_ALLOCATE( tau(1:gr%mesh%np, 1:st%d%nspin))

    rho = M_ZERO
    call density_calc(st, gr, rho)
    call states_calc_quantities(gr%der, st, kinetic_energy_density = tau)

    pressure = M_ZERO
    do is = 1, st%d%spin_channels
      lrho = M_ZERO
      call dderivatives_lapl(gr%der, rho(:, is), lrho)

      pressure(:) = pressure(:) + &
        tau(:, is)/M_THREE - lrho(:)/M_FOUR
    end do

    do ii = 1, gr%mesh%np
      dens = sum(rho(ii,1:st%d%spin_channels))

      p_tf = M_TWO/M_FIVE*(M_THREE*M_PI**2)**(M_TWO/M_THREE)* &
        dens**(M_FIVE/M_THREE)

      ! add XC pressure
      pressure(ii) = pressure(ii) + (dens*hm%vxc(ii,1) - hm%energy%exchange - hm%energy%correlation)

      pressure(ii) = pressure(ii)/p_tf
      pressure(ii) = M_HALF*(M_ONE + pressure(ii)/sqrt(M_ONE + pressure(ii)**2))
    end do

    POP_SUB(calc_electronic_pressure)
  end subroutine calc_electronic_pressure


  ! ---------------------------------------------------------
  subroutine output_berkeleygw_init(nst, bgw)
    integer,            intent(in)  :: nst
    type(output_bgw_t), intent(out) :: bgw
  
    PUSH_SUB(output_berkeleygw_init)
  
    ! conditions to die: spinors, not 3D, parallel in states or k-points (if spin-polarized), non-local functionals, SIC
    ! nlcc, not a kgrid, st%smear%method == SMEAR_FIXED_OCC, single precision

#ifndef HAVE_BERKELEYGW
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1)
#endif
  
    !%Variable BerkeleyGW_NumberBands
    !%Type integer
    !%Default all states
    !%Section Output::BerkeleyGW
    !%Description
    !% Wavefunctions for bands up to this number will be output.
    !%End
    call parse_integer(datasets_check('BerkeleyGW_NumberBands'), nst, bgw%nbands)
  
    !%Variable BerkeleyGW_Vxc_diag_nmin
    !%Type integer
    !%Default 1
    !%Section Output::BerkeleyGW
    !%Description
    !% Lowest band for which to write diagonal exchange-correlation matrix elements.
    !%End
    call parse_integer(datasets_check('BerkeleyGW_Vxc_diag_nmin'), 1, bgw%vxc_diag_nmin)
    
    !%Variable BerkeleyGW_Vxc_diag_nmax
    !%Type integer
    !%Default nst
    !%Section Output::BerkeleyGW
    !%Description
    !% Highest band for which to write diagonal exchange-correlation matrix elements.
    !%End
    call parse_integer(datasets_check('BerkeleyGW_Vxc_diag_nmax'), nst, bgw%vxc_diag_nmax)
    
    !%Variable BerkeleyGW_Vxc_offdiag_nmin
    !%Type integer
    !%Default 1
    !%Section Output::BerkeleyGW
    !%Description
    !% Lowest band for which to write off-diagonal exchange-correlation matrix elements.
    !% If < 1, off-diagonals will be skipped.
    !%End
    call parse_integer(datasets_check('BerkeleyGW_Vxc_offdiag_nmin'), 1, bgw%vxc_offdiag_nmin)
    
    !%Variable BerkeleyGW_Vxc_offdiag_nmax
    !%Type integer
    !%Default nst
    !%Section Output::BerkeleyGW
    !%Description
    !% Highest band for which to write off-diagonal exchange-correlation matrix elements.
    !% If < 1, off-diagonals will be skipped.
    !%End
    call parse_integer(datasets_check('BerkeleyGW_Vxc_offdiag_nmax'), nst, bgw%vxc_offdiag_nmax)
    
    !%Variable BerkeleyGW_Complex
    !%Type logical
    !%Default false
    !%Section Output::BerkeleyGW
    !%Description
    !% Even when wavefunctions, density, and XC potential could be real in reciprocal space,
    !% they will be output as complex.
    !%End
    call parse_logical(datasets_check('BerkeleyGW_Complex'), .false., bgw%complex)
    
    !%Variable BerkeleyGW_WFN_filename
    !%Type string
    !%Default WFN
    !%Section Output::BerkeleyGW
    !%Description
    !% Filename for the wavefunctions.
    !%End
    call parse_string(datasets_check('BerkeleyGW_WFN_filename'), 'WFN', bgw%wfn_filename)
  
    POP_SUB(output_berkeleygw_init)
  end subroutine output_berkeleygw_init


  ! ---------------------------------------------------------
  subroutine output_berkeleygw(bgw, dir, st, gr, xc, geo)
    type(output_bgw_t), intent(in) :: bgw
    character(len=*),   intent(in) :: dir
    type(states_t),     intent(in) :: st
    type(grid_t),       intent(in) :: gr
    type(xc_t),         intent(in) :: xc
    type(geometry_t),   intent(in) :: geo

    integer :: ik, is, ikk, ist, itran, iunit, iatom, mtrx(3, 3, 48), ig, ix, iy, iz
    integer, pointer :: ifmin(:,:), ifmax(:,:), atyp(:), ngk(:)
    character*3 :: sheader
    FLOAT :: adot(3,3), bdot(3,3), recvol, tnp(3, 48)
    FLOAT, pointer :: energies(:,:,:), occupations(:,:,:), apos(:,:), vxc(:,:)
    CMPLX, pointer :: field_g(:,:)
    type(cube_t) :: dcube
    type(cube_function_t) :: cf
    type(fourier_shell_t) :: shell

    PUSH_SUB(output_berkeleygw)

#ifdef HAVE_BERKELEYGW

    SAFE_ALLOCATE(vxc(gr%mesh%np, st%d%nspin))
    vxc(:,:) = M_ZERO
    ! we should not include core rho here. that is why we do not just use hm%vxc
    call xc_get_vxc(gr%der, xc, st, st%rho, st%d%ispin, -minval(st%eigenval(st%nst, :)), st%qtot, vxc = vxc)

    message(1) = "BerkeleyGW output: vxc.dat"
    call messages_info(1)

    if(states_are_real(st)) then
      call dbgw_vxc_dat(bgw, dir, st, gr, vxc)
    else
      call zbgw_vxc_dat(bgw, dir, st, gr, vxc)
    endif

    call cube_init(dcube, gr%mesh%idx%ll, gr%sb, fft_type=FFT_COMPLEX)
    call cube_function_null(cf)
    call dcube_function_alloc_rs(dcube, cf)
    call dcube_function_alloc_fs(dcube, cf)
    call fourier_shell_init(shell, dcube, gr%mesh)
    SAFE_ALLOCATE(field_g(shell%ngvectors, st%d%nspin))

    adot(:,:) = matmul(gr%sb%rlattice, gr%sb%rlattice)
    bdot(:,:) = matmul(gr%sb%klattice, gr%sb%klattice)
    recvol = (M_TWO * M_PI)**3 / gr%sb%rcell_volume

    ! symmetry is not analyzed by Octopus for finite systems, but we only need it for periodic ones
    do itran = 1, symmetries_number(gr%sb%symm)
      mtrx(:,:, itran) = symm_op_rotation_matrix(gr%sb%symm%ops(itran))
      tnp(:, itran) = symm_op_translation_vector(gr%sb%symm%ops(itran))
    enddo

    SAFE_ALLOCATE(ifmin(gr%sb%kpoints%reduced%npoints, st%d%nspin))
    SAFE_ALLOCATE(ifmax(gr%sb%kpoints%reduced%npoints, st%d%nspin))
    SAFE_ALLOCATE(energies(st%nst, gr%sb%kpoints%reduced%npoints, st%d%nspin))
    SAFE_ALLOCATE(occupations(st%nst, gr%sb%kpoints%reduced%npoints, st%d%nspin))
    SAFE_ALLOCATE(ngk(gr%sb%kpoints%reduced%npoints))
    ifmin(:,:) = 1
    ifmax(:,:) = 0
    ngk(:) = shell%ngvectors
    do ik = 1, st%d%nik
      is = states_dim_get_spin_index(st%d, ik)
      ikk = states_dim_get_kpoint_index(st%d, ik)
      energies(1:st%nst, ikk, is) = st%eigenval(1:st%nst,ik) * M_TWO
      occupations(1:st%nst, ikk, is) = st%occ(1:st%nst, ik) / st%smear%el_per_state
      do ist = 1, st%nst
        if(st%eigenval(ist, ik) > st%smear%e_fermi) then
          ifmax(ikk, is) = ist - 1
          exit
        endif
      enddo
    enddo

    SAFE_ALLOCATE(atyp(geo%natoms))
    SAFE_ALLOCATE(apos(3, geo%natoms))
    do iatom = 1, geo%natoms
      atyp(iatom) = species_index(geo%atom(iatom)%spec)
      apos(1:3, iatom) = geo%atom(iatom)%x(1:3)
    enddo

    message(1) = "BerkeleyGW output: VXC"
    call messages_info(1)

    sheader = 'VXC'
    if(mpi_grp_is_root(mpi_world)) call bgw_write_header(sheader, iunit)
    call dbgw_write_FS(iunit, vxc, field_g, shell, st%d%nspin, gr, dcube, cf)
    if(mpi_grp_is_root(mpi_world)) call io_close(iunit)
    SAFE_DEALLOCATE_P(vxc)

    message(1) = "BerkeleyGW output: RHO"
    call messages_info(1)

    sheader = 'RHO'
    if(mpi_grp_is_root(mpi_world)) call bgw_write_header(sheader, iunit)
    call dbgw_write_FS(iunit, st%rho, field_g, shell, st%d%nspin, gr, dcube, cf)
    if(mpi_grp_is_root(mpi_world)) call io_close(iunit)

    call fourier_shell_end(shell)
    call dcube_function_free_fs(dcube, cf)
    call dcube_function_free_rs(dcube, cf)

    ! wfns

    SAFE_DEALLOCATE_P(field_g)
    SAFE_DEALLOCATE_P(ifmin)
    SAFE_DEALLOCATE_P(ifmax)
    SAFE_DEALLOCATE_P(ngk)
    SAFE_DEALLOCATE_P(energies)
    SAFE_DEALLOCATE_P(occupations)
    SAFE_DEALLOCATE_P(atyp)
    SAFE_DEALLOCATE_P(apos)

#else
    message(1) = "Cannot do BerkeleyGW output: the library was not linked."
    call messages_fatal(1)
#endif

    POP_SUB(output_berkeleygw)

  contains
    
    subroutine bgw_write_header(sheader, iunit)
      character(len=3), intent(inout) :: sheader
      integer,          intent(out)   :: iunit
      
      PUSH_SUB(output_berkeleygw.bgw_write_header)

      if(mpi_grp_is_root(mpi_world)) then
        iunit = io_open(trim(dir) // sheader, form = 'unformatted', action = 'write')
        call write_binary_header(iunit, sheader, iflavor = 2, ns = st%d%nspin, ng = shell%ngvectors, &
          ntran = symmetries_number(gr%sb%symm), cell_symmetry = 0, nat = geo%natoms, &
          nk = gr%sb%kpoints%reduced%npoints, nbands = st%nst, ngkmax = shell%ngvectors, ecutrho = shell%ekin_cutoff / 2,  &
          ecutwfc = shell%ekin_cutoff / 2, kmax = gr%mesh%idx%ll, kgrid = gr%sb%kpoints%nik_axis, kshift = gr%sb%kpoints%shifts, &
          celvol = gr%sb%rcell_volume, alat = M_ONE, avec = gr%sb%rlattice, adot = adot, recvol = recvol, &
          blat = M_ONE, bvec = gr%sb%klattice, bdot = bdot, mtrx = mtrx, tnp = tnp, atyp = atyp, &
          apos = apos, ngk = ngk, kw = gr%sb%kpoints%reduced%weight, kpt = gr%sb%kpoints%reduced%red_point, &
          ifmin = ifmin, ifmax = ifmax, energies = energies, occupations = occupations, warn = .false.)

        call write_binary_gvectors(iunit, shell%ngvectors, shell%ngvectors, shell%red_gvec)
      endif

      POP_SUB(output_berkeleygw.bgw_write_header)
    end subroutine bgw_write_header

  end subroutine output_berkeleygw

#include "output_etsf_inc.F90"

#include "output_states_inc.F90"
#include "output_h_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_linear_response_inc.F90"
#include "output_berkeleygw_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "output_linear_response_inc.F90"
#include "output_berkeleygw_inc.F90"

end module output_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
