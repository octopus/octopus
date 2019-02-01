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

module output_oct_m
  use basins_oct_m
  use batch_oct_m
  use comm_oct_m 
  use cube_function_oct_m
  use cube_oct_m
  use current_oct_m
  use density_oct_m
  use derivatives_oct_m
  use dos_oct_m
  use elf_oct_m
  use fft_oct_m
  use fourier_shell_oct_m
  use fourier_space_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use io_oct_m
  use io_function_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use lda_u_io_oct_m
  use linear_response_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use magnetic_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use orbitalset_oct_m
  use output_me_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use periodic_copy_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use string_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use submesh_oct_m
  use symm_op_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use vtk_oct_m
  use vdw_ts_oct_m
  use young_oct_m
  use xc_oct_m

  implicit none

  private
  public ::              &
    output_t,            &
    output_init,         &
    output_end,          &
    output_states,       &
    output_hamiltonian,  &
    output_all,          &
    output_current_flow, &
    output_kick,         &
    output_scalar_pot

  type output_t
    !> General output variables:
    integer(8) :: what                !< what to output
    integer(8) :: whatBZ              !< what to output - for k-point resolved output
    integer(8) :: what_lda_u          !< what to output for the LDA+U part
    integer(8) :: how                 !< how to output

    type(output_me_t) :: me        !< this handles the output of matrix elements

    !> These variables fine-tune the output for some of the possible output options:
    integer :: output_interval     !< output every iter
    integer :: restart_write_interval
    logical :: duringscf
    character(len=80) :: wfs_list  !< If output_wfs, this list decides which wavefunctions to print.
    character(len=MAX_PATH_LEN) :: iter_dir  !< The folder name, if information will be output while iterating.

    type(mesh_plane_t) :: plane    !< This is to calculate the current flow across a plane
    type(mesh_line_t)  :: line     !< or through a line (in 2D)

  end type output_t

  integer(8), parameter, public ::              &
    OPTION__OUTPUT__J_FLOW          =     32768
  
contains

  subroutine output_init(outp, sb, nst, ks)
    type(output_t),       intent(out)   :: outp
    type(simul_box_t),    intent(in)    :: sb
    integer,              intent(in)    :: nst
    type(v_ks_t),         intent(inout) :: ks

    type(block_t) :: blk
    FLOAT :: norm
    character(len=80) :: nst_string, default
    integer :: what_no_how, what_no_how_u

    PUSH_SUB(output_init)

    !%Variable Output
    !%Type flag
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. The output files are written at the end of the run into the output directory for the
    !% relevant kind of run (<i>e.g.</i> <tt>static</tt> for <tt>CalculationMode = gs</tt>).
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% (available for <tt>CalculationMode</tt> = <tt>gs</tt>, <tt>unocc</tt>,  <tt>td</tt>, and <tt>opt_control</tt>)
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir</tt>.
    !% For linear-response run modes, the derivatives of many quantities can be printed, as listed in
    !% the options below. Indices in the filename are labelled as follows:
    !% <tt>sp</tt> = spin (or spinor component), <tt>k</tt> = <i>k</i>-point, <tt>st</tt> = state/band.
    !% There is no tag for directions, given as a letter. The perturbation direction is always
    !% the last direction for linear-response quantities, and a following +/- indicates the sign of the frequency.
    !% Example: <tt>density + potential</tt>
    !%Option potential  bit(0)
    !% Outputs Kohn-Sham potential, separated by parts. File names are <tt>v0</tt> for 
    !% the local part of the ionic potential, <tt>vc</tt> for the classical potential (if it exists),
    !% <tt>vh</tt> for the Hartree potential, <tt>vks</tt> for the local part of the Kohn-Sham potential, and
    !% <tt>vxc-</tt> for the exchange-correlation potentials. For <tt>vks</tt> and <tt>vxc</tt>,
    !% a suffix for spin is added in the spin-polarized case.
    !%Option density bit(1)
    !% Outputs density. The output file is called <tt>density-</tt>, or <tt>lr_density-</tt> in linear response.
    !%Option wfs bit(2)
    !% Outputs wavefunctions. Which wavefunctions are to be printed is specified
    !% by the variable <tt>OutputWfsNumber</tt> -- see below. The output file is called
    !% <tt>wf-</tt>, or <tt>lr_wf-</tt> in linear response.
    !%Option wfs_sqmod bit(3)
    !% Outputs modulus squared of the wavefunctions. 
    !% The output file is called <tt>sqm-wf-</tt>. For linear response, the filename is <tt>sqm_lr_wf-</tt>.
    !%Option geometry bit(4)
    !% Outputs file containing the coordinates of the atoms treated within quantum mechanics.
    !% If <tt>OutputFormat = xyz</tt>, the file is called <tt>geometry.xyz</tt>; a
    !% file <tt>crystal.xyz</tt> is written with a supercell geometry if the system is periodic;
    !% if point charges were defined in the PDB file (see <tt>PDBCoordinates</tt>), they will be output
    !% in the file <tt>geometry_classical.xyz</tt>.
    !% If <tt>OutputFormat = xcrysden</tt>, a file called <tt>geometry.xsf</tt> is written.
    !%Option current bit(5)
    !% Outputs the total current density. The output file is called <tt>current-</tt>.
    !% For linear response, the filename is <tt>lr_current-</tt>.
    !%Option ELF bit(6)
    !% Outputs electron localization function (ELF). The output file is called <tt>elf-</tt>,
    !% or <tt>lr_elf-</tt> in linear response, in which case the associated function D is also written,
    !% as <tt>lr_elf_D-</tt>. Only in 2D and 3D.
    !%Option ELF_basins bit(7)
    !% Outputs basins of attraction of the ELF. The output file is called
    !% <tt>elf_rs_basins.info</tt>. Only in 2D and 3D.
    !%Option Bader bit(9)
    !% Outputs Laplacian of the density which shows lone pairs, bonded charge concentrations
    !% and regions subject to electrophilic or nucleophilic attack.
    !% See RF Bader, <i>Atoms in Molecules: A Quantum Theory</i> (Oxford Univ. Press, Oxford, 1990).
    !%Option el_pressure bit(10)
    !% Outputs electronic pressure. See Tao, Vignale, and Tokatly, <i>Phys Rev Lett</i> <b>100</b>, 206405 (2008).
    !%Option matrix_elements bit(11)
    !% Outputs a series of matrix elements of the Kohn-Sham states. What is output can
    !% be controlled by the <tt>OutputMatrixElements</tt> variable.
    !%Option pol_density bit(12)
    !% Outputs dipole-moment density <tt>dipole_density-</tt>, or polarizability density <tt>alpha_density-</tt>
    !% in linear response. If <tt>ResponseMethod = finite_differences</tt>, the hyperpolarizability density
    !% <tt>beta_density-</tt> is also printed.
    !%Option mesh_r bit(13)
    !% Outputs values of the coordinates over the grid. Files
    !% will be called <tt>mesh_r-</tt> followed by the direction.
    !%Option kinetic_energy_density bit(14)
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
    !%Option dos bit(16)
    !% Outputs density of states. See <tt>DOSEnergyMax</tt>, <tt>DOSEnergyMin</tt>, <tt>DOSEnergyPoints</tt>,
    !% and <tt>DOSGamma</tt>.
    !%Option tpa bit(17)
    !% Outputs transition-potential approximation (TPA) matrix elements, using <math>\vec{q}</math>-vector specified
    !% by <tt>MomentumTransfer</tt>.
    !%Option forces bit(18)
    !% Outputs file <tt>forces.xsf</tt> containing structure and forces on the atoms as 
    !% a vector associated with each atom, which can be visualized with XCrySDen.
    !%Option xc_density bit(20)
    !% Outputs the XC density, which is the charge density that
    !% generates the XC potential. (This is <math>-1/4\pi</math> times
    !% the Laplacian of the XC potential). The files are called <tt>nxc</tt>.
    !%Option PES_wfs bit(21)
    !% Outputs the photoelectron wavefunctions. The file name is <tt>pes_wfs-</tt>  
    !% plus the orbital number.
    !%Option PES_density bit(22)
    !% Outputs the photolectron density. Output file is <tt>pes_dens-</tt> plus spin species if
    !% spin-polarized calculation is performed. 
    !%Option PES bit(23)
    !% Outputs the time-dependent photoelectron spectrum.
    !%Option delta_perturbation bit(25)
    !% Outputs the "kick", or time-delta perturbation applied to compute optical response in real time.
    !%Option external_td_potential bit(26)
    !% Outputs the (scalar) time-dependent potential.
    !%Option potential_gradient bit(31)
    !% Prints the gradient of the potential.
    !%Option energy_density bit(32)
    !% Outputs the total energy density to a file called
    !% <tt>energy_density</tt>.
    !%Option heat_current bit(33)
    !% Outputs the total heat current density. The output file is
    !% called <tt>heat_current-</tt>.
    !%End
    call parse_variable('Output', 0, outp%what)

    ! cannot calculate the ELF in 1D
    if(bitand(outp%what, OPTION__OUTPUT__ELF) /= 0 .or. bitand(outp%what, OPTION__OUTPUT__ELF_BASINS) /= 0) then
       if(sb%dim /= 2 .and. sb%dim /= 3) then
         outp%what = bitand(outp%what, not(OPTION__OUTPUT__ELF + OPTION__OUTPUT__ELF_BASINS))
         write(message(1), '(a)') 'Cannot calculate ELF except in 2D and 3D.'
         call messages_warning(1)
       end if
    end if

    if(.not.varinfo_valid_option('Output', outp%what, is_flag=.true.)) then
      call messages_input_error('Output')
    end if

    if(bitand(outp%what, OPTION__OUTPUT__ENERGY_DENSITY) /= 0) call messages_experimental("'Output = energy_density'")
    if(bitand(outp%what, OPTION__OUTPUT__HEAT_CURRENT) /= 0) call messages_experimental("'Output = heat_current'")
    
    if(bitand(outp%what, OPTION__OUTPUT__WFS) /= 0  .or.  bitand(outp%what, OPTION__OUTPUT__WFS_SQMOD) /= 0 ) then

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
      call parse_variable('OutputWfsNumber', default, outp%wfs_list)
    end if

    if(parse_block('CurrentThroughPlane', blk) == 0) then
      outp%what = ior(outp%what, OPTION__OUTPUT__J_FLOW)

      !%Variable CurrentThroughPlane
      !%Type block
      !%Section Output
      !%Description
      !% The code can calculate current
      !% traversing a user-defined portion of a plane, as specified by this block.
      !% A small plain-text file <tt>current-flow</tt> will be written containing this information.
      !% Only available for 1D, 2D, or 3D.
      !% In the format below, <tt>origin</tt> is a point in the plane.
      !% <tt>u</tt> and <tt>v</tt> are the (dimensionless) vectors defining the plane;
      !% they will be normalized. <tt>spacing</tt> is the fineness of the mesh
      !% on the plane. Integers <tt>nu</tt> and <tt>mu</tt> are the length and
      !% width of the portion of the plane, in units of <tt>spacing</tt>.
      !% Thus, the grid points included in the plane are
      !% <tt>x_ij = origin + i*spacing*u + j*spacing*v</tt>,
      !% for <tt>nu <= i <= mu</tt> and <tt>nv <= j <= mv</tt>.
      !% Analogously, in the 2D case, the current flow is calculated through a line;
      !% in the 1D case, the current flow is calculated through a point. Note that the spacing
      !% can differ from the one used in the main calculation; an interpolation will be performed.
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
          write(message(1), '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        end if
        outp%plane%u(1:3) = outp%plane%u(1:3) / norm

        norm = sqrt(sum(outp%plane%v(1:3)**2))
        if(norm < M_EPSILON) then
          write(message(1), '(a)') 'v-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        end if
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
          write(message(1), '(a)') 'u-vector for CurrentThroughPlane cannot have norm zero.'
          call messages_fatal(1)
        end if
        outp%line%u(1:2) = outp%line%u(1:2) / norm

        outp%line%n(1) = -outp%line%u(2)
        outp%line%n(2) =  outp%line%u(1)

      case(1)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)

      case default

        call messages_not_implemented("CurrentThroughPlane for 4D or higher")

      end select
      call parse_block_end(blk)
    end if

    if(bitand(outp%what, OPTION__OUTPUT__MATRIX_ELEMENTS) /= 0) then
      call output_me_init(outp%me, sb)
    end if

    !%Variable OutputLDA_U
    !%Type flag
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print, related to LDA+U. 
    !% The output files are written at the end of the run into the output directory for the
    !% relevant kind of run (<i>e.g.</i> <tt>static</tt> for <tt>CalculationMode = gs</tt>).
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% (available for <tt>CalculationMode</tt> = <tt>gs</tt>, <tt>unocc</tt>,  <tt>td</tt>, and <tt>opt_control</tt>)
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir/effectiveU</tt>.
    !% For linear-response run modes, the derivatives of many quantities can be printed, as listed in
    !% the options below. Indices in the filename are labelled as follows:
    !% <tt>sp</tt> = spin (or spinor component), <tt>k</tt> = <i>k</i>-point, <tt>st</tt> = state/band.
    !% There is no tag for directions, given as a letter. The perturbation direction is always
    !% the last direction for linear-response quantities, and a following +/- indicates the sign of the frequency.
    !% Example: <tt>occ_matrices + effectiveU</tt>
    !%Option occ_matrices  bit(0)
    !% Outputs the occupation matrices of LDA+U
    !%Option effectiveU bit(1)
    !% Outputs the value of the effectiveU for each atoms 
    !%Option magnetization bit(2)
    !% Outputs file containing structure and magnetization of the localized subspace 
    !% on the atoms as a vector associated with each atom, which can be visualized.
    !% For the moment, it only works if a +U is added on one type of orbital per atom. 
    !%Option local_orbitals bit(3)
    !% Outputs the localized orbitals that form the correlated subspace
    !%Option kanamoriU bit(4)
    !% Outputs the Kanamori interaction parameters U, U`, and J.
    !% These parameters are not determined self-consistently, but are taken from the 
    !% occupation matrices and Coulomb integrals comming from a standard +U calculation.
    !%End
    call parse_variable('OutputLDA_U', 0_8, outp%what_lda_u)

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
    !%End
    call parse_variable('OutputInterval', 50, outp%output_interval)
    call messages_obsolete_variable("OutputEvery", "OutputInterval/RestartWriteInterval")
    if(outp%output_interval < 0) then
      message(1) = "OutputInterval must be >= 0."
      call messages_fatal(1)
    end if

    !%Variable OutputDuringSCF
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% During <tt>gs</tt> and <tt>unocc</tt> runs, if this variable is set to yes, 
    !% output will be written after every <tt>OutputInterval</tt> iterations.
    !%End
    call parse_variable('OutputDuringSCF', .false., outp%duringscf) 

    !%Variable RestartWriteInterval
    !%Type integer
    !%Default 50
    !%Section Execution::IO
    !%Description
    !% Restart data is written when the iteration number is a multiple
    !% of the <tt>RestartWriteInterval</tt> variable. For
    !% time-dependent runs this includes the update of the output
    !% controlled by the <tt>TDOutput</tt> variable. (Other output is
    !% controlled by <tt>OutputInterval</tt>.)
    !%End
    call parse_variable('RestartWriteInterval', 50, outp%restart_write_interval)
    if(outp%restart_write_interval <= 0) then
      message(1) = "RestartWriteInterval must be > 0."
      call messages_fatal(1)
    end if

    ! these kinds of Output do not have a how
    what_no_how = OPTION__OUTPUT__MATRIX_ELEMENTS + OPTION__OUTPUT__DOS + &
      OPTION__OUTPUT__TPA + OPTION__OUTPUT__J_FLOW
    what_no_how_u = OPTION__OUTPUTLDA_U__OCC_MATRICES + OPTION__OUTPUTLDA_U__EFFECTIVEU + &
      OPTION__OUTPUTLDA_U__MAGNETIZATION + OPTION__OUTPUTLDA_U__KANAMORIU

    if(bitand(outp%what, OPTION__OUTPUT__CURRENT) /= 0 .or. bitand(outp%what, OPTION__OUTPUT__HEAT_CURRENT) /= 0) then
      call v_ks_calculate_current(ks, .true.)
    else
      call v_ks_calculate_current(ks, .false.)
    end if

   
    !%Variable Output_KPT
    !%Type flag
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what to print. The output files are written at the end of the run into the output directory for the
    !% relevant kind of run (<i>e.g.</i> <tt>static</tt> for <tt>CalculationMode = gs</tt>).
    !% Time-dependent simulations print only per iteration, including always the last. The frequency of output per iteration
    !% (available for <tt>CalculationMode</tt> = <tt>gs</tt>, <tt>unocc</tt>,  <tt>td</tt>, and <tt>opt_control</tt>)
    !% is set by <tt>OutputInterval</tt> and the directory is set by <tt>OutputIterDir</tt>.
    !% For linear-response run modes, the derivatives of many quantities can be printed, as listed in
    !% the options below. Indices in the filename are labelled as follows:
    !% <tt>sp</tt> = spin (or spinor component), <tt>k</tt> = <i>k</i>-point, <tt>st</tt> = state/band.
    !% There is no tag for directions, given as a letter. The perturbation direction is always
    !% the last direction for linear-response quantities, and a following +/- indicates the sign of the frequency.
    !% Example: <tt>current_kpt</tt>
    !%Option current_kpt  bit(0)
    !% Outputs the current density resolved in momentum space. The output file is called <tt>current_kpt-</tt>.
    !%Option density_kpt bit(1)
    !% Outputs the electronic density resolved in momentum space. 
    !%End
    call parse_variable('Output_KPT', 0_8, outp%whatBZ)

    if(.not.varinfo_valid_option('Output_KPT', outp%whatBZ, is_flag=.true.)) then
      call messages_input_error('Output_KPT')
    end if

    if(bitand(outp%whatBZ, OPTION__OUTPUT_KPT__CURRENT_KPT) /= 0) then
     call v_ks_calculate_current(ks, .true.) 
    end if

    !%Variable OutputIterDir
    !%Default "output_iter"
    !%Type string
    !%Section Output
    !%Description
    !% The name of the directory where <tt>Octopus</tt> stores information
    !% such as the density, forces, etc. requested by variable <tt>Output</tt>
    !% in the format specified by <tt>OutputFormat</tt>.
    !% This information is written while iterating <tt>CalculationMode = gs</tt>, <tt>unocc</tt>, or <tt>td</tt>,
    !% according to <tt>OutputInterval</tt>, and has nothing to do with the restart information.
    !%End
    call parse_variable('OutputIterDir', "output_iter", outp%iter_dir)
    if(outp%what + outp%whatBZ + outp%what_lda_u /= 0 .and. outp%output_interval > 0) then
      call io_mkdir(outp%iter_dir)
    end if
    call add_last_slash(outp%iter_dir)

    ! we are using a what that has a how.
    if(bitand(outp%what, not(what_no_how)) /= 0 .or. outp%whatBZ /= 0 .or. bitand(outp%what_lda_u, not(what_no_how_u)) /= 0) then
      call io_function_read_how(sb, outp%how)
    else
      outp%how = 0
    end if


    POP_SUB(output_init)
  end subroutine output_init

  ! ---------------------------------------------------------

  subroutine output_end(outp)
    type(output_t), intent(inout) :: outp

    PUSH_SUB(output_end)
    
    POP_SUB(output_end)

  end subroutine output_end

  ! ---------------------------------------------------------
  subroutine output_all(outp, gr, geo, st, hm, ks, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    type(v_ks_t),         intent(in)    :: ks
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir

    integer :: idir, ierr
    character(len=80) :: fname
    type(profile_t), save :: prof
    
    PUSH_SUB(output_all)
    call profiling_in(prof, "OUTPUT_ALL")

    if(outp%what+outp%whatBZ+outp%what_lda_u /= 0) then
      message(1) = "Info: Writing output to " // trim(dir)
      call messages_info(1)
      call io_mkdir(dir)
    end if

    if(bitand(outp%what, OPTION__OUTPUT__MESH_R) /= 0) then
      do idir = 1, gr%mesh%sb%dim
        write(fname, '(a,a)') 'mesh_r-', index2axis(idir)
        call dio_function_output(outp%how, dir, fname, gr%mesh, gr%mesh%x(:,idir), &
          units_out%length, ierr, geo = geo)
      end do
    end if
    
    call output_states(st, gr, geo, hm, dir, outp)
    call output_hamiltonian(hm, st, gr%der, dir, outp, geo, gr, st%st_kpt_mpi_grp)
    call output_localization_funct(st, hm, gr, dir, outp, geo)
    call output_current_flow(gr, st, dir, outp)

    if(bitand(outp%what, OPTION__OUTPUT__GEOMETRY) /= 0) then
      if(bitand(outp%how, OPTION__OUTPUTFORMAT__XCRYSDEN) /= 0) then        
        call write_xsf_geometry_file(dir, "geometry", geo, gr%mesh)
      end if
      if(bitand(outp%how, OPTION__OUTPUTFORMAT__XYZ) /= 0) then
        call geometry_write_xyz(geo, trim(dir)//'/geometry')
        if(simul_box_is_periodic(gr%sb))  call periodic_write_crystal(gr%sb, geo, dir)
      end if
      if(bitand(outp%how, OPTION__OUTPUTFORMAT__VTK) /= 0) then
        call vtk_output_geometry(trim(dir)//'/geometry', geo)
      end if     
    end if

    if(bitand(outp%what, OPTION__OUTPUT__FORCES) /= 0) then
      if(bitand(outp%how, OPTION__OUTPUTFORMAT__BILD) /= 0) then
        call write_bild_forces_file(dir, "forces", geo, gr%mesh)
      else
        call write_xsf_geometry_file(dir, "forces", geo, gr%mesh, write_forces = .true.)
      end if
    end if

    if(bitand(outp%what, OPTION__OUTPUT__MATRIX_ELEMENTS) /= 0) then
      call output_me(outp%me, dir, st, gr, geo, hm)
    end if

    call output_energy_density(hm, ks, st, gr%der, dir, outp, geo, gr, st%st_kpt_mpi_grp)

    if(hm%lda_u_level /= DFT_U_NONE) then
      if(iand(outp%what_lda_u, OPTION__OUTPUTLDA_U__OCC_MATRICES) /= 0)&
        call lda_u_write_occupation_matrices(dir, hm%lda_u, st)

      if(iand(outp%what_lda_u, OPTION__OUTPUTLDA_U__EFFECTIVEU) /= 0)&
        call lda_u_write_effectiveU(dir, hm%lda_u)

      if(iand(outp%what_lda_u, OPTION__OUTPUTLDA_U__MAGNETIZATION) /= 0)&
        call lda_u_write_magnetization(dir, hm%lda_u, geo, gr%mesh, st)

      if(iand(outp%what_lda_u, OPTION__OUTPUTLDA_U__LOCAL_ORBITALS) /= 0)&
        call output_dftu_orbitals(dir, hm%lda_u, outp, st, gr%mesh, geo, associated(hm%hm_base%phase))

      if(iand(outp%what_lda_u, OPTION__OUTPUTLDA_U__KANAMORIU) /= 0)&
        call lda_u_write_kanamoriU(dir, st, hm%lda_u)
    end if

    call profiling_out(prof)
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
    
    mpi_grp = st%dom_st_kpt_mpi_grp

    ! if SPIN_POLARIZED, the ELF contains one extra channel: the total ELF
    imax = st%d%nspin
    if(st%d%ispin == SPIN_POLARIZED) imax = 3

    SAFE_ALLOCATE(f_loc(1:gr%mesh%np, 1:imax))

    ! First the ELF in real space
    if(bitand(outp%what, OPTION__OUTPUT__ELF) /= 0 .or. bitand(outp%what, OPTION__OUTPUT__ELF_BASINS) /= 0) then
      ASSERT(gr%mesh%sb%dim /= 1)

      call elf_calc(st, gr, f_loc)
      
      ! output ELF in real space
      if(bitand(outp%what, OPTION__OUTPUT__ELF) /= 0) then
        write(fname, '(a)') 'elf_rs'
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,imax), unit_one, ierr, geo = geo, grp = mpi_grp)
        ! this quantity is dimensionless

        if(st%d%ispin /= UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,i1)') 'elf_rs-sp', is
            call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
              f_loc(:, is), unit_one, ierr, geo = geo, grp = mpi_grp)
            ! this quantity is dimensionless
          end do
        end if
      end if

      if(bitand(outp%what, OPTION__OUTPUT__ELF_BASINS) /= 0) &
        call out_basins(f_loc(:,1), "elf_rs_basins")
    end if

    ! Now Bader analysis
    if(bitand(outp%what, OPTION__OUTPUT__BADER) /= 0) then
      do is = 1, st%d%nspin
        call dderivatives_lapl(gr%der, st%rho(:,is), f_loc(:,is))
        write(fname, '(a,i1)') 'bader-sp', is
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,is), units_out%length**(-2 - gr%sb%dim), ierr, &
          geo = geo, grp = mpi_grp)

        write(fname, '(a,i1)') 'bader_basins-sp', is
        call out_basins(f_loc(:,1), fname)
      end do
    end if

    ! Now the pressure
    if(bitand(outp%what, OPTION__OUTPUT__EL_PRESSURE) /= 0) then
      call calc_electronic_pressure(st, hm, gr, f_loc(:,1))
      call dio_function_output(outp%how, dir, "el_pressure", gr%mesh, &
        f_loc(:,1), unit_one, ierr, geo = geo, grp = mpi_grp)
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
      call basins_analyze(basins, gr%mesh, ff(:), st%rho, CNST(0.01))

      call dio_function_output(outp%how, dir, trim(filename), gr%mesh, &
        real(basins%map, REAL_PRECISION), unit_one, ierr, geo = geo, grp = mpi_grp)
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
    call states_calc_quantities(gr%der, st, .false., kinetic_energy_density = tau)

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
  subroutine output_energy_density(hm, ks, st, der, dir, outp, geo, gr, grp)
    type(hamiltonian_t),       intent(in)    :: hm
    type(v_ks_t),              intent(in)    :: ks
    type(states_t),            intent(inout) :: st
    type(derivatives_t),       intent(inout) :: der
    character(len=*),          intent(in)    :: dir
    type(output_t),            intent(in)    :: outp
    type(geometry_t),          intent(in)    :: geo
    type(grid_t),              intent(in)    :: gr
    type(mpi_grp_t), optional, intent(in)    :: grp !< the group that shares the same data, must contain the domains group

    integer :: is, ispin, ierr, ip
    character(len=MAX_PATH_LEN) :: fname
    type(unit_t) :: fn_unit
    FLOAT, allocatable :: energy_density(:, :)
    FLOAT, allocatable :: ex_density(:)
    FLOAT, allocatable :: ec_density(:)

    PUSH_SUB(output_hamiltonian)
   
    if(bitand(outp%what, OPTION__OUTPUT__ENERGY_DENSITY) /= 0) then
      fn_unit = units_out%energy*units_out%length**(-gr%mesh%sb%dim)
      SAFE_ALLOCATE(energy_density(1:gr%mesh%np, 1:st%d%nspin))

      ! the kinetic energy density
      call states_calc_quantities(gr%der, st, .true., kinetic_energy_density = energy_density)

      ! the external potential energy density
      forall(ip = 1:gr%fine%mesh%np, is = 1:st%d%nspin) energy_density(ip, is) = energy_density(ip, is) + st%rho(ip, is)*hm%ep%vpsl(ip)

      ! the hartree energy density
      forall(ip = 1:gr%fine%mesh%np, is = 1:st%d%nspin) energy_density(ip, is) = energy_density(ip, is) + CNST(0.5)*st%rho(ip, is)*hm%vhartree(ip)

      ! the XC energy density
      SAFE_ALLOCATE(ex_density(1:gr%mesh%np))
      SAFE_ALLOCATE(ec_density(1:gr%mesh%np))

      ASSERT(.not. gr%have_fine_mesh)

      call xc_get_vxc(gr%fine%der, ks%xc, st, st%rho, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, ex_density = ex_density, ec_density = ec_density)
      forall(ip = 1:gr%fine%mesh%np, is = 1:st%d%nspin) energy_density(ip, is) = energy_density(ip, is) + ex_density(ip) + ec_density(ip)
      
      SAFE_DEALLOCATE_A(ex_density)
      SAFE_DEALLOCATE_A(ec_density)
      
      select case(st%d%ispin)
      case(UNPOLARIZED)
        write(fname, '(a)') 'energy_density'
        call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
          energy_density(:,1), unit_one, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
      case(SPIN_POLARIZED, SPINORS)
        do is = 1, 2
          write(fname, '(a,i1)') 'energy_density-sp', is
          call dio_function_output(outp%how, dir, trim(fname), gr%mesh, &
            energy_density(:, is), unit_one, ierr, geo = geo, grp = st%dom_st_kpt_mpi_grp)
        end do
      end select
      SAFE_DEALLOCATE_A(energy_density)
    end if
 
    POP_SUB(output_hamiltonian)
  end subroutine output_energy_density

   ! ---------------------------------------------------------
  subroutine output_dftu_orbitals(dir, this, outp, st, mesh, geo, has_phase)
    character(len=*),  intent(in) :: dir
    type(lda_u_t),     intent(in) :: this
    type(output_t),    intent(in) :: outp
    type(states_t),    intent(in) :: st
    type(mesh_t),      intent(in) :: mesh
    type(geometry_t),  intent(in) :: geo
    logical,           intent(in) :: has_phase

    integer :: ios, im, ik, idim, ierr
    CMPLX, allocatable :: tmp(:)
    FLOAT, allocatable :: dtmp(:)
    type(orbitalset_t), pointer :: os
    type(unit_t) :: fn_unit
    character(len=32) :: fname

    PUSH_SUB(output_dftu_orbitals)

    fn_unit = sqrt(units_out%length**(-mesh%sb%dim))

    if(.not.(has_phase .and. .not.this%basis%submeshforperiodic &
           .and.simul_box_is_periodic(mesh%sb)).and. .not. this%basisfromstates) then
      if(states_are_real(st)) then
        SAFE_ALLOCATE(dtmp(1:mesh%np))
      else
        SAFE_ALLOCATE(tmp(1:mesh%np))
      end if
    end if

    do ios = 1, this%norbsets
      os => this%orbsets(ios)
      do ik = st%d%kpt%start, st%d%kpt%end
        do im = 1, this%orbsets(ios)%norbs
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i1,a,i3.3,a,i4.4,a,i1)') 'orb', im, '-os', ios, '-k', ik, '-sp', idim
              else
                write(fname, '(a,i1,a,i3.3,a,i4.4)') 'orb', im, '-os', ios, '-k', ik
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i1,a,i3.3,a,i1)') 'orb', im, '-os', ios, '-sp', idim
              else
                write(fname, '(a,i1,a,i3.3)') 'orb', im, '-os', ios
              end if
            end if
            if(has_phase) then
              if(simul_box_is_periodic(mesh%sb) .and. .not. this%basis%submeshforperiodic) then
               call zio_function_output(outp%how, dir, fname, mesh, &
                  os%eorb_mesh(1:mesh%np,im,idim,ik), fn_unit, ierr, geo = geo)
              else
               tmp = M_Z0
               call submesh_add_to_mesh(os%sphere, os%eorb_submesh(1:os%sphere%np,idim,im,ik), tmp)
               call zio_function_output(outp%how, dir, fname, mesh, tmp, fn_unit, ierr, geo = geo)
              end if
            else
              if(this%basisfromstates) then
                if (states_are_real(st)) then
                  call dio_function_output(outp%how, dir, fname, mesh, &
                      os%dorb(1:mesh%np,idim,im), fn_unit, ierr, geo = geo)
                else
                  call zio_function_output(outp%how, dir, fname, mesh, &
                      os%zorb(1:mesh%np,idim,im), fn_unit, ierr, geo = geo)
                end if
              else
                if (states_are_real(st)) then
                  dtmp = M_Z0
                  call submesh_add_to_mesh(os%sphere, os%dorb(1:os%sphere%np,idim,im), dtmp)
                  call dio_function_output(outp%how, dir, fname, mesh, dtmp, fn_unit, ierr, geo = geo)
                else
                  tmp = M_Z0
                  call submesh_add_to_mesh(os%sphere, os%zorb(1:os%sphere%np,idim,im), tmp)
                  call zio_function_output(outp%how, dir, fname, mesh, tmp, fn_unit, ierr, geo = geo)
                end if
              end if
            end if
          end do
        end do
      end do
    end do

    if(.not.(has_phase .and. .not.this%basis%submeshforperiodic &
               .and.simul_box_is_periodic(mesh%sb)).and. .not. this%basisfromstates) then
      SAFE_DEALLOCATE_A(tmp)
      SAFE_DEALLOCATE_A(dtmp)
    end if

    POP_SUB(output_dftu_orbitals)
  end subroutine output_dftu_orbitals

  ! ---------------------------------------------------------

#include "output_states_inc.F90"

#include "output_h_inc.F90"

end module output_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
