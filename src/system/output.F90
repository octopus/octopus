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
  use batch_oct_m
  use comm_oct_m 
  use cube_function_oct_m
  use cube_oct_m
  use current_oct_m
  use density_oct_m
  use derivatives_oct_m
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
  use loct_oct_m
  use loct_math_oct_m
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
    output_kick,         &
    output_scalar_pot

  type output_t
    !> General output variables:
    integer(8) :: what                !< what to output
    integer(8) :: whatBZ              !< what to output - for k-point resolved output
    integer(8) :: how                 !< how to output

    type(output_me_t) :: me        !< this handles the output of matrix elements

    !> These variables fine-tune the output for some of the possible output options:
    integer :: output_interval     !< output every iter
    integer :: restart_write_interval
    logical :: duringscf
    character(len=80) :: wfs_list  !< If output_wfs, this list decides which wavefunctions to print.
    character(len=MAX_PATH_LEN) :: iter_dir  !< The folder name, if information will be output while iterating.
  end type output_t

contains

  subroutine output_init(outp, sb, nst, ks)
    type(output_t),       intent(out)   :: outp
    type(simul_box_t),    intent(in)    :: sb
    integer,              intent(in)    :: nst
    type(v_ks_t),         intent(inout) :: ks

    type(block_t) :: blk
    FLOAT :: norm
    character(len=80) :: nst_string, default
    integer :: what_no_how

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
    !%Option forces bit(18)
    !% Outputs file <tt>forces.xsf</tt> containing structure and forces on the atoms as 
    !% a vector associated with each atom, which can be visualized with XCrySDen.
    !%Option delta_perturbation bit(25)
    !% Outputs the "kick", or time-delta perturbation applied to compute optical response in real time.
    !%Option external_td_potential bit(26)
    !% Outputs the (scalar) time-dependent potential.
    !%Option potential_gradient bit(31)
    !% Prints the gradient of the potential.
    !%End
    call parse_variable('Output', 0, outp%what)

    if(.not.varinfo_valid_option('Output', outp%what, is_flag=.true.)) then
      call messages_input_error('Output')
    end if

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

    if(bitand(outp%what, OPTION__OUTPUT__MATRIX_ELEMENTS) /= 0) then
      call output_me_init(outp%me, sb)
    end if

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
    what_no_how = OPTION__OUTPUT__MATRIX_ELEMENTS
   
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
    if(outp%what + outp%whatBZ /= 0 .and. outp%output_interval > 0) then
      call io_mkdir(outp%iter_dir)
    end if
    call add_last_slash(outp%iter_dir)

    ! we are using a what that has a how.
    if(bitand(outp%what, not(what_no_how)) /= 0 .or. outp%whatBZ /= 0) then
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

    if(outp%what /= 0) then
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

    call profiling_out(prof)
    POP_SUB(output_all)
  end subroutine output_all

#include "output_states_inc.F90"

#include "output_h_inc.F90"

end module output_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
