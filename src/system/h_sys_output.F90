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
!! $Id: h_sys_output.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

module h_sys_output_m
  use basins_m
  use cube_function_m
  use datasets_m
  use density_m
  use derivatives_m
  use elf_m
#if defined(HAVE_ETSF_IO)
  use etsf_io
#endif
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
  use young_m

  implicit none

  private
  public ::                    &
    h_sys_output_t,            &
    h_sys_output_init,         &
    h_sys_output_states,       &
    h_sys_output_modelmb,      &
    h_sys_output_hamiltonian,  &
    h_sys_output_all,          &
    h_sys_output_current_flow, &
    dh_sys_output_lr,          &
    zh_sys_output_lr

  type h_sys_output_t
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

  end type h_sys_output_t

  integer, parameter, public ::           &
    output_potential       =        1,    &
    output_density         =        2,    &
    output_wfs             =        4,    &
    output_wfs_sqmod       =        8,    &
    output_geometry        =       16,    &
    output_current         =       32,    &
    output_ELF             =       64,    &
    output_ELF_basins      =      128,    &
    output_ELF_FS          =      256,    &
    output_Bader           =      512,    &
    output_el_pressure     =     1024,    &
    output_matrix_elements =     2048,    &
    output_pol_density     =     4096,    &
    output_r               =     8192,    &
    output_ked             =    16384,    &
    output_j_flow          =    32768,    &
    output_dos             =    65536,    &
    output_tpa             =   131072,    &
    output_density_matrix  =   262144,    &
    output_modelmb         =   524288,    &
    output_forces          =  1048576

contains

  subroutine h_sys_output_init(sb, nst, outp)
    type(simul_box_t),    intent(in)  :: sb
    integer,              intent(in)  :: nst
    type(h_sys_output_t), intent(out) :: outp

    type(block_t) :: blk
    FLOAT :: norm
    character(len=80) :: nst_string, default

    PUSH_SUB(h_sys_output_init)

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
    !%End
    call parse_integer(datasets_check('Output'), 0, outp%what)

    if(iand(outp%what, output_elf_fs) .ne. 0) then
      call messages_experimental("ELF in Fourier space")
    endif

    ! cannot calculate the ELF in 1D
    if(iand(outp%what, output_elf) .ne. 0 .or. iand(outp%what, output_elf_basins) .ne. 0 &
       .or. iand(outp%what, output_elf_fs) .ne. 0) then
       if(sb%dim .ne. 2 .and. sb%dim .ne. 3) then
         outp%what = iand(outp%what, not(output_elf + output_elf_basins + output_elf_fs))
         write(message(1), '(a)') 'Cannot calculate ELF except in 2D and 3D.'
         call write_warning(1)
       endif
    endif

    if(.not.varinfo_valid_option('Output', outp%what, is_flag=.true.)) then
      call input_error('Output')
    end if

    if(iand(outp%what, output_modelmb) .ne. 0 .or. iand(outp%what, output_density_matrix) .ne. 0) then
      call messages_experimental("Model many-body and density matrix")
    endif

    if(iand(outp%what, output_modelmb) .ne. 0) then
      write(message(1),'(a)') 'Model many-body quantities will be output, according to the presence of'
      write(message(2),'(a)') '  wfs, density, or density_matrix in Output.'
      call write_info(2)
    end if

    if(iand(outp%what, output_density_matrix) .ne. 0) then
      write(message(1),'(a)') 'Info: The density matrix will be calculated, traced'
      write(message(2),'(a)') 'over the second dimension, diagonalized, and output.'
      call write_info(2)
      if(iand(outp%what, output_modelmb) .eq. 0) then
        write(message(1),'(a)') 'Note that density matrix only works for model MB calculations for the moment.'
        call write_info(1)
      end if
      ! NOTES:
      !   could be made into block to be able to specify which dimensions to trace
      !   in principle all combinations are interesting, but this means we need to
      !   be able to output density matrices for multiple particles or multiple
      !   dimensions. The current 1D 1-particle case is simple.
    end if

    if(iand(outp%what, output_wfs) .ne. 0  .or.  iand(outp%what, output_wfs_sqmod) .ne. 0 ) then

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
      outp%what = ior(outp%what, output_j_flow)

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
          call write_fatal(1)
        endif
        outp%plane%u(1:3) = outp%plane%u(1:3) / norm

        norm = sqrt(sum(outp%plane%v(1:3)**2))
        if(norm < M_EPSILON) then
          write(1, '(a)') 'v-vector for CurrentThroughPlane cannot have norm zero.'
          call write_fatal(1)
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
          call write_fatal(1)
        endif
        outp%line%u(1:2) = outp%line%u(1:2) / norm

        outp%line%n(1) = -outp%line%u(2)
        outp%line%n(2) =  outp%line%u(1)

      case(1)

        call parse_block_float(blk, 0, 0, outp%line%origin(1), units_inp%length)

      end select
    end if

    if(iand(outp%what, output_matrix_elements) .ne. 0) then
      call output_me_init(outp%me, sb)
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

    if(outp%what .ne. 0 .and. outp%what .ne. output_matrix_elements) call io_function_read_how(sb, outp%how)

    POP_SUB(h_sys_output_init)
  end subroutine h_sys_output_init


  ! ---------------------------------------------------------
  subroutine h_sys_output_all(outp, gr, geo, st, hm, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    type(h_sys_output_t), intent(in)    :: outp
    character(len=*),     intent(in)    :: dir

    integer :: idir
    FLOAT   :: offset(1:MAX_DIM)
    
    PUSH_SUB(h_sys_output_all)
    
    call h_sys_output_states(st, gr, geo, dir, outp)
    call h_sys_output_hamiltonian(hm, gr%mesh, dir, outp, geo)
    call h_sys_output_localization_funct(st, hm, gr, dir, outp, geo)
    call h_sys_output_current_flow(gr, st, dir, outp, geo)

    if(iand(outp%what, output_geometry) .ne. 0) then
      if(iand(outp%how, output_xcrysden) .ne. 0) then
        
        offset = M_ZERO
        ! The corner of the cell is always (0,0,0) to XCrySDen
        ! so the offset is applied to the atomic coordinates.
        ! Offset in periodic directions:
        offset(1:3) = -matmul(gr%mesh%sb%rlattice_primitive(1:3,1:3), gr%mesh%sb%lsize(1:3))
        ! Offset in aperiodic directions:
        do idir = gr%mesh%sb%periodic_dim + 1, 3
          offset(idir) = -(gr%mesh%idx%ll(idir) - 1)/2*gr%mesh%spacing(idir)
        end do

        call write_xsf_geometry_file(dir, "geometry", geo, gr%sb, offset = offset)
      endif
      if(iand(outp%how, output_xyz) .ne. 0) then
        call atom_write_xyz(dir, "geometry", geo, gr%sb%dim)
        if(simul_box_is_periodic(gr%sb)) &
          call periodic_write_crystal(gr%sb, geo, dir)
      endif
    end if

    if(iand(outp%what, output_forces) .ne. 0) then
      call write_xsf_geometry_file(dir, "forces", geo, gr%sb, write_forces = .true.)
    endif

    if(iand(outp%what, output_matrix_elements) .ne. 0) then
      call output_me(outp%me, dir, st, gr, geo, hm)
    end if

    if (iand(outp%how, output_etsf) .ne. 0) then
      call h_sys_output_etsf(st, gr, geo, dir, outp)
    end if
    
    POP_SUB(h_sys_output_all)
  end subroutine h_sys_output_all


  ! ---------------------------------------------------------
  subroutine h_sys_output_localization_funct(st, hm, gr, dir, outp, geo)
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: dir
    type(h_sys_output_t),   intent(in)    :: outp
    type(geometry_t),       intent(in)    :: geo

    FLOAT, allocatable :: f_loc(:,:)
    character(len=256) :: fname
    integer :: is, ierr, imax
    type(mpi_grp_t) :: mpi_grp

    PUSH_SUB(h_sys_output_localization_funct)
    
    mpi_grp = gr%mesh%mpi_grp
    if(st%parallel_in_states) mpi_grp = st%mpi_grp
    if(st%d%kpt%parallel) mpi_grp = st%d%kpt%mpi_grp

    ! if SPIN_POLARIZED, the ELF contains one extra channel: the total ELF
    imax = st%d%nspin
    if(st%d%ispin == SPIN_POLARIZED) imax = 3

    SAFE_ALLOCATE(f_loc(1:gr%mesh%np, 1:imax))

    ! First the ELF in real space
    if(iand(outp%what, output_elf) .ne. 0 .or. iand(outp%what, output_elf_basins) .ne. 0) then
      ASSERT(gr%mesh%sb%dim .ne. 1)

      call elf_calc(st, gr, f_loc)
      
      ! output ELF in real space
      if(iand(outp%what, output_elf) .ne. 0) then
        write(fname, '(a)') 'elf_rs'
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,imax), unit_one, ierr, is_tmp = .false., geo = geo)
        ! this quantity is dimensionless

        if(st%d%ispin .ne. UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,i1)') 'elf_rs-sp', is
            call doutput_function(outp%how, dir, trim(fname), gr%mesh, &
              f_loc(:, is), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
            ! this quantity is dimensionless
          end do
        end if
      end if

      if(iand(outp%what, output_elf_basins) .ne. 0) &
        call out_basins(f_loc(:,1), "elf_rs_basins")
    end if

    ! Second, ELF in Fourier space.
    if(iand(outp%what, output_elf_fs) .ne. 0) then
      call elf_calc_fs(st, gr, f_loc)
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'elf_fs-sp', is
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,is), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
        ! this quantity is dimensionless
      end do
    end if

    ! Now Bader analysis
    if(iand(outp%what, output_bader) .ne. 0) then
      do is = 1, st%d%nspin
        call dderivatives_lapl(gr%der, st%rho(:,is), f_loc(:,is))
        write(fname, '(a,i1)') 'bader-sp', is
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, &
          f_loc(:,is), units_out%length**(-2 - gr%sb%dim), ierr, is_tmp = .false., &
          geo = geo, grp = mpi_grp)

        write(fname, '(a,i1)') 'bader_basins-sp', is
        call out_basins(f_loc(:,1), fname)
      end do
    end if

    ! Now the pressure
    if(iand(outp%what, output_el_pressure) .ne. 0) then
      call h_sys_calc_electronic_pressure(st, hm, gr, f_loc(:,1))
      call doutput_function(outp%how, dir, "el_pressure", gr%mesh, &
        f_loc(:,1), unit_one, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
      ! this quantity is dimensionless
    end if

    SAFE_DEALLOCATE_A(f_loc)

    POP_SUB(h_sys_output_localization_funct)

  contains
    ! ---------------------------------------------------------
    subroutine out_basins(ff, filename)
      FLOAT,            intent(in)    :: ff(:)
      character(len=*), intent(in)    :: filename

      character(len=256) :: fname
      type(basins_t)     :: basins
      integer            :: iunit

      PUSH_SUB(h_sys_output_localization_funct.out_basins)

      call basins_init(basins, gr%mesh)
      call basins_analyze(basins, gr%mesh, st%d%nspin, ff(:), st%rho, CNST(0.01))

      call doutput_function(outp%how, dir, trim(filename), gr%mesh, &
        real(basins%map, REAL_PRECISION), unit_one, ierr, is_tmp = .false., geo = geo)
      ! this quantity is dimensionless

      write(fname,'(4a)') trim(dir), '/', trim(filename), '.info'
      iunit = io_open(file=trim(fname), action = 'write')
      call basins_write(basins, gr%mesh, iunit)
      call io_close(iunit)

      call basins_end(basins)

      POP_SUB(h_sys_output_localization_funct.out_basins)
    end subroutine out_basins

  end subroutine h_sys_output_localization_funct

  
  ! ---------------------------------------------------------
  subroutine h_sys_calc_electronic_pressure(st, hm, gr, pressure)
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(grid_t),           intent(inout) :: gr
    FLOAT,                  intent(out)   :: pressure(:)

    FLOAT, allocatable :: rho(:,:), lrho(:), tau(:,:)
    FLOAT   :: p_tf, dens
    integer :: is, ii

    PUSH_SUB(h_sys_calc_electronic_pressure)

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

    POP_SUB(h_sys_calc_electronic_pressure)
  end subroutine h_sys_calc_electronic_pressure


#include "output_etsf_inc.F90"

#include "output_states_inc.F90"
#include "output_h_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_linear_response_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "output_linear_response_inc.F90"

end module h_sys_output_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
