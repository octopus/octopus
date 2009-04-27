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
!! $Id: out.F90 3613 2007-11-29 16:47:41Z xavier $

#include "global.h"

module h_sys_output_m
  use datasets_m
  use density_matrix_m
  use derivatives_m
  use basins_m
  use cube_function_m
  use elf_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_function_m
  use io_m
  use linear_response_m
  use loct_m
  use loct_math_m
  use loct_parser_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_dim_m
  use solids_m
  use mesh_m
  use units_m
  use varinfo_m
#if defined(HAVE_ETSF_IO)
  use etsf_io
#endif
  implicit none

  private
  public ::                    &
    h_sys_output_t,            &
    h_sys_output_init,         &
    h_sys_output_states,       &
    h_sys_output_hamiltonian,  &
    h_sys_output_all,          &
    h_sys_output_current_flow, &
    dh_sys_output_lr,          &
    zh_sys_output_lr

  type h_sys_output_t
    ! General output variables:
    integer :: what                ! what to output
    integer :: how                 ! how to output

    ! This variables fine-tune the output for some of the possible output options:
    integer :: iter                ! output every iter
    logical :: duringscf

    character(len=80) :: wfs_list  ! If output_wfs, this list decides which wavefunctions to print.
    integer :: ksmultipoles        ! If output_ksdipole, this number sets up which matrix elements will
                                   ! be printed: e.g. if ksmultipoles = 3, the dipole, quadrupole and 
                                   ! octopole matrix elements (between Kohn-Sham or single-particle orbitals).

    type(mesh_plane_t) :: plane    ! This is to calculate the current flow across a plane
    type(mesh_line_t)  :: line     ! or though a line (in 2D)

  end type h_sys_output_t

  integer, parameter, public ::               &
    output_potential      =     1,    &
    output_density        =     2,    &
    output_wfs            =     4,    &
    output_wfs_sqmod      =     8,    &
    output_geometry       =    16,    &
    output_current        =    32,    &
    output_ELF            =    64,    &
    output_ELF_basins     =   128,    &
    output_ELF_FS         =   256,    &
    output_Bader          =   512,    &
    output_el_pressure    =  1024,    &
    output_ksdipole       =  2048,    &
    output_pol_density    =  4096,    &
    output_r              =  8192,    &
    output_ked            = 16384,    &
    output_j_flow         = 32768,    &
    output_dos            = 65536,    &
    output_tpa            =131072,    &
    output_density_matrix =262144

contains

  subroutine h_sys_output_init(sb, outp)
    type(simul_box_t),    intent(in)  :: sb
    type(h_sys_output_t), intent(out) :: outp

    type(block_t) :: blk

    !%Variable Output
    !%Type flag
    !%Default no
    !%Section Output
    !%Description
    !% Specifies what to print. The output files go into the "static" directory, except when
    !% running a time-dependent simulation, when the directory "td.XXXXXXX" is used.
    !% Example: "density + potential"
    !%Option potential 1
    !% Prints out Kohn-Sham potential, separated by parts. File names would be "v0" for 
    !% the local part, "vc" for the classical potential (if it exists), "vh" for the
    !% Hartree potential, and "vxc-x" for each of the exchange and correlation potentials
    !% of a give spin channel, where "x" stands for the spin channel.
    !%Option density 2
    !% Prints out the density. The output file is called "density-i", where "i" stands for 
    !% the spin channel. In a linear-response run mode, instead the linear-response density
    !% is printed, in a file "linear/freq_XXX/lr_density-i-j" where j is the perturbation
    !% direction. 
    !%Option wfs 4
    !% Prints out wave-functions. Which wavefunctions are to be printed is specified
    !% by the variable <tt>OutputWfsNumber</tt> -- see below. The output file is called
    !% "wf-k-p-i", where k stands for the <i>k</i> number, p for the state, and
    !% i for the spin channel.
    !%Option wfs_sqmod 8
    !% Prints out squared modulus of the wave-functions. 
    !% The output file is called "sqm-wf-k-p-i",
    !% where k stands for the <i>k</i> number, p for the state,
    !% and i for the spin channel.
    !%Option geometry 16
    !% Outputs a XYZ file called "geometry.xyz" containing the coordinates of the atoms
    !% treated within Quantum Mechanics. If point charges were defined
    !% in the PDB file (see <tt>PDBCoordinates</tt>), they will be output
    !% in the file "geometry_classical.xyz".
    !%Option current 32
    !% Prints out the paramagnetic current density. The output file is called "current-i-(x,y,z)", 
    !% where "i" stands for the spin channel and x, y, z indicates the vector component
    !% of the current.
    !%Option ELF 64
    !% Prints out the electron localization function, ELF. The output file is called
    !% "elf-i", where i stands for the spin channel.
    !%Option ELF_basins 128
    !% Prints out the basins of attraction of the ELF. The output file is called
    !% "elf_rs_basins.info".
    !%Option ELF_FS 256
    !% Prints the electron localization function in Fourier space. The output file is called
    !% "elf_FS-i", where i stands for the spin channel. (EXPERIMENTAL)
    !%Option Bader 512
    !% Prints the Laplacian of the density which shows lone pairs, bonded charge concentrations
    !% and regions subject to electrophilic or nucleophilic attack.
    !% See Bader, RF Atoms in Molecules: A Quantum Theory (Oxford, Oxford, 1990)
    !%Option el_pressure 1024
    !% Prints the electronic pressure. See Tao, Vignale, and Tokatly, Phys Rev Lett 100, 206405
    !%Option ksdipole 2048
    !% Prints out the multipole matrix elements between Kohn-Sham states (or just the single-
    !% particle states, in independent-electrons mode). Note that despite the name
    !% ("ksdipole"), the program may print higher-order multipoles (the order can be
    !% set with the variable OutputMatrixElementsL).
    !%Option pol_density 4096
    !% Prints out the density of dipole moment. For pol and pol_lr modules, 
    !% prints the density of polarizability, in the linear directory.
    !%Option mesh_r 8192
    !% Prints out the values of the coordinates over the grid. Files
    !% will be in the 'exec/' directory.
    !%Option kinetic_energy_density 16384
    !% Prints out the kinetic energy density, defined as:
    !%
    !% <math>\tau_\sigma(\vec{r}) = \sum_{i=1}^{N_\sigma} 
    !%  \vert \nabla \phi_{i\sigma}(\vec{r}) \vert^2\,. </math>
    !%
    !% The index <math>\sigma</math> is the spin index for the spin-polarized case,
    !% or if you are using spinors. For spin-unpolarized calculations, you
    !% get the total kinetic-energy density. The previous expression assumes full 
    !% or null occupations. If fractional occupation numbers, each term in the sum
    !% is weighted by the occupation. Also, if we are working with an infinite 
    !% system, all k-points are summed up, with their corresponding weights. The
    !% files will be called "tau-1" and "tau-2", if the spin-resolved kinetic
    !% energy density is produced (runs in spin-polarized and spinors mode), or
    !% only "tau" if the run is in spin-unpolarized mode.
    !%Option dos 65536
    !% Prints out the density of states.
    !%Option tpa 131072
    !% Prints out the transition potential approximation (tpa) matrix elements
    !%Option density_matrix 262144
    !% Calculates, and outputs, the reduced density
    !% matrix. For the moment the trace is made over the second dimension, and
    !% the code is limited to 2D. The idea is to model N particles in 1D as a N
    !% dimensional non-interacting problem, then to trace out N-1 coordinates.
    !% 
    !% WARNING: NOT TESTED YET.
    !%End
    call loct_parse_int(datasets_check('Output'), 0, outp%what)

    ! cannot calculate the ELF in 1D
    if(sb%dim == 1) outp%what = iand(outp%what, not(output_elf + output_elf_basins))

    if(.not.varinfo_valid_option('Output', outp%what, is_flag=.true.)) then
      call input_error('Output')
    end if

    if(iand(outp%what, output_density_matrix).ne.0) then
      write(message(1),'(a)') 'Info: The density matrix will be calculated, traced'
      write(message(2),'(a)') 'over the second dimension, diagonalized, and output.'
      call write_info(2)
      ! NOTES:
      !   could be made into block to be able to specify which dimensions to trace
      !   in principle all combinations are interesting, but this means we need to
      !   be able to output density matrices for multiple particles or multiple
      !   dimensions. The current 1D 1-particle case is simple.
    end if

    if(iand(outp%what, output_wfs).ne.0  .or.  iand(outp%what, output_wfs_sqmod).ne.0 ) then
      !%Variable OutputWfsNumber
      !%Type string
      !%Default "1-1024"
      !%Section Output
      !%Description
      !% Which wavefunctions to print, in list form: i.e., "1-5" to print the first
      !% five states, "2,3" to print the second and the third state, etc.
      !%End
      call loct_parse_string(datasets_check('OutputWfsNumber'), "1-1024", outp%wfs_list)
    end if

    if(iand(outp%what, output_ksdipole).ne.0) then
      !%Variable OutputMatrixElementsL
      !%Type integer
      !%Default 1
      !%Section Output
      !%Description
      !% If Output contains ksdipole, then this variable decides which matrix
      !% elements are printed out: e.g., if OutputMatrixElementsL = 1, then the
      !% program will plot three files, matrix_elements.x (x=1,2,3), containing
      !% respectively the (1,-1), (1,0) and (1,1) multipole matrix elements
      !% between Kohn-Sham states.
      !%End
      call loct_parse_int(datasets_check('OutputMatrixElementsL'), 1, outp%ksmultipoles)
    end if

    if(loct_parse_block(datasets_check('CurrentThroughPlane'), blk) == 0) then
      outp%what = ior(outp%what, output_j_flow)

      select case(sb%dim)
      case(3)

        !%Variable CurrentThroughPlane
        !%Type block
        !%Section States
        !%Description
        !% At the end of the ground state calculation, the code may calculate
        !% the steady current that the obtained ground state electronic state
        !% transverses a user-defined portion of a plane....
        !%
        !% In the 2D case, the current flow should be calculated through a line.
        !%
        !% Example (3D):
        !%
        !% <tt>%CurrentThroughPlane
        !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 0.0  # origin
        !% <br>&nbsp;&nbsp; 0.0 | 1.0 | 0.0  # u
        !% <br>&nbsp;&nbsp; 0.0 | 0.0 | 1.0  # v
        !% <br>&nbsp;&nbsp; 0.2              # spacing
        !% <br>&nbsp;&nbsp; 0 | 50           # nu | mu
        !% <br>&nbsp;&nbsp; -50 | 50         # nm | mv
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
        
        call loct_parse_block_float(blk, 0, 0, outp%plane%origin(1))
        call loct_parse_block_float(blk, 0, 1, outp%plane%origin(2))
        call loct_parse_block_float(blk, 0, 2, outp%plane%origin(3))
        call loct_parse_block_float(blk, 1, 0, outp%plane%u(1))
        call loct_parse_block_float(blk, 1, 1, outp%plane%u(2))
        call loct_parse_block_float(blk, 1, 2, outp%plane%u(3))
        call loct_parse_block_float(blk, 2, 0, outp%plane%v(1))
        call loct_parse_block_float(blk, 2, 1, outp%plane%v(2))
        call loct_parse_block_float(blk, 2, 2, outp%plane%v(3))
        call loct_parse_block_float(blk, 3, 0, outp%plane%spacing)
        call loct_parse_block_int(blk, 4, 0, outp%plane%nu)
        call loct_parse_block_int(blk, 4, 1, outp%plane%mu)
        call loct_parse_block_int(blk, 5, 0, outp%plane%nv)
        call loct_parse_block_int(blk, 5, 1, outp%plane%mv)

        outp%plane%n(1) = outp%plane%u(2)*outp%plane%v(3) - outp%plane%u(3)*outp%plane%v(2)
        outp%plane%n(2) = outp%plane%u(3)*outp%plane%v(1) - outp%plane%u(1)*outp%plane%v(3)
        outp%plane%n(3) = outp%plane%u(1)*outp%plane%v(2) - outp%plane%u(2)*outp%plane%v(1)

      case(2)

        call loct_parse_block_float(blk, 0, 0, outp%line%origin(1))
        call loct_parse_block_float(blk, 0, 1, outp%line%origin(2))
        call loct_parse_block_float(blk, 1, 0, outp%line%u(1))
        call loct_parse_block_float(blk, 1, 1, outp%line%u(2))
        call loct_parse_block_float(blk, 2, 0, outp%line%spacing)
        call loct_parse_block_int(blk, 3, 0, outp%line%nu)
        call loct_parse_block_int(blk, 3, 1, outp%line%mu)

        outp%line%n(1) = -outp%line%u(2)
        outp%line%n(2) =  outp%line%u(1)

      case(1)

        call loct_parse_block_float(blk, 0, 0, outp%line%origin(1))

      end select
    end if

    !%Variable OutputEvery
    !%Type integer
    !%Default 1000
    !%Section Time Dependent::TD Output
    !%Description
    !% The output is saved when the iteration number is a multiple of the
    !% <tt>OutputEvery</tt> variable.
    !%End
    call loct_parse_int(datasets_check('OutputEvery'), 1000, outp%iter)

    call loct_parse_logical(datasets_check('OutputDuringSCF'), .false., outp%duringscf)

    if(outp%what.ne.0) call io_function_read_how(sb, outp%how)

  end subroutine h_sys_output_init


  ! ---------------------------------------------------------
  subroutine h_sys_output_all(outp, gr, geo, st, hm, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    type(h_sys_output_t), intent(in)    :: outp
    character(len=*),     intent(in)    :: dir
    
    call h_sys_output_states(st, gr, geo, dir, outp)
    call h_sys_output_hamiltonian(hm, gr%mesh, gr%sb, dir, outp, geo)
    call h_sys_output_localization_funct(st, hm, gr, dir, outp, geo)
    call h_sys_output_current_flow(gr, st, dir, outp, geo)

    if(iand(outp%what, output_geometry).ne.0) then
      call atom_write_xyz(dir, "geometry", geo, gr%mesh%sb%dim)
      if(simul_box_is_periodic(gr%sb)) &
        call periodic_write_crystal(gr%sb, geo, dir)
    end if

#if defined(HAVE_ETSF_IO)
    if (outp%how == output_etsf) then
      call h_sys_output_etsf(st, gr, geo, dir, outp)
    end if
#endif
    
  end subroutine h_sys_output_all

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
    
    mpi_grp = gr%mesh%mpi_grp
    if(st%parallel_in_states) mpi_grp = st%mpi_grp
    if(st%d%kpt%parallel) mpi_grp = st%d%kpt%mpi_grp

    ! if SPIN_POLARIZED, the ELF contains one extra channel: the total ELF
    imax = st%d%nspin
    if(st%d%ispin == SPIN_POLARIZED) imax = 3

    SAFE_ALLOCATE(f_loc(1:gr%mesh%np, 1:imax))

    ! First the ELF in real space
    if(iand(outp%what, output_elf).ne.0 .or. iand(outp%what, output_elf_basins).ne.0) then
      ASSERT(gr%mesh%sb%dim.ne.1)

      call elf_calc(st, gr, f_loc)
      
      ! output ELF in real space
      if(iand(outp%what, output_elf).ne.0) then
        write(fname, '(a)') 'elf_rs'
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
          f_loc(:,imax), M_ONE, ierr, is_tmp = .false., geo = geo)

        if(st%d%ispin.ne.UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,a,i1)') 'elf_rs', '-', is
            call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
              f_loc(:, is), M_ONE, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
          end do
        end if
      end if

      if(iand(outp%what, output_elf_basins).ne.0) &
        call out_basins(f_loc(:,1), "elf_rs_basins")
    end if

    ! Second, ELF in Fourier space.
    if(iand(outp%what, output_elf_fs).ne.0) then
      call elf_calc_fs(st, gr, f_loc)
      do is = 1, st%d%nspin
        write(fname, '(a,a,i1)') 'elf_fs', '-', is
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
          f_loc(:,is), M_ONE, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
      end do
    end if

    ! Now Bader analysis
    if(iand(outp%what, output_bader).ne.0) then
      do is = 1, st%d%nspin
        call dderivatives_lapl(gr%der, st%rho(:,is), f_loc(:,is))
        write(fname, '(a,a,i1)') 'bader', '-', is
        call doutput_function(outp%how, dir, trim(fname), gr%mesh, gr%sb, &
          f_loc(:,is), M_ONE, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)

        write(fname, '(a,a,i1)') 'bader_basins', '-', is
        call out_basins(f_loc(:,1), fname)
      end do
    end if

    ! Now the pressure
    if(iand(outp%what, output_el_pressure).ne.0) then
      call h_sys_calc_electronic_pressure(st, hm, gr, f_loc(:,1))
      call doutput_function(outp%how, dir, "el_pressure", gr%mesh, gr%sb, &
        f_loc(:,1), M_ONE, ierr, is_tmp = .false., geo = geo, grp = mpi_grp)
    end if

    SAFE_DEALLOCATE_A(f_loc)

  contains
    ! ---------------------------------------------------------
    subroutine out_basins(ff, filename)
      FLOAT,            intent(in)    :: ff(:)
      character(len=*), intent(in)    :: filename

      character(len=256) :: fname
      type(basins_t)     :: basins
      integer            :: iunit

      call basins_init(basins, gr%mesh)
      call basins_analyze(basins, gr%mesh, st%d%nspin, ff(:), st%rho, CNST(0.01))

      call doutput_function(outp%how, dir, trim(filename), gr%mesh, gr%sb, &
        real(basins%map, REAL_PRECISION), M_ONE, ierr, is_tmp = .false., geo = geo)
        
      write(fname,'(4a)') trim(dir), '/', trim(filename), '.info'
      iunit = io_open(file=trim(fname), action = 'write')
      call basins_write(basins, gr%mesh, iunit)
      call io_close(iunit)

      call basins_end(basins)

    end subroutine out_basins

  end subroutine h_sys_output_localization_funct

  
  subroutine h_sys_calc_electronic_pressure(st, hm, gr, pressure)
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(grid_t),           intent(inout) :: gr
    FLOAT,                  intent(out)   :: pressure(:)

    FLOAT, allocatable :: rho(:,:), lrho(:), tau(:,:)
    FLOAT   :: p_tf, dens
    integer :: is, ii

    SAFE_ALLOCATE( rho(1:gr%mesh%np_part, 1:st%d%nspin))
    SAFE_ALLOCATE(lrho(1:gr%mesh%np))
    SAFE_ALLOCATE( tau(1:gr%mesh%np, 1:st%d%nspin))

    rho = M_ZERO
    call states_calc_dens(st, gr%mesh%np, rho)
    call states_calc_tau_jp_gn(gr, st, tau=tau)

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

      ! add xc pressure
      pressure(ii) = pressure(ii) + (dens*hm%vxc(ii,1) - hm%ex - hm%ec)

      pressure(ii) = pressure(ii)/p_tf
      pressure(ii) = M_HALF*(M_ONE + pressure(ii)/sqrt(M_ONE + pressure(ii)**2))
    end do

  end subroutine h_sys_calc_electronic_pressure


#include "output_etsf.F90"

#include "output_states.F90"
#include "output_h.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_linear_response.F90"

#include "undef.F90"
#include "real.F90"
#include "output_linear_response.F90"

end module h_sys_output_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
