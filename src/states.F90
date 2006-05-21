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

module states_m
  use global_m
  use varinfo_m
  use string_m
  use messages_m
  use datasets_m
  use profiling_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use math_m
  use mesh_m
  use grid_m
  use simul_box_m
  use functions_m
  use mesh_function_m
  use cube_function_m
  use output_m
  use geometry_m
  use crystal_m
  use mpi_m
  use multicomm_m

  implicit none

  private
  public ::                         &
    states_t,                       &
    states_dim_t,                   &
    states_pair_t,                  &
    states_init,                    &
    states_read_user_def_orbitals,  &
    states_densities_init,          &
    states_allocate_wfns,           &
    states_deallocate_wfns,         &
    states_null,                    &
    states_end,                     &
    states_copy,                    &
    states_generate_random,         &
    states_fermi,                   &
    states_calculate_multipoles,    &
    states_eigenvalues_sum,         &
    states_write_eigenvalues,       &
    states_write_bands,             &
    states_spin_channel,            &
    states_calc_projection,         &
    states_magnetization_dens,      &
    states_magnetic_moment,         &
    states_local_magnetic_moments,  &
    states_calc_physical_current,   &
    states_calc_dens,               &
    states_output,                  &
    states_calc_elf,                &
    states_calc_elf_fs,             &
    kpoints_write_info


  public ::                         &
    dstates_gram_schmidt,           &
    zstates_gram_schmidt,           &
    dstates_dotp,                   &
    zstates_dotp,                   &
    dstates_nrm2,                   &
    zstates_nrm2,                   &
    dstates_normalize_orbital,      &
    zstates_normalize_orbital,      &
    dstates_residue,                &
    zstates_residue,                &
    dstates_mpdotp,                 &
    zstates_mpdotp,                 &
    dstates_angular_momentum,       &
    zstates_angular_momentum,       &
    states_distribute_nodes

  type states_t

    type(states_dim_t), pointer :: d
    integer :: nst                  ! Number of states in each irreducible subspace

    ! pointers to the wavefunctions 
    FLOAT, pointer :: dpsi(:,:,:,:) ! dpsi(sys%NP_PART, st%d%dim, st%nst, st%d%nik)
    CMPLX, pointer :: zpsi(:,:,:,:) ! zpsi(sys%NP_PART, st%d%dim, st%nst, st%d%nik)

    ! used for the user defined wavefunctions (they are stored as formula strings)
    character(len=1024), pointer :: user_def_states(:,:,:) ! (st%d%dim, st%nst, st%d%nik)

    ! the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)      ! rho(gr%m%np_part, st%d%nspin)
    FLOAT, pointer :: j(:,:,:)      !   j(gr%m%np_part, gr%sb%dim, st%d%nspin)

    FLOAT, pointer :: rho_core(:)   ! core charge for nl core corrections

    FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
    logical        :: fixed_occ     ! should the occupation numbers be fixed?
    FLOAT, pointer :: occ(:,:)      ! the occupation numbers
    FLOAT, pointer :: mag(:, :, :)

    FLOAT :: qtot                   ! (-) The total charge in the system (used in Fermi)
    FLOAT :: val_charge             ! valence charge

    FLOAT :: el_temp                ! electronic temperature for the Fermi function
    FLOAT :: ef                     ! the fermi energy

    ! This is stuff needed for the parallelization in states
    logical :: parallel_in_states   ! am I parallel in states?
    type(mpi_grp_t) :: mpi_grp   ! the MPI group related to the parallelization in states

    integer :: st_start, st_end     ! needed for some parallel parts
    integer, pointer :: node(:)     ! To which node belongs each state.
  end type states_t

  type states_dim_t
    integer :: wfs_type            ! real (M_REAL) or complex (M_CMPLX) wavefunctions
    integer :: dim                  ! Dimension of the state (one or two for spinors)
    integer :: nik                  ! Number of irreducible subspaces
    integer :: nik_axis(MAX_DIM)    ! Number of kpoints per axis
    integer :: ispin                ! spin mode (unpolarized, spin polarized, spinors)
    integer :: nspin                ! dimension of rho (1, 2 or 4)
    integer :: spin_channels        ! 1 or 2, wether spin is or not considered.
    logical :: cdft                 ! Are we using Current-DFT or not?
    FLOAT, pointer :: kpoints(:,:)  ! obviously the kpoints
    FLOAT, pointer :: kweights(:)   ! weights for the kpoint integrations
  end type states_dim_t

  type states_pair_t
    integer :: i
    integer :: a
    integer :: sigma
  end type states_pair_t

  ! Parameters...
  integer, public, parameter ::     &
    UNPOLARIZED    = 1,             &
    SPIN_POLARIZED = 2,             &
    SPINORS        = 3

  integer, public, parameter ::     &
    M_REAL  = 1,             &
    M_CMPLX = 2

  interface assignment (=)
    module procedure states_copy
  end interface

  interface dstates_gram_schmidt
    module procedure dstates_gram_schmidt1, dstates_gram_schmidt2
  end interface

  interface zstates_gram_schmidt
    module procedure zstates_gram_schmidt1, zstates_gram_schmidt2
  end interface

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_t), intent(out) :: st

    nullify(st%dpsi, st%zpsi, st%rho, st%j, st%rho_core, st%eigenval, st%occ, st%mag)
    nullify(st%d)
    ALLOCATE(st%d, 1)
    nullify(st%d%kpoints, st%d%kweights)

    ! By default, calculations use real wave-functions
    st%d%wfs_type = M_REAL

  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr

    FLOAT :: excess_charge, r
    integer :: nempty, i, j
    integer(POINTER_SIZE) :: blk

    call push_sub('states.states_init')

    call states_null(st)

    !%Variable SpinComponents
    !%Type integer
    !%Default unpolarized
    !%Section States
    !%Description
    !% The calculations may be done in three different ways: spin-restricted (TD)DFT (i.e., doubly
    !% occupied "closed shells"), spin-unsrestricted or "spin-polarized" (TD)DFT (i.e. we have two
    !% electronic systes, one with spin up and one with spin down), or making use of two-component
    !% spinors.
    !%Option unpolarized 1
    !% Spin-restricted calculations.
    !%Option polarized 2
    !%Option spin_polarized 2
    !% Spin unrestricted, also know as spin-DFT, SDFT. This mode will double the number of wave
    !% functions will double the number of wave-functions necessary for a spin-unpolarised
    !% calculation.
    !%Option non_collinear 3
    !%Option spinors 3
    !% The spin-orbitals are two-component spinors. This effectively allows the spin-density to
    !% arrange non-collinearly - i.e. the magnetization vector is allowed to take different
    !% directions in different points.
    !%End
    call loct_parse_int(check_inp('SpinComponents'), UNPOLARIZED, st%d%ispin)
    if(.not.varinfo_valid_option('SpinComponents', st%d%ispin)) call input_error('SpinComponents')
    call messages_print_var_option(stdout, 'SpinComponents', st%d%ispin)
    ! Use of Spinors requires complex wave-functions
    if (st%d%ispin == SPINORS) st%d%wfs_type = M_CMPLX


    !%Variable ExcessCharge
    !%Type float
    !%Default 0.0
    !%Section States
    !%Description
    !% The net charge of the system. A negative value means that we are adding 
    !% electrons, while a positive value means we are taking electrons
    !% from the system.
    !%End
    call loct_parse_float(check_inp('ExcessCharge'), M_ZERO, excess_charge)


    !%Variable ExtraStates
    !%Type integer
    !%Default 0
    !%Section States
    !%Description
    !% The number of states is in principle calculated considering the minimum
    !% numbers of states necessary to hold the electrons present in the system.
    !% The number of electrons is
    !% in turn calculated considering the nature of the species supplied in the
    !% <tt>Species</tt> block, and the value of the <tt>ExcessCharge</tt> variable.
    !% However, one may command <tt>octopus</tt> to put more states, which is necessary if one wants to
    !% use fractional occupational numbers, either fixed from the origin through
    !% the <tt>Occupations</tt> block or by prescribing
    !% an electronic temperature with <tt>ElectronicTemperature</tt>.
    !%
    !% Note that this number is unrelated to <tt>CalculationMode == unocc</tt>.
    !%End
    call loct_parse_int(check_inp('ExtraStates'), 0, nempty)
    if (nempty < 0) then
      write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
      message(2) = '(0 <= ExtraStates)'
      call write_fatal(2)
    end if

    call geometry_val_charge(gr%geo, st%val_charge)
    st%qtot = -(st%val_charge + excess_charge)

    select case(st%d%ispin)
    case(UNPOLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nspin = 1
      st%d%spin_channels = 1
    case(SPIN_POLARIZED)
      st%d%dim = 1
      st%nst = int(st%qtot/2)
      if(st%nst*2 < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nik = st%d%nik*2
      st%d%nspin = 2
      st%d%spin_channels = 2
    case(SPINORS)
      st%d%dim = 2
      st%nst = int(st%qtot)
      if(st%nst < st%qtot) st%nst = st%nst + 1
      st%nst = st%nst + nempty
      st%d%nspin = 4
      st%d%spin_channels = 2
    end select

    ! current
    call loct_parse_logical(check_inp('CurrentDFT'), .false., st%d%cdft)
    if (st%d%cdft) then
      ! Use of CDFT requires complex wave-functions
      st%d%wfs_type = M_CMPLX

      if(st%d%ispin == SPINORS) then
        message(1) = "Sorry, Current DFT not working yet for spinors"
        call write_fatal(1)
      end if
      message(1) = "Info: Using Current DFT"
      call write_info(1)
    end if

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, gr%geo)

    ! Periodic systems require complex wave-functions
    if(simul_box_is_periodic(gr%sb)) st%d%wfs_type = M_CMPLX


    ! we now allocate some arrays
    ALLOCATE(st%occ     (st%nst, st%d%nik), st%nst*st%d%nik)
    ALLOCATE(st%eigenval(st%nst, st%d%nik), st%nst*st%d%nik)
    ! allocate space for formula strings that define user defined states
    ALLOCATE(st%user_def_states(st%d%dim, st%nst, st%d%nik), st%d%dim*st%nst*st%d%nik)
    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%mag(st%nst, st%d%nik, 2), st%nst*st%d%nik*2)
    end if

    ! initially we mark all 'formulas' as undefined
    st%user_def_states(:, :, :) = 'undefined'

    !%Variable Occupations
    !%Type block
    !%Section States
    !%Description
    !% The occupation numbers of the orbitals can be fixed through the use of this
    !% variable. For example:
    !%
    !% <tt>%Occupations
    !% <br>&nbsp;&nbsp;2.0 | 2.0 | 2.0 | 2.0 | 2.0
    !% <br>%</tt>
    !%
    !% would fix the occupations of the five states to <i>2.0</i>. There must be
    !% as many columns as states in the calculation. If <tt>SpinComponents == polarized</tt>
    !% this block should contain two lines, one for each spin channel.
    !% This variable is very useful when dealing with highly symmetric small systems
    !% (like an open shell atom), for it allows us to fix the occupation numbers
    !% of degenerate states in order to help <tt>octopus</tt> to converge. This is to
    !% be used in conjuction with <tt>ExtraStates</tt>. For example, to calculate the
    !% carbon atom, one would do:
    !%
    !% <tt>ExtraStates = 2
    !% <br>%Occupations
    !% <br>&nbsp;&nbsp;2 | 2/3 | 2/3 | 2/3
    !% <br>%</tt>
    !%End
    occ_fix: if(loct_parse_block(check_inp('Occupations'), blk)==0) then
      ! read in occupations
      st%fixed_occ = .true.

      do i = 1, st%d%nik
        do j = 1, st%nst
          call loct_parse_block_float(blk, i-1, j-1, st%occ(j, i))
        end do
      end do
      call loct_parse_block_end(blk)
    else
      st%fixed_occ = .false.

      ! first guess for occupation...paramagnetic configuration
      if(st%d%ispin == UNPOLARIZED) then
        r = M_TWO
      else
        r = M_ONE
      end if
      st%occ  = M_ZERO
      st%qtot = M_ZERO

      do j = 1, st%nst
        do i = 1, st%d%nik
          st%occ(j, i) = min(r, -(st%val_charge + excess_charge) - st%qtot)
          st%qtot = st%qtot + st%occ(j, i)

        end do
      end do

      !%Variable ElectronicTemperature
      !%Type float
      !%Default 0.0
      !%Section States
      !%Description
      !% If <tt>Occupations</tt> is not set, <tt>ElectronicTemperature</tt> is the
      !% temperature in the Fermi-Dirac function used to distribute the electrons
      !% among the existing states.
      !%End
      call loct_parse_float(check_inp('ElectronicTemperature'), M_ZERO, st%el_temp)
      st%el_temp = st%el_temp * units_inp%energy%factor
    end if occ_fix

    st%st_start = 1; st%st_end = st%nst
    ALLOCATE(st%node(st%nst), st%nst)
    st%node(:) = 0

    call mpi_grp_init(st%mpi_grp, -1)
    st%parallel_in_states = .false.

    nullify(st%dpsi, st%zpsi)

    call pop_sub()
  end subroutine states_init

  ! ---------------------------------------------------------
  ! Allocates the KS wavefunctions defined within an states_t
  ! structure.
  subroutine states_allocate_wfns(st, m, wfs_type)
    type(states_t), intent(inout) :: st
    type(mesh_t),    intent(in)    :: m
    integer, optional, intent(in) :: wfs_type

    integer :: n

    call push_sub('states.states_allocate_wfns')

    if (present(wfs_type)) then
      ASSERT(wfs_type == M_REAL .or. wfs_type == M_CMPLX)
      st%d%wfs_type = wfs_type
    end if

    n = m%np_part * st%d%dim * (st%st_end-st%st_start+1) * st%d%nik
    if (st%d%wfs_type == M_REAL) then
      ALLOCATE(st%dpsi(m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), n)
      st%dpsi = M_ZERO
    else
      ALLOCATE(st%zpsi(m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), n)
      st%zpsi = M_Z0
    end if

    call pop_sub()
  end subroutine states_allocate_wfns

  ! ---------------------------------------------------------
  ! Deallocates the KS wavefunctions defined within an states_t
  ! structure.
  subroutine states_deallocate_wfns(st)
    type(states_t), intent(inout) :: st

    call push_sub('states.states_deallocate_wfns')

    if (st%d%wfs_type == M_REAL) then
      deallocate(st%dpsi); nullify(st%dpsi)
    else
      deallocate(st%zpsi); nullify(st%zpsi)
    end if
    
    call pop_sub()
  end subroutine states_deallocate_wfns

  ! ---------------------------------------------------------
  ! the routine reads formulas for user defined wavefunctions 
  ! from the input file and fills the respective orbitals
  subroutine states_read_user_def_orbitals(mesh, st)
    type(mesh_t),      intent(in) :: mesh
    type(states_t), intent(inout) :: st    

    integer(POINTER_SIZE) :: blk
    integer :: ip, id, is, ik, nstates
    integer :: ib, idim, inst, inik
    FLOAT   :: x(MAX_DIM), r, psi_re, psi_im

    call push_sub('td.read_user_def_states')

    !%Variable UserDefinedStates
    !%Type block
    !%Section States
    !%Description
    !% Instead of using the ground state as initial state for
    !% time propagations it might be interesting in some cases 
    !% to specify alternative states. Similar to user defined
    !% potentials this block allows to specify formulas for
    !% the orbitals at t=0.
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedStates
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
    !% <br>%</tt>
    !%
    !% The first column specifies the component of the spinor, 
    !% the second column the number of the state and the third 
    !% contains kpoint and spin quantum numbers. Finally column 
    !% four contains a formula for the correspondig orbital.
    !% 
    !% Octopus reads first the ground state orbitals from
    !% the restart_gs directory. Only the states that are
    !% specified in the above block will be overwritten with
    !% the given analytical expression for the orbital.
    !%
    !%End
    if(loct_parse_block(check_inp('UserDefinedStates'), blk) == 0) then

      ! find out how many lines (i.e. states) the block has
      nstates = loct_parse_block_n(blk)

      ! read all lines
      do ib = 1, nstates
        call loct_parse_block_int(blk, ib-1, 0, idim)
        call loct_parse_block_int(blk, ib-1, 1, inst)
        call loct_parse_block_int(blk, ib-1, 2, inik)

        ! read formula strings and convert to C strings
        do id = 1, st%d%dim
          do is = 1, st%nst
            do ik = 1, st%d%nik   

              ! does the block entry match and is this node responsible?
              if(.not.(id.eq.idim .and. is.eq.inst .and. ik.eq.inik    &
                .and. st%st_start.le.is .and. st%st_end.ge.is) ) cycle

              ! parse formula string
              call loct_parse_block_string(                            &
                blk, ib-1, 3, st%user_def_states(id, is, ik))
              ! convert to C string
              call conv_to_C_string(st%user_def_states(id, is, ik))

              ! fill states with user defined formulas
              do ip = 1, mesh%np
                x = mesh%x(ip, :)
                r = sqrt(sum(x(:)**2))

                ! parse user defined expressions
                call loct_parse_expression(psi_re, psi_im,             &
                  x(1), x(2), x(3), r, st%user_def_states(id, is, ik))
                ! fill state
                st%zpsi(ip, id, is, ik) = psi_re + M_zI*psi_im
              end do

              ! normalize orbital
              call zstates_normalize_orbital(mesh, st%d%dim, st%zpsi(:,:, is, ik))

            end do
          end do
        end do

      end do
      call loct_parse_block_end(blk)

    else
      message(1) = '"UserDefinesStates" has to be specified as block.'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine states_read_user_def_orbitals


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr

    ! allocate arrays for charge and current densities
    ALLOCATE(st%rho(NP_PART, st%d%nspin), NP_PART*st%d%nspin)
    ALLOCATE(st%j(NP_PART, NDIM, st%d%nspin), NP_PART*NDIM*st%d%nspin)
    st%rho = M_ZERO
    st%j   = M_ZERO
    if(gr%geo%nlcc) ALLOCATE(st%rho_core(gr%m%np), gr%m%np)

  end subroutine states_densities_init

  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(in)  :: stin
    type(states_t), intent(out) :: stout

    integer :: i

    call states_null(stout)

    stout%d%wfs_type = stin%d%wfs_type
    stout%d%dim = stin%d%dim
    stout%d%nik = stin%d%nik
    stout%d%nik_axis(:) = stout%d%nik_axis(:)
    stout%d%ispin = stin%d%ispin
    stout%d%nspin = stin%d%nspin
    stout%d%spin_channels = stin%d%spin_channels
    stout%d%cdft = stin%d%cdft
    stout%nst           = stin%nst
    stout%qtot = stin%qtot
    stout%el_temp = stin%el_temp
    stout%ef = stin%ef
    stout%st_start = stin%st_start
    stout%st_end = stin%st_end
    if(associated(stin%dpsi)) then
      i = size(stin%dpsi, 1)*stin%d%dim*(stin%st_end-stin%st_start+1)*stin%d%nik
      ALLOCATE(stout%dpsi(size(stin%dpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik), i)
      stout%dpsi = stin%dpsi
    end if
    if(associated(stin%zpsi)) then
      i = size(stin%zpsi, 1)*stin%d%dim*(stin%st_end-stin%st_start+1)*stin%d%nik
      ALLOCATE(stout%zpsi(size(stin%zpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik), i)
      stout%zpsi = stin%zpsi
    end if
    if(associated(stin%rho)) then
      i = size(stin%rho, 1)*size(stin%rho, 2)
      ALLOCATE(stout%rho(size(stin%rho, 1), size(stin%rho, 2)), i)
      stout%rho = stin%rho
    end if
    if(associated(stin%j)) then
      i = size(stin%j, 1)*size(stin%j, 2)*size(stin%j, 3)
      ALLOCATE(stout%j(size(stin%j, 1), size(stin%j, 2), size(stin%j, 3)), i)
      stout%j = stin%j
    end if
    if(associated(stin%rho_core)) then
      i = size(stin%rho_core, 1)
      ALLOCATE(stout%rho_core(size(stin%rho_core, 1)), i)
      stout%rho_core = stin%rho_core
    end if
    if(associated(stin%eigenval)) then
      i = (stin%st_end-stin%st_start)*stin%d%nik
      ALLOCATE(stout%eigenval(stin%st_start:stin%st_end, stin%d%nik), i)
      stout%eigenval = stin%eigenval
    end if
    stout%fixed_occ = stin%fixed_occ
    if(associated(stin%occ)) then
      i = size(stin%occ, 1)*size(stin%occ, 2)
      ALLOCATE(stout%occ(size(stin%occ, 1), size(stin%occ, 2)), i)
      stout%occ = stin%occ
    end if
    if(associated(stin%mag)) then
      i = size(stin%mag, 1)*size(stin%mag, 2)*size(stin%mag, 3)
      ALLOCATE(stout%mag(size(stin%mag, 1), size(stin%mag, 2), size(stin%mag, 3)), i)
      stout%mag = stin%mag
    end if
    if(associated(stin%d%kpoints)) then
      i = size(stin%d%kpoints, 1)*size(stin%d%kpoints, 2)
      ALLOCATE(stout%d%kpoints(size(stin%d%kpoints, 1), size(stin%d%kpoints, 2)), i)
      stout%d%kpoints = stin%d%kpoints
    end if
    if(associated(stin%d%kweights)) then
      i = size(stin%d%kpoints, 1)
      ALLOCATE(stout%d%kweights(size(stin%d%kpoints, 1)), i)
      stout%d%kweights = stin%d%kweights
    end if
    if(associated(stin%node)) then
      i = size(stin%node)
      ALLOCATE(stout%node(size(stin%node)), i)
      stout%node = stin%node
    end if
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_t), intent(inout) :: st

    call push_sub('states.states_end')

    if(associated(st%rho)) then
      deallocate(st%rho, st%occ, st%eigenval, st%node)
      nullify   (st%rho, st%occ, st%eigenval, st%node)
    end if

    if(associated(st%j)) then
      deallocate(st%j)
      nullify(st%j)
    end if

    if(associated(st%rho_core)) then
      deallocate(st%rho_core)
      nullify(st%rho_core)
    end if

    if(st%d%ispin==SPINORS .and. associated(st%mag)) then
      deallocate(st%mag); nullify(st%mag)
    end if

    if(associated(st%dpsi)) then
      deallocate(st%dpsi); nullify(st%dpsi)
    end if

    if(associated(st%zpsi)) then
      deallocate(st%zpsi); nullify(st%zpsi)
    end if

    if(associated(st%d%kpoints)) then
      deallocate(st%d%kpoints); nullify(st%d%kpoints)
    end if

    if(associated(st%d%kweights)) then
      deallocate(st%d%kweights); nullify(st%d%kweights)
    end if

    if(associated(st%user_def_states)) then
      deallocate(st%user_def_states); nullify(st%user_def_states)
    end if

    call pop_sub()
  end subroutine states_end

  ! ---------------------------------------------------------
  ! Calculates the new density out the wavefunctions and
  ! occupations...
  subroutine states_calc_dens(st, np, rho)
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: np
    FLOAT,          intent(out) :: rho(:,:)

    integer :: i, ik, p, sp
    CMPLX   :: c

#ifdef HAVE_MPI
    FLOAT,  allocatable :: reduce_rho(:)
#endif

    call push_sub('states.states_calc_dens')

    if(st%d%ispin == SPIN_POLARIZED) then
      sp = 2
    else
      sp = 1
    end if

    rho = M_ZERO
    do ik = 1, st%d%nik, sp
      do p  = st%st_start, st%st_end
        do i = 1, np

          if (st%d%wfs_type == M_REAL) then
            rho(i, 1) = rho(i, 1) + st%d%kweights(ik)  *st%occ(p, ik)*abs(st%dpsi(i, 1, p, ik))**2
          else
            rho(i, 1) = rho(i, 1) + st%d%kweights(ik)  *st%occ(p, ik)*abs(st%zpsi(i, 1, p, ik))**2
          end if
          select case(st%d%ispin)

          case(SPIN_POLARIZED)
            if (st%d%wfs_type == M_REAL) then
              rho(i, 2) = rho(i, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1)*abs(st%dpsi(i, 1, p, ik+1))**2
            else
              rho(i, 2) = rho(i, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1)*abs(st%zpsi(i, 1, p, ik+1))**2
            end if
          case(SPINORS) ! in this case wave-functions are always complex
            rho(i, 2) = rho(i, 2) + st%d%kweights(ik)  *st%occ(p, ik)*abs(st%zpsi(i, 2, p, ik))**2

            c = st%d%kweights(ik)*st%occ(p, ik) * st%zpsi(i, 1, p, ik) * conjg(st%zpsi(i, 2, p, ik))
            rho(i, 3) = rho(i, 3) + real(c, PRECISION)
            rho(i, 4) = rho(i, 4) + aimag(c)
          end select

        end do
      end do
    end do

#ifdef HAVE_MPI
    ! reduce density
    if(st%parallel_in_states) then
      ALLOCATE(reduce_rho(1:np), np)
      do i = 1, st%d%nspin
        call MPI_Allreduce(rho(1, i), reduce_rho(1), np, &
           MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        rho(1:np, i) = reduce_rho(1:np)
      end do
      deallocate(reduce_rho)
    end if
#endif

    call pop_sub()
  end subroutine states_calc_dens

  ! ---------------------------------------------------------
  ! generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, m, ist_start_, ist_end_)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: m
    integer, optional, intent(in) :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end

    call push_sub('states.states_generate_random')

    ist_start = 1
    if(present(ist_start_)) ist_start = ist_start_
    ist_end = st%nst
    if(present(ist_end_)) ist_end = ist_end_

    do ik = 1, st%d%nik
      do ist = ist_start, ist_end
        do id = 1, st%d%dim
          if (st%d%wfs_type == M_REAL) then
            call dmf_random(m, st%dpsi(:, id, ist, ik))
          else
            call zmf_random(m, st%zpsi(:, id, ist, ik))
          end if
        end do
        st%eigenval(ist, ik) = M_ZERO
      end do
      ! Do not orthonormalize if we have parallelization in states.
      if(.not.st%parallel_in_states) then
        if (st%d%wfs_type == M_REAL) then
          call dstates_gram_schmidt(ist_end, m, st%d%dim, st%dpsi(:,:,1:ist_end,ik), start = ist_start)
        else
          call zstates_gram_schmidt(ist_end, m, st%d%dim, st%zpsi(:,:,1:ist_end,ik), start = ist_start)
        end if
      end if
    end do

    call pop_sub()
  end subroutine states_generate_random


  ! ---------------------------------------------------------
  subroutine states_fermi(st, m)
    type(states_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: m

    ! Local variables
    integer :: ie, ik, iter
    integer, parameter :: nitmax = 200
    FLOAT :: drange, t, emin, emax, sumq
    FLOAT, parameter :: tol = CNST(1.0e-10)
    logical :: conv

    call push_sub('states.fermi')

    if(st%fixed_occ) then ! nothing to do
      ! Calculate magnetizations...
      if(st%d%ispin == SPINORS) then
        do ik = 1, st%d%nik
          do ie = 1, st%nst
            if (st%d%wfs_type == M_REAL) then
              st%mag(ie, ik, 1) = dmf_nrm2(m, st%dpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
              st%mag(ie, ik, 2) = dmf_nrm2(m, st%dpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
            else
              st%mag(ie, ik, 1) = zmf_nrm2(m, st%zpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
              st%mag(ie, ik, 2) = zmf_nrm2(m, st%zpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
            end if
          end do
        end do
      end if
      call pop_sub()
      return
    end if

    ! Initializations
    emin = minval(st%eigenval)
    emax = maxval(st%eigenval)

    if(st%d%ispin == SPINORS) then
      sumq = real(st%nst, PRECISION)
    else
      sumq = M_TWO*st%nst
    end if

    t = max(st%el_temp, CNST(1.0e-6))
    st%ef = emax

    conv = .true.
    if (abs(sumq - st%qtot) > tol) conv = .false.
    if (conv) then ! all orbitals are full; nothing to be done
      st%occ = M_TWO/st%d%spin_channels!st%d%nspin
      ! Calculate magnetizations...
      if(st%d%ispin == SPINORS) then
        do ik = 1, st%d%nik
          do ie = 1, st%nst
            if (st%d%wfs_type == M_REAL) then
              st%mag(ie, ik, 1) = dmf_nrm2(m, st%dpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
              st%mag(ie, ik, 2) = dmf_nrm2(m, st%dpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
            else
              st%mag(ie, ik, 1) = zmf_nrm2(m, st%zpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
              st%mag(ie, ik, 2) = zmf_nrm2(m, st%zpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
            end if
          end do
        end do
      end if
      call pop_sub()
      return
    end if

    if (sumq < st%qtot) then ! not enough states
      message(1) = 'Fermi: Not enough states'
      write(message(2),'(6x,a,f12.6,a,f12.6)')'(total charge = ', st%qtot, &
        ' max charge = ', sumq
      call write_fatal(2)
    end if

    drange = t*sqrt(-log(tol*CNST(.01)))

    emin = emin - drange
    emax = emax + drange

    do iter = 1, nitmax
      st%ef = M_HALF*(emin + emax)
      sumq  = M_ZERO

      do ik = 1, st%d%nik
        do ie =1, st%nst
          sumq = sumq + st%d%kweights(ik)/st%d%spin_channels * & !st%d%nspin * &
            stepf((st%eigenval(ie, ik) - st%ef)/t)
        end do
      end do

      conv = .true.
      if(abs(sumq - st%qtot) > tol) conv = .false.
      if(conv) exit

      if(sumq <= st%qtot ) emin = st%ef
      if(sumq >= st%qtot ) emax = st%ef
    end do

    if(iter == nitmax) then
      message(1) = 'Fermi: did not converge'
      call write_fatal(1)
    end if

    do ik = 1, st%d%nik
      do ie = 1, st%nst
        st%occ(ie, ik) = stepf((st%eigenval(ie, ik) - st%ef)/t)/st%d%spin_channels!st%d%nspin
      end do
    end do

    ! Calculate magnetizations...
    if(st%d%ispin == SPINORS) then
      do ik = 1, st%d%nik
        do ie = 1, st%nst
          if (st%d%wfs_type == M_REAL) then
            st%mag(ie, ik, 1) = dmf_nrm2(m, st%dpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = dmf_nrm2(m, st%dpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
          else
            st%mag(ie, ik, 1) = zmf_nrm2(m, st%zpsi(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = zmf_nrm2(m, st%zpsi(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
          end if
        end do
      end do
    end if

    call pop_sub()
  end subroutine states_fermi


  ! -----------------------------------------------------------------------------
  ! This routine calculates the multipoles of the *electronic* charge
  ! distribution, defined in the following way:
  ! multipole(1, is) is the electronic charge (defined to be positive; integral
  !   of the electronic density.
  ! multipole(2:4, is) contains the dipole: integral of the electronic density
  !   times x, y or z.
  ! multipole(5:9, is) contains the quadrupole, defined in the usual way using
  !   the spherical harmonics: multipole(5, is) = Integral [ rho * Y_{2,-2} ],
  !   multipole(6, is) = Integral [ rho * Y_{2, -1} ].
  ! And so on.
  ! -----------------------------------------------------------------------------
  subroutine states_calculate_multipoles(gr, st, lmax, multipole)
    type(grid_t),   intent(in)  :: gr
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: lmax
    FLOAT,          intent(out) :: multipole(:, :) ! multipole((lmax + 1)**2, st%d%nspin)

    integer :: i, is, l, lm, add_lm
    FLOAT :: x(MAX_DIM), r, ylm
    FLOAT, allocatable :: f(:)

    call push_sub('states.states_calculate_multipoles')

    ALLOCATE(f(NP), NP)
    f = M_ZERO

    do is = 1, st%d%nspin

      f(1:NP) = st%rho(1:NP, is)
      multipole(1, is) = dmf_integrate(gr%m, f)

      if(lmax>0) then
        do i = 1, 3
          f(1:NP) = st%rho(1:NP, is)*gr%m%x(1:NP, i)
          multipole(i+1, is) = dmf_integrate(gr%m, f)
        end do
      end if

      if(lmax>1) then
        add_lm = 5
        do l = 2, lmax
          do lm = -l, l
            do i = 1, NP
              call mesh_r(gr%m, i, r, x=x)
              ylm = loct_ylm(x(1), x(2), x(3), l, lm)
              f(i) = st%rho(i, is) * ylm * r**l
            end do
            multipole(add_lm, is) = dmf_integrate(gr%m, f)
            add_lm = add_lm + 1
          end do
        end do
      end if

    end do

    deallocate(f)
    call pop_sub()
  end subroutine states_calculate_multipoles


  ! ---------------------------------------------------------
  ! function to calculate the eigenvalues sum using occupations as weights
  function states_eigenvalues_sum(st, x) result(e)
    type(states_t), intent(in)  :: st
    FLOAT                       :: e
    FLOAT, optional, intent(in) :: x(:, :)

    integer :: ik
#ifdef HAVE_MPI
    FLOAT :: s
#endif

    e = M_ZERO
    do ik = 1, st%d%nik
      if(present(x)) then
        e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
          x(st%st_start:st%st_end, ik))
      else
        e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
          st%eigenval(st%st_start:st%st_end, ik))
      end if
    end do

#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(e, s, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      e = s
    end if
#endif

  end function states_eigenvalues_sum


  ! ---------------------------------------------------------
  subroutine states_write_eigenvalues(iunit, nst, st, sb, error)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    FLOAT,             intent(in), optional :: error(nst, st%d%nik)

    integer ik, j, ns, is
    FLOAT :: o, oplus, ominus
    character(len=80) tmp_str(MAX_DIM), cspin

    call push_sub('states.states_write_eigenvalues')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_out%energy%abbrev) // ']'
    call write_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
      call write_info(1, iunit)
    end if

    if(.not.mpi_grp_is_root(mpi_world)) return

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
          st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
          st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
          st%d%kpoints(3, ik)*units_out%length%factor, ')'
        call write_info(1, iunit)
      end if

      if(present(error)) then
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x,a10)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', 'Error'
      else
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x)')       &
          '#st',' Spin',' Eigenvalue', 'Occupation '
      end if
      call write_info(1, iunit)

      do j = 1, nst
        do is = 0, ns-1
          if(j > st%nst) then
            o = M_ZERO
            if(st%d%ispin == SPINORS) oplus = M_ZERO; ominus = M_ZERO
          else
            o = st%occ(j, ik+is)
            if(st%d%ispin == SPINORS) then
              oplus  = st%mag(j, ik+is, 1)
              ominus = st%mag(j, ik+is, 2)
            end if
          end if

          if(is.eq.0) cspin = 'up'
          if(is.eq.1) cspin = 'dn'
          if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

          write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
          if(simul_box_is_periodic(sb)) then
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,f5.2,a1,f5.2)') &
                (st%eigenval(j, ik)-st%ef)/units_out%energy%factor, oplus, '/', ominus
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                (st%eigenval(j, ik+is))/units_out%energy%factor, o
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik), ')'
            end if
          else
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,f5.2,a1,f5.2)') &
                st%eigenval(j, ik)/units_out%energy%factor, oplus, '/', ominus
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                st%eigenval(j, ik+is)/units_out%energy%factor, o
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(j, ik), ')'
            end if
          end if
          if(present(error)) then
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))//trim(tmp_str(3))
          else
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
          end if
          call write_info(1, iunit)
        end do
      end do
    end do

    call pop_sub()
  end subroutine states_write_eigenvalues


  ! ---------------------------------------------------------
  subroutine states_write_bands(iunit, nst, st, sb)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: i, ik, j, mode, ns
    FLOAT :: factor(MAX_DIM)

    integer, parameter :: &
      GNUPLOT = 1024,     &
      XMGRACE = 2048

    if(.not.mpi_grp_is_root(mpi_world)) return

    !%Variable OutputBandsMode
    !%Type integer
    !%Default 1024
    !%Section Output
    !%Description
    !% Chose if the band file is to be written in gnuplot or xmgrace friendly format
    !%Option gnuplot 1024
    !% gnuplot format
    !%Option xmgrace 2048
    !% xmgrace format
    !%End
    call loct_parse_int(check_inp('OutputBandsMode'), GNUPLOT, mode)
    if(.not.varinfo_valid_option('OutputBandsMode', mode)) call input_error('OutputBandsMode')
    call messages_print_var_option(stdout, "OutputBandsMode", mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do i =1,3
      factor(i) = M_ONE
      if (sb%klat(i,i) /= M_ZERO) factor(i) = sb%klat(i,i)
    end do

    select case(mode)
    case(GNUPLOT)
      ! output bands in gnuplot format
      do j = 1, nst
        do ik = 1, st%d%nik, ns
          write(iunit, '(1x,3(f10.4))', advance='no')  st%d%kpoints(:,ik)/factor(:)
          write(iunit, '(3x,f12.6))',   advance='yes') st%eigenval(j, ik)/units_out%energy%factor
        end do
        write(iunit, '(a)')' '
      end do
    case(XMGRACE)
      ! output bands in xmgrace format, i.e.:
      ! k_x, k_y, k_z, e_1, e_2, ..., e_n
      do ik = 1, st%d%nik, ns
        write(iunit, '(1x,3(f10.4))', advance='no')   st%d%kpoints(:,ik)/factor(:)
        write(iunit, '(3x,20f12.6))', advance='yes') (st%eigenval(j, ik)/units_out%energy%factor, j=1,nst)
      end do
    end select

  end subroutine states_write_bands


  ! -------------------------------------------------------
  integer function states_spin_channel(ispin, ik, dim)
    integer, intent(in) :: ispin, ik, dim

    ASSERT(ispin >= UNPOLARIZED .or. ispin <= SPINORS)
    ASSERT(ik > 0)
    ASSERT(dim==1 .or. dim==2)
    ASSERT(.not.(ispin.ne.3 .and. dim==2))

    select case(ispin)
    case(1); states_spin_channel = 1
    case(2); states_spin_channel = mod(ik+1, 2)+1
    case(3); states_spin_channel = dim
    case default; states_spin_channel = -1
    end select

  end function states_spin_channel


  ! ---------------------------------------------------------
  ! This subroutine calculates:
  ! p(uist, ist, ik) = < phi0(uist, k) | phi(ist, ik) (t) >
  ! ---------------------------------------------------------
  subroutine states_calc_projection(m, st, gs_st, p)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    type(states_t), intent(in)  :: gs_st
    CMPLX,          intent(out) :: p(st%nst, gs_st%nst, st%d%nik)

    integer :: uist, ist, ik

    call push_sub('states.calc_projection')

    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do uist = 1, gs_st%nst
          p(ist, uist, ik) = zstates_dotp(m, st%d%dim, st%zpsi(:, :, ist, ik), gs_st%zpsi(:, :, uist, ik))
        end do
      end do
    end do

    call pop_sub()
  end subroutine states_calc_projection


  ! ---------------------------------------------------------
  subroutine states_magnetization_dens(m, st, rho, md)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    FLOAT,             intent(in)  :: rho(:,:) ! (np, st%d%nspin)
    FLOAT,             intent(out) :: md(:,:)   ! (np, 3)

    call push_sub('states.states_magnetization_dens')

    select case (st%d%ispin)
    case (UNPOLARIZED)
      md = M_ZERO
    case (SPIN_POLARIZED)
      md = M_ZERO
      md(1:m%np, 3) = rho(1:m%np, 1) - rho(1:m%np, 2)
    case (SPINORS)
      md(1:m%np, 1) =  M_TWO*rho(1:m%np, 3)
      md(1:m%np, 2) = -M_TWO*rho(1:m%np, 4)
      md(1:m%np, 3) = rho(1:m%np, 1) - rho(1:m%np, 2)
    end select

    call pop_sub()
  end subroutine states_magnetization_dens


  ! ---------------------------------------------------------
  subroutine states_magnetic_moment(m, st, rho, mm)
    type(mesh_t),   intent(in)  :: m
    type(states_t), intent(in)  :: st
    FLOAT,          intent(in)  :: rho(:,:) ! (m%np_part, st%d%nspin)
    FLOAT,          intent(out) :: mm(3)

    FLOAT, allocatable :: md(:,:)

    call push_sub('states.states_magnetic_moment')

    ALLOCATE(md(m%np, 3), m%np*3)
    call states_magnetization_dens(m, st, rho, md)
    mm(1) = dmf_integrate(m, md(:, 1))
    mm(2) = dmf_integrate(m, md(:, 2))
    mm(3) = dmf_integrate(m, md(:, 3))
    deallocate(md)

    call pop_sub()
  end subroutine states_magnetic_moment


  ! ---------------------------------------------------------
  subroutine states_local_magnetic_moments(m, st, geo, rho, r, lmm)
    type(mesh_t),     intent(in)  :: m
    type(states_t),   intent(in)  :: st
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(in)  :: rho(:,:) ! (m%np_part, st%d%nspin)
    FLOAT,            intent(in)  :: r
    FLOAT,            intent(out) :: lmm(3, geo%natoms)

    integer :: ia, i
    FLOAT :: ri
    FLOAT, allocatable :: md(:, :), aux(:, :)

    call push_sub('states.states_local_magnetic_moments')

    ALLOCATE(md (m%np, 3), m%np*3)
    ALLOCATE(aux(m%np, 3), m%np*3)
    call states_magnetization_dens(m, st, rho, md)
    lmm = M_ZERO
    do ia = 1, geo%natoms
      aux = M_ZERO
      do i = 1, m%np
        call mesh_r(m, i, ri, a=geo%atom(ia)%x)
        if (ri > r) cycle
        aux(i, :) = md(i, :)
      end do
      lmm(1, ia) = dmf_integrate(m, aux(:, 1))
      lmm(2, ia) = dmf_integrate(m, aux(:, 2))
      lmm(3, ia) = dmf_integrate(m, aux(:, 3))
    end do
    deallocate(md, aux)

    call pop_sub()
  end subroutine states_local_magnetic_moments


  ! ---------------------------------------------------------
  ! This routine (obviously) assumes complex wave-functions
  subroutine calc_paramagnetic_current(gr, st, jp)
    type(grid_t),   intent(inout) :: gr
    type(states_t), intent(inout) :: st
    FLOAT,          intent(out)   :: jp(:,:,:)  ! (NP, NDIM, st%d%nspin)

    integer :: ik, p, sp, k
    CMPLX, allocatable :: grad(:,:)
#ifdef  HAVE_MPI
    FLOAT, allocatable :: red(:,:,:)
#endif

    call push_sub('states.calc_paramagnetic_current')

    ASSERT(st%d%wfs_type == M_CMPLX)

    if(st%d%ispin == SPIN_POLARIZED) then
      sp = 2
    else
      sp = 1
    end if

    jp = M_ZERO
    ALLOCATE(grad(NP_PART, NDIM), NP_PART*NDIM)

    do ik = 1, st%d%nik, sp
      do p  = st%st_start, st%st_end
        call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 1, p, ik), grad)

        ! spin-up density
        do k = 1, NDIM
          jp(:, k, 1) = jp(:, k, 1) + st%d%kweights(ik)*st%occ(p, ik)       &
            * aimag(conjg(st%zpsi(:, 1, p, ik)) * grad(:, k))
        end do

        ! spin-down density
        if(st%d%ispin == SPIN_POLARIZED) then
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 1, p, ik+1), grad)

          do k = 1, NDIM
            jp(:, k, 2) = jp(:, k, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1) &
              * aimag(conjg(st%zpsi(:, 1, p, ik+1)) * grad(:, k))
          end do

          ! WARNING: the next lines DO NOT work properly
        else if(st%d%ispin == SPINORS) then ! off-diagonal densities
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 2, p, ik), grad)

          do k = 1, NDIM
            jp(:, k, 2) = jp(:, k, 2) + st%d%kweights(ik)*st%occ(p, ik)     &
              * aimag(conjg(st%zpsi(:, 2, p, ik)) * grad(:, k))
          end do
        end if

      end do
    end do
    deallocate(grad)

#if defined(HAVE_MPI)
    ALLOCATE(red(NP_PART, NDIM, st%d%nspin), NP_PART*NDIM*st%d%nspin)
    call MPI_Allreduce(jp(1, 1, 1), red(1, 1, 1), NP*NDIM*st%d%nspin,       &
      MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    jp = red
    deallocate(red)
#endif

    call pop_sub()
  end subroutine calc_paramagnetic_current


  ! ---------------------------------------------------------
  subroutine states_calc_physical_current(gr, st, j)
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: st
    FLOAT,            intent(out)   :: j(:,:,:)   ! j(NP, NDIM, st%d%nspin)

    call push_sub('states.states_calc_physical_current')

    ! Paramagnetic contribution to the physical current
    call calc_paramagnetic_current(gr, st, j)

    ! TODO
    ! Diamagnetic contribution to the physical current

    call pop_sub()
  end subroutine states_calc_physical_current


  ! ---------------------------------------------------------
  subroutine states_distribute_nodes(st, mc)
    type(states_t),    intent(inout) :: st
    type(multicomm_t), intent(in)    :: mc

#ifdef HAVE_MPI
    integer :: sn, sn1, r, j, k, i, ii, st_start, st_end
#endif

    ! defaults
    st%node(:)  = 0
    st%st_start = 1
    st%st_end   = st%nst
    st%parallel_in_states = .false.
    call mpi_grp_init(st%mpi_grp, -1)

#if defined(HAVE_MPI)
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then

      st%parallel_in_states = .true.
      call mpi_grp_init(st%mpi_grp, mc%group_comm(P_STRATEGY_STATES))

      if(st%nst < st%mpi_grp%size) then
        message(1) = "Have more processors than necessary"
        write(message(2),'(i4,a,i4,a)') st%mpi_grp%size, " processors and ", st%nst, " states."
        call write_fatal(2)
      end if

      do k = 0, st%mpi_grp%size-1
        i = st%nst / st%mpi_grp%size
        ii = st%nst - i*st%mpi_grp%size
        if(ii > 0 .and. k < ii) then
          i = i + 1
          st_start = k*i + 1
          st_end = st_start + i - 1
        else
          st_end = st%nst - (st%mpi_grp%size - k - 1)*i
          st_start = st_end - i + 1
        end if
        write(message(1),'(a,i4,a,i4,a,i4)') 'Info: Nodes in states-group ', k,&
                                             ' will manage states', st_start, " - ", st_end
        call write_info(1)
        if(st%mpi_grp%rank .eq. k) then
          st%st_start = st_start
          st%st_end = st_end
        endif
      end do

      sn  = st%nst/st%mpi_grp%size
      sn1 = sn + 1
      r  = mod(st%nst, st%mpi_grp%size)
      do j = 1, r
        st%node((j-1)*sn1+1:j*sn1) = j - 1
      end do
      k = sn1*r
      call MPI_Barrier(st%mpi_grp%comm, mpi_err)
      do j = 1, st%mpi_grp%size - r
        st%node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
      end do
    end if
#endif

  end subroutine states_distribute_nodes


  ! ---------------------------------------------------------
  subroutine states_output(st, gr, dir, outp)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    character(len=*), intent(in)    :: dir
    type(output_t),   intent(in)    :: outp

    integer :: ik, ist, idim, is, id, ierr
    character(len=80) :: fname
    FLOAT :: u
    FLOAT, allocatable :: dtmp(:), elf(:,:)

    call push_sub('states_inc.states_output')
    
    u = M_ONE/units_out%length%factor**NDIM

    if(iand(outp%what, output_density).ne.0) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp%how, dir, fname, gr%m, gr%sb, st%rho(:, is), u, ierr)
      end do
    end if

    if(iand(outp%what, output_pol_density).ne.0) then
      ALLOCATE(dtmp(NP), NP)
      do idim=1, NDIM
        do is = 1, st%d%nspin
          dtmp(1:NP)=st%rho(:,is)*gr%m%x(:,idim)
          write(fname, '(a,i1,a,i1)') 'dipole_density-', is, '-',idim
          call doutput_function(outp%how, dir, fname, gr%m, gr%sb, dtmp(:), u, ierr)
        end do
      end do
    end if


    if(iand(outp%what, output_current).ne.0) then
      ! calculate current first
      call calc_paramagnetic_current(gr, st, st%j)
      do is = 1, st%d%nspin
        do id = 1, NDIM
          write(fname, '(a,i1,a,a)') 'current-', is, '-', index2axis(id)
          call doutput_function(outp%how, dir, fname, gr%m, gr%sb, st%j(:, id, is), u, ierr)
        end do
      end do
    end if


    if(iand(outp%what, output_wfs).ne.0) then
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = 1, st%d%nik
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i3.3,a,i1)') 'wf-', ik, '-', ist, '-', idim
              if (st%d%wfs_type == M_REAL) then
                call doutput_function(outp%how, dir, fname, gr%m, gr%sb, &
                     st%dpsi(1:, idim, ist, ik), sqrt(u), ierr)
              else
                call zoutput_function(outp%how, dir, fname, gr%m, gr%sb, &
                     st%zpsi(1:, idim, ist, ik), sqrt(u), ierr)
              end if
            end do
          end do
        end if
      end do
    end if

    if(iand(outp%what, output_wfs_sqmod).ne.0) then
      ALLOCATE(dtmp(NP_PART), NP_PART)
      do ist = 1, st%nst
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = 1, st%d%nik
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
              if (st%d%wfs_type == M_REAL) then
                dtmp = abs(st%dpsi(:, idim, ist, ik))**2
              else
                dtmp = abs(st%zpsi(:, idim, ist, ik))**2
              end if
              call doutput_function (outp%how, dir, fname, gr%m, gr%sb, dtmp, u, ierr)
            end do
          end do
        end if
      end do
      deallocate(dtmp)
    end if

    if(NDIM .eq. 3) then ! If the dimensions is not three, the ELF calculation will not work.
      if(  iand(outp%what, output_elf).ne.0  ) then ! First, ELF in real space.
        ALLOCATE(elf(1:gr%m%np,1:st%d%nspin),gr%m%np*st%d%nspin)
        call states_calc_elf(st, gr, elf)
        do is = 1, st%d%nspin
          write(fname, '(a,a,i1)') 'elf_rs', '-', is
          call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, elf(:,is), M_ONE, ierr)
        end do
        deallocate(elf)
      end if

      if(  iand(outp%what, output_elf_fs).ne.0  ) then ! Second, ELF in Fourier space.
        ALLOCATE(elf(1:gr%m%np,1:st%d%nspin),gr%m%np*st%d%nspin)
        call states_calc_elf_fs(st, gr, elf)
        do is = 1, st%d%nspin
          write(fname, '(a,a,i1)') 'elf_fs', '-', is
          call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, elf(:,is), M_ONE, ierr)
        end do
        deallocate(elf)
      end if
    end if

    call pop_sub()
  end subroutine states_output


  ! ---------------------------------------------------------
  ! (time-dependent) electron localization function, (TD)ELF.
  ! ---------------------------------------------------------
  subroutine states_calc_elf(st,gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,             intent(inout):: elf(:,:)
    FLOAT,   optional, intent(inout):: de(:,:)

    FLOAT :: f, d, s
    integer :: i, is, ik, ist, idim
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: r(:), gradr(:,:), j(:,:)

    FLOAT, parameter :: dmin = CNST(1e-10)
#if defined(HAVE_MPI)
    FLOAT, allocatable :: reduce_elf(:)
#endif

    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, NDIM), NP*NDIM)
      ALLOCATE(    j(NP, NDIM), NP*NDIM)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO
      
      elf(1:NP,is) = M_ZERO

      ALLOCATE(psi_fs(NP_PART),  NP_PART)
      ALLOCATE(gpsi  (NP, NDIM), NP*NDIM)
      do ik = is, st%d%nik, st%d%nspin
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim

            if (st%d%wfs_type == M_REAL) then
              psi_fs(:) = cmplx(st%dpsi(:, idim, ist, ik), KIND=PRECISION)
            else
              psi_fs(:) = st%zpsi(:, idim, ist, ik)
            end if

            call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, NDIM
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              elf(i,is) = elf(i,is) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                   sum(abs(gpsi(i, 1:NDIM))**2)
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ALLOCATE(reduce_elf(1:NP), NP)
        call MPI_Allreduce(elf(1, is), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        elf(1:NP, is) = reduce_elf(1:NP)
        call MPI_Allreduce(r(1), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        r(1:NP) = reduce_elf(1:NP)
        do i = 1, NDIM
          call MPI_Allreduce(gradr(1, i), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
          gradr(1:NP, i) = reduce_elf(1:NP)
          call MPI_Allreduce(j(1, i), reduce_elf(1), NP, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
          j(1:NP, i) = reduce_elf(1:NP)
        end do
        deallocate(reduce_elf)
      end if
#endif

      do i = 1, NP
        if(r(i) >= dmin) then
          elf(i,is) = elf(i,is) - (M_FOURTH*sum(gradr(i, 1:NDIM)**2) + sum(j(i, 1:NDIM)**2))/(s*r(i))
        end if
      end do
    
      if(present(de)) de(1:NP,is)=elf(1:NP,is)

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          elf(i,is) = M_ONE/(M_ONE + (elf(i,is)/d)**2)
        else
          elf(i,is) = M_ZERO
        end if
      end do

      deallocate(r, gradr, j)

    end do do_is

  end subroutine states_calc_elf

  ! ---------------------------------------------------------
  ! ELF function in Fourier space. Not tested.
  ! ---------------------------------------------------------
  subroutine states_calc_elf_fs(st, gr, elf, de)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    FLOAT,             intent(inout):: elf(:,:)
    FLOAT,   optional, intent(inout):: de(:,:)

    FLOAT :: f, d, s
    integer :: i, is, ik, ist, idim
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: r(:), gradr(:,:), j(:,:)
    type(dcf_t) :: dcf_tmp
    type(zcf_t) :: zcf_tmp

    FLOAT, parameter :: dmin = CNST(1e-10)
    
    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if
 
    if (st%d%wfs_type == M_REAL) then
      call dcf_new(gr%m%l, dcf_tmp)
      call dcf_fft_init(dcf_tmp, gr%sb)
    else
      call zcf_new(gr%m%l, zcf_tmp)
      call zcf_fft_init(zcf_tmp, gr%sb)
    end if

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, NDIM), NP*NDIM)
      ALLOCATE(    j(NP, NDIM), NP*NDIM)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO

      elf(1:NP,is) = M_ZERO

      ALLOCATE(psi_fs(NP_PART),  NP_PART)
      ALLOCATE(gpsi  (NP, NDIM), NP*NDIM)
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            
            if (st%d%wfs_type == M_REAL) then
              call dmf2mf_RS2FS(gr%m, st%dpsi(:, idim, ist, ik), psi_fs(:), dcf_tmp)
              call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)
              do i = 1, NDIM
                gpsi(:,i) = gpsi(:,i) * gr%m%h(i)**2 * real(dcf_tmp%n(i), PRECISION) / (M_TWO*M_PI)
              end do
            else
              call zmf2mf_RS2FS(gr%m, st%zpsi(:, idim, ist, ik), psi_fs(:), zcf_tmp)
              call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)
              do i = 1, NDIM
                gpsi(:,i) = gpsi(:,i) * gr%m%h(i)**2 * real(zcf_tmp%n(i), PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, NDIM
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              if(r(i) >= dmin) then
                elf(i,is) = elf(i,is) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                     sum(abs(gpsi(i, 1:NDIM))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, NP
        if(r(i) >= dmin) then
          elf(i,is) = elf(i,is) - (M_FOURTH*sum(gradr(i, 1:NDIM)**2) + sum(j(i, 1:NDIM)**2))/(s*r(i))
        end if
      end do
    
      if(present(de)) de(1:NP,is)=elf(1:NP,is)

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          elf(i,is) = M_ONE/(M_ONE + (elf(i,is)/d)**2)
        else
          elf(i,is) = M_ZERO
        end if
      end do

      deallocate(r, gradr, j)

    end do do_is

    if (st%d%wfs_type == M_REAL) then
      call dcf_free(dcf_tmp)
    else
      call zcf_free(zcf_tmp)
    end if

  contains

    subroutine dmf2mf_RS2FS(m, fin, fout, c)
      type(mesh_t),  intent(in)    :: m
      FLOAT,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(dcf_t),   intent(inout) :: c
    
      call dcf_alloc_RS(c)
      call dcf_alloc_FS(c)
      call dmf2cf(m, fin, c)
      call dcf_RS2FS(c)
      call dcf_FS2mf(m, c, fout)
      call dcf_free_RS(c)
      call dcf_free_FS(c)
    end subroutine dmf2mf_RS2FS

    subroutine zmf2mf_RS2FS(m, fin, fout, c)
      type(mesh_t),  intent(in)    :: m
      CMPLX,         intent(in)    :: fin(:)
      CMPLX,         intent(out)   :: fout(:)
      type(zcf_t),   intent(inout) :: c
    
      call zcf_alloc_RS(c)
      call zcf_alloc_FS(c)
      call zmf2cf(m, fin, c)
      call zcf_RS2FS(c)
      call zcf_FS2mf(m, c, fout)
      call zcf_free_RS(c)
      call zcf_free_FS(c)
    end subroutine zmf2mf_RS2FS

  end subroutine states_calc_elf_fs


#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states_m
