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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module states_m
  use crystal_m
  use datasets_m
  use functions_m
  use geometry_m
  use grid_m
  use io_m
  use lib_basic_alg_m
  use lib_oct_m
  use lib_oct_parser_m
  use math_m
  use mesh_function_m
  use mesh_m
  use mpi_m
  use multicomm_m
  use profiling_m
  use simul_box_m
  use string_m
  use varinfo_m

  implicit none

  private
  public ::                         &
    states_t,                       &
    states_dim_t,                   &
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
    states_write_dos,               &
    states_write_bands,             &
    states_write_fermi_energy,      &
    states_degeneracy_matrix,       &
    states_spin_channel,            &
    states_calc_projection,         &
    states_calc_dens,               &
    kpoints_write_info,             &
    wfs_are_complex,                &
    wfs_are_real,                   &
    states_dump,                    &
    assignment(=)


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
    dstates_calc_momentum,          &
    zstates_calc_momentum,          &
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
    FLOAT, pointer :: momentum(:, :, :)

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
    integer :: wfs_type             ! real (M_REAL) or complex (M_CMPLX) wavefunctions
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


  ! Parameters...
  integer, public, parameter ::     &
    UNPOLARIZED    = 1,             &
    SPIN_POLARIZED = 2,             &
    SPINORS        = 3

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
    type(states_t), intent(inout) :: st
    call push_sub('states.states_null')

    nullify(st%dpsi, st%zpsi, st%rho, st%j, st%rho_core, st%eigenval, st%occ, st%mag)

    nullify(st%d)
    ALLOCATE(st%d, 1)
    nullify(st%d%kpoints, st%d%kweights)

    ! By default, calculations use real wave-functions
    st%d%wfs_type = M_REAL

    call pop_sub()
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo


    FLOAT :: excess_charge, r
    integer :: nempty, i, j
    C_POINTER :: blk

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

    call geometry_val_charge(geo, st%val_charge)
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
    call states_choose_kpoints(st%d, gr%sb, geo)

    ! Periodic systems require complex wave-functions
    if(simul_box_is_periodic(gr%sb)) st%d%wfs_type = M_CMPLX


    ! we now allocate some arrays
    ALLOCATE(st%occ     (st%nst, st%d%nik),      st%nst*st%d%nik)
    ALLOCATE(st%eigenval(st%nst, st%d%nik),      st%nst*st%d%nik)
    ALLOCATE(st%momentum(3, st%nst, st%d%nik), 3*st%nst*st%d%nik)
    ! allocate space for formula strings that define user defined states
    ALLOCATE(st%user_def_states(st%d%dim, st%nst, st%d%nik), st%d%dim*st%nst*st%d%nik)
    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%mag(st%nst, st%d%nik, 2), st%nst*st%d%nik*2)
    end if

    ! initially we mark all 'formulas' as undefined
    st%user_def_states(1:st%d%dim, 1:st%nst, 1:st%d%nik) = 'undefined'

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
    st%node(1:st%nst) = 0

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

    C_POINTER :: blk
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

      call messages_print_stress(stdout, trim('Substitution of orbitals'))

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

              write(message(1), '(a,3i5)') 'Substituting state of orbital with k, ist, dim = ', ik, is, id
              write(message(2), '(2a)') '  with the expression:'
              write(message(3), '(2a)') '  ',trim(st%user_def_states(id, is, ik))
              call write_info(3)

              ! convert to C string
              call conv_to_C_string(st%user_def_states(id, is, ik))

              ! fill states with user defined formulas
              do ip = 1, mesh%np
                x = mesh%x(ip, :)
                r = sqrt(sum(x(:)**2))

                ! parse user defined expressions
                call loct_parse_expression(psi_re, psi_im,             &
                  x(1), x(2), x(3), r, M_ZERO, st%user_def_states(id, is, ik))
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
      call messages_print_stress(stdout)

    else
      message(1) = '"UserDefinesStates" has to be specified as block.'
      call write_fatal(1)
    end if

    call pop_sub()
  end subroutine states_read_user_def_orbitals


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr, geo)
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo
    call push_sub('states.states_densities_init')

    ! allocate arrays for charge and current densities
    ALLOCATE(st%rho(NP_PART, st%d%nspin), NP_PART*st%d%nspin)
    ALLOCATE(st%j(NP_PART, NDIM, st%d%nspin), NP_PART*NDIM*st%d%nspin)
    st%rho = M_ZERO
    st%j   = M_ZERO
    if(geo%nlcc) ALLOCATE(st%rho_core(gr%m%np), gr%m%np)

    call pop_sub()
  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_t), intent(inout) :: stout
    type(states_t), intent(in)  :: stin

    integer :: i

    call states_null(stout)

    stout%d%wfs_type = stin%d%wfs_type
    stout%d%dim = stin%d%dim
    stout%d%nik = stin%d%nik
    stout%d%nik_axis(1:MAX_DIM) = stout%d%nik_axis(1:MAX_DIM)
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
      i = size(stin%d%kweights, 1)
      ALLOCATE(stout%d%kweights(size(stin%d%kweights, 1)), i)
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
      deallocate(st%rho, st%occ, st%eigenval, st%momentum, st%node)
      nullify   (st%rho, st%occ, st%eigenval, st%momentum, st%node)
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
  subroutine states_write_bands(dir, nst, st, sb)
    character(len=*),  intent(in) :: dir    
    integer,           intent(in) :: nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: i, ik, j, ns, is
    integer, allocatable :: iunit(:)
    FLOAT   :: factor(MAX_DIM)
    logical :: grace_mode, gnuplot_mode
    character(len=80) :: filename    

    call push_sub('states.states_write_bands')

    if(.not.mpi_grp_is_root(mpi_world)) return

    !%Variable OutputBandsGnuplotMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in gnuplot friendly format
    !%End
    call loct_parse_logical(check_inp('OutputBandsGnuplotMode'), .true., gnuplot_mode)

    !%Variable OutputBandsGraceMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in grace friendly format
    !%End
    call loct_parse_logical(check_inp('OutputBandsGraceMode'), .false., grace_mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ALLOCATE(iunit(0:ns-1), ns)

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do i = 1, MAX_DIM
      factor(i) = M_ONE
      if (sb%klat(i,i) /= M_ZERO) factor(i) = sb%klat(i,i)      
    end do

    if (gnuplot_mode) then
      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i1.1,a)') 'bands-gp-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-gp.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in gnuplot format
      do j = 1, nst
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            write(iunit(is), '(1x,6f14.8,3x,f14.8)')            &
              st%d%kpoints(1:MAX_DIM, ik+is),                   & ! unscaled
              st%d%kpoints(1:MAX_DIM, ik+is)/factor(1:MAX_DIM), & ! scaled
              st%eigenval(j, ik+is)/units_out%energy%factor
          end do
        end do
        do is = 0, ns-1
          write(iunit(is), '(a)') ''
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end if

    if (grace_mode) then
      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i1.1,a)') 'bands-grace-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-grace.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a, i6)') '# kx ky kz (unscaled), kx ky kz (scaled), bands:', nst
      end do

      ! output bands in xmgrace format, i.e.:
      ! k_x, k_y, k_z, e_1, e_2, ..., e_n
      do ik = 1, st%d%nik, ns
        do is = 0, ns-1
          write(iunit(is), '(1x,6f14.8,3x,16384f14.8)')         &
            st%d%kpoints(1:MAX_DIM, ik+is),                     & ! unscaled
            st%d%kpoints(1:MAX_DIM, ik+is)/factor(1:MAX_DIM),   & ! scaled
            (st%eigenval(j, ik+is)/units_out%energy%factor, j = 1, nst)
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do        
    end if

    deallocate(iunit)

    call pop_sub()
  end subroutine states_write_bands


  ! ---------------------------------------------------------
  subroutine states_write_dos(dir, st)
    character(len=*), intent(in) :: dir
    type(states_t),   intent(in) :: st

    integer :: ie, ik, ist, epoints, is, ns
    integer, allocatable :: iunit(:)
    FLOAT   :: emin, emax, de, gamma, energy
    FLOAT   :: evalmax, evalmin, tdos, eextend
    FLOAT, allocatable :: dos(:,:,:)
    character(len=64)  :: filename

    call push_sub('states.states_write_dos')

    evalmin = minval(st%eigenval/units_out%energy%factor)
    evalmax = maxval(st%eigenval/units_out%energy%factor)
    ! we extend the energy mesh by this amount
    eextend  = (evalmax - evalmin) / M_FOUR

    !%Variable DOSEnergyMin
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Lower bound for the energy mesh of the DOS
    !%End
    call loct_parse_float(check_inp('DOSEnergyMin'), evalmin - eextend, emin)

    !%Variable DOSEnergyMax
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS
    !%End
    call loct_parse_float(check_inp('DOSEnergyMax'), evalmax + eextend, emax)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points octopus should use for 
    !% the DOS energy grid
    !%End
    call loct_parse_int(check_inp('DOSEnergyPoints'), 500, epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used to sum 
    !% up the DOS sum
    !%End
    call loct_parse_float(check_inp('DOSGamma'), &
      CNST(0.008)/units_out%energy%factor, gamma)

    ! spacing for energy mesh
    de = (emax - emin) / (epoints - 1)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! space for state dependent DOS
    ALLOCATE(dos(epoints, st%nst, 0:ns-1), epoints*st%nst*ns)
    ALLOCATE(iunit(0:ns-1), ns)    

    ! compute band/spin resolved density of states
    do ist = 1, st%nst

      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i4.4,a,i1.1,a)') 'dos-', ist, '-', is+1,'.dat'
        else
          write(filename, '(a,i4.4,a)') 'dos-', ist, '.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, band resolved DOS'
      end do
      
      do ie = 1, epoints
        energy = emin + (ie - 1) * de
        dos(ie, ist, :) = M_ZERO
        ! sum up Lorentzians
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            dos(ie, ist, is) = dos(ie, ist, is) + st%d%kweights(ik+is) * M_ONE/M_Pi * &
              gamma / ( (energy - st%eigenval(ist, ik+is)/units_out%energy%factor)**2 + gamma**2 )
          end do
        end do
        do is = 0, ns-1
          write(message(1), '(2f12.6)') energy, dos(ie, ist, is)
          call write_info(1, iunit(is))
        end do
      end do

      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end do

    ! for spin polarized calculations also output spin resolved tdos
    if(st%d%nspin .gt. 1) then    
      do is = 0, ns-1
        write(filename, '(a,i1.1,a)') 'total-dos-', is+1,'.dat'
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, total DOS (spin resolved)'

        do ie = 1, epoints
          energy = emin + (ie - 1) * de
          tdos = M_ZERO
          do ist = 1, st%nst
            tdos = tdos + dos(ie, ist, is)
          end do
          write(message(1), '(2f12.6)') energy, tdos
          call write_info(1, iunit(is))
        end do

        call io_close(iunit(is))
      end do
    end if


    iunit(0) = io_open(trim(dir)//'/'//'total-dos.dat', action='write')    

    ! compute total density of states
    do ie = 1, epoints
      energy = emin + (ie - 1) * de
      tdos = M_ZERO
      do ist = 1, st%nst
        do is = 0, ns-1
          tdos = tdos + dos(ie, ist, is)
        end do
      end do
      write(message(1), '(2f12.6)') energy, tdos
      call write_info(1, iunit(0))
    end do

    call io_close(iunit(0))

    deallocate(iunit, dos)

    call pop_sub()
  end subroutine states_write_dos


  ! ---------------------------------------------------------
  subroutine states_write_fermi_energy(dir, st, m, sb)
    character(len=*),  intent(in) :: dir
    type(states_t), intent(inout) :: st
    type(mesh_t),      intent(in) :: m
    type(simul_box_t), intent(in) :: sb

    integer :: iunit, i
    FLOAT :: scale, maxdos
    FLOAT :: factor(MAX_DIM)

    call push_sub('states.states_write_fermi_energy')

    call states_fermi(st, m)
    
    iunit = io_open(trim(dir)//'/'//'bands-efermi.dat', action='write')    

    scale = units_out%energy%factor

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do i = 1, MAX_DIM
      factor(i) = M_ONE
      if (sb%klat(i,i) /= M_ZERO) factor(i) = sb%klat(i,i)
    end do

    ! write fermi energy in a format that can be used together 
    ! with bands.dat
    write(message(1), '(a)') '# Fermi energy in a format compatible with bands-gp.dat'

    write(message(2), '(7f12.6)')          &
      minval(st%d%kpoints(1,:)),           &
      minval(st%d%kpoints(2,:)),           &
      minval(st%d%kpoints(3,:)),           &
      minval(st%d%kpoints(1,:)/factor(1)), &
      minval(st%d%kpoints(2,:)/factor(2)), &
      minval(st%d%kpoints(3,:)/factor(3)), &
      st%ef/scale

    ! Gamma point
    write(message(3), '(7f12.6)')          &
      (M_ZERO, i = 1, 6),                  &
      st%ef/scale
    
    write(message(4), '(7f12.6)')          &
      maxval(st%d%kpoints(1,:)),           &
      maxval(st%d%kpoints(2,:)),           &
      maxval(st%d%kpoints(3,:)),           &
      maxval(st%d%kpoints(1,:)/factor(1)), &
      maxval(st%d%kpoints(2,:)/factor(2)), &
      maxval(st%d%kpoints(3,:)/factor(3)), &
      st%ef/scale

    call write_info(4, iunit)
    call io_close(iunit)

    ! now we write the same information so that it can be used 
    ! together with total-dos.dat
    iunit = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')    

    write(message(1), '(a)') '# Fermi energy in a format compatible with total-dos.dat'    
    
    ! this is the maximum that tdos can reach
    maxdos = sum(st%d%kweights) * st%nst

    write(message(2), '(4f12.6)') st%ef/scale, M_ZERO
    write(message(3), '(4f12.6)') st%ef/scale, maxdos

    call write_info(3, iunit)
    call io_close(iunit)

    call pop_sub()
  end subroutine states_write_fermi_energy


  ! -------------------------------------------------------
  subroutine states_degeneracy_matrix(st)
    type(states_t), intent(in) :: st

    integer :: is, js, inst, inik, dsize, iunit
    integer, allocatable :: eindex(:,:), sindex(:)
    integer, allocatable :: degeneracy_matrix(:, :)
    FLOAT,   allocatable :: eigenval_sorted(:)
    FLOAT :: degen_thres, evis, evjs

    call push_sub('states.states_degeneracy_matrix')

    ALLOCATE(eigenval_sorted(st%nst*st%d%nik),   st%nst*st%d%nik)
    ALLOCATE(         sindex(st%nst*st%d%nik),   st%nst*st%d%nik)
    ALLOCATE(      eindex(2, st%nst*st%d%nik), 2*st%nst*st%d%nik)
    dsize = st%nst*st%d%nik * st%nst*st%d%nik
    ALLOCATE(degeneracy_matrix(st%nst*st%d%nik, st%nst*st%d%nik), dsize)
    
    ! convert double index "inst, inik" to single index "is"
    ! and keep mapping array
    is = 1
    do inst = 1, st%nst
      do inik = 1, st%d%nik
        eigenval_sorted(is) = st%eigenval(inst, inik)        
        eindex(1, is) = inst
        eindex(2, is) = inik
        is = is + 1
      end do
    end do
    
    ! sort eigenvalues
    call sort(eigenval_sorted, sindex)

    !%Variable DegeneracyMatrixThreshold
    !%Type float
    !%Default 1e-5
    !%Section States
    !%Description
    !% A state j with energy E_j will be considered degenerate with a state
    !% with energy E_i, if  E_i - threshold < E_j < E_i + threshold.
    !%End
    call loct_parse_float(check_inp('DegeneracyMatrixThreshold'), CNST(1e-5), degen_thres)    

    ! setup degeneracy matrix. the matrix summarizes the degeneracy relations 
    ! among the states
    degeneracy_matrix = 0

    do is = 1, st%nst*st%d%nik
      do js = 1, st%nst*st%d%nik

        ! a state is always degenerate to itself
        if ( is.eq.js ) cycle

        evis = st%eigenval(eindex(1, sindex(is)), eindex(2, sindex(is)))
        evjs = st%eigenval(eindex(1, sindex(js)), eindex(2, sindex(js)))

        ! is evjs in the "evis plus minus threshold" bracket?
        if( (evjs.gt.evis - degen_thres).and.(evjs.lt.evis + degen_thres) ) then
          ! mark forward scattering states with +1 and backward scattering
          ! states with -1
          degeneracy_matrix(is, js) = &
            sign(M_ONE, st%momentum(1, eindex(1, sindex(js)), eindex(2, sindex(js))))
        end if

      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then

      ! write matrix to "tmp/restart_gs" directory
      iunit = io_open('tmp/restart_gs/degeneracy_matrix', action='write', is_tmp = .true.)    

      write(iunit, '(a)') '# index  kx ky kz  eigenvalue  degeneracy matrix'

      do is = 1, st%nst*st%d%nik
        write(iunit, '(i6,4e24.16,32767i3)') is, st%d%kpoints(:, eindex(2, sindex(is))), &
          eigenval_sorted(is), (degeneracy_matrix(is, js), js = 1, st%nst*st%d%nik)
      end do
    
      call io_close(iunit)
    
      ! write index vectors to "tmp/restart_gs" directory
      iunit = io_open('tmp/restart_gs/index_vectors', action='write', is_tmp = .true.)    
      
      write(iunit, '(a)') '# index  sindex  eindex1 eindex2'

      do is = 1, st%nst*st%d%nik
        write(iunit,'(4i6)') is, sindex(is), eindex(1, sindex(is)), eindex(2, sindex(is))
      end do
    
      call io_close(iunit)
    end if

    deallocate(eigenval_sorted, sindex, eindex)
    deallocate(degeneracy_matrix)
    
    call pop_sub()
  end subroutine states_degeneracy_matrix


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
  logical function wfs_are_complex(st) result (wac)
    type(states_t),    intent(in) :: st
    wac = (st%d%wfs_type == M_CMPLX)
  end function wfs_are_complex


  ! ---------------------------------------------------------
  logical function wfs_are_real(st) result (war)
    type(states_t),    intent(in) :: st
    war = (st%d%wfs_type == M_REAL)
  end function wfs_are_real


  ! ---------------------------------------------------------
  subroutine states_dump(st, iunit)
    type(states_t), intent(in) :: st
     integer,       intent(in) :: iunit

     call push_sub('states.states_dump')
     
     write(iunit, '(a20,1i10)')  'nst=                ', st%nst
     write(iunit, '(a20,1i10)')  'dim=                ', st%d%dim
     write(iunit, '(a20,1i10)')  'nik=                ', st%d%nik

     call pop_sub()
  end subroutine states_dump

#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states_m
