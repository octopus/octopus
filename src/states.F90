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

module states
  use global
  use varinfo
  use messages
  use syslabels
  use lib_oct
  use lib_oct_parser
  use io
  use lib_basic_alg
  use lib_adv_alg
  use math
  use mesh
  use grid
  use simul_box
  use functions
  use mesh_function
  use cube_function
  use output
  use geometry
  use crystal
  use mpi_mod
  use multicomm_mod

  implicit none

  private
  public ::                         &
    states_type,                    &
    states_dim_type,                &
    states_init,                    &
    states_densities_init,          &
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
    kpoints_write_info,             &
    dstates_calc_dens,              &
    zstates_calc_dens,              &
    dstates_gram_schmidt,           &
    zstates_gram_schmidt,           &
    dstates_dotp,                   &
    zstates_dotp,                   &
    dstates_nrm2,                   &
    zstates_nrm2,                   &
    dstates_residue,                &
    zstates_residue,                &
    dstates_output,                 &
    zstates_output,                 &
    dstates_mpdotp,                 &
    zstates_mpdotp,                 &
    dstates_calc_angular,           &
    zstates_calc_angular,           &
    states_distribute_nodes

  type states_type

    type(states_dim_type), pointer :: d
    integer :: nst                  ! Number of states in each irreducible subspace

    ! pointers to the wavefunctions
    FLOAT, pointer :: dpsi(:,:,:,:)
    CMPLX, pointer :: zpsi(:,:,:,:)

    ! the densities and currents (after all we are doing DFT :)
    FLOAT, pointer :: rho(:,:)
    FLOAT, pointer :: j(:,:,:)

    FLOAT, pointer :: rho_core(:)   ! core charge for nl core corrections

    FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
    logical           :: fixed_occ  ! should the occupation numbers be fixed?
    FLOAT, pointer :: occ(:,:)      ! the occupation numbers
    FLOAT, pointer :: mag(:, :, :)

    FLOAT :: qtot                   ! (-) The total charge in the system (used in Fermi)
    FLOAT :: val_charge             ! valence charge

    FLOAT :: el_temp                ! electronic temperature for the Fermi function
    FLOAT :: ef                     ! the fermi energy

    ! This is stuff needed for the parallelization in states
    logical :: parallel_in_states   ! am I parallel in states?
    integer :: numprocs             ! how many nodes are in this group
    integer :: rank                 ! who am I in this group?
    integer :: comm                 ! my communicator

    integer :: st_start, st_end     ! needed for some parallel parts
    integer, pointer :: node(:)     ! To which node belongs each state.
  end type states_type

  type states_dim_type
    integer :: dim                  ! Dimension of the state (one or two for spinors)
    integer :: nik                  ! Number of irreducible subspaces
    integer :: nik_axis(3)          ! Number of kpoints per axis
    integer :: ispin                ! spin mode (unpolarized, spin polarized, spinors)
    integer :: nspin                ! dimension of rho (1, 2 or 4)
    integer :: spin_channels        ! 1 or 2, wether spin is or not considered.
    logical :: cdft                 ! Are we using Current-DFT or not?
    FLOAT, pointer :: kpoints(:,:)  ! obviously the kpoints
    FLOAT, pointer :: kweights(:)   ! weights for the kpoint integrations
  end type states_dim_type

  ! Parameters...
  integer, public, parameter ::     &
    UNPOLARIZED    = 1,             &
    SPIN_POLARIZED = 2,             &
    SPINORS        = 3

  interface assignment (=)
    module procedure states_copy
  end interface

contains

  ! ---------------------------------------------------------
  subroutine states_null(st)
    type(states_type), intent(out) :: st

    nullify(st%dpsi, st%zpsi, st%rho, st%j, st%rho_core, st%eigenval, st%occ, st%mag)
    nullify(st%d); allocate(st%d)
    nullify(st%d%kpoints, st%d%kweights)
  end subroutine states_null


  ! ---------------------------------------------------------
  subroutine states_init(st, gr)
    type(states_type),    intent(inout) :: st
    type(grid_type),      intent(in)    :: gr

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
#if !defined(COMPLEX_WFNS)
    if(st%d%ispin == SPINORS) then
      message(1) = "Cannot use spinors with an executable compiled for real wavefunctions."
      call write_fatal(1)
    end if
#endif

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
#ifdef COMPLEX_WFNS
    if (st%d%cdft .and. st%d%ispin == SPINORS) then
      message(1) = "Sorry, Current DFT not working yet for spinors"
      call write_fatal(1)
    elseif (st%d%cdft) then
      message(1) = "Info: Using Current DFT"
      call write_info(1)
    end if
#else
    if (st%d%cdft) then
      message(1) = "Cannot use Current DFT with an executable compiled"
      message(2) = "for real wavefunctions."
      call write_fatal(2)
    end if
#endif

    ! For non-periodic systems this should just return the Gamma point
    call states_choose_kpoints(st%d, gr%sb, gr%geo)

    ! we now allocate some arrays
    allocate(st%occ     (st%nst, st%d%nik))
    allocate(st%eigenval(st%nst, st%d%nik))
    if(st%d%ispin == SPINORS) then
      allocate(st%mag(st%nst, st%d%nik, 2))
    end if

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

      ! first guest for occupation...paramagnetic configuration
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
    allocate(st%node(st%nst))
    st%node(:) = 0

    nullify(st%dpsi, st%zpsi)

    call pop_sub()
  end subroutine states_init


  ! ---------------------------------------------------------
  subroutine states_densities_init(st, gr)
    type(states_type),    intent(inout) :: st
    type(grid_type),      intent(in)    :: gr

    ! allocate arrays for charge and current densities
    allocate(st%rho(gr%m%np, st%d%nspin))
    if (st%d%cdft) then
      allocate(st%j(NP, NDIM, st%d%nspin))
      st%j = M_ZERO
    end if
    if(gr%geo%nlcc) allocate(st%rho_core(gr%m%np))

  end subroutine states_densities_init


  ! ---------------------------------------------------------
  subroutine states_copy(stout, stin)
    type(states_type), intent(in)  :: stin
    type(states_type), intent(out) :: stout

    call states_null(stout)

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
      allocate(stout%dpsi(size(stin%dpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik))
      stout%dpsi = stin%dpsi
    end if
    if(associated(stin%zpsi)) then
      allocate(stout%zpsi(size(stin%zpsi, 1), stin%d%dim, stin%st_start:stin%st_end, stin%d%nik))
      stout%zpsi = stin%zpsi
    end if
    if(associated(stin%rho)) then
      allocate(stout%rho(size(stin%rho, 1), size(stin%rho, 2)))
      stout%rho = stin%rho
    end if
    if(associated(stin%j)) then
      allocate(stout%j(size(stin%j, 1), size(stin%j, 2), size(stin%j, 3)))
      stout%j = stin%j
    end if
    if(associated(stin%rho_core)) then
      allocate(stout%rho_core(size(stin%rho_core, 1)))
      stout%rho_core = stin%rho_core
    end if
    if(associated(stin%eigenval)) then
      allocate(stout%eigenval(stin%st_start:stin%st_end, stin%d%nik))
      stout%eigenval = stin%eigenval
    end if
    stout%fixed_occ = stin%fixed_occ
    if(associated(stin%occ)) then
      allocate(stout%occ(size(stin%occ, 1), size(stin%occ, 2)))
      stout%occ = stin%occ
    end if
    if(associated(stin%mag)) then
      allocate(stout%mag(size(stin%mag, 1), size(stin%mag, 2), size(stin%mag, 3)))
      stout%mag = stin%mag
    end if
    if(associated(stin%d%kpoints)) then
      allocate(stout%d%kpoints(size(stin%d%kpoints, 1), size(stin%d%kpoints, 2)))
      stout%d%kpoints = stin%d%kpoints
    end if
    if(associated(stin%d%kweights)) then
      allocate(stout%d%kweights(size(stin%d%kpoints, 1)))
      stout%d%kweights = stin%d%kweights
    end if
    if(associated(stin%node)) then
      allocate(stout%node(size(stin%node)))
      stout%node = stin%node
    end if
  end subroutine states_copy


  ! ---------------------------------------------------------
  subroutine states_end(st)
    type(states_type), intent(inout) :: st

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

    call pop_sub()
  end subroutine states_end


  ! ---------------------------------------------------------
  ! generate a hydrogen s-wavefunction around a random point
  subroutine states_generate_random(st, m, ist_start_, ist_end_)
    type(states_type), intent(inout) :: st
    type(mesh_type),   intent(in)    :: m
    integer, optional, intent(in)    :: ist_start_, ist_end_

    integer :: ist, ik, id, ist_start, ist_end

    call push_sub('states.states_generate_random')

    ist_start = 1
    if(present(ist_start_)) ist_start = ist_start_
    ist_end = st%nst
    if(present(ist_end_)) ist_end = ist_end_

    do ik = 1, st%d%nik
      do ist = ist_start, ist_end
        do id = 1, st%d%dim
          call X(mf_random)(m, st%X(psi)(:, id, ist, ik))
        end do
        st%eigenval(ist, ik) = M_ZERO
      end do
      call X(states_gram_schmidt)(ist_end, m, st%d%dim, st%X(psi)(:,:,1:ist_end,ik), start = ist_start)
    end do

    call pop_sub()
  end subroutine states_generate_random


  ! ---------------------------------------------------------
  subroutine states_fermi(st, m)
    type(states_type), intent(inout) :: st
    type(mesh_type),   intent(in)    :: m

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
            st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
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
            st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
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
          st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
          st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
        end do
      end do
    end if

    call pop_sub(); return
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
    type(grid_type),   intent(in)  :: gr
    type(states_type), intent(in)  :: st
    integer,           intent(in)  :: lmax
    FLOAT,             intent(out) :: multipole(:, :) ! multipole((lmax + 1)**2, st%d%nspin)

    integer :: i, is, l, lm, add_lm
    FLOAT :: x(3), r, ylm
    FLOAT, allocatable :: f(:)

    call push_sub('states.states_calculate_multipoles')

    allocate(f(gr%m%np)); f = M_ZERO

    do is = 1, st%d%nspin

      f(:) = st%rho(:, is)
      multipole(1, is) = dmf_integrate(gr%m, f)

      if(lmax>0) then
        do i = 1, 3
          f(:) = st%rho(1:gr%m%np, is)*gr%m%x(1:gr%m%np, i)
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
  function states_eigenvalues_sum(st) result(e)
    type(states_type), intent(in) :: st
    FLOAT                         :: e

    integer :: ik

    e = M_ZERO
    do ik = 1, st%d%nik
      e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
        st%eigenval(st%st_start:st%st_end, ik))
    end do

  end function states_eigenvalues_sum


  ! ---------------------------------------------------------
  subroutine states_write_eigenvalues(iunit, nst, st, sb, error)
    integer,           intent(in) :: iunit, nst
    type(states_type), intent(in) :: st
    type(simul_box_type), intent(in) :: sb
    FLOAT,             intent(in), optional :: error(nst, st%d%nik)

    integer ik, j, ns, is
    FLOAT :: o, oplus, ominus

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_out%energy%abbrev) // ']'
    call write_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
      call write_info(1, iunit)
    end if

#ifdef HAVE_MPI
    if(mpiv%node == 0) then
#endif

      do ik = 1, st%d%nik, ns
        if(st%d%nik > ns) then
          write(iunit, '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
            st%d%kpoints(1, ik)*units_out%length%factor, ',',           &
            st%d%kpoints(2, ik)*units_out%length%factor, ',',           &
            st%d%kpoints(3, ik)*units_out%length%factor, ')'
        end if

        do is = 1, ns
          write(iunit, '(a4)', advance='no') '#st'
          if(present(error)) then
            write(iunit, '(1x,a12,3x,a12,2x,a10,i3,a1)', advance='no') &
              ' Eigenvalue', 'Occupation ', 'Error (', is, ')'
          else
            write(iunit, '(1x,a12,3x,a12,2x)', advance='no') &
              ' Eigenvalue', 'Occupation '
          end if
        end do
        write(iunit, '(1x)', advance='yes')

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

            write(iunit, '(i4)', advance='no') j
            if(simul_box_is_periodic(sb)) then
              if(st%d%ispin == SPINORS) then
                write(iunit, '(1x,f12.6,3x,f5.2,a1,f5.2)', advance='no') &
                  (st%eigenval(j, ik)-st%ef)/units_out%energy%factor, oplus, '/', ominus
                if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
              else
                write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                  (st%eigenval(j, ik+is))/units_out%energy%factor, o
                if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
              end if
            else
              if(st%d%ispin == SPINORS) then
                write(iunit, '(1x,f12.6,3x,f5.2,a1,f5.2)', advance='no') &
                  st%eigenval(j, ik)/units_out%energy%factor, oplus, '/', ominus
                if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
              else
                write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                  st%eigenval(j, ik+is)/units_out%energy%factor, o
                if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
              end if
            end if
          end do
          write(iunit, '(1x)', advance='yes')
        end do
      end do

#ifdef HAVE_MPI
    end if
#endif

  end subroutine states_write_eigenvalues


  ! ---------------------------------------------------------
  subroutine states_write_bands(iunit, nst, st, sb)
    integer,           intent(in) :: iunit, nst
    type(states_type), intent(in) :: st
    type(simul_box_type), intent(in) :: sb

    integer :: i, ik, j, mode, ns
    FLOAT :: factor(3)

    integer, parameter :: GNUPLOT = 1, &
      XMGRACE = 2

    call loct_parse_int(check_inp('BandsOutputMode'), GNUPLOT, mode)
    if(mode /= GNUPLOT .and. mode /= XMGRACE) then
      message(1) = "Input: BandsOutputMode must be 1 (gnuplot) or 2 (xmgrace)"
      call write_fatal(1)
    end if

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

#ifdef HAVE_MPI
    if(mpiv%node == 0) then
#endif

      ! define the scaling factor to output k_i/G_i, instead of k_i
      do i =1,3
        factor(i) = M_ONE
        if (sb%klat(i,i) /= M_ZERO) factor(i) = sb%klat(i,i)
      end do

      select case(mode)
      case(1)
        ! output bands in gnuplot format
        do j = 1, nst
          do ik = 1, st%d%nik, ns
            write(iunit, '(1x,3(f10.4))', advance='no') &
              st%d%kpoints(:,ik)/factor(:)
            write(iunit, '(3x,f12.6))', advance='yes') st%eigenval(j, ik)/units_out%energy%factor
          end do
          write(iunit, '(a)')' '
        end do
      case(2)
        ! output bands in xmgrace format, i.e.:
        ! k_x, k_y, k_z, e_1, e_2, ..., e_n
        do ik = 1, st%d%nik, ns
          write(iunit, '(1x,3(f10.4))', advance='no') &
            st%d%kpoints(:,ik)/factor(:)
          write(iunit, '(3x,20f12.6))', advance='yes') (st%eigenval(j, ik)/units_out%energy%factor, j=1,nst)
        end do
      end select

#ifdef HAVE_MPI
    end if
#endif

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
    type(mesh_type),   intent(in)  :: m
    type(states_type), intent(in)  :: st
    type(states_type), intent(in)  :: gs_st
    CMPLX,             intent(out) :: p(st%nst, gs_st%nst, st%d%nik)

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
  subroutine states_magnetization_dens(st, np, rho, m)
    type(states_type), intent(in)  :: st
    integer,           intent(in)  :: np
    FLOAT,             intent(in)  :: rho(np, st%d%nspin)
    FLOAT,             intent(out) :: m(np, 3)

    call push_sub('states.states_magnetization_dens')

    select case (st%d%ispin)
    case (UNPOLARIZED)
      m = M_ZERO
    case (SPIN_POLARIZED)
      m = M_ZERO
      m(:, 3) = rho(:, 1) - rho(:, 2)
    case (SPINORS)
      m(:, 1) =  M_TWO*rho(:, 3)
      m(:, 2) = -M_TWO*rho(:, 4)
      m(:, 3) = rho(:, 1) - rho(:, 2)
    end select

    call pop_sub()
  end subroutine states_magnetization_dens


  ! ---------------------------------------------------------
  subroutine states_magnetic_moment(m, st, rho, mm)
    type(mesh_type),   intent(in)  :: m
    type(states_type), intent(in)  :: st
    FLOAT,             intent(in)  :: rho(m%np, st%d%nspin)
    FLOAT,             intent(out) :: mm(3)

    FLOAT, allocatable :: md(:,:)

    call push_sub('states.states_magnetic_moment')

    allocate(md(m%np, 3))
    call states_magnetization_dens(st, m%np, rho, md)
    mm(1) = dmf_integrate(m, md(:, 1))
    mm(2) = dmf_integrate(m, md(:, 2))
    mm(3) = dmf_integrate(m, md(:, 3))
    deallocate(md)

    call pop_sub()
  end subroutine states_magnetic_moment


  ! ---------------------------------------------------------
  subroutine states_local_magnetic_moments(m, st, geo, rho, r, lmm)
    type(mesh_type),     intent(in)  :: m
    type(states_type),   intent(in)  :: st
    type(geometry_type), intent(in)  :: geo
    FLOAT,               intent(in)  :: rho(m%np, st%d%nspin)
    FLOAT,               intent(in)  :: r
    FLOAT,               intent(out) :: lmm(3, geo%natoms)

    integer :: ia, i
    FLOAT :: ri
    FLOAT, allocatable :: md(:,:)

    call push_sub('states.states_local_magnetic_moments')

    allocate(md(m%np, 3))
    call states_magnetization_dens(st, m%np, rho, md)
    lmm = M_ZERO
    do ia = 1, geo%natoms
      do i = 1, m%np
        call mesh_r(m, i, ri, a=geo%atom(ia)%x)
        if (ri > r) cycle
        lmm(:, ia) = lmm(:, ia) + md(i, :)*m%vol_pp(i)
      end do
    end do
    deallocate(md)

    call pop_sub()
  end subroutine states_local_magnetic_moments


  ! ---------------------------------------------------------
  ! This routine (obviously) assumes complex wave-functions
  subroutine calc_paramagnetic_current(gr, st, jp)
    type(grid_type),   intent(inout) :: gr
    type(states_type), intent(in)    :: st
    FLOAT,             intent(out)   :: jp(:,:,:)  ! (NP, NDIM, st%d%nspin)

    integer :: ik, p, sp, k
    CMPLX, allocatable :: grad(:,:)
#if defined(HAVE_MPI)
    integer :: ierr
    FLOAT, allocatable :: red(:,:,:)
#endif

    call push_sub('states.calc_paramagnetic_current')

    if(st%d%ispin == SPIN_POLARIZED) then
      sp = 2
    else
      sp = 1
    end if

    jp = M_ZERO
    allocate(grad(NP, NDIM))

    do ik = 1, st%d%nik, sp
      do p  = st%st_start, st%st_end
        call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, 1, p, ik), grad)

        ! spin-up density
        do k = 1, NDIM
          jp(:, k, 1) = jp(:, k, 1) + st%d%kweights(ik)*st%occ(p, ik)  &
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
            jp(:, k, 2) = jp(:, k, 2) + st%d%kweights(ik)*st%occ(p, ik) &
              * aimag(conjg(st%zpsi(:, 2, p, ik)) * grad(:, k))
          end do
        end if

      end do
    end do
    deallocate(grad)

#if defined(HAVE_MPI)
    allocate(red(NP, NDIM, st%d%nspin))
    call MPI_ALLREDUCE(jp(1, 1, 1), red(1, 1, 1), NP*NDIM*st%d%nspin, &
      MPI_FLOAT, MPI_SUM, st%comm, ierr)
    jp = red
    deallocate(red)
#endif

    call pop_sub()
  end subroutine calc_paramagnetic_current


  ! ---------------------------------------------------------
  subroutine states_calc_physical_current(gr, st, j)
    type(grid_type),     intent(inout) :: gr
    type(states_type),   intent(in)    :: st
    FLOAT,               intent(out)   :: j(:,:,:)   ! j(NP, NDIM, st%d%nspin)

    call push_sub('states.states_calc_physical_current')

    ! Paramagnetic contribution to the physical current
    call calc_paramagnetic_current(gr, st, j)

    ! TODO
    ! Diamagnetic contribution to the physical current

    call pop_sub()
  end subroutine states_calc_physical_current


  ! ---------------------------------------------------------
  subroutine states_distribute_nodes(st, mc)
    type(states_type),    intent(inout) :: st
    type(multicomm_type), intent(in)    :: mc

#if defined(HAVE_MPI)
    integer :: sn, sn1, r, j, k, ierr, i, ii
#endif

    ! defaults
    st%node(:)  = 0
    st%numprocs = mc%n_node ! 1 
    st%rank     = 0
    st%st_start = 1
    st%st_end   = st%nst
    st%parallel_in_states = .false.

#if defined(HAVE_MPI)
    if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then

      st%parallel_in_states = .true.
      st%comm = mc%group_comm(P_STRATEGY_STATES)
      call MPI_Comm_size(st%comm, st%numprocs, ierr)
      call MPI_Comm_rank(st%comm, st%rank, ierr)

      if(st%nst < st%numprocs) then
        message(1) = "Have more processors than necessary"
        write(message(2),'(i4,a,i4,a)') st%numprocs, " processors and ", st%nst, " states."
        call write_fatal(2)
      end if

      i = st%nst / st%numprocs
      ii = st%nst - i*st%numprocs
      if(ii > 0 .and. st%rank < ii) then
        i = i + 1
        st%st_start = st%rank*i + 1
        st%st_end = st%st_start + i - 1
      else
        st%st_end = st%nst - (st%numprocs - st%rank - 1)*i
        st%st_start = st%st_end - i + 1
      end if
      call MPI_Barrier(st%comm, i)
      write(stdout, '(a,i4,a,i4,a,i4)') "Info: Node ", st%rank, " will propagate state ", &
        st%st_start, " - ", st%st_end
      call MPI_Barrier(st%comm, i)

      sn  = st%nst/st%numprocs
      sn1 = sn + 1
      r  = mod(st%nst, st%numprocs)
      do j = 1, r
        st%node((j-1)*sn1+1:j*sn1) = j - 1
      end do
      k = sn1*r
      call MPI_BARRIER(st%comm, ierr)
      do j = 1, st%numprocs - r
        st%node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
      end do
    end if
#endif

  end subroutine states_distribute_nodes

#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
