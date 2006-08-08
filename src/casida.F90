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

module casida_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use units_m
  use lib_oct_m
  use io_m
  use lib_adv_alg_m
  use math_m, only : sort
  use mesh_m
  use functions_m
  use mesh_function_m
  use poisson_m
  use states_m
  use xc_m
  use system_m
  use hamiltonian_m
  use restart_m
  use mpi_m
  use multicomm_m
  use output_m
  
  implicit none

  private
  public :: casida_run

  integer, parameter ::    &
    CASIDA_EPS_DIFF   = 1, &
    CASIDA_PETERSILKA = 2, &
    CASIDA_CASIDA     = 3

  type casida_t
    integer :: type          ! CASIDA_EPS_DIFF | CASIDA_PETERSILKA | CASIDA_CASIDA

    integer, pointer  :: n_occ(:)       ! number of occupied states
    integer, pointer  :: n_unocc(:)     ! number of unoccupied states
    character(len=80) :: wfn_list

    integer           :: n_pairs        ! number of pairs to take into acount
    type(states_pair_t), pointer :: pair(:)

    FLOAT,   pointer  :: mat(:,:)       ! general purpose matrix
    FLOAT,   pointer  :: w(:)           ! The excitation energies.
    FLOAT,   pointer  :: tm(:, :)       ! The transition matrix elements (between the many-particle states)
    FLOAT,   pointer  :: f(:)           ! The (dipole) strengths
    FLOAT,   pointer  :: s(:)           ! The diagonal part of the S-matrix

    logical           :: parallel_in_eh_pairs
    type(mpi_grp_t)   :: mpi_grp
  end type casida_t

contains

  ! ---------------------------------------------------------
  subroutine casida_run(sys, h, fromScratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    logical,             intent(inout) :: fromScratch

    type(casida_t) :: cas
    integer :: i, ierr, kpoints, dim, nst, n_filled, n_partially_filled, n_half_filled
    character(len=80) :: trandens

    call push_sub('casida.casida_run')

    message(1) = 'Info: Starting linear response calculation.'
    call write_info(1)

    call restart_look(trim(tmpdir)//'restart_gs', sys%gr%m, kpoints, dim, nst, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    if(sys%st%d%ispin == SPINORS) then
      message(1) = 'Linear response TDDFT ("Casida" mode) is not implemented for spinors-DFT.'
      call write_fatal(1)
    end if

    sys%st%nst    = nst
    sys%st%st_end = nst
    deallocate(sys%st%eigenval, sys%st%occ)

    call states_allocate_wfns(sys%st, sys%gr%m)

    ALLOCATE(sys%st%eigenval(sys%st%nst, sys%st%d%nik), sys%st%nst*sys%st%d%nik)
    ALLOCATE(sys%st%occ(sys%st%nst, sys%st%d%nik), sys%st%nst*sys%st%d%nik)

    if(sys%st%d%ispin == SPINORS) then
      ALLOCATE(sys%st%mag(sys%st%nst, sys%st%d%nik, 2), sys%st%nst*sys%st%d%nik*2)
      sys%st%mag = M_ZERO
    end if
    sys%st%eigenval = huge(PRECISION)
    sys%st%occ      = M_ZERO

    call restart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    ALLOCATE(cas%n_occ(sys%st%d%nspin), sys%st%d%nspin)
    ALLOCATE(cas%n_unocc(sys%st%d%nspin), sys%st%d%nspin)

    cas%n_occ(:) = 0
    do i = 1, sys%st%d%nspin
      call occupied_states(sys%st, i, n_filled, n_partially_filled, n_half_filled)
      cas%n_occ(i) = n_filled + n_partially_filled + n_half_filled
      cas%n_unocc(i) = sys%st%nst - cas%n_occ(i)
    end do

    select case(sys%st%d%ispin)
    case(UNPOLARIZED)
      write(message(1),'(a,i4,a)') "Info: Found",cas%n_occ(1)," occupied states."
      write(message(2),'(a,i4,a)') "Info: Found",cas%n_unocc(1)," unoccupied states."
      call write_info(2)
    case(SPIN_POLARIZED)
      write(message(1),'(a,i4,a)') "Info: Found",cas%n_occ(1)," occupied states with spin up."
      write(message(2),'(a,i4,a)') "Info: Found",cas%n_unocc(1)," unoccupied states with spin up."
      write(message(3),'(a,i4,a)') "Info: Found",cas%n_occ(2)," occupied states with spin down."
      write(message(4),'(a,i4,a)') "Info: Found",cas%n_unocc(2)," unoccupied states with spin down."
      call write_info(4)
    end select


    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call system_h_setup(sys, h)

    !%Variable LinearResponseKohnShamStates
    !%Type string
    !%Section Linear Response::Casida
    !%Description
    !% The calculation of the excitation spectrum of a system in the frequency-domain
    !% formulation of linear-response time-dependent density functional theory (TDDFT)
    !% implies the use of a basis set of occupied/unoccupied Kohn-Sham orbitals. This
    !% basis set should, in principle, include all pairs formed by all occupied states,
    !% and an infinite number of unoccupied states. In practice, one has to truncate this
    !% basis set, selecting a number of occupied and unoccupied states that will form the
    !% pairs. These states are specified with this variable. If there are, say, 10 occupied
    !% states, and one sets this variable to the value "10-18", this means that occupied
    !% states from 10 to 15, and unoccupied states from 16 to 18 will be considered.
    !%
    !% This variable is a string in list form, i.e. expressions such as "1,2-5,8-15" are
    !% valid. You should include a non-null number of unoccupied states and a non-null number
    !% of occupied states.
    !%End
    call loct_parse_string(check_inp('LinearResponseKohnShamStates'), "1-1024", cas%wfn_list)
    write(message(1),'(a,a)') "Info: States that form the basis: ",trim(cas%wfn_list)
    Call write_info(1)

    !%Variable LinearResponseTransitionDensities
    !%Type string
    !%Section Linear Response::Casida
    !%Description
    !% Specifies which transition densities are to be calculated and written down. The
    !% transition density for the many-body state n will be written to a file called
    !% linear/rho0n.
    !% 
    !% By default, no transition density is calculated. 
    !%
    !% This variable is a string in list form, i.e. expressions such as "1,2-5,8-15" are
    !% valid.
    !%End
    call loct_parse_string(check_inp('LinearResponseTransitionDensities'), "0", trandens)

    ! Initialize structure
    call casida_type_init(cas, sys%gr%sb%dim, sys%st%d%nspin, sys%mc)

    if(fromScratch) call loct_rm(trim(tmpdir)//'restart_casida')

    ! First, print the differences between KS eigenvalues (first approximation to the
    ! excitation energies, or rather, to the DOS.
    message(1) = "Info: Approximating resonance energies through KS eigenvalue differences"
    call write_info(1)
    cas%type = CASIDA_EPS_DIFF
    call casida_work(sys, h, cas)
    call casida_write(cas, 'eps-diff')

    ! Then, calculate the excitation energies by making use of the Petersilka approximation
    message(1) = "Info: Calculating resonance energies a la Petersilka"
    call write_info(1)
    cas%type = CASIDA_PETERSILKA
    call casida_work(sys, h, cas)
    call casida_write(cas, 'petersilka')

    ! And finally, solve the full Casida problem.
    message(1) = "Info: Calculating resonance energies a la Casida"
    call write_info(1)
    cas%type = CASIDA_CASIDA
    call casida_work(sys, h, cas)
    call casida_write(cas, 'casida')

    ! Calculate and write the transition matrix
    if (sys%st%d%wfs_type == M_REAL) then
      call dget_transition_densities(cas, sys, trandens)
    else
      call zget_transition_densities(cas, sys, trandens)
    end if

    call casida_type_end(cas)

    call pop_sub()
  end subroutine casida_run

  ! ---------------------------------------------------------
  ! allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, dim, nspin, mc)
    type(casida_t), intent(inout) :: cas
    integer, intent(in) :: dim, nspin
    type(multicomm_t), intent(in) :: mc

    integer :: i, a, j, k

    call push_sub('casida.casida_type_init')

    ! count pairs
    cas%n_pairs = 0
    do k = 1, nspin
      do a = cas%n_occ(k)+1, cas%n_occ(k) + cas%n_unocc(k)
        if(loct_isinstringlist(a, cas%wfn_list)) then
          do i = 1, cas%n_occ(k)
            if(loct_isinstringlist(i, cas%wfn_list)) then
              cas%n_pairs = cas%n_pairs + 1
            end if
          end do
        end if
      end do
    end do

    if(cas%n_pairs < 1) then
      message(1) = "Error: Maybe there are no unoccupied states?"
      call write_fatal(1)
    end if

    ! allocate stuff
    ALLOCATE(cas%pair(cas%n_pairs), cas%n_pairs)
    ALLOCATE(cas%mat(cas%n_pairs, cas%n_pairs), cas%n_pairs*cas%n_pairs)
    ALLOCATE(cas%tm(cas%n_pairs, dim), cas%n_pairs*dim)
    ALLOCATE(cas%f(cas%n_pairs), cas%n_pairs)
    ALLOCATE(cas%s(cas%n_pairs), cas%n_pairs)
    ALLOCATE(cas%w(cas%n_pairs), cas%n_pairs)

    ! create pairs
    j = 1
    do k = 1, nspin
      do a = cas%n_occ(k)+1, cas%n_occ(k) + cas%n_unocc(k)
        if(loct_isinstringlist(a, cas%wfn_list)) then
          do i = 1, cas%n_occ(k)
            if(loct_isinstringlist(i, cas%wfn_list)) then
              cas%pair(j)%i = i
              cas%pair(j)%a = a
              cas%pair(j)%sigma = k
              j = j + 1
            end if
          end do
        end if
      end do
    end do

    ! now let us take care of initializing the parallel stuff
    cas%parallel_in_eh_pairs = multicomm_strategy_is_parallel(mc, P_STRATEGY_OTHER)
    if(cas%parallel_in_eh_pairs) then
      call mpi_grp_init(cas%mpi_grp, mc%group_comm(P_STRATEGY_OTHER))
    else
      call mpi_grp_init(cas%mpi_grp, -1)
    end if

    call pop_sub()
  end subroutine casida_type_init


  ! ---------------------------------------------------------
  subroutine casida_type_end(cas)
    type(casida_t), intent(inout) :: cas
    call push_sub('casida.casida_type_end')

    ASSERT(associated(cas%pair))
    deallocate(cas%pair, cas%mat, cas%tm, cas%s, cas%f, cas%w)
    nullify   (cas%pair, cas%mat, cas%tm, cas%s, cas%f, cas%w)

    call pop_sub()
  end subroutine casida_type_end


  ! ---------------------------------------------------------
  ! this subroutine calculates electronic excitation energies using
  ! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine casida_work(sys, h, cas)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(in)    :: h
    type(casida_t),         intent(inout) :: cas

    logical, allocatable :: saved_K(:, :)         ! which matrix elements have been loaded
    type(states_t), pointer :: st
    type(mesh_t),   pointer :: m

    FLOAT, allocatable :: rho(:, :), fxc(:,:,:), pot(:)
    integer :: is, j_old, b_old, mu_old

    call push_sub('casida.casida_work')

    ! sanity checks
    ASSERT(cas%type>=CASIDA_EPS_DIFF.and.cas%type<=CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    m  => sys%gr%m

    ! initialize stuff
    ALLOCATE(saved_K(cas%n_pairs, cas%n_pairs), cas%n_pairs*cas%n_pairs)
    cas%mat = M_ZERO
    saved_K = .false.
    cas%tm  = M_ZERO
    cas%f   = M_ZERO
    cas%w   = M_ZERO
    cas%s   = M_ZERO

    ! load saved matrix elements
    call load_saved()

    ! This is to be allocated here, and is used inside K_term.
    ALLOCATE(pot(m%np), m%np)
    j_old = -1; b_old = -1; mu_old = -1

    ! We calculate here the kernel, since it will be needed later.
    ALLOCATE(rho(m%np, st%d%nspin), m%np*st%d%nspin)
    ALLOCATE(fxc(m%np, st%d%nspin, st%d%nspin), m%np*st%d%nspin*st%d%nspin)
    rho = M_ZERO; fxc = M_ZERO

    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(1:m%np, is) = st%rho(1:m%np, is) + st%rho_core(1:m%np)/st%d%spin_channels
      end do
    else
      rho(1:m%np, 1:st%d%nspin) = st%rho(1:m%np, 1:st%d%nspin)
    end if
    call xc_get_fxc(sys%ks%xc, m, rho, st%d%ispin, fxc)

    select case(cas%type)
    case(CASIDA_EPS_DIFF)
      call solve_petersilka()
    case(CASIDA_PETERSILKA)
      call solve_petersilka()
    case(CASIDA_CASIDA)
      call solve_casida()
    end select

    ! clean up
    deallocate(fxc, rho, pot, saved_K)

    call pop_sub()
  contains

    ! ---------------------------------------------------------
    subroutine solve_petersilka
      integer :: ia, iunit, k
      FLOAT   :: f
      FLOAT, allocatable :: deltav(:), x(:)

      call push_sub('casida.solve_petersilka')

      ! initialize progress bar
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, cas%n_pairs)

      ! file to save matrix elements
      iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append', is_tmp=.true.)

      do ia = 1, cas%n_pairs
        cas%w(ia) = st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                    st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma)

        if(cas%type == CASIDA_PETERSILKA) then
          if(saved_K(ia, ia)) then
            f = cas%mat(ia, ia)
          else
            f = K_term(cas%pair(ia), cas%pair(ia))
            write(iunit, *) ia, ia, f
          end if
          cas%w(ia) = cas%w(ia) + M_TWO*f
        end if

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ia, cas%n_pairs)
      end do

      ALLOCATE(x(cas%n_pairs), cas%n_pairs)
      ALLOCATE(deltav(m%np), m%np)
      do k = 1, m%sb%dim
        
        !WARNING: should x always be real?
        x = dks_matrix_elements(cas, st, m, deltav)

        !FIXME: Calculate the oscillator strengths and matrix elements a la Petersilka
      end do
      deallocate(x, deltav)

      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      ! close restart file
      call io_close(iunit)
      call pop_sub()
    end subroutine solve_petersilka


    ! ---------------------------------------------------------
    subroutine solve_casida()
      FLOAT :: temp
      integer :: ia, jb, k
      integer :: max, actual, iunit, counter
      FLOAT, allocatable :: deltav(:)

      FLOAT, allocatable :: dx(:), tmp(:,:)
      CMPLX, allocatable :: zx(:)
      type(states_pair_t), pointer :: p, q
#ifdef HAVE_MPI
      FLOAT, allocatable :: mpi_mat(:,:)
#endif

      call push_sub('casida.solve_casida')

      max = (cas%n_pairs*(1 + cas%n_pairs)/2)/cas%mpi_grp%size
      counter = 0
      actual = 0
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, max)

      if(.not.mpi_grp_is_root(mpi_world)) cas%mat = M_ZERO

      ! calculate the matrix elements of (v + fxc)
      do jb = 1, cas%n_pairs
        actual = actual + 1
        if(mod(actual, cas%mpi_grp%size) .ne. cas%mpi_grp%rank) cycle
        do ia = jb, cas%n_pairs
          counter = counter + 1
          ! if not loaded, then calculate matrix element
          if(.not.saved_K(ia, jb)) then
            cas%mat(ia, jb) = K_term(cas%pair(ia), cas%pair(jb))
          end if
          if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
        end do
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter, max)
      end do

      ! sum all matrix elements
#ifdef HAVE_MPI
      if(cas%parallel_in_eh_pairs) then
        ALLOCATE(mpi_mat(cas%n_pairs, cas%n_pairs), cas%n_pairs*cas%n_pairs)
        call MPI_Allreduce(cas%mat(1,1), mpi_mat(1,1), cas%n_pairs**2, &
          MPI_FLOAT, MPI_SUM, cas%mpi_grp%comm, mpi_err)
        cas%mat = mpi_mat
        deallocate(mpi_mat)
      end if
#endif
      !if(mpi_grp_is_root(cas%mpi_grp)) print *, "mat =", cas%mat

      ! all processors with the exception of the first are done
      if (mpi_grp_is_root(cas%mpi_grp)) then

        ! complete progress bar
        if(mpi_grp_is_root(mpi_world)) write(stdout, '(1x)')

        ! complete the matrix and output the restart file
        iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append', is_tmp=.true.)
        do ia = 1, cas%n_pairs
          p => cas%pair(ia)
          temp = st%eigenval(p%a, p%sigma) - st%eigenval(p%i, p%sigma)

          do jb = ia, cas%n_pairs
            q => cas%pair(jb)
            if(.not.saved_K(ia, jb)) write(iunit, *) ia, jb, cas%mat(ia, jb)

            if(sys%st%d%ispin == UNPOLARIZED) then
              cas%mat(ia, jb)  = M_FOUR * sqrt(temp) * cas%mat(ia, jb) * &
                sqrt(st%eigenval(q%a, 1) - st%eigenval(q%i, 1))
            else if(sys%st%d%ispin == SPIN_POLARIZED) then
              cas%mat(ia, jb)  = M_TWO * sqrt(temp) * cas%mat(ia, jb) * &
                sqrt(st%eigenval(q%a, q%sigma) - st%eigenval(q%i, q%sigma))
            end if

            if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
          end do
          cas%mat(ia, ia) = temp**2 + cas%mat(ia, ia)
        end do
        call io_close(iunit)
        ALLOCATE(tmp(1:cas%n_pairs,1:cas%n_pairs),cas%n_pairs*cas%n_pairs)
        tmp(1:cas%n_pairs,1:cas%n_pairs) = cas%mat(1:cas%n_pairs,1:cas%n_pairs)
        ! now we diagonalize the matrix
        call lalg_eigensolve(cas%n_pairs, tmp, cas%mat, cas%w)
        DEALLOCATE(tmp)

        do ia = 1, cas%n_pairs
          if(cas%w(ia) < M_ZERO) then
            write(message(1),'(a,i4,a)') 'For whatever reason, the excitation energy',ia,' is negative.'
            write(message(2),'(a)')      'This should not happen.'
            call write_warning(2)
            cas%w(ia) = M_ZERO
          else
            cas%w(ia) = sqrt(cas%w(ia))
          end if
        end do

        ! And let us now get the S matrix...
        do ia = 1, cas%n_pairs
          if(sys%st%d%ispin == UNPOLARIZED) then
            cas%s(ia) = M_HALF / ( st%eigenval(cas%pair(ia)%a, 1) - st%eigenval(cas%pair(ia)%i, 1) )
          elseif(sys%st%d%ispin == SPIN_POLARIZED) then
            cas%s(ia) = M_ONE / ( st%eigenval(cas%pair(ia)%a, cas%pair(ia)%sigma) - &
                                  st%eigenval(cas%pair(ia)%i, cas%pair(ia)%sigma) )
          end if
        end do

        ALLOCATE(deltav(m%np), m%np)
        if (st%d%wfs_type == M_REAL) then
          ALLOCATE(dx(cas%n_pairs), cas%n_pairs)
          do k = 1, m%sb%dim
            deltav(1:m%np) = m%x(1:m%np, k)
            ! let us get now the x vector.
            dx = dks_matrix_elements(cas, st, m, deltav)
            ! And now we are able to get the transition matrix elements between many-electron states.
            do ia = 1, cas%n_pairs
              cas%tm(ia, k) = dtransition_matrix_element(cas, ia, dx)
            end do
          end do
          deallocate(deltav, dx)
        else
          ALLOCATE(zx(cas%n_pairs), cas%n_pairs)
          do k = 1, m%sb%dim
            deltav(1:m%np) = m%x(1:m%np, k)
            ! let us get now the x vector.
            zx = zks_matrix_elements(cas, st, m, deltav)
            ! And now we are able to get the transition matrix elements between many-electron states.
            do ia = 1, cas%n_pairs
              cas%tm(ia, k) = ztransition_matrix_element(cas, ia, zx)
            end do
          end do
          deallocate(deltav, zx)
        end if


        ! And the oscillatory strengths.
        do ia = 1, cas%n_pairs
          cas%f(ia) = (M_TWO/m%sb%dim) * cas%w(ia) * sum( (abs(cas%tm(ia, :)))**2 )
        end do

      end if

#if defined(HAVE_MPI)
      if(cas%parallel_in_eh_pairs) then
        call MPI_Barrier(cas%mpi_grp%comm, mpi_err)
      end if
#endif

      call pop_sub()
    end subroutine solve_casida


    ! return the matrix element of <i(p),a(p)|v + fxc|j(q),b(q)>
    function K_term(p, q)
      FLOAT :: K_term
      type(states_pair_t), intent(in) :: p, q

      integer :: i, j, sigma, a, b, mu
      FLOAT, allocatable :: rho_i(:), rho_j(:)

      i = p%i; a = p%a; sigma = p%sigma
      j = q%i; b = q%a; mu = q%sigma

      ALLOCATE(rho_i(m%np), m%np)
      ALLOCATE(rho_j(m%np), m%np)

      if (st%d%wfs_type == M_REAL) then
        rho_i(:) =  st%dpsi(1:m%np, 1, i, sigma) * st%dpsi(1:m%np, 1, a, sigma)
        rho_j(:) =  st%dpsi(1:m%np, 1, j, mu) * st%dpsi(1:m%np, 1, b, mu)
      else
        rho_i(:) =  st%zpsi(1:m%np, 1, i, sigma) * conjg(st%zpsi(1:m%np, 1, a, sigma))
        rho_j(:) =  conjg(st%zpsi(1:m%np, 1, j, mu)) * st%zpsi(1:m%np, 1, b, mu)
      end if

      !  first the Hartree part (only works for real wfs...)
      if( j.ne.j_old  .or.   b.ne.b_old   .or.  mu.ne.mu_old) then
        pot = M_ZERO
        if( (.not.h%ip_app) ) call dpoisson_solve(sys%gr, pot, rho_j)
      end if

      K_term = dmf_dotp(m, rho_i(:), pot(:))
      rho(1:m%np, 1) = rho_i(1:m%np) * rho_j(1:m%np) * fxc(1:m%np, sigma, mu)
      K_term = K_term + dmf_integrate(m, rho(:, 1))

      j_old = j; b_old = b; mu_old = mu

      deallocate(rho_i, rho_j)
    end function K_term

    ! ---------------------------------------------------------
    subroutine load_saved
      integer :: iunit, err
      integer :: ia, jb
      FLOAT   :: val

      call push_sub('casida.load_saved')

      iunit = io_open(trim(tmpdir)//'restart_casida', action='read', status='old', die=.false., is_tmp=.true.)
      err = min(iunit, 0)

      do
        read(iunit, fmt=*, iostat=err) ia, jb, val
        if(err.ne.0) exit
        if((ia > 0.and.ia <= cas%n_pairs) .and. (jb > 0.and.jb <= cas%n_pairs)) then
          cas%mat(ia, jb) = val
          saved_K(ia, jb) = .true.
          cas%mat(jb, ia) = val
          saved_K(jb, ia) = .true.
        end if
      end do

      if(iunit > 0) call io_close(iunit)
      call pop_sub()
    end subroutine load_saved

  end subroutine casida_work


  ! ---------------------------------------------------------
  subroutine casida_write(cas, filename)
    type(casida_t), intent(in) :: cas
    character(len=*),  intent(in) :: filename

    character(len=5) :: str
    integer :: iunit, ia, jb, dim
    FLOAT   :: temp
    integer, allocatable :: ind(:)
    FLOAT, allocatable :: w(:)

    if(.not.mpi_grp_is_root(cas%mpi_grp)) return

    call push_sub('casida.casida_write')

    dim = size(cas%tm, 2)

    ALLOCATE(w(cas%n_pairs), cas%n_pairs)
    ALLOCATE(ind(cas%n_pairs), cas%n_pairs)
    w = cas%w
    call sort(w, ind)

    ! output excitation energies and oscillator strengths
    call io_mkdir('linear')
    iunit = io_open('linear/'//trim(filename), action='write')

    if(cas%type == CASIDA_EPS_DIFF) write(iunit, '(2a4)', advance='no') 'From', ' To '

    select case(dim)
    case(1); write(iunit, '(3(a15,1x))') 'E' , '<x>', '<f>'
    case(2); write(iunit, '(4(a15,1x))') 'E' , '<x>', '<y>', '<f>'
    case(3); write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
    end select
    do ia = 1, cas%n_pairs
      if((cas%type==CASIDA_EPS_DIFF).or.(cas%type==CASIDA_PETERSILKA)) then
        write(iunit, '(2i4)', advance='no') cas%pair(ind(ia))%i, cas%pair(ind(ia))%a
      end if
      write(iunit, '(5(es15.8,1x))') cas%w(ind(ia)) / units_out%energy%factor, &
        cas%tm(ind(ia), 1:dim) / units_out%length%factor, cas%f(ind(ia))
    end do
    call io_close(iunit)

    ! output eigenvectors in casida approach

    if(cas%type.ne.CASIDA_CASIDA) return

    call io_mkdir('linear/excitations')
    do ia = 1, cas%n_pairs
      write(str,'(i5.5)') ia
      iunit = io_open('linear/excitations/'//trim(str), action='write')
      ! First, a little header
      write(iunit,'(a,es14.5)') '# Energy ['// trim(units_out%energy%abbrev) // '] = ', &
                                cas%w(ind(ia)) / units_out%energy%factor
        write(iunit,'(a,es14.5)') '# <X> ['//trim(units_out%length%abbrev)// '] = ', &
                                  cas%tm(ind(ia),1) / units_out%length%factor
      if(dim > 1) &
        write(iunit,'(a,es14.5)') '# <Y> ['//trim(units_out%length%abbrev)// '] = ', &
                                  cas%tm(ind(ia),2) / units_out%length%factor
      if(dim > 2) &
        write(iunit,'(a,es14.5)') '# <Z> ['//trim(units_out%length%abbrev)// '] = ', &
                                  cas%tm(ind(ia),3) / units_out%length%factor

      temp = M_ONE
      ! I do not know what this does, or what is for.
      !if( maxval(cas%mat(:, ind(ia))) < abs(minval(cas%mat(:, ind(ia)))) ) temp = -temp

      do jb = 1, cas%n_pairs
        write(iunit,*) cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%sigma, temp * cas%mat(jb, ind(ia))
      end do
      call io_close(iunit)
    end do

    deallocate(w, ind)
    call pop_sub()
  end subroutine casida_write


#include "undef.F90"
#include "real.F90"
#include "casida_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "casida_inc.F90"

end module casida_m
