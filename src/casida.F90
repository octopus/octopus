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

module casida
  use global
  use messages
  use syslabels
  use lib_oct_parser
  use units
  use lib_oct
  use io
  use lib_adv_alg
  use math, only : sort
  use mesh
  use functions
  use mesh_function
  use poisson
  use states
  use xc
  use system
  use hamiltonian
  use restart
  use mpi_mod
  use multicomm_mod

  implicit none

  integer, parameter ::    &
    CASIDA_EPS_DIFF   = 1, &
    CASIDA_PETERSILKA = 2, &
    CASIDA_CASIDA     = 3

  type casida_type
    integer :: type          ! CASIDA_EPS_DIFF | CASIDA_PETERSILKA | CASIDA_CASIDA

    integer, pointer  :: n_occ(:)       ! number of occupied states
    integer, pointer  :: n_unocc(:)     ! number of unoccupied states
    character(len=80) :: wfn_list

    integer          :: n_pairs         ! number of pairs to take into acount
    integer, pointer :: pair_i(:)       ! holds the separated indices of compund index ia
    integer, pointer :: pair_a(:)
    integer, pointer :: pair_sigma(:)
    FLOAT,   pointer :: mat(:,:)        ! general purpose matrix
    FLOAT,   pointer :: w(:)            ! The excitation energies.
    FLOAT,   pointer :: tm(:, :)        ! The transition matrix elements (between the many-particle states)
    FLOAT,   pointer :: f(:)            ! The (dipole) strengths
    FLOAT,   pointer :: s(:)            ! The diagonal part of the S-matrix

    logical            :: parallel_in_eh_pairs
    type(mpi_grp_type) :: mpi_grp
  end type casida_type

contains

  ! ---------------------------------------------------------
  integer function casida_run(sys, h, fromScratch) result(ierr)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(casida_type) ::  cas
    integer :: i, err, kpoints, dim, nst, ist

    call push_sub('casida.casida_run')

    ierr = 0
    message(1) = 'Info: Starting linear response calculation.'
    call write_info(1)

    call restart_look (trim(tmpdir)//'restart_gs', sys%gr%m, kpoints, dim, nst, err)
    if(err.ne.0) then
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
    allocate(sys%st%X(psi) (sys%gr%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
    allocate(sys%st%eigenval(sys%st%nst, sys%st%d%nik), sys%st%occ(sys%st%nst, sys%st%d%nik))
    if(sys%st%d%ispin == SPINORS) then
      allocate(sys%st%mag(sys%st%nst, sys%st%d%nik, 2))
      sys%st%mag = M_ZERO
    end if
    sys%st%eigenval = huge(PRECISION)
    sys%st%occ      = M_ZERO

    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr%m, err)
    if(err.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    allocate(cas%n_occ(sys%st%d%nspin), cas%n_unocc(sys%st%d%nspin))
    cas%n_occ(:) = 0
    do i = 1, sys%st%d%nspin
      do ist = 1, sys%st%nst
        if(sys%st%occ(ist, i) > M_ZERO) cas%n_occ(i) = cas%n_occ(i) + 1
      end do
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
    call X(system_h_setup) (sys, h)

    ! which states to take into account
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

    ! Initialize structure
    call casida_type_init(cas, sys%gr%sb%dim, sys%st%d%nspin, sys%mc)

    if(fromScratch) call loct_rm(trim(tmpdir)//'restart_casida')

    ! First, print the differences between KS eigenvalues (first approximation to the
    ! excitation energies, or rather, to the DOS.
    message(1) = "Info: Approximating resonance energies through KS eigenvalue differences"
    call write_info(1)
    cas%type = CASIDA_EPS_DIFF
    call casida_work(sys, cas)
    call casida_write(cas, 'eps-diff')

    ! Then, calculate the excitation energies by making use of the Petersilka approximation
    message(1) = "Info: Calculating resonance energies a la Petersilka"
    call write_info(1)
    cas%type = CASIDA_PETERSILKA
    call casida_work(sys, cas)
    call casida_write(cas, 'petersilka')

    ! And finally, solve the full Casida problem.
    message(1) = "Info: Calculating resonance energies a la Casida"
    call write_info(1)
    cas%type = CASIDA_CASIDA
    call casida_work(sys, cas)
    call casida_write(cas, 'casida')

    call casida_type_end(cas)

    call pop_sub()
  end function casida_run


  ! allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas, dim, nspin, mc)
    type(casida_type), intent(inout) :: cas
    integer, intent(in) :: dim, nspin
    type(multicomm_type), intent(in) :: mc

    integer :: i, a, j, k
#if defined(HAVE_MPI)
    integer :: mpi_err
#endif

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
    allocate(cas%pair_i(cas%n_pairs), cas%pair_a(cas%n_pairs), cas%pair_sigma(cas%n_pairs))
    allocate(cas%mat(cas%n_pairs, cas%n_pairs))
    allocate(cas%tm(cas%n_pairs, dim), cas%f(cas%n_pairs), cas%s(cas%n_pairs), cas%w(cas%n_pairs))

    ! create pairs
    j = 1
    do k = 1, nspin
      do a = cas%n_occ(k)+1, cas%n_occ(k) + cas%n_unocc(k)
        if(loct_isinstringlist(a, cas%wfn_list)) then
          do i = 1, cas%n_occ(k)
            if(loct_isinstringlist(i, cas%wfn_list)) then
              cas%pair_i(j) = i
              cas%pair_a(j) = a
              cas%pair_sigma(j) = k
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
    type(casida_type), intent(inout) :: cas
    call push_sub('casida.casida_type_end')

    ASSERT(associated(cas%pair_i))

    deallocate(cas%pair_i, cas%pair_a, cas%mat, cas%tm, cas%s, cas%f, cas%w)
    nullify   (cas%pair_i, cas%pair_a, cas%mat, cas%tm, cas%s, cas%f, cas%w)

    call pop_sub()
  end subroutine casida_type_end


  ! this subroutine calculates electronic excitation energies using
  ! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine casida_work(sys, cas)
    type(system_type), target, intent(inout) :: sys
    type(casida_type),         intent(inout) :: cas

    logical, allocatable :: saved_K(:, :)         ! which matrix elements have been loaded
    type(states_type), pointer :: st
    type(mesh_type),   pointer :: m

    FLOAT, allocatable :: rho(:, :), fxc(:,:,:), pot(:)
    integer :: is, j_old, b_old, mu_old

    call push_sub('casida.casida_work')

    ! sanity checks
    ASSERT(cas%type>=CASIDA_EPS_DIFF.and.cas%type<=CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    m  => sys%gr%m

    ! initialize stuff
    allocate(saved_K(cas%n_pairs, cas%n_pairs))
    cas%mat = M_ZERO
    saved_K = .false.
    cas%tm  = R_TOTYPE(M_ZERO)
    cas%f   = M_ZERO
    cas%w   = M_ZERO
    cas%s   = M_ZERO

    ! load saved matrix elements
    call load_saved()

    ! This is to be allocated here, and is used inside K_term.
    allocate(pot(m%np)); j_old = -1; b_old = -1; mu_old = -1

    ! We calculate here the kernel, since it will be needed later.
    allocate(rho(m%np, st%d%nspin), fxc(m%np, st%d%nspin, st%d%nspin))
    rho = M_ZERO; fxc = M_ZERO
    if(associated(st%rho_core)) then
      do is = 1, st%d%spin_channels
        rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
      end do
    else
      rho = st%rho
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
    deallocate(rho, fxc, pot, saved_K)

    call pop_sub()
  contains

    ! ---------------------------------------------------------
    subroutine solve_petersilka
      integer :: ia, a, i, iunit, k, sigma
      FLOAT   :: f
      FLOAT, allocatable :: deltav(:), x(:)

      call push_sub('casida.solve_petersilka')

      ! initialize progress bar
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, cas%n_pairs-1)

      ! file to save matrix elements
      iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append')

      do ia = 1, cas%n_pairs
        a = cas%pair_a(ia)
        i = cas%pair_i(ia)
        sigma = cas%pair_sigma(ia)
        cas%w(ia) = st%eigenval(a, sigma) - st%eigenval(i, sigma)

        if(cas%type == CASIDA_PETERSILKA) then
          if(saved_K(ia, ia)) then
            f = cas%mat(ia, ia)
          else
            f = K_term(i, a, sigma, i, a, sigma)
            write(iunit, *) ia, ia, f
          end if
          cas%w(ia) = cas%w(ia) + M_TWO*f
        end if

        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ia-1, cas%n_pairs-1)
      end do

      allocate(x(cas%n_pairs), deltav(m%np))
      do k = 1, m%sb%dim
        x = ks_matrix_elements(cas, st, m, deltav)
        !FIXME: Calculate the oscillator strengths and matrix elements a la Petersilka
      end do
      deallocate(x, deltav)

      ! complete progress bar
      if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")

      ! close restart file
      call io_close(iunit)
      call pop_sub()
    end subroutine solve_petersilka


    ! ---------------------------------------------------------
    subroutine solve_casida()
      FLOAT :: temp
      integer :: ia, jb, i, j, a, b, k, sigma, mu
      integer :: max, actual, iunit, counter
      FLOAT, allocatable :: deltav(:)
      R_TYPE, allocatable :: x(:)
#ifdef HAVE_MPI
      integer :: mpi_err
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
          i = cas%pair_i(ia)
          a = cas%pair_a(ia)
          sigma = cas%pair_sigma(ia)

          ! if not loaded, then calculate matrix element
          if(.not.saved_K(ia, jb)) then
            j = cas%pair_i(jb)
            b = cas%pair_a(jb)
            mu = cas%pair_sigma(jb)
            cas%mat(ia, jb) = K_term(i, a, sigma, j, b, mu)
          end if
          if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric

        end do
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter-1, max)
      end do

      ! sum all matrix elements
#ifdef HAVE_MPI
      if(cas%parallel_in_eh_pairs) then
        allocate(mpi_mat(cas%n_pairs, cas%n_pairs))
        call MPI_ALLREDUCE(cas%mat(1,1), mpi_mat(1,1), cas%n_pairs**2, &
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
        iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append')
        do ia = 1, cas%n_pairs
          i = cas%pair_i(ia)
          a = cas%pair_a(ia)
          sigma = cas%pair_sigma(ia)
          temp = st%eigenval(a, sigma) - st%eigenval(i, sigma)

          do jb = ia, cas%n_pairs
            j = cas%pair_i(jb)
            b = cas%pair_a(jb)
            mu = cas%pair_sigma(jb)

            if(.not.saved_K(ia, jb)) write(iunit, *) ia, jb, cas%mat(ia, jb)

            if(sys%st%d%ispin == UNPOLARIZED) then
              cas%mat(ia, jb)  = M_FOUR * sqrt(temp) * cas%mat(ia, jb) * &
                sqrt(st%eigenval(b, 1) - st%eigenval(j, 1))
            else if(sys%st%d%ispin == SPIN_POLARIZED) then
              cas%mat(ia, jb)  = M_TWO * sqrt(temp) * cas%mat(ia, jb) * &
                sqrt(st%eigenval(b, mu) - st%eigenval(j, mu))
            end if

            if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
          end do
          cas%mat(ia, ia) = temp**2 + cas%mat(ia, ia)
        end do
        call io_close(iunit)

        ! now we diagonalize the matrix
        call lalg_eigensolve(cas%n_pairs, cas%mat, cas%mat, cas%w)
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
          i = cas%pair_i(ia)
          a = cas%pair_a(ia)
          sigma = cas%pair_sigma(ia)
          if(sys%st%d%ispin == UNPOLARIZED) then
            cas%s(ia) = M_HALF / ( st%eigenval(cas%pair_a(ia), 1) - st%eigenval(cas%pair_i(ia), 1) )
          elseif(sys%st%d%ispin == SPIN_POLARIZED) then
            cas%s(ia) = M_ONE / ( st%eigenval(cas%pair_a(ia), sigma) - st%eigenval(cas%pair_i(ia), sigma) )
          end if
        end do

        allocate(deltav(m%np), x(cas%n_pairs))
        do k = 1, m%sb%dim
          deltav(1:m%np) = m%x(1:m%np, k)
          ! let us get now the x vector.
          x = ks_matrix_elements(cas, st, m, deltav)
          ! And now we are able to get the transition matrix elements between many-electron states.
          do ia = 1, cas%n_pairs
            cas%tm(ia, k) = transition_matrix_element(cas, ia, x)
          end do
        end do
        deallocate(deltav, x)

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


    ! return the matrix element of <i,a|v + fxc|j,b>
    function K_term(i, a, sigma, j, b, mu)
      FLOAT :: K_term
      integer, intent(in) :: i, j, sigma, a, b, mu

      FLOAT, allocatable :: rho_i(:), rho_j(:)

      allocate(rho_i(m%np), rho_j(m%np))

      rho_i(:) =  st%X(psi) (1:m%np, 1, i, sigma) * R_CONJ(st%X(psi) (1:m%np, 1, a, sigma))
      rho_j(:) =  R_CONJ(st%X(psi) (1:m%np, 1, j, mu)) * st%X(psi) (1:m%np, 1, b, mu)

      !  first the Hartree part (only works for real wfs...)
      if( j.ne.j_old  .or.   b.ne.b_old   .or.  mu.ne.mu_old) then
        pot = M_ZERO
        call dpoisson_solve(sys%gr, pot, rho_j)
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

      iunit = io_open(trim(tmpdir)//'restart_casida', action='read', status='old', die=.false.)
      err = min(iunit, 0)

      do while(err .eq. 0)
        read(iunit, fmt=*, iostat=err) ia, jb, val
        if(err.eq.0 .and. (ia > 0.and.ia <= cas%n_pairs) .and. (jb > 0.and.jb <= cas%n_pairs)) then
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
  function ks_matrix_elements(cas, st, m, dv) result(x)
    type(casida_type), intent(in) :: cas
    type(states_type), intent(in) :: st
    type(mesh_type),   intent(in) :: m
    FLOAT, intent(in)   :: dv(:)
    FLOAT :: x(cas%n_pairs)

    R_TYPE, allocatable :: f(:)
    integer :: k, ia, i, a, sigma

    allocate(f(m%np))
    do ia = 1, cas%n_pairs
      i     = cas%pair_i(ia)
      a     = cas%pair_a(ia)
      sigma = cas%pair_sigma(ia)
      do k = 1, m%np
        f(k) = dv(k) * R_CONJ(st%X(psi) (k, 1, i, sigma)) * st%X(psi) (k, 1, a, sigma)
      end do
      x(ia) = X(mf_integrate)(m, f)
    end do

    deallocate(f)
  end function ks_matrix_elements


  ! ---------------------------------------------------------
  R_TYPE function transition_matrix_element(cas, ia, x) result(z)
    type(casida_type), intent(in) :: cas
    integer,           intent(in) :: ia
    R_TYPE,            intent(in) :: x(:)

    integer :: jb

    z = R_TOTYPE(M_ZERO)
    if(cas%w(ia) > M_ZERO) then
      do jb = 1, cas%n_pairs
        z = z + x(jb) * (M_ONE/sqrt(cas%s(jb))) * cas%mat(jb, ia)
      end do
      z = (M_ONE/sqrt(cas%w(ia))) * z
    end if

  end function transition_matrix_element


  ! ---------------------------------------------------------
  ! FIXME It needs to re-order also the oscillator strengths, and the transition matrix elements.
  subroutine sort_energies(cas)
    type(casida_type), intent(inout) :: cas
    integer :: i
    integer, allocatable :: ind(:), itmp(:)
    FLOAT, allocatable :: tmp(:)
    R_TYPE, allocatable :: xtmp(:)
    integer :: dim
    allocate(ind(cas%n_pairs), itmp(cas%n_pairs), tmp(cas%n_pairs), xtmp(cas%n_pairs))

    call push_sub('casida.sort_energies')

    call sort(cas%w(1:cas%n_pairs), ind)
    dim = size(cas%tm, 2)
    do i = 1, dim
      xtmp(:) = cas%tm(:, i); cas%tm(:, i) = xtmp(ind(:))
    end do
    tmp(:) = cas%f(:); cas%f(:) = tmp(ind(:))
    itmp(:) = cas%pair_i(:) ; cas%pair_i(:) = itmp(ind(:))
    itmp(:) = cas%pair_a(:) ; cas%pair_a(:) = itmp(ind(:))
    deallocate(ind, itmp, tmp, xtmp)
    call pop_sub()
  end subroutine sort_energies


  ! ---------------------------------------------------------
  subroutine casida_write(cas, filename)
    type(casida_type), intent(in) :: cas
    character(len=*),  intent(in) :: filename

    integer :: iunit, ia, jb, dim
    FLOAT   :: temp

    type(casida_type) :: casp

    if(.not.mpi_grp_is_root(cas%mpi_grp)) return

    call push_sub('casida.casida_write')

    casp%type      = cas%type
    casp%n_occ     = cas%n_occ
    casp%n_unocc   = cas%n_unocc
    casp%wfn_list  = cas%wfn_list
    casp%n_pairs   = cas%n_pairs
    allocate(casp%pair_i(casp%n_pairs), casp%pair_a(casp%n_pairs), casp%mat(casp%n_pairs, casp%n_pairs))
    dim = size(cas%tm, 2)
    allocate(casp%tm(casp%n_pairs, dim), casp%s(casp%n_pairs), casp%f(casp%n_pairs), casp%w(casp%n_pairs))
    casp%tm        = cas%tm
    casp%s         = cas%s
    casp%f         = cas%f
    casp%w         = cas%w
    casp%pair_i    = cas%pair_i
    casp%pair_a    = cas%pair_a
    casp%mat       = cas%mat
    call sort_energies(casp)

    ! output excitation energies and oscillator strengths
    call io_mkdir('linear')
    iunit = io_open('linear/'//trim(filename), action='write')

    if(casp%type == CASIDA_EPS_DIFF) write(iunit, '(2a4)', advance='no') 'From', ' To '

    select case(dim)
    case(1); write(iunit, '(3(a15,1x))') 'E' , '<x>', '<f>'
    case(2); write(iunit, '(4(a15,1x))') 'E' , '<x>', '<y>', '<f>'
    case(3); write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
    end select
    do ia = 1, casp%n_pairs
      if((casp%type==CASIDA_EPS_DIFF).or.(casp%type==CASIDA_PETERSILKA)) then
        write(iunit, '(2i4)', advance='no') casp%pair_i(ia), casp%pair_a(ia)
      end if
      write(iunit, '(5(es15.8,1x))') casp%w(ia) / units_out%energy%factor, &
        casp%w(ia)*abs(casp%tm(ia, 1:dim))**2, casp%f(ia)
    end do
    call io_close(iunit)

    ! output eigenvectors in casida approach
    if(casp%type.ne.CASIDA_CASIDA) return

    iunit = io_open('linear/'//trim(filename)//'.vec', action='write')
    write(iunit, '(a14)', advance = 'no') ' value '
    do ia = 1, casp%n_pairs
      write(iunit, '(3x,i4,a1,i4,2x)', advance='no') casp%pair_i(ia), ' - ', casp%pair_a(ia)
    end do
    write(iunit, '(1x)')

    do ia = 1, casp%n_pairs
      write(iunit, '(es14.6)', advance='no') casp%w(ia) / units_out%energy%factor
      temp = M_ONE
      if(maxval(casp%mat(:, ia)) < abs(minval(casp%mat(:, ia)))) temp = -temp
      do jb = 1, casp%n_pairs
        write(iunit, '(es14.6)', advance='no') temp*casp%mat(jb, ia)
      end do
      write(iunit, '(1x)')
    end do

    call io_close(iunit)
    call casida_type_end(casp)
    call pop_sub()
  end subroutine casida_write

end module casida
