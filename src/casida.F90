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

  implicit none

  integer, parameter ::       &
       CASIDA_EPS_DIFF   = 1, &
       CASIDA_PETERSILKA = 2, &
       CASIDA_CASIDA     = 3

  type casida_type
     integer :: type          ! CASIDA_EPS_DIFF | CASIDA_PETERSILKA | CASIDA_CASIDA

     integer :: n_occ         ! number of occupied states
     integer :: n_unocc       ! number of unoccupied states
     integer :: wfn_flags(32) ! flags determining which wfs to take into account

     integer          :: n_pairs         ! number of pairs to take into acount
     integer, pointer :: pair_i(:)       ! holds the separated indices of compund index ia
     integer, pointer :: pair_a(:)
     FLOAT,   pointer :: mat(:,:)        ! general purpose matrix
     FLOAT,   pointer :: energies(:,:)   ! excitation energies and intensities
  end type casida_type

contains

  ! ---------------------------------------------------------
  integer function casida_run(sys, h, fromScratch) result(ierr)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(casida_type) ::  cas
    integer :: i, n, err, kpoints, dim, nst, nocc
    character(len=100) :: ch
    logical :: l

    call push_sub('casida.casida_run')

    ierr = 0
    message(1) = 'Info: Starting linear response calculation.'
    call write_info(1)

    call restart_look (trim(tmpdir)//'restart_gs', sys%gr%m, kpoints, dim, nst, nocc, err)
    if(err.ne.0) then
       message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
       call write_fatal(1)
    end if

    sys%st%nst    = nst  
    sys%st%st_end = nst
    cas%n_occ     = nocc
    cas%n_unocc   = nst - nocc
    write(message(1),'(a,i4,a)') "Info: Found",cas%n_occ," occupied states."
    write(message(2),'(a,i4,a)') "Info: Found",cas%n_unocc," unoccupied states."
    call write_info(2)

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

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    ! which states to take into account
    !%Variable LinearResponseKohnShamStates
    !%Type string
    !%Section 5 External Utilities
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
    !%
    !% FIXME: This variable should go into another section, since it is no longer an
    !% external utility. But it does not go into any other section, so it will wait
    !% until we redo the manual.
    !%End
    call loct_parse_string(check_inp('LinearResponseKohnShamStates'), "1-1024", ch)
    call loct_wfs_list(ch, cas%wfn_flags)
    write(message(1),'(a,a)') "Info: States that form the basis: ",trim(ch)
    Call write_info(1)

    ! Initialize structure
    call casida_type_init(cas)

    if(fromScratch) call loct_rm(trim(tmpdir)//'restart_casida')

    ! First, print the differences between KS eigenvalues (first approximation to the
    ! excitation energies, or rather, to the DOS.
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

    ! if I do not change these values, unocc will get crazy
    sys%st%nst    = sys%st%nst - cas%n_unocc
    sys%st%st_end = sys%st%nst
    call pop_sub()
  end function casida_run


  ! allocates stuff, and constructs the arrays pair_i and pair_j
  subroutine casida_type_init(cas)
    type(casida_type), intent(inout) :: cas

    integer :: i, is, a, as, j

    ! count pairs
    cas%n_pairs = 0
    do a = cas%n_occ+1, cas%n_occ + cas%n_unocc
       as = cas%wfn_flags((a-1)/32 + 1)
       if(iand(as, 2**(modulo(a-1, 32))).ne.0) then
          message(1) = " "
          do i = 1, cas%n_occ
             is = cas%wfn_flags((i-1)/32 + 1)
             if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
                cas%n_pairs = cas%n_pairs + 1
             end if
          end do
       end if
    end do

    if(cas%n_pairs < 1) then
       message(1) = "Error: Maybe there are no unoccupied states?"
       call write_fatal(1)
    end if

    ! allocate stuff
    allocate(cas%pair_i(cas%n_pairs), cas%pair_a(cas%n_pairs))
    allocate(cas%mat(cas%n_pairs, cas%n_pairs), cas%energies(cas%n_pairs, 4))

    ! create pairs
    j = 1
    do a = cas%n_occ+1, cas%n_occ + cas%n_unocc
       as = cas%wfn_flags((a-1)/32 + 1)
       if(iand(as, 2**(modulo(a-1, 32))).ne.0) then
          do i = 1, cas%n_occ
             is = cas%wfn_flags((i-1)/32 + 1)
             if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
                cas%pair_i(j) = i
                cas%pair_a(j) = a
                j = j + 1
             end if
          end do
       end if
    end do

  end subroutine casida_type_init


  ! ---------------------------------------------------------
  subroutine casida_type_end(cas)
    type(casida_type), intent(inout) :: cas

    ASSERT(associated(cas%pair_i))

    deallocate(cas%pair_i, cas%pair_a, cas%mat, cas%energies)
    nullify   (cas%pair_i, cas%pair_a, cas%mat, cas%energies)
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
    integer :: is, j_old, b_old

    ! sanity checks
    ASSERT(cas%type>=CASIDA_EPS_DIFF.and.cas%type<=CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    m  => sys%gr%m

    ! initialize stuff
    allocate(saved_K(cas%n_pairs, cas%n_pairs))
    cas%mat      = M_ZERO
    saved_K      = .false.
    cas%energies = M_ZERO

    ! load saved matrix elements
    call load_saved()

    ! This is to be allocated here, and is used inside K_term.
    allocate(pot(m%np)); j_old = -1; b_old = -1

    ! We calculate here the kernel, since it will be needed later.
    allocate(rho(m%np, st%d%nspin), fxc(m%np, st%d%nspin, st%d%nspin))
    if(associated(st%rho_core)) then
       do is = 1, st%d%spin_channels
          rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
       enddo
    else
       rho = st%rho
    endif
    call xc_get_fxc(sys%ks%xc, m, rho, st%d%ispin, fxc)

    select case(cas%type)
      case(CASIDA_CASIDA)
         call solve_casida()
      case(CASIDA_PETERSILKA)
         call solve_petersilka()
    end select

    ! clean up
    deallocate(rho, fxc, pot, saved_K)

  contains

    ! ---------------------------------------------------------
    subroutine solve_petersilka
      integer :: ia, a, i, iunit
      FLOAT   :: f

      ! initialize progress bar
      call loct_progress_bar(-1, cas%n_pairs-1)

      ! file to save matrix elements
      iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append')

      do ia = 1, cas%n_pairs
         a = cas%pair_a(ia)
         i = cas%pair_i(ia)
         cas%energies(ia, 1) = st%eigenval(a, 1) - st%eigenval(i, 1)

         if(cas%type == CASIDA_PETERSILKA) then
            if(saved_K(ia, ia)) then
               f = cas%mat(ia, ia)
            else
               f = K_term(i, a, i, a)
               write(iunit, *) ia, ia, f
            end if
            cas%energies(ia, 1) = cas%energies(ia, 1) + M_TWO*f
         end if

         ! oscilator strengths?
         call dipole_matrix_elem(i, a, cas%energies(ia, 2:4))
         cas%energies(ia, 2:4) = M_TWO * (cas%energies(ia, 2:4))**2 * &
              (st%eigenval(a, 1) - st%eigenval(i, 1))

         call loct_progress_bar(ia-1, cas%n_pairs-1)
      end do

      ! complete progress bar
      write(*, "(1x)")

      ! close restart file
      call io_close(iunit)

    end subroutine solve_petersilka


    ! ---------------------------------------------------------
    subroutine solve_casida()
      FLOAT, allocatable :: os(:,:)
      FLOAT :: temp
      integer :: ia, jb, i, j, a, b
      integer :: max, actual, iunit, counter
#ifdef HAVE_MPI
      integer :: ierr
      FLOAT, allocatable :: mpi_mat(:,:)
#endif

      max = cas%n_pairs*(1 + cas%n_pairs)/2 - 1
      counter = 0
      actual = 0
      if (mpiv%node == 0) call loct_progress_bar(-1, max)

      ! calculate the matrix elements of (v + fxc)
      do jb = 1, cas%n_pairs
         actual = actual + 1
         if(mod(actual, mpiv%numprocs) .ne. mpiv%node) cycle

         do ia = jb, cas%n_pairs
            counter = counter + 1
            i = cas%pair_i(ia)
            a = cas%pair_a(ia)

            ! if not loaded, then calculate matrix element
            if(.not.saved_K(ia, jb)) then
               j = cas%pair_i(jb)
               b = cas%pair_a(jb)

               cas%mat(ia, jb) = K_term(i, a, j, b)
            endif
            if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric

         end do
         if (mpiv%node == 0) call loct_progress_bar(counter-1, max)
      end do

      ! sum all matrix elements
#ifdef HAVE_MPI
      allocate(mpi_mat(cas%n_pairs, cas%n_pairs))
      call MPI_ALLREDUCE(cas%mat(1,1), mpi_mat(1,1), cas%n_pairs**2, &
           MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
      cas%mat = mpi_mat
      deallocate(mpi_mat)
#endif

      ! all processors with the exception of the first are done
      if (mpiv%node.ne.0) return

      ! complete progress bar
      write(stdout, '(1x)')

      ! complete the matrix and output the restart file
      iunit = io_open(trim(tmpdir)//'restart_casida', action='write', position='append')
      do ia = 1, cas%n_pairs
         i = cas%pair_i(ia)
         a = cas%pair_a(ia)
         temp = st%eigenval(a, 1) - st%eigenval(i, 1)

         do jb = ia, cas%n_pairs
            j = cas%pair_i(jb)
            b = cas%pair_a(jb)

            if(.not.saved_K(ia, jb)) write(iunit, *) ia, jb, cas%mat(ia, jb)

            cas%mat(ia, jb)  = M_FOUR * sqrt(temp) * cas%mat(ia, jb) * &
                 sqrt(st%eigenval(b, 1) - st%eigenval(j, 1))

            if(jb /= ia) cas%mat(jb, ia) = cas%mat(ia, jb) ! the matrix is symmetric
         end do
         cas%mat(ia, ia) = temp**2 + cas%mat(ia, ia)
      end do
      call io_close(iunit)

      ! now we diagonalize the matrix
      call lalg_eigensolve(cas%n_pairs, cas%mat, cas%mat, cas%energies(:, 1))
      cas%energies(:, 1) = sqrt(cas%energies(:, 1))

      ! let us get now the oscillator strengths
      allocate(os(cas%n_pairs, 3))
      do ia = 1, cas%n_pairs
         i = cas%pair_i(ia)
         a = cas%pair_a(ia)
         call dipole_matrix_elem(i, a, os(ia,:))
      end do

      do ia = 1, cas%n_pairs
         do j = 1, 3
            cas%energies(ia, 1+j) = M_TWO * (sum(os(:,j)*cas%mat(:,ia)        &
                 *sqrt(st%eigenval(cas%pair_a(:), 1) - st%eigenval(cas%pair_i(:), 1)) ))**2
         end do
      end do

    end subroutine solve_casida


    ! return the matrix element of <i,a|v + fxc|j,b>
    function K_term(i, a, j, b)
      FLOAT :: K_term
      integer, intent(in) :: i, j, a, b

      FLOAT, allocatable :: rho_i(:), rho_j(:)

      allocate(rho_i(m%np), rho_j(m%np))

      rho_i(:) =  st%X(psi) (1:m%np, 1, i, 1) * R_CONJ(st%X(psi) (1:m%np, 1, a, 1))
      rho_j(:) =  R_CONJ(st%X(psi) (1:m%np, 1, j, 1)) * st%X(psi) (1:m%np, 1, b, 1)

      !  first the Hartree part (only works for real wfs...)
      if( j.ne.j_old  .or.   b.ne.b_old) then
         pot = M_ZERO
         call dpoisson_solve(sys%gr, pot, rho_j)
      endif

      K_term = dmf_dotp(m, rho_i(:), pot(:))
      rho(1:m%np, 1) = rho_i(1:m%np) * rho_j(1:m%np) * fxc(1:m%np, 1, 1)
      K_term = K_term + dmf_integrate(m, rho(:, 1))


      j_old = j; b_old = b

      deallocate(rho_i, rho_j)
    end function K_term


    ! ---------------------------------------------------------
    subroutine load_saved
      integer :: iunit, err
      integer :: ia, jb
      FLOAT   :: val

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

    end subroutine load_saved




    ! ---------------------------------------------------------
    subroutine dipole_matrix_elem(i, j, s)
      integer, intent(in) :: i, j
      FLOAT, intent(out) :: s(3)
      R_TYPE, allocatable :: f(:, :)

      integer :: k

      allocate(f(m%np, 3))
      do k = 1, m%np
         f(k, 1:3) = m%x(k, 1:3) * R_CONJ(st%X(psi) (k, 1, i, 1)) * st%X(psi) (k, 1, j, 1)
      enddo
      s(1) = X(mf_integrate)(m, f(:, 1))
      s(2) = X(mf_integrate)(m, f(:, 2))
      s(3) = X(mf_integrate)(m, f(:, 3))
      deallocate(f)

    end subroutine dipole_matrix_elem

  end subroutine casida_work

  ! ---------------------------------------------------------
  subroutine sort_energies(cas)
    type(casida_type), intent(inout) :: cas
    integer :: i
    integer, allocatable :: ind(:), itmp(:)
    FLOAT, allocatable :: tmp(:)
    allocate(ind(cas%n_pairs), itmp(cas%n_pairs), tmp(cas%n_pairs))

    call sort(cas%energies(:, 1), ind)
    do i = 2, 4
       tmp(:) = cas%energies(:, i) ; cas%energies(:, i) = tmp(ind(:))
    enddo
    itmp(:) = cas%pair_i(:) ; cas%pair_i(:) = itmp(ind(:))
    itmp(:) = cas%pair_a(:) ; cas%pair_a(:) = itmp(ind(:))
    deallocate(ind, itmp, tmp)
  end subroutine sort_energies



  ! ---------------------------------------------------------
  subroutine casida_write(cas, filename)
    type(casida_type), intent(in) :: cas
    character(len=*),  intent(in) :: filename

    integer :: iunit, ia, jb
    FLOAT   :: temp

    type(casida_type) :: casp

    casp%type      = cas%type
    casp%n_occ     = cas%n_occ
    casp%n_unocc   = cas%n_unocc
    casp%wfn_flags = cas%wfn_flags
    casp%n_pairs   = cas%n_pairs
    allocate(casp%pair_i(casp%n_pairs), casp%pair_a(casp%n_pairs), &
             casp%mat(casp%n_pairs, casp%n_pairs), casp%energies(casp%n_pairs, 4))
    casp%pair_i    = cas%pair_i
    casp%pair_a    = cas%pair_a
    casp%mat       = cas%mat
    casp%energies  = cas%energies

    call sort_energies(casp)

    ! output excitation energies and oscillator strengths
    call io_mkdir('linear')
    iunit = io_open('linear/'//trim(filename), action='write')

    if(casp%type == CASIDA_EPS_DIFF) write(iunit, '(2a4)', advance='no') 'From', ' To '

    write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
    do ia = 1, casp%n_pairs
       if((casp%type==CASIDA_EPS_DIFF).or.(casp%type==CASIDA_PETERSILKA)) then
          write(iunit, '(2i4)', advance='no') casp%pair_i(ia), casp%pair_a(ia)
       end if
       write(iunit, '(5(e15.8,1x))') casp%energies(ia,1) / units_out%energy%factor, &
            casp%energies(ia, 2:4), M_TWOTHIRD*sum(casp%energies(ia, 2:4))
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
       write(iunit, '(es14.6)', advance='no') casp%energies(ia, 1) / units_out%energy%factor
       temp = M_ONE
       if(maxval(casp%mat(:, ia)) < abs(minval(casp%mat(:, ia)))) temp = -temp
       do jb = 1, casp%n_pairs
          write(iunit, '(es14.6)', advance='no') temp*casp%mat(jb, ia)
       end do
       write(iunit, '(1x)')
    end do

    call io_close(iunit)
    call casida_type_end(casp)
  end subroutine casida_write

end module casida
