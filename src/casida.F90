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

#include "global.h"

module casida
  use global
  use units
  use lib_oct
  use io
  use lib_adv_alg
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

  integer, parameter :: &
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

  private
  public :: casida_run

contains

  ! ---------------------------------------------------------
  integer function casida_run(sys, h, fromScratch) result(ierr)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(casida_type) ::  cas
    integer :: i, n, err
    character(len=100) :: ch
    logical :: l


    message(1) = 'Info: Starting linear response calculation.'
    call write_info(1)

    ierr = 0
    call init_()

    call X(restart_read) ('tmp/restart_unocc', sys%st, sys%m, err)
    if(err.ne.0) then
      message(1) = "Could not read wave-functions from 'tmp/restart_unocc'"
      call write_warning(1)
      
      ierr = 1
      call end_()
      return
    end if

    ! now we count the number of occupied states
    cas%n_occ = 0
    n         = 0
    do i = 1, sys%st%nst
      if(sys%st%occ(i, 1) == M_ZERO) then
        n = n + 1
      else
        cas%n_occ = cas%n_occ + 1
      end if
    end do
    if(n.ne.cas%n_unocc) then
      message(1) = "Inconsistency between variable 'NumberUnoccStates' and file"
      message(2) = "  'tmp/restart_unocc/occs'"
      call write_fatal(2)
    else
      write(message(1),'(a,i4,a)') "Info: Found",cas%n_occ," occupied states."
      write(message(2),'(a,i4,a)') "Info: Found",cas%n_unocc," unoccupied states."
      call write_info(2)
    endif
    

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    ! which states to take into account
    call loct_parse_string("ExciteStates", "1-1024", ch)
    call loct_wfs_list(ch, cas%wfn_flags)
    write(message(1),'(a,a)') "Info: States that form the basis: ",trim(ch)
    Call Write_info(1)

    ! Initialize structure
    call casida_type_init(cas)

    if(fromScratch) call loct_rm('tmp/restart_casida')

    ! calculate resonances
    call loct_parse_logical("CasEigenvalues", .true., l)
    if(l) then
      message(1) = "Info: Eigenvalue differences"
      call write_info(1)

      cas%type = CASIDA_EPS_DIFF
      call casida_work(sys, h, cas)
      call casida_write(cas, 'eps-diff')
    end if
  

    call loct_parse_logical("CasPetersilka", .true., l)
    if(l) then
      message(1) = "Info: Calculating resonance energies a la Petersilka"
      call write_info(1)

      cas%type = CASIDA_PETERSILKA
      call casida_work(sys, h, cas)
      call casida_write(cas, 'petersilka')
    end if

    call loct_parse_logical("LinCasida", .true., l)
    if(l) then
      message(1) = "Info: Calculating resonance energies a la Casida"
      call write_info(1)
      
      cas%type = CASIDA_CASIDA
      call casida_work(sys, h, cas)
      call casida_write(cas, 'casida')
    end if

    call casida_type_end(cas)
    call end_()
  contains

    ! ---------------------------------------------------------
    subroutine init_()

      call push_sub('casida_run')

      call loct_parse_int("NumberUnoccStates", 5, cas%n_unocc)
      if(cas%n_unocc <= 0) then
        message(1) = "Input: NumberUnoccStates must be > 0"
        call write_fatal(1)
      end if

      ! fix states: THIS IS NOT NICE
      sys%st%nst = sys%st%nst + cas%n_unocc
      sys%st%st_end = sys%st%nst

      deallocate(sys%st%eigenval, sys%st%occ)
      allocate(sys%st%X(psi) (sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
      allocate(sys%st%eigenval(sys%st%nst, sys%st%d%nik), sys%st%occ(sys%st%nst, sys%st%d%nik))
      if(sys%st%d%ispin == SPINORS) then
        allocate(sys%st%mag(sys%st%nst, sys%st%d%nik, 2))
        sys%st%mag = M_ZERO
      end if
      sys%st%eigenval = huge(PRECISION)
      sys%st%occ      = M_ZERO

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      ! if I do not change these values, unocc will get crazy
      sys%st%nst    = sys%st%nst - cas%n_unocc
      sys%st%st_end = sys%st%nst

      call pop_sub()
    end subroutine end_
    
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
  subroutine casida_work(sys, h, cas)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(casida_type),      intent(inout) :: cas

    logical, allocatable :: saved_K(:, :)         ! which matrix elements have been loaded
    type(states_type), pointer :: st
    type(mesh_type),   pointer :: m

    ! sanity checks
    ASSERT(cas%type>=CASIDA_EPS_DIFF.and.cas%type<=CASIDA_CASIDA)

    ! some shortcuts
    st => sys%st
    m  => sys%m

    ! initialize stuff
    allocate(saved_K(cas%n_pairs, cas%n_pairs))
    cas%mat      = M_ZERO
    saved_K      = .false.
    cas%energies = M_ZERO

    ! load saved matrix elements
    call load_saved()
    
    if(cas%type == CASIDA_CASIDA) then 
      call solve_casida()              ! solve casida matrix

    else if (mpiv%node == 0) then      ! this is not yet parallel
      call solve_petersilka()          ! eigenvalues or petersilka formula
      call sort_energies()             ! energies may be out of order
    end if
    
    ! clean up
    deallocate(saved_K)

  contains

    ! ---------------------------------------------------------
    subroutine solve_petersilka
      integer :: ia, a, i, iunit
      FLOAT   :: f

      ! initialize progress bar
      call loct_progress_bar(-1, cas%n_pairs-1)

      ! file to save matrix elements
      iunit = io_open('tmp/restart_casida', action='write', position='append')
      
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
      integer :: max, actual, iunit
#ifdef HAVE_MPI
      integer :: ierr
      FLOAT, allocatable :: mpi_mat(:,:)
#endif

      max = cas%n_pairs*(1 + cas%n_pairs)/2 - 1
      actual = 0
      if (mpiv%node == 0) call loct_progress_bar(-1, max)

      ! calculate the matrix elements of (v + fxc)
      do ia = 1, cas%n_pairs
        i = cas%pair_i(ia)
        a = cas%pair_a(ia)

        do jb = ia, cas%n_pairs
          actual = actual + 1
          if(mod(actual, mpiv%numprocs) .ne. mpiv%node) cycle

          ! if not loaded, then calculate matrix element
          if(.not.saved_K(ia, jb)) then
            j = cas%pair_i(jb)
            b = cas%pair_a(jb)

            cas%mat(ia, jb) = K_term(i, a, j, b)
          end if

          if (mpiv%node == 0) call loct_progress_bar(actual-1, max)
        end do
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
      iunit = io_open('tmp/restart_casida', action='write', position='append')
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
    
      integer :: is, ik
      FLOAT :: ffxc
      FLOAT, allocatable :: rho(:, :), fxc(:,:,:), rho_i(:), rho_j(:), pot(:)

      allocate(rho_i(m%np), rho_j(m%np), pot(m%np))
    
      rho_i(:) =  st%X(psi) (1:m%np, 1, i, 1) * R_CONJ(st%X(psi) (1:m%np, 1, a, 1))
      rho_j(:) =  R_CONJ(st%X(psi) (1:m%np, 1, j, 1)) * st%X(psi) (1:m%np, 1, b, 1)
    
      !  first the Hartree part (only works for real wfs...)
      pot = M_ZERO
      call dpoisson_solve(m, sys%f_der, pot, rho_j)
      K_term = dmf_dotp(m, rho_i(:), pot(:))
      deallocate(pot)
      
      ! now we have fxc
      allocate(rho(m%np, st%d%nspin), fxc(m%np, st%d%nspin, st%d%nspin))
      if(associated(st%rho_core)) then
        do is = 1, st%d%spin_channels
          rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
        enddo
      else
        rho = st%rho
      endif
      call xc_get_fxc(sys%ks%xc, m, rho, st%d%ispin, fxc)

      K_term = K_term + sum(rho_i(:)*rho_j(:)*fxc(:,1,1)*m%vol_pp(:))

      deallocate(rho_i, rho_j, rho, fxc)
    end function K_term


    ! ---------------------------------------------------------
    subroutine load_saved
      integer :: iunit, err
      integer :: ia, jb
      FLOAT   :: val
      
      iunit = io_open('tmp/restart_casida', action='read', status='old', die=.false.)
      err = min(iunit, 0)

      do while(err .eq. 0)
        read(iunit, fmt=*, iostat=err) ia, jb, val
        if(err.eq.0 .and. (ia > 0.and.ia <= cas%n_pairs) .and. (jb > 0.and.jb <= cas%n_pairs)) then
          cas%mat(ia, jb) = val
          saved_K(ia, jb) = .true.
        end if
      end do

      if(iunit > 0) call io_close(iunit)

    end subroutine load_saved


    ! ---------------------------------------------------------
    subroutine sort_energies
      FLOAT :: tmp(4), emin
      integer ia, jb, min, itmp

      ! stupid algorith, but who cares
      do ia = 1, cas%n_pairs
        min = ia
        emin = cas%energies(ia, 1)
        do jb = ia + 1, cas%n_pairs
          if(cas%energies(jb, 1) < emin) then
            emin = cas%energies(jb, 1)
            min = jb
          end if
        end do
        if(min .ne. ia) then
          tmp = cas%energies(ia, :)
          cas%energies(ia, :) = cas%energies(min, :)
          cas%energies(min, :) = tmp
          
          itmp = cas%pair_i(ia); cas%pair_i(ia) = cas%pair_i(min); cas%pair_i(min) = itmp
          itmp = cas%pair_a(ia); cas%pair_a(ia) = cas%pair_a(min); cas%pair_a(min) = itmp
        end if

      end do
    end subroutine sort_energies


    ! ---------------------------------------------------------
    subroutine dipole_matrix_elem(i, j, s)
      integer, intent(in) :: i, j
      FLOAT, intent(out) :: s(3)

      FLOAT :: x(3)
      integer :: k

      s = M_ZERO
      do k = 1, m%np
        call mesh_xyz(m, k, x)
        s = s + x * R_CONJ(st%X(psi) (k, 1, i, 1)) * st%X(psi) (k, 1, j, 1) * m%vol_pp(k)
      end do
      
    end subroutine dipole_matrix_elem

  end subroutine casida_work


  ! ---------------------------------------------------------
  subroutine casida_write(cas, filename)
    type(casida_type), intent(in) :: cas
    character(len=*),  intent(in) :: filename

    integer :: iunit, ia, jb
    FLOAT   :: temp

    ! output excitation energies and oscillator strengths
    call io_mkdir('linear')
    iunit = io_open('linear/'//trim(filename))

    if(cas%type == CASIDA_EPS_DIFF) write(iunit, '(2a4)', advance='no') 'From', ' To '

    write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
    do ia = 1, cas%n_pairs
      if((cas%type==CASIDA_EPS_DIFF).or.(cas%type==CASIDA_PETERSILKA)) then
        write(iunit, '(2i4)', advance='no') cas%pair_i(ia), cas%pair_a(ia)
      end if
      write(iunit, '(5(e15.8,1x))') cas%energies(ia,1) / units_out%energy%factor, &
         cas%energies(ia, 2:4), M_TWOTHIRD*sum(cas%energies(ia, 2:4))
    end do
    call io_close(iunit)

    ! output eigenvectors in casida approach
    if(cas%type.ne.CASIDA_CASIDA) return
      
    iunit = io_open('linear/'//trim(filename)//'.vec')
    write(iunit, '(a14)', advance = 'no') ' value '
    do ia = 1, cas%n_pairs
      write(iunit, '(3x,i4,a1,i4,2x)', advance='no') cas%pair_i(ia), ' - ', cas%pair_a(ia)
    end do
    write(iunit, '(1x)')
      
    do ia = 1, cas%n_pairs
      write(iunit, '(es14.6)', advance='no') cas%energies(ia, 1) / units_out%energy%factor
      temp = M_ONE
      if(maxval(cas%mat(:, ia)) < abs(minval(cas%mat(:, ia)))) temp = -temp
      do jb = 1, cas%n_pairs
        write(iunit, '(es14.6)', advance='no') temp*cas%mat(jb, ia)
      end do
      write(iunit, '(1x)')
    end do
      
    call io_close(iunit)

  end subroutine casida_write

end module casida
