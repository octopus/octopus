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

#define RESTART_FILE "tmp/restart_linear"

module linear
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

  implicit none

contains
  ! this subroutine calculates electronic excitation energies using
  ! the matrix formulation of M. Petersilka, or of M. Casida
  subroutine LR_el_excitations(type, st, m, f_der, n_occ, n_unocc, flags, dir, fname, from_scratch)
    integer,           intent(in)    :: type         ! 1 => eps diff; 2 => petersilka; 3 => casida
    type(states_type), intent(IN)    :: st           ! the wave-functions
    type(mesh_type),   intent(inout) :: m            ! the mesh
    type(f_der_type),  intent(in)    :: f_der
    integer,           intent(in)    :: n_occ, n_unocc ! number of occupied and unoccupied states
    integer,           intent(in)    :: flags(32)    ! which electron-hole pairs to take into account
    character(len=*),  intent(in)    :: dir, fname   ! directory and filename to write
    logical,           intent(in)    :: from_scratch ! should I start from scratch, or use information from restart file

    ! variables global to the subroutine
    integer, allocatable :: pair_i(:), pair_a(:)  ! holds separated indices of compund index ia
    FLOAT,   allocatable :: mat(:,:)              ! general purpose matrix
    FLOAT,   allocatable :: energies(:,:)         ! excitation energies and intensities
    logical, allocatable :: saved_K(:, :)         ! which matrix elements have been loaded
    integer :: n_pairs                            ! number of particle-hole pairs to take into account

    ! sanity checks
    ASSERT(type>=0.and.type<=2)

    ! output
    if (mpiv%node == 0) call loct_mkdir(trim(dir))

    ! get occupied/unoccupied pairs
    call get_pairs()

    ! initialize stuff
    allocate(mat(n_pairs, n_pairs), saved_K(n_pairs, n_pairs), energies(n_pairs, 4))
    mat      = M_ZERO
    saved_K  = .false.
    energies = M_ZERO

    ! load saved matrix elements
    if(.not.from_scratch) call load_saved()
    
    if(type == 2) then 
      call solve_casida()              ! solve casida matrix
    else if (mpiv%node == 0) then      ! this is not yet parallel
      call solve_petersilka()          ! eigenvalues or petersilka formula
      call sort_energies()             ! energies may be out of order
    end if
    
    ! write excitation energies and oscillator strengths to file
    if(mpiv%node == 0) call write_LR()

    ! clean up
    call free_pairs()
    deallocate(mat, saved_K)

  contains

    subroutine solve_petersilka
      integer :: ia, a, i, iunit
      FLOAT   :: f

      ! initialize progress bar
      call loct_progress_bar(-1, n_pairs-1)

      ! file to save matrix elements
      call io_assign(iunit)
      open(iunit, file=RESTART_FILE, action='write', status='unknown', position='append')
      
      do ia = 1, n_pairs
        a = pair_a(ia)
        i = pair_i(ia)
        energies(ia, 1) = st%eigenval(a, 1) - st%eigenval(i, 1)
        if(type == 1) then
          if(saved_K(ia, ia)) then
            f = mat(ia, ia)
          else
            f = K_term(i, a, i, a)
            write(iunit, *) ia, ia, f
          end if
          energies(ia, 1) = energies(ia, 1) + M_TWO*f
        end if

        ! oscilator strengths?
        call dipole_matrix_elem(i, a, energies(ia, 2:4))
        energies(ia, 2:4) = M_TWO * (energies(ia, 2:4))**2 * &
           (st%eigenval(a, 1) - st%eigenval(i, 1))

        call loct_progress_bar(ia-1, n_pairs-1)
      end do

      ! complete progress bar
      write(*, "(1x)")

      ! close restart file
      call io_close(iunit)

    end subroutine solve_petersilka
    
    subroutine solve_casida()
      FLOAT, allocatable :: os(:,:)
      FLOAT :: temp
      integer :: ia, jb, i, j, a, b
      integer :: max, actual, iunit
#ifdef HAVE_MPI
      integer :: ierr
      FLOAT, allocatable :: mpi_mat(:,:)
#endif

      max = n_pairs*(1 + n_pairs)/2 - 1
      actual = 0
      if (mpiv%node == 0) call loct_progress_bar(-1, max)

      ! calculate the matrix elements of (v + fxc)
      do ia = 1, n_pairs
        i = pair_i(ia)
        a = pair_a(ia)

        do jb = ia, n_pairs
          actual = actual + 1
          if(mod(actual, mpiv%numprocs) .ne. mpiv%node) cycle

          ! if not loaded, then calculate matrix element
          if(.not.saved_K(ia, jb)) then
            j = pair_i(jb)
            b = pair_a(jb)

            mat(ia, jb) = K_term(i, a, j, b)
          end if

          if (mpiv%node == 0) call loct_progress_bar(actual-1, max)
        end do
      end do

      ! sum all matrix elements
#ifdef HAVE_MPI
      allocate(mpi_mat(n_pairs, n_pairs))
      call MPI_ALLREDUCE(mat(1,1), mpi_mat(1,1), n_pairs**2, &
         MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)      
      mat = mpi_mat
      deallocate(mpi_mat)
#endif
      
      ! all processors with the exception of the first are done
      if (mpiv%node.ne.0) return

      ! complete progress bar
      write(stdout, '(1x)')

      ! complete the matrix and output the restart file
      call io_assign(iunit)
      open(iunit, file=RESTART_FILE, action='write', status='unknown', position='append')
      do ia = 1, n_pairs
        i = pair_i(ia)
        a = pair_a(ia)
        temp = st%eigenval(a, 1) - st%eigenval(i, 1)
        
        do jb = ia, n_pairs
          j = pair_i(jb)
          b = pair_a(jb)
          
          if(.not.saved_K(ia, jb)) write(iunit, *) ia, jb, mat(ia, jb)

          mat(ia, jb)  = M_FOUR * sqrt(temp) * mat(ia, jb) * sqrt(st%eigenval(b, 1) - st%eigenval(j, 1))
          
          if(jb /= ia) mat(jb, ia) = mat(ia, jb) ! the matrix is symmetric
        end do
        mat(ia, ia) = temp**2 + mat(ia, ia)
      end do
      call io_close(iunit)

      ! now we diagonalize the matrix
      call lalg_eigensolve(n_pairs, mat, mat, energies(:, 1))
      energies(:, 1) = sqrt(energies(:, 1))
      
      ! let us get now the oscillator strengths
      allocate(os(n_pairs, 3))
      do ia = 1, n_pairs
        i = pair_i(ia)
        a = pair_a(ia)
        call dipole_matrix_elem(i, a, os(ia,:))
      end do
      
      do ia = 1, n_pairs
        do j = 1, 3
          energies(ia, 1+j) = M_TWO * (sum(os(:,j)*mat(:,ia)        &
             *sqrt(st%eigenval(pair_a(:), 1) - st%eigenval(pair_i(:), 1)) ))**2 
        end do
      end do
      
    end subroutine solve_casida

    ! return the matrix element of <i,a|v + fxc|j,b>
    function K_term(i, a, j, b)
      FLOAT :: K_term
      integer, intent(in) :: i, j, a, b
    
      integer :: ik
      FLOAT :: fxc
      FLOAT, allocatable :: rho_i(:), rho_j(:), pot(:)

      allocate(rho_i(m%np), rho_j(m%np), pot(m%np))
    
      rho_i(:) =  st%X(psi) (1:m%np, 1, i, 1) * st%X(psi) (1:m%np, 1, a, 1)
      rho_j(:) =  st%X(psi) (1:m%np, 1, j, 1) * st%X(psi) (1:m%np, 1, b, 1)
    
      !  first the Hartree part (only works for real wfs...)
      pot = M_ZERO
      call poisson_solve(m, f_der, pot, rho_j)
      K_term = dmf_dotp(m, rho_i(:), pot(:))
      
      ! now we have fxc
      do ik = 1, m%np
        call fxc_LDA(st%rho(ik, 1), fxc)
        K_term = K_term + rho_i(ik)*rho_j(ik)*fxc*m%vol_pp(ik)
      end do
      
      deallocate(rho_i, rho_j, pot)
    end function K_term

    ! reads flags, and constructs the arrays pair_i and pair_j
    subroutine get_pairs()
      integer :: i, is, a, as, j

      ! count pairs
      n_pairs = 0
      do a = n_occ+1, n_occ + n_unocc
        as = flags((a-1)/32 + 1)
        if(iand(as, 2**(modulo(a-1, 32))).ne.0) then
          do i = 1, n_occ
            is = flags((i-1)/32 + 1)
            if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
              n_pairs = n_pairs + 1
            end if
          end do
        end if
      end do

      if(n_pairs < 1) then
        message(1) = "Error: Maybe there are no unoccupied states?"
        call write_fatal(1)
      end if
      
      ! allocate stuff
      allocate(pair_i(n_pairs), pair_a(n_pairs))
      
      ! create pairs
      j = 1
      do a = n_occ+1, n_occ + n_unocc
        as = flags((a-1)/32 + 1)
        if(iand(as, 2**(modulo(a-1, 32))).ne.0) then
          do i = 1, n_occ
            is = flags((i-1)/32 + 1)
            if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
              pair_i(j) = i
              pair_a(j) = a
              j = j + 1
            end if
          end do
        end if
      end do

    end subroutine get_pairs

    subroutine free_pairs
      deallocate(pair_i, pair_a)
    end subroutine free_pairs

    subroutine load_saved
      integer :: iunit, err
      integer :: ia, jb
      FLOAT   :: val
      
      call io_assign(iunit)
      open(unit=iunit, file=RESTART_FILE, iostat=err, action='read', status='old')
      do while(err .eq. 0)
        read(iunit, fmt=*, iostat=err) ia, jb, val
        if(err.eq.0 .and. (ia > 0.and.ia <= n_pairs) .and. (jb > 0.and.jb <= n_pairs)) then
          mat    (ia, jb) = val
          saved_K(ia, jb) = .true.
        end if
      end do
      call io_close(iunit)

    end subroutine load_saved

    subroutine sort_energies
      FLOAT :: tmp(4), emin
      integer ia, jb, min, itmp

      ! stupid algorith, but who cares
      do ia = 1, n_pairs
        min = ia
        emin = energies(ia, 1)
        do jb = ia + 1, n_pairs
          if(energies(jb, 1) < emin) then
            emin = energies(jb, 1)
            min = jb
          end if
        end do
        if(min .ne. ia) then
          tmp = energies(ia, :)
          energies(ia, :) = energies(min, :)
          energies(min, :) = tmp
          
          itmp = pair_i(ia); pair_i(ia) = pair_i(min); pair_i(min) = itmp
          itmp = pair_a(ia); pair_a(ia) = pair_a(min); pair_a(min) = itmp
        end if

      end do
    end subroutine sort_energies

    subroutine write_LR()
      integer :: iunit, ia, jb
      FLOAT   :: temp

      ! output excitation energies and oscillator strengths
      call io_assign(iunit)
      open(iunit, file=trim(dir) // "/" // trim(fname), status='unknown')

      if(type == 0) write(iunit, '(2a4)', advance='no') 'From', ' To '
      write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
      do ia = 1, n_pairs
        if(type == 0.or.type == 1) then
          write(iunit, '(2i4)', advance='no') pair_i(ia), pair_a(ia)
        end if
        write(iunit, '(5(e15.8,1x))') energies(ia,1) / units_out%energy%factor, &
           energies(ia, 2:4), M_TWOTHIRD*sum(energies(ia, 2:4))
      end do
      call io_close(iunit)

      ! output eigenvectors in casida approach
      if(type.ne.2) return
      
      call io_assign(iunit)
      open(iunit, file=trim(dir) // "/" // trim(fname) // ".vec", status='unknown')
      write(iunit, '(a14)', advance = 'no') ' value '
      do ia = 1, n_pairs
        write(iunit, '(3x,i4,a1,i4,2x)', advance='no') pair_i(ia), ' - ', pair_a(ia)
      end do
      write(iunit, '(1x)')
      
      do ia = 1, n_pairs
        write(iunit, '(es14.6)', advance='no') energies(ia, 1) / units_out%energy%factor
        temp = M_ONE
        if(maxval(mat(:, ia)) < abs(minval(mat(:, ia)))) temp = -temp
        do jb = 1, n_pairs
          write(iunit, '(es14.6)', advance='no') temp*mat(jb, ia)
        end do
        write(iunit, '(1x)')
      end do
      
      call io_close(iunit)

    end subroutine write_LR

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

  end subroutine LR_el_excitations


  ! WARNING This should be very improved...
  subroutine fxc_LDA(n, fxc)
    FLOAT, intent(in) :: n
    FLOAT, intent(out) :: fxc

    FLOAT, parameter :: &
         OPF  = CNST(1.5), C014 = CNST(0.014),                                 &
         TRD  = M_ONE/M_THREE, FTRD = M_FOUR*TRD, ALP = M_TWO*TRD,             &
         TFTM = CNST(0.519842099789746380 ),  A0 = CNST(0.521061761197848080), &
         CRS  = CNST(0.6203504908994000870), CXP = -M_THREE*ALP / (M_PI*A0),   &
         CXF  = CNST(1.259921049894873190)

    FLOAT, parameter :: &
         C0311 = CNST(0.03110), C0014 = CNST(0.00140), &
         C0538 = CNST(0.05380), C0096 = CNST(0.00960), C096  = CNST(0.0960 ), &
         C0622 = CNST(0.06220), C004  = CNST(0.0040 ), C0232 = CNST(0.02320), &
         C1686 = CNST(0.16860), C1P398= CNST(1.39810), C2611 = CNST(0.26110), &
         C2846 = CNST(0.28460), C1P053= CNST(1.05290), C3334 = CNST(0.33340)

    !    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
    FLOAT, parameter :: &
                                      CON2 = CNST(0.0080 )/M_THREE, &
        CON3 = CNST(0.35020)/M_THREE, CON4 = CNST(0.05040)/M_THREE, &
        CON5 = CNST(0.00280)/M_THREE, CON6 = CNST(0.19250)/M_THREE, &
        CON7 = CNST(0.02060)/M_THREE, CON8 = CNST(9.78670)/M_SIX,   &
        CON9 = CNST(1.0444 )/M_THREE, CON10= CNST(7.37030)/M_SIX,   &
        CON11= CNST(1.33360)/M_THREE

    FLOAT :: rs, sqrs, rslog, te, be, dte, dbe, ecp

    ! calculate rs
    if(n < CNST(1e-30)) then
      fxc = M_ZERO
      return
    end if
    rs = CRS / n**TRD

    ! first the exchange part
    fxc = - CXP/(rs*rs)

    ! now PZ correlation
    if (rs .gt. M_ONE) then
      sqrs = sqrt(rs)
      te   = M_ONE + CON10*sqrs  + CON11*rs
      dte  = CON10/(M_TWO*sqrs)  + CON11
      be   = M_ONE + C1P053*sqrs + C3334*rs
      dbe  = C1P053/(M_TWO*sqrs) + C3334
      ecp  = -(C2846/be)
      fxc  = fxc + ecp/be**2 * (dte*be - M_TWO*dbe*te)
    else
      rslog = log(rs)
      fxc  = fxc + C0622/rs + CON2*(M_ONE + rslog) - CON4
    end if

    fxc = - fxc*rs/(n*M_THREE) ! missing factor d rs/d n
    fxc =   fxc / M_TWO        ! Rydbergs -> Hartree

  end subroutine fxc_LDA
end module linear

#undef RESTART_FILE
