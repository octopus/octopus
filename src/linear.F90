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

#include "config_F90.h"

module linear
  use global
  use io
  use mesh
  use hartree
  use states

  implicit none

contains
  subroutine calc_petersilka(type, st, m, hart, n_occ, n_unocc, flags, dir, fname)
    type(states_type), intent(IN) :: st
    type(mesh_type), intent(inout) :: m
    type(hartree_type), intent(inout) :: hart
    integer, intent(IN) :: type, n_occ, n_unocc, flags(32)
    character(len=*), intent(IN) :: dir, fname

    integer, allocatable :: pair_i(:), pair_a(:)
    real(r8), allocatable :: energies(:,:)
    integer :: iunit, n_pairs, i, a, ia

    ! output
    call oct_mkdir(C_string(trim(dir)))

    ! get occupied/unoccupied pairs
    call fix_pairs()
    
    if(type == 0.or.type == 1) then ! eigenvalues or petersilka formula
      call oct_progress_bar(-1, n_pairs) ! initialize bar

      do ia = 1, n_pairs
        a = pair_a(ia)
        i = pair_i(ia)
        energies(ia, 1) = st%eigenval(a, 1) - st%eigenval(i, 1)
        if(type == 1) then 
          energies(ia, 1) = energies(ia, 1) + 2._r8*K_term(i, a, i, a)
        end if

        ! oscilator strengths?
        call matrix_elem(i, a, energies(ia, 2:4))
        energies(ia, 2:4) = 2._r8 * (energies(ia, 2:4))**2 * &
             (st%eigenval(a, 1) - st%eigenval(i, 1))

        call oct_progress_bar(ia-1, n_pairs-1)
      end do

    else if(type == 2) then ! solve full matrix
      call solve_matrix()
    else
      call sort_energies()
    end if
    write(*, "(1x)")

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname), status='unknown')

    if(type == 0) write(iunit, '(2a4)', advance='no') 'From', ' To '
    write(iunit, '(5(a15,1x))') 'E' , '<x>', '<y>', '<z>', '<f>'
    do ia = 1, n_pairs
      if(type == 0.or.type == 1) then
        write(iunit, '(2i4)', advance='no') pair_i(ia), pair_a(ia)
      end if
      write(iunit, '(5(e15.8,1x))') energies(ia,1) / units_out%energy%factor, &
           energies(ia, 2:4), 2._r8/3._r8*sum(energies(ia, 2:4))
    end do
    call io_close(iunit)

    ! clean up
    deallocate(pair_i, pair_a, energies)

  contains
    
    subroutine solve_matrix()
      real(r8), allocatable :: mat(:,:), w(:), work(:), os(:,:)
      real(r8) :: temp
      integer :: ia, jb, i, j, a, b, lwork, info
      integer :: max, actual, iunit

      allocate(mat(n_pairs, n_pairs))
      mat = 0._r8

      max = n_pairs*(1 + n_pairs)/2 - 1
      actual = 0
      call oct_progress_bar(-1, max)
      do ia = 1, n_pairs
        i = pair_i(ia)
        a = pair_a(ia)
        do jb = ia, n_pairs
          j = pair_i(jb)
          b = pair_a(jb)
          mat(ia, jb) = 4._r8 * K_term(i, a, j, b) &
                        * sqrt(st%eigenval(b, 1) - st%eigenval(j, 1))

          if(jb /= ia) mat(jb, ia) = mat(ia, jb) ! the matrix is symmetric

          actual = actual + 1
          call oct_progress_bar(actual, max)
        end do

        temp = st%eigenval(a, 1) - st%eigenval(i, 1)
        mat(ia, :)  = sqrt(temp)*mat(ia, :)
        mat(ia, ia) = temp**2 + mat(ia, ia)
      end do
      write(stdout, '(1x)')

      ! now we diagonalise the matrix using LAPACK
      lwork = 3*n_pairs - 1
      allocate(work(lwork), w(n_pairs))
      call dsyev ('v', 'u', n_pairs, mat, n_pairs, w, work, lwork, info)
      if(info.ne.0) then
        write(message(1),'(a,i5)') 'LAPACK "zheev/dsyev" returned error code ', info
        call write_fatal(1)
      endif

      energies(:, 1) = sqrt(w(:))
      deallocate(work, w)

      ! let us get now the oscillator strengths
      allocate(os(n_pairs, 3))
      do ia = 1, n_pairs
        i = pair_i(ia)
        a = pair_a(ia)
        call matrix_elem(i, a, os(ia,:))
      end do

      do ia = 1, n_pairs
        do j = 1, 3
          energies(ia, 1+j) = 2._r8 * (sum(os(:,j)*mat(:,ia)        &
               *sqrt(st%eigenval(pair_a(:), 1) - st%eigenval(pair_i(:), 1)) ))**2 
        end do
      end do
      
      ! output eigenvectors
      call io_assign(iunit)
      open(iunit, file=trim(dir)+"/"+trim(fname)+".vec", status='unknown')
      write(iunit, '(a14)', advance = 'no') ' value '
      do ia = 1, n_pairs
        write(iunit, '(3x,i4,a1,i4,2x)', advance='no') pair_i(ia), ' - ', pair_a(ia)
      end do
      write(iunit, '(1x)')

      do ia = 1, n_pairs
        write(iunit, '(es14.6)', advance='no') energies(ia, 1) / units_out%energy%factor
        temp = 1._r8
        if(maxval(mat(:, ia)) < abs(minval(mat(:, ia)))) temp = -temp
        do j = 1, n_pairs
          write(iunit, '(es14.6)', advance='no') temp*mat(j, ia)
        end do
        write(iunit, '(1x)')
      end do

      call io_close(iunit)

      deallocate(mat)
    end subroutine solve_matrix

    function K_term(i, a, j, b)
      real(r8) :: K_term
      integer, intent(in) :: i, j, a, b
    
      integer :: ik
      real(r8) :: fxc
      real(r8), allocatable :: rho_i(:,:), rho_j(:,:), pot(:)
      allocate(rho_i(m%np, 1), rho_j(m%np, 1), pot(m%np))
    
      rho_i(:, 1) =  st%R_FUNC(psi) (1:m%np, 1, i, 1) * st%R_FUNC(psi) (1:m%np, 1, a, 1)
      rho_j(:, 1) =  st%R_FUNC(psi) (1:m%np, 1, j, 1) * st%R_FUNC(psi) (1:m%np, 1, b, 1)
    
     !  first the Hartree part (only works for real wfs...)
      pot = 0._r8
      call hartree_solve(hart, m, pot, rho_j(:, 1:1))
      K_term = sum(rho_i(:,1)*pot(:))*m%vol_pp
      
      ! now we have fxc
      do ik = 1, m%np
        call fxc_LDA(st%rho(ik, 1), fxc)
        K_term = K_term + rho_i(ik, 1)*rho_j(ik, 1)*fxc*m%vol_pp
      end do
      
      deallocate(rho_i, rho_j, pot)
    end function K_term

    subroutine fix_pairs()
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
      allocate(energies(n_pairs, 4)) ! excitations + intensities
      energies = 0._r8
      
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

    end subroutine fix_pairs

    subroutine sort_energies
      real(r8) :: tmp(4), emin
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

    subroutine matrix_elem(i, j, s)
      integer, intent(in) :: i, j
      real(r8), intent(out) :: s(3)

      real(r8) :: x(3)
      integer :: k

      s = 0._r8
      do k = 1, m%np
        call mesh_xyz(m, k, x)
        s = s + x * R_CONJ(st%R_FUNC(psi) (k, 1, i, 1)) * st%R_FUNC(psi) (k, 1, j, 1)
      end do

      s = s*m%vol_pp
      
    end subroutine matrix_elem

  end subroutine calc_petersilka


  ! WARNING This should be very improved...
  subroutine fxc_LDA(n, fxc)
    real(r8), intent(in) :: n
    real(r8), intent(out) :: fxc

    real(r8), parameter :: &
         ZERO=0.0_r8, ONE=1.0_r8, HALF=.5_r8, OPF=1.5_r8, C014=0.014_r8, &
         TRD = ONE/3.0_r8, FTRD = 4.0_r8*TRD, ALP = 2.0_r8 * TRD, &
         TFTM = 0.519842099789746380_r8, A0   = 0.521061761197848080_r8, &
         CRS  = 0.6203504908994000870_r8, CXP  = (- 3.0_r8) * ALP / (M_PI*A0), &
         CXF  = 1.259921049894873190_r8

    real(r8), parameter :: &
         C0311=0.03110_r8, C0014=0.00140_r8, &
         C0538=0.05380_r8, C0096=0.00960_r8, C096=0.0960_r8, &
         C0622=0.06220_r8, C004=0.0040_r8, C0232=0.02320_r8, &
         C1686=0.16860_r8, C1P398=1.39810_r8, C2611=0.26110_r8, &
         C2846=0.28460_r8, C1P053=1.05290_r8, C3334=0.33340_r8

    !    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
    real(r8), parameter :: &
         CON2=0.0080_r8/3, CON3=0.35020_r8/3, &
         CON4=0.05040_r8/3, CON5=0.00280_r8/3, CON6=0.19250_r8/3, &
         CON7=0.02060_r8/3, CON8=9.78670_r8/6, CON9=1.0444_r8/3, &
         CON10=7.37030_r8/6, CON11=1.33360_r8/3

    real(r8) :: rs, sqrs, rslog, te, be, dte, dbe, exp, ecp

    fxc = 0._r8

    ! calculate rs
    if(n < 1e-30_r8) then
      fxc = 0._r8
      return
    end if
    rs = CRS / n**TRD

    ! first the exchange part
    fxc = - CXP/(rs*rs)

    ! now PZ correlation
    if (rs .gt. ONE) then
      sqrs = sqrt(rs)
      te   = ONE + CON10*sqrs  + CON11*rs
      dte  = CON10/(2._r8*sqrs)  + CON11
      be   = ONE + C1P053*sqrs + C3334*rs
      dbe  = C1P053/(2._r8*sqrs) + C3334
      ecp  = -(C2846/be)
      fxc  = fxc + ecp/be**2 * (dte*be - 2._r8*dbe*te)
    else
      rslog = log(rs)
      fxc  = fxc + C0622/rs + CON2*(1._r8 + rslog) - CON4
    end if

    fxc = - fxc*rs/(n*3._r8) ! missing factor d rs/d n
    fxc =   fxc / 2._r8      ! Rydbergs -> Hartree

  end subroutine fxc_LDA
end module linear

