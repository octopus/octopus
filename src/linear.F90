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
  subroutine calc_matrix_elem(st, m, n_occ, n_unocc, flags, dir, fname)
    type(states_type), intent(IN) :: st
    type(mesh_type), intent(IN) :: m
    integer, intent(IN) :: n_occ, n_unocc, flags(32)
    character(len=*), intent(IN) :: dir, fname

    integer :: iunit(3), nst, i, is, j, js, k, file
    real(r8) :: a, s(3), x(3)
    character(len=1) :: name(3)

    ! do not bother with errors
    call oct_mkdir(C_string(trim(dir)))

    name(1) = 'x'; name(2) = 'y'; name(3) = 'z';
    do file = 1, 3
      call io_assign(iunit(file))
      open(iunit(file), file=trim(dir)+"/"+trim(fname)+"."+name(file), status='unknown')

      write(iunit(file), '(a10,1x)', advance='no') 'occ/unocc'
      do i = 1, n_occ
        is = flags((i-1)/32 + 1)
        if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
          write(iunit(file), '(i17,1x)', advance='no') i
        end if
      end do
      write(iunit(file), '(1x)')
    end do

    do j = n_occ+1, n_occ + n_unocc
      js = flags((j-1)/32 + 1)
      if(iand(js, 2**(modulo(j-1, 32))).ne.0) then
        do file = 1, 3
          write(iunit(file), '(i10,1x)', advance='no') j
        end do

        do i = 1, n_occ
          is = flags((i-1)/32 + 1)
          if(iand(is, 2**(modulo(i-1, 32))).ne.0) then

            s = R_TOTYPE(0._r8)
            do k = 1, m%np
              call mesh_xyz(m, k, x)
              s = s + x * R_CONJ(st%R_FUNC(psi) (k, 1, i, 1)) * st%R_FUNC(psi) (k, 1, j, 1)
            end do

            do file = 1, 3
              write(iunit(file), '(e17.10,1x)', advance='no') R_ABS(s(file))**2*m%vol_pp
            end do
          end if
        end do

        do file = 1, 3
          write(iunit(file), '(1x)')
        end do

      end if
    end do

    do file = 1, 3
      close(iunit(file))
    end do

  end subroutine calc_matrix_elem

  subroutine calc_petersilka(st, m, n_occ, n_unocc, flags, dir, fname)
    type(states_type), intent(IN) :: st
    type(mesh_type), intent(inout) :: m
    integer, intent(IN) :: n_occ, n_unocc, flags(32)
    character(len=*), intent(IN) :: dir, fname

    type(hartree_type) :: hart
    integer :: iunit, i, is, j, js, ik
    real(r8) :: k, fxc
    real(r8), allocatable :: rho(:,:), pot(:)

    ! initialize Hartree potential
    call hartree_init(hart, m)
    allocate(rho(m%np, 1), pot(m%np))

    ! do not bother with errors
    call oct_mkdir(C_string(trim(dir)))

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname), status='unknown')

    write(iunit, '(a10,1x)', advance='no') 'occ/unocc'
    do i = 1, n_occ
      is = flags((i-1)/32 + 1)
      if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
        write(iunit, '(i17,1x)', advance='no') i
      end if
    end do
    write(iunit, '(1x)')

    do j = n_occ+1, n_occ + n_unocc
      js = flags((j-1)/32 + 1)
      if(iand(js, 2**(modulo(j-1, 32))).ne.0) then
        write(iunit, '(i10,1x)', advance='no') j

        do i = 1, n_occ
          is = flags((i-1)/32 + 1)
          if(iand(is, 2**(modulo(i-1, 32))).ne.0) then
            print *, i, j
            k = st%eigenval(j, 1) - st%eigenval(i, 1) 

            rho(:, 1) =  R_CONJ(st%R_FUNC(psi) (1:m%np, 1, i, 1)) * st%R_FUNC(psi) (1:m%np, 1, j, 1)

            ! first the Hartree part (only works for real wfs...)
            call hartree_solve(hart, m, pot, rho)
            k = k + sum(rho(:,1)*pot(:))*m%vol_pp

            ! now we have fxc
            do ik = 1, m%np
              call fxc_LDA(st%rho(ik, 1), fxc)
              k = k - rho(ik, 1)**2 * fxc * m%vol_pp
            end do
            
            write(iunit, '(e17.10,1x)', advance='no') k / units_out%energy%factor
          end if
        end do

        write(iunit, '(1x)')
      end if
    end do

    deallocate(rho, pot)
    call hartree_end(hart)

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

    real(r8) :: rs, sqrs, rslog, te, be, dte, dbe, exp, vcp, ecp

    fxc = 0._r8

    ! calculate rs
    if(n < 1e-20_r8) then
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
      be   = ONE + C1P053*sqrs + C3334*rs
      dte  = CON10 /(2._r8*sqrs) + CON11
      dbe  = C1P053/(2._r8*sqrs) + C3334
      ecp  = -(C2846/be)
      vcp  = ecp*te/be
      fxc  = fxc + ecp/be*(dbe + (dte*be - dbe*te)/be)
    else
      rslog = log(rs)
      fxc  = fxc + C0622/rs + CON2*(1 + rslog) - CON4
    end if

    fxc = - fxc*n*rs/3._r8 ! missing factor d rs/d n
    fxc =   fxc / 2._r8    ! Rydbergs -> Hartree

  end subroutine fxc_LDA
end module linear

