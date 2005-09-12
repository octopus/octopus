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

module td_exp
  use global
  use messages
  use lib_oct
  use lib_oct_parser
  use blas
  use lib_basic_alg
  use math
  use mesh
#ifdef HAVE_FFT
  use cube_function
#endif
  use td_exp_split

  implicit none

  private
  public :: td_exp_type, &
       td_exp_init,      &
       td_exp_end,       &
       td_exp_dt


  type td_exp_type
     integer  :: exp_method      ! which method is used to apply the exponential

     FLOAT :: lanczos_tol     ! tolerance for the lanczos method
     integer  :: exp_order       ! order to which the propagator is expanded

     type(zcf) :: cf             ! auxiliary cube for split operator methods
  end type td_exp_type

  integer, parameter :: &
       SPLIT_OPERATOR     = 0, &
       SUZUKI_TROTTER     = 1, &
       LANCZOS_EXPANSION  = 2, &
       FOURTH_ORDER       = 3, &
       CHEBYSHEV          = 4
contains

  subroutine td_exp_init(gr, te)
    type(grid_type),   intent(in)  :: gr
    type(td_exp_type), intent(out) :: te

    call loct_parse_int(check_inp('TDExponentialMethod'), FOURTH_ORDER, te%exp_method)
    select case(te%exp_method)
    case(FOURTH_ORDER)
       message(1) = 'Info: Exponential method: 4th order expansion.'

    case(CHEBYSHEV)
       message(1) = 'Info: Exponential method: Chebyshev.'

    case(LANCZOS_EXPANSION)
       call loct_parse_float(check_inp('TDLanczosTol'), CNST(1e-5), te%lanczos_tol)
       if (te%lanczos_tol <= M_ZERO) then
          write(message(1),'(a,f14.6,a)') "Input: '", te%lanczos_tol, "' is not a valid TDLanczosTol"
          message(2) = '(0 < TDLanczosTol)'
          call write_fatal(2)
       end if

       message(1) = 'Info: Exponential method: Lanczos subspace approximation.'

#if defined(HAVE_FFT)
    case(SPLIT_OPERATOR)
       message(1) = 'Info: Exponential method: Split-Operator.'
    case(SUZUKI_TROTTER)
       message(1) = 'Info: Exponential method: Suzuki-Trotter.'
#endif

    case default
       write(message(1), '(a,i6,a)') "Input: '", te%exp_method, "' is not a valid TDEvolutionMethod"
       message(2) = '(1 <= TDExponentialMethod <= 4)'
       call write_fatal(2)
    end select
    call write_info(1)

    if(te%exp_method==FOURTH_ORDER.or.te%exp_method==CHEBYSHEV.or.te%exp_method==LANCZOS_EXPANSION) then
       call loct_parse_int(check_inp('TDExpOrder'), 4, te%exp_order)
       if (te%exp_order < 2) then
          write(message(1), '(a,i6,a)') "Input: '", te%exp_order, "' is not a valid TDExpOrder"
          message(2) = '(2 <= TDExpOrder)'
          call write_fatal(2)
       end if
#if defined(HAVE_FFT)
    else if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
       call zcf_new(gr%m%l, te%cf)
       call zcf_fft_init(te%cf, gr%sb)
#endif
    end if

  end subroutine td_exp_init

  subroutine td_exp_end(te)
    type(td_exp_type), intent(inout) :: te

    if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
       call zcf_free(te%cf)
    end if
  end subroutine td_exp_end

  subroutine td_exp_dt(te, gr, h, zpsi, ik, timestep, t, order, vmagnus)
    type(td_exp_type),      intent(inout)   :: te
    type(grid_type),        intent(inout)   :: gr
    type(hamiltonian_type), intent(in)      :: h
    integer,                intent(in)      :: ik
    CMPLX,                  intent(inout)   :: zpsi(:, :)
    FLOAT,                  intent(in)      :: timestep, t
    integer,      optional, intent(out)     :: order ! For the methods that rely on Hamiltonian-vector
    ! multiplication, the number of these.
    FLOAT,        optional, intent(in)      :: vmagnus(NP, h%d%nspin, 2)

    logical :: apply_magnus

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    call push_sub('td_exp.td_dtexp')

    select case(te%exp_method)
    case(FOURTH_ORDER);      call fourth
    case(LANCZOS_EXPANSION)
       if(h%ab .eq. IMAGINARY_ABSORBING) then
          message(1) = 'Lanczos expansion exponential method not supported for non-hermtian'
          message(2) = 'Hamiltonians (e.g., with imaginary potential absorbing boundaries).'
          call write_fatal(2)
       endif
       call lanczos
#if defined(HAVE_FFT)
    case(SPLIT_OPERATOR);    call split
    case(SUZUKI_TROTTER);    call suzuki
#endif
    case(CHEBYSHEV);         call cheby
    end select

    call pop_sub()
  contains

    subroutine operate(psi, oppsi)
      CMPLX, intent(inout) :: psi(:, :)
      CMPLX, intent(inout) :: oppsi(:, :)

      if(apply_magnus) then
         call zmagnus(h, gr, psi, oppsi, ik, vmagnus)
      else
         call zHpsi(h, gr, psi, oppsi, ik, t)
      endif

    end subroutine operate

    subroutine fourth
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:), hzpsi1(:,:)
      integer :: i

      call push_sub('td_exp.fourth')

      allocate(zpsi1(NP, h%d%dim), hzpsi1(NP, h%d%dim))
      zfact = M_z1
      call lalg_copy(NP, h%d%dim, zpsi(:,:), zpsi1(:,:))
      do i = 1, te%exp_order
         zfact = zfact*(-M_zI*timestep)/i
         call operate(zpsi1, hzpsi1)
         call lalg_axpy(NP, h%d%dim, zfact, hzpsi1(:,:), zpsi(:,:))
         if(i .ne. te%exp_order) call lalg_copy(NP, h%d%dim, hzpsi1(:,:), zpsi1(:,:))
      end do
      deallocate(zpsi1, hzpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine fourth

    subroutine cheby
      ! Calculates the exponential of the hamiltonian through a expansion in Chebyshev polynomials.
      ! For that purposes it uses the closed form of the coefficients[1] and Clenshaw-Gordons[2]
      ! recursive algorithm.
      ! [1] H. Tal-Ezer and R. Kosloff, J. Chem. Phys 81, 3967 (1984).
      ! [2] C. W. Clenshaw, MTAC 9, 118 (1955).
      ! Since I don't have access to MTAC, I copied next Maple algorithm from Dr. F. G. Lether's
      ! (University of Georgia) homepage: (www.math.uga.edu/~fglether):
      ! {twot := t + t; u0 := 0; u1 := 0;
      !  for k from n to 0 by -1 do
      !    u2 := u1; u1 := u0;
      !    u0 := twot*u1 - u2 + c[k];
      !  od;
      !  ChebySum := 0.5*(u0 - u2);}
      integer :: j
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:,:)

      call push_sub('td_exp.cheby')

      allocate(zpsi1(NP, h%d%dim, 0:2))
      zpsi1 = M_z0
      do j = te%exp_order-1, 0, -1
         call lalg_copy(NP, h%d%dim, zpsi1(:,:, 1), zpsi1(:,:, 2))
         call lalg_copy(NP, h%d%dim, zpsi1(:,:, 0), zpsi1(:,:, 1))

         call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
         zfact = 2*(-M_zI)**j*loct_bessel(j, h%spectral_half_span*timestep)

         call lalg_axpy(NP, h%d%dim, cmplx(-h%spectral_middle_point, M_ZERO, PRECISION), &
              zpsi1(:,:, 1), zpsi1(:,:, 0))
         call lalg_scal(NP, h%d%dim, cmplx(M_TWO/h%spectral_half_span, M_ZERO, PRECISION), zpsi1(:,:, 0))
         call lalg_axpy(NP, h%d%dim, zfact, zpsi(:,:), zpsi1(:,:, 0))
         call lalg_axpy(NP, h%d%dim, cmplx(-M_ONE, M_ZERO, PRECISION), zpsi1(:,:, 2),  zpsi1(:,:, 0))
      end do

      zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
      call lalg_scal(NP, h%d%dim, exp(-M_zI*h%spectral_middle_point*timestep), zpsi(:,:))
      deallocate(zpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine cheby

    subroutine lanczos
      integer ::  korder, n, k, iflag, lwsp, ns, i, j, iexph
      integer, allocatable :: ipiv(:)
      CMPLX, allocatable :: hm(:,:), v(:,:,:), expo(:,:), wsp(:)
      FLOAT :: beta, res, tol !, nrm

      call push_sub('td_exp.lanczos')

      tol    = te%lanczos_tol
      allocate(v(NP, h%d%dim, te%exp_order+1), &
           hm(te%exp_order+1, te%exp_order+1), &
           expo(te%exp_order+1, te%exp_order+1))
      hm = M_z0; expo = M_z0

      lwsp = 4*(te%exp_order)**2+7
      allocate(wsp(lwsp), ipiv(te%exp_order+1))

      ! Normalize input vector, and put it into v(:, :, 1)
      beta = zstates_nrm2(gr%m, h%d%dim, zpsi)
      v(:, :, 1) = zpsi(:, :)/beta

      ! This is the Lanczos loop...
      do n = 1, te%exp_order
         call operate(v(:,:, n), v(:, :, n+1))
         korder = n

         do k = max(1, n-1), n
            hm(k, n) = zstates_dotp(gr%m, h%d%dim, v(:, :, k), v(:, :, n+1))
            call lalg_axpy(NP, h%d%dim, -hm(k, n), v(:, :, k), v(:, :, n+1))
         enddo
         hm(n+1, n) = zstates_nrm2(gr%m, h%d%dim, v(:, :, n+1))
         call lalg_scal(NP, h%d%dim, M_z1/hm(n+1, n), v(:, :, n+1))
         call zgpadm(6, n, timestep, -M_zI*hm(1:n, 1:n), n, wsp, lwsp, ipiv(1:n), iexph, ns, iflag)
         k = 0
         do i = 1, n
            do j = 1, n
               expo(j, i) = wsp(iexph + k); k = k + 1
            enddo
         enddo

         res = abs(hm(n+1, n)*abs(expo(1, n)))
         if(abs(hm(n+1, n)) < CNST(1.0e-12)) exit ! (very unlikely) happy breakdown!!! Yuppi!!!
         if(n>2 .and. res<tol) exit
      enddo

      if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
         write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
         call write_warning(1)
      endif

      call zgpadm(6, korder+1, timestep, -M_zI*hm(1:korder+1, 1:korder+1), korder+1, wsp, &
           lwsp, ipiv(1:korder+1), iexph, ns, iflag)
      k = 0
      do i = 1, korder+1
         do j = 1, korder+1
            expo(j, i) = wsp( iexph + k); k = k + 1
         enddo
      enddo
      ! zpsi = nrm * V * expo(1:korder, 1) = nrm * V * expo * V^(T) * zpsi
      call lalg_gemv(NP, h%d%dim, korder+1, M_z1*beta, v, expo(1:korder+1, 1), M_z0, zpsi)

      if(present(order)) order = korder
      deallocate(v, hm, expo, ipiv, wsp)
      call pop_sub()
    end subroutine lanczos


#if defined(HAVE_FFT)
    subroutine split
      call push_sub('td_exp.split')

      if(h%gauge == 2) then
         message(1) = 'Split operator does not work well if velocity gauge is used.'
         call write_fatal(1)
      endif

      call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*timestep/M_TWO)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*timestep/M_TWO, .true.)
      call zexp_kinetic(gr, h, zpsi, te%cf, -M_zI*timestep)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*timestep/M_TWO, .false.)
      call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*timestep/M_TWO)

      if(present(order)) order = 0
      call pop_sub()
    end subroutine split

    subroutine suzuki
      FLOAT :: dt(5), p, pp(5)
      integer :: k

      call push_sub('td_exp.suzuki')

      if(h%gauge == 2) then
         message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
         call write_fatal(1)
      endif

      p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
      pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
      dt(1:5) = pp(1:5)*timestep

      do k = 1, 5
         call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
         if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*dt(k)/M_TWO, .true.)
         call zexp_kinetic(gr, h, zpsi, te%cf, -M_zI*dt(k))
         if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*dt(k)/M_TWO, .false.)
         call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
      end do

      if(present(order)) order = 0
      call pop_sub()
    end subroutine suzuki
#endif

  end subroutine td_exp_dt

end module td_exp
