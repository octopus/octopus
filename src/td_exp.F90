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

module td_exp
  use global
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

  type td_exp_type
    integer  :: exp_method      ! which method is used to apply the exponential

    FLOAT :: lanczos_tol     ! tolerance for the lanczos method
    integer  :: exp_order       ! order to which the propagator is expanded

    type(zcf) :: cf             ! auxiliary cube for split operator methods
  end type td_exp_type

  integer, parameter, private :: SPLIT_OPERATOR     = 0, &
                                 SUZUKI_TROTTER     = 1, &
                                 LANCZOS_EXPANSION  = 2, &
                                 FOURTH_ORDER       = 3, &
                                 CHEBYSHEV          = 4
contains

  subroutine td_exp_init(m, te)
    type(mesh_type),   intent(IN) :: m
    type(td_exp_type), intent(out) :: te

    call loct_parse_int("TDExponentialMethod", FOURTH_ORDER, te%exp_method)
    select case(te%exp_method)
    case(FOURTH_ORDER)
      message(1) = 'Info: Exponential method: 4th order expansion.'

    case(CHEBYSHEV)
      message(1) = 'Info: Exponential method: Chebyshev.'

    case(LANCZOS_EXPANSION)
      call loct_parse_float("TDLanczosTol", CNST(5e-4), te%lanczos_tol)
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
      call loct_parse_int("TDExpOrder", 4, te%exp_order)
      if (te%exp_order < 2) then
        write(message(1), '(a,i6,a)') "Input: '", te%exp_order, "' is not a valid TDExpOrder"
        message(2) = '(2 <= TDExpOrder)'
        call write_fatal(2)
      end if
#if defined(HAVE_FFT)
    else if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
      call zcf_new(m%l, te%cf)
      call zcf_fft_init(te%cf)
#endif
    end if

  end subroutine td_exp_init

  subroutine td_exp_end(te)
    type(td_exp_type), intent(inout) :: te
    
    if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
      call zcf_free(te%cf)
    end if
  end subroutine td_exp_end

  subroutine td_exp_dt(te, m, h, zpsi, ik, timestep, t, order, vmagnus)
    type(td_exp_type),      intent(inout)   :: te
    type(mesh_type),        intent(IN)      :: m
    type(hamiltonian_type), intent(IN)      :: h
    integer,                intent(in)      :: ik
    CMPLX,                  intent(inout)   :: zpsi(m%np, h%d%dim)
    FLOAT,                  intent(in)      :: timestep, t
    integer,      optional, intent(out)     :: order ! For the methods that rely on Hamiltonian-vector
                                                     ! multiplication, the number of these.
    FLOAT,        optional, intent(IN)      :: vmagnus(m%np, h%d%nspin, 2)

    logical :: apply_magnus

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    call push_sub('td_dtexp')

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
      CMPLX, intent(IN)    :: psi(m%np, h%d%dim)
      CMPLX, intent(inout) :: oppsi(m%np, h%d%dim)

      if(apply_magnus) then
        call zmagnus(h, m, psi, oppsi, ik, vmagnus)
      else
        call zHpsi(h, m, psi, oppsi, ik, t)
      endif

    end subroutine operate

    subroutine fourth
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:), hzpsi1(:,:)
      integer :: i
      
      call push_sub('fourth')
      
      allocate(zpsi1(m%np, h%d%dim), hzpsi1(m%np, h%d%dim))
      zfact = M_z1
      call lalg_copy(m%np, h%d%dim, zpsi(:,:), zpsi1(:,:))
      do i = 1, te%exp_order
        zfact = zfact*(-M_zI*timestep)/i
        call operate(zpsi1, hzpsi1)
        call lalg_axpy(m%np, h%d%dim, zfact, hzpsi1(:,:), zpsi(:,:))
        if(i .ne. te%exp_order) call lalg_copy(m%np, h%d%dim, hzpsi1(:,:), zpsi1(:,:))
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
      
      call push_sub('cheby')
      
      allocate(zpsi1(m%np, h%d%dim, 0:2))
      zpsi1 = M_z0
      do j = te%exp_order-1, 0, -1
        call lalg_copy(m%np, h%d%dim, zpsi1(:,:, 1), zpsi1(:,:, 2))
        call lalg_copy(m%np, h%d%dim, zpsi1(:,:, 0), zpsi1(:,:, 1))

        call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
        zfact = 2*(-M_zI)**j*loct_bessel(j, h%spectral_half_span*timestep)

        call lalg_axpy(m%np, h%d%dim, cmplx(-h%spectral_middle_point, M_ZERO, PRECISION), &
             zpsi1(:,:, 1), zpsi1(:,:, 0))
        call lalg_scal(m%np, h%d%dim, cmplx(M_TWO/h%spectral_half_span, M_ZERO, PRECISION), zpsi1(:,:, 0))
        call lalg_axpy(m%np, h%d%dim, zfact, zpsi(:,:), zpsi1(:,:, 0))
        call lalg_axpy(m%np, h%d%dim, cmplx(-M_ONE, M_ZERO, PRECISION), zpsi1(:,:, 2),  zpsi1(:,:, 0))
      end do

      zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
      call lalg_scal(m%np, h%d%dim, exp(-M_zI*h%spectral_middle_point*timestep), zpsi(:,:))
      deallocate(zpsi1)
      
      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine cheby
    
    subroutine lanczos
      integer ::  korder, n, np
      CMPLX, allocatable :: hm(:,:), v(:,:,:), expo(:,:)
      FLOAT :: alpha, beta, res, tol, nrm
      
      call push_sub('lanczos')
      
      np = m%np*h%d%dim
      korder = te%exp_order
      tol = te%lanczos_tol
      allocate(v(m%np, h%d%dim, korder), &
           hm(korder, korder),         &
           expo(korder, korder))
      
      ! Normalize input vector, and put it into v(:, :, 1)
      nrm = zstates_nrm2(m, h%d%dim, zpsi(:,:))
      v(:, :, 1) = zpsi(:, :)/nrm
      
      ! Operate on v(:, :, 1) and place it onto w.
      call operate(v(:, :, 1), zpsi)
      alpha = zstates_dotp(m, h%d%dim, v(:, :, 1), zpsi)
      call lalg_axpy(m%np, h%d%dim, -M_z1*alpha, v(:,:, 1), zpsi(:,:))
      
      hm = M_z0; hm(1, 1) = alpha
      beta = zstates_nrm2(m, h%d%dim, zpsi(:,:))

      do n = 1, korder - 1
        v(:,:, n + 1) = zpsi(:,:)/beta
        hm(n+1, n) = beta
        call operate(v(:,:, n+1), zpsi)

        hm(n    , n + 1) = zstates_dotp(m, h%d%dim, v(:,:, n)  , zpsi)
        hm(n + 1, n + 1) = zstates_dotp(m, h%d%dim, v(:,:, n+1), zpsi)
        call lalg_gemv(np, 2, -M_z1, v(:,:, n), hm(n:, n + 1), M_z1, zpsi(:, 1))
        call zmatexp(n+1, hm(1:n+1, 1:n+1), expo(1:n+1, 1:n+1), -M_zI*timestep, method = 2)

        res = abs(beta*abs(expo(1, n+1)))
        beta = zstates_nrm2(m, h%d%dim, zpsi)

        if(beta < CNST(1.0e-12)) exit
        if(n>1 .and. res<tol) exit
      enddo
      korder = min(korder, n + 1)
      if(res > tol .and. beta > CNST(1.0e-12)) then
        write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
        call write_warning(1)
      endif
      
      ! zpsi = nrm * V * expo(1:korder, 1) = nrm * V * expo * V^(T) * zpsi
      call blas_gemv('N', np, korder, M_z1*nrm, v(1,1,1), np, expo(1,1), 1, M_z0, zpsi(1,1), 1)
      
      if(present(order)) order = korder
      deallocate(v, hm, expo)
      call pop_sub()
    end subroutine lanczos

#if defined(HAVE_FFT)
    subroutine split
      call push_sub('split')
      
      if(h%gauge == 2) then
        message(1) = 'Split operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      endif
      
      call zexp_vlpsi (m, h, zpsi, ik, t, -M_zI*timestep/M_TWO)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (m, h, zpsi, ik, -M_zI*timestep/M_TWO, .true.)
      call zexp_kinetic(m, h, zpsi, ik, te%cf, -M_zI*timestep)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (m, h, zpsi, ik, -M_zI*timestep/M_TWO, .false.)
      call zexp_vlpsi (m, h, zpsi, ik, t, -M_zI*timestep/M_TWO)
      
      if(present(order)) order = 0
      call pop_sub()
    end subroutine split
    
    subroutine suzuki
      FLOAT :: dt(5), p, pp(5)
      integer :: k
      
      call push_sub('suzuki')
      
      if(h%gauge == 2) then
        message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      endif
      
      p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
      pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
      dt(1:5) = pp(1:5)*timestep
      
      do k = 1, 5
        call zexp_vlpsi (m, h, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
        if (h%ep%nvnl > 0) call zexp_vnlpsi (m, h, zpsi, ik, -M_zI*dt(k)/M_TWO, .true.)
        call zexp_kinetic(m, h, zpsi, ik, te%cf, -M_zI*dt(k))
        if (h%ep%nvnl > 0) call zexp_vnlpsi (m, h, zpsi, ik, -M_zI*dt(k)/M_TWO, .false.)
        call zexp_vlpsi (m, h, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
      end do
      
      if(present(order)) order = 0
      call pop_sub()
    end subroutine suzuki
#endif
    
  end subroutine td_exp_dt

end module td_exp
