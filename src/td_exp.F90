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

module td_exp_m
  use global_m
  use messages_m
  use lib_oct_m
  use lib_oct_parser_m
  use datasets_m
  use blas_m
  use lib_basic_alg_m
  use math_m
  use mesh_m
#ifdef HAVE_FFT
  use cube_function_m
#endif
  use td_exp_split_m
  use varinfo_m

  implicit none

  private
  public ::                 &
    td_exp_t,               &
    td_exp_init,            &
    td_exp_end,             &
    td_exp_dt

  integer, parameter ::     &
    SPLIT_OPERATOR     = 0, &
    SUZUKI_TROTTER     = 1, &
    LANCZOS_EXPANSION  = 2, &
    FOURTH_ORDER       = 3, &
    CHEBYSHEV          = 4

  type td_exp_t
    integer   :: exp_method  ! which method is used to apply the exponential
    FLOAT     :: lanczos_tol ! tolerance for the lanczos method
    integer   :: exp_order   ! order to which the propagator is expanded
    type(zcf_t) :: cf          ! auxiliary cube for split operator methods
  end type td_exp_t


contains

  ! ---------------------------------------------------------
  subroutine td_exp_init(gr, te)
    type(grid_t),   intent(in)  :: gr
    type(td_exp_t), intent(out) :: te

    !%Variable TDExponentialMethod
    !%Type integer
    !%Default standard
    !%Section Time Dependent::Propagation
    !%Description
    !% Method used to numerically calculate the exponential of the Hamiltonian,
    !% a core part of the full algorithm used to approximate the evolution
    !% operator, specified through the variable <tt>TDEvolutionMethod</tt>.
    !% In the case of using the Magnus method, described below, the action of the exponential
    !% of the Magnus operator is also calculated through the algorithm specified
    !% by this variable.
    !%Option split 0
    !% It is important to distinguish between applying the split operator method
    !% to calculate the exponential of the Hamiltonian at a given time -- which
    !% is what this variable is referring to -- from the split operator method
    !% as an algorithm to approximate the full evolution operator <math>U(t+\delta t, t)</math>,
    !% and which will be described below as one of the possibilities
    !% of the variable <tt>TDEvolutionMethod</tt>.
    !% The equation that describes the split operator scheme is well known:
    !%
    !% <MATH>\exp_{\rm SO} (-i \delta t H) = \exp (-i \delta t/2 V) \exp (-i \delta t T) \exp (-i \delta t/2 V).</MATH>
    !%
    !% Note that this is a "kinetic referenced SO", since the kinetic term is sandwiched in the
    !% middle. This is so because in <tt>octopus</tt>, the states spend most of its time in real-space; doing
    !% it "potential referenced" would imply 4 FFTs instead of 2.
    !% This split-operator technique may be used in combination with, for example,
    !% the exponential midpoint rule as a means to approximate the evolution operator.
    !% In that case, the potential operator <i>V</i> that appears in the equation would be
    !% calculated at time <math>t+\delta t/2</math>, that is, in the middle of the time-step.
    !% However, note that if the split-operator method is invoked as a means to approximate
    !% the evolution operator (<tt>TDEvolutionMethod = 0</tt>), a different procedure is taken -- it
    !% will be described below --, and in fact the variable <tt>TDExponentialMethod</tt> has no
    !% effect at all.
    !%Option suzuki_trotter 1
    !% This is a higher-order SO based algorithm. See O. Sugino and Y. Miyamoto,
    !% Phys. Rev. B <b>59</b>, 2579 (1999). Allows for larger time-steps,
    !% but requires five times more time than the normal SO.
    !%
    !% The considerations made above for the SO algorithm about the distinction
    !% between using the method as a means to approximate <math>U(t+\delta t)</math> or as a
    !% means to approximate the exponential also apply here. Setting <tt>TDEvolutionMethod = 1</tt>
    !% enforces the use of the ST as an algorithm to approximate the full evolution operator,
    !% which is slightly different (see below).
    !%Option lanczos 2
    !% Allows for larger time-steps.
    !% However, the larger the time-step, the longer the computational time per time-step. 
    !% In certain cases, if the time-step is too large, the code will emit a warning
    !% whenever it considers that the evolution may not be properly proceeding --
    !% the Lanczos process did not converge. The method consists in a Krylov
    !% subspace approximation of the action of the exponential
    !% (see M. Hochbruck and C. Lubich, SIAM J. Numer. Anal. <b>34</b>, 1911 (1997) for details). 
    !% Two more variables control the performance of the method: the maximum dimension
    !% of this subspace (controlled by variable <tt>TDExpOrder</tt>), and
    !% the stopping criterium (controlled by variable <tt>TDLanczosTol</tt>).
    !% The smaller the stopping criterium, the more precisely the exponential
    !% is calculated, but also the larger the dimension of the Arnoldi
    !% subspace. If the maximum dimension allowed by <tt>TDExpOrder</tt> is not
    !% enough to meet the criterium, the above-mentioned warning is emitted.
    !%Option standard 3
    !% This method amounts to a straightforward application of the definition of
    !% the exponential of an operator, in terms of it Taylor expansion.
    !%
    !% <MATH>\exp_{\rm STD} (-i\delta t H) = \sum_{i=0}^{k} {(-i\delta t)^i\over{i!}} H^i.</MATH>
    !%
    !% The order <i>k</i> is determined by variable <i>TDExpOder</i>.
    !% Some numerical considerations (by Jeff Giansiracusa and George F. Bertsch;
    !% see http://www.phys.washington.edu/~bertsch/num3.ps)
    !% suggest the 4th order as especially suitable and stable.
    !%Option chebyshev 4
    !% In principle, the Chebyshev expansion
    !% of the exponential represents it more accurately than the canonical or standard expansion. 
    !% As in the latter case, <tt>TDExpOrder</tt> determines the order of the expansion.
    !%
    !% There exists a closed analytical form for the coefficients of the exponential in terms
    !% of Chebyshev polynomials:
    !%
    !% <MATH>\exp_{\rm CHEB} \left( -i\delta t H \right) = \sum_{k=0}^{\infty} (2-\delta_{k0})(-i)^{k}J_k(\delta t) T_k(H),</MATH>
    !%
    !% where <math>J_k</math> are the Bessel functions of the first kind, and H has te be previously
    !% scaled to <math>[-1,1]</math>.
    !% See H. Tal-Ezer and R. Kosloff, J. Chem. Phys. <b>81</b>,
    !% 3967 (1984); R. Kosloff, Annu. Rev. Phys. Chem. <b>45</b>, 145 (1994);
    !% C. W. Clenshaw, MTAC <b>9</b>, 118 (1955).
    !%End
    call loct_parse_int(check_inp('TDExponentialMethod'), FOURTH_ORDER, te%exp_method)

    select case(te%exp_method)
    case(FOURTH_ORDER)
    case(CHEBYSHEV)
    case(LANCZOS_EXPANSION)
      !%Variable TDLanczosTol
      !%Type float
      !%Default 1e-5
      !%Section Time Dependent::Propagation
      !%Description
      !% An internal tolerance variable for the Lanczos method. The smaller, the more
      !% precisely the exponential is calculated, and also the bigger the dimension
      !% of the Krylov subspace needed to perform the algorithm. One should carefully
      !% make sure that this value is not too big, or else the evolution will be
      !% wrong.
      !%End
      call loct_parse_float(check_inp('TDLanczosTol'), CNST(1e-5), te%lanczos_tol)
      if (te%lanczos_tol <= M_ZERO) call input_error('TDLanczosTol')

#if defined(HAVE_FFT)
    case(SPLIT_OPERATOR)
    case(SUZUKI_TROTTER)
#endif

    case default
      call input_error('TDExponentialMethod')
    end select
    call messages_print_var_option(stdout, 'TDExponentialMethod', te%exp_method)

    if(te%exp_method==FOURTH_ORDER.or.te%exp_method==CHEBYSHEV.or.te%exp_method==LANCZOS_EXPANSION) then
      !%Variable TDExpOrder
      !%Type integer
      !%Default 4
      !%Section Time Dependent::Propagation
      !%Description
      !% For <tt>TDExponentialMethod</tt> equal <tt>standard</tt> or <tt>chebyshev</tt>, the order to which
      !% the exponential is expanded. For the Lanczos approximation, it is the maximum
      !% Lanczos-subspace dimension.
      !%End
      call loct_parse_int(check_inp('TDExpOrder'), 4, te%exp_order)
      if (te%exp_order < 2) call input_error('TDExpOrder')

#if defined(HAVE_FFT)
    else if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
      call zcf_new(gr%m%l, te%cf)
      call zcf_fft_init(te%cf, gr%sb)
#endif
    end if

  end subroutine td_exp_init


  ! ---------------------------------------------------------
  subroutine td_exp_end(te)
    type(td_exp_t), intent(inout) :: te

    if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
      call zcf_free(te%cf)
    end if
  end subroutine td_exp_end


  ! ---------------------------------------------------------
  subroutine td_exp_dt(te, gr, h, zpsi, ik, timestep, t, order, vmagnus)
    type(td_exp_t),      intent(inout) :: te
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: h
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    FLOAT,               intent(in)    :: timestep, t
    integer, optional,   intent(out)   :: order ! For the methods that rely on Hamiltonian-vector
    ! multiplication, the number of these.
    FLOAT,     optional, intent(in)    :: vmagnus(NP, h%d%nspin, 2)

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
      end if
      call lanczos
#if defined(HAVE_FFT)
    case(SPLIT_OPERATOR);    call split
    case(SUZUKI_TROTTER);    call suzuki
#endif
    case(CHEBYSHEV);         call cheby
    end select

    call pop_sub()
  contains

    ! ---------------------------------------------------------
    subroutine operate(psi, oppsi)
      CMPLX, intent(inout) :: psi(:, :)
      CMPLX, intent(inout) :: oppsi(:, :)

      if(apply_magnus) then
        call zmagnus(h, gr, psi, oppsi, ik, vmagnus)
      else
        call zHpsi(h, gr, psi, oppsi, ik, t)
      end if

    end subroutine operate


    ! ---------------------------------------------------------
    subroutine fourth()
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:), hzpsi1(:,:)
      integer :: i, idim

      call push_sub('td_exp.fourth')

      ALLOCATE(zpsi1 (NP_PART, h%d%dim), NP_PART*h%d%dim)
      ALLOCATE(hzpsi1(NP,      h%d%dim), NP     *h%d%dim)
      zfact = M_z1
      do idim = 1, h%d%dim
        call lalg_copy(NP, zpsi(:, idim), zpsi1(:, idim))
      end do
      do i = 1, te%exp_order
        zfact = zfact*(-M_zI*timestep)/i
        call operate(zpsi1, hzpsi1)
        do idim = 1, h%d%dim
          call lalg_axpy(NP, zfact, hzpsi1(:, idim), zpsi(:, idim))
        end do
        if(i .ne. te%exp_order) then
          do idim = 1, h%d%dim
            call lalg_copy(NP, hzpsi1(:, idim), zpsi1(:, idim))
          end do
        end if
      end do
      deallocate(zpsi1, hzpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine fourth


    ! ---------------------------------------------------------
    subroutine cheby()
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
      integer :: j, idim
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:,:)

      call push_sub('td_exp.cheby')

      ALLOCATE(zpsi1(NP_PART, h%d%dim, 0:2), NP_PART*h%d%dim*3)
      zpsi1 = M_z0
      do j = te%exp_order-1, 0, -1
        do idim = 1, h%d%dim
          call lalg_copy(NP, zpsi1(:, idim, 1), zpsi1(:, idim, 2))
          call lalg_copy(NP, zpsi1(:, idim, 0), zpsi1(:, idim, 1))
        end do

        call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
        zfact = 2*(-M_zI)**j*loct_bessel(j, h%spectral_half_span*timestep)

        do idim = 1, h%d%dim
          call lalg_axpy(NP, cmplx(-h%spectral_middle_point, M_ZERO, PRECISION), &
               zpsi1(:, idim, 1), zpsi1(:, idim, 0))
          call lalg_scal(NP, cmplx(M_TWO/h%spectral_half_span, M_ZERO, PRECISION), zpsi1(:, idim, 0))
          call lalg_axpy(NP, zfact, zpsi(:, idim), zpsi1(:, idim, 0))
          call lalg_axpy(NP, cmplx(-M_ONE, M_ZERO, PRECISION), zpsi1(:, idim, 2),  zpsi1(:, idim, 0))
        end do
      end do

      zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
      do idim = 1, h%d%dim
        call lalg_scal(NP, exp(-M_zI*h%spectral_middle_point*timestep), zpsi(:, idim))
      end do
      deallocate(zpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine cheby


    ! ---------------------------------------------------------
    subroutine lanczos
      integer ::  korder, n, k, iflag, lwsp, ns, i, j, iexph, idim
      integer, allocatable :: ipiv(:)
      CMPLX, allocatable :: hm(:,:), v(:,:,:), expo(:,:), wsp(:)
      FLOAT :: beta, res, tol !, nrm

      call push_sub('td_exp.lanczos')

      tol    = te%lanczos_tol
      ALLOCATE( v(NP_PART, h%d%dim, te%exp_order+1), NP_PART*h%d%dim*(te%exp_order+1))
      ALLOCATE(  hm(te%exp_order+1, te%exp_order+1), (te%exp_order+1)*(te%exp_order+1))
      ALLOCATE(expo(te%exp_order+1, te%exp_order+1), (te%exp_order+1)*(te%exp_order+1))
      v = M_z0; hm = M_z0; expo = M_z0

      lwsp = 4*(te%exp_order+1)**2+7
      ALLOCATE(wsp(lwsp), lwsp)
      ALLOCATE(ipiv(te%exp_order+1), (te%exp_order+1))

      ! Normalize input vector, and put it into v(:, :, 1)
      beta = zstates_nrm2(gr%m, h%d%dim, zpsi)
      v(1:NP, :, 1) = zpsi(1:NP, :)/beta

      ! This is the Lanczos loop...
      do n = 1, te%exp_order
        call operate(v(:,:, n), v(:, :, n+1))
        korder = n

        do k = max(1, n-1), n
          hm(k, n) = zstates_dotp(gr%m, h%d%dim, v(:, :, k), v(:, :, n+1))
          do idim = 1, h%d%dim
            call lalg_axpy(NP, -hm(k, n), v(:, idim, k), v(:, idim, n+1))
          end do
        end do
        hm(n+1, n) = zstates_nrm2(gr%m, h%d%dim, v(:, :, n+1))
        do idim = 1, h%d%dim
          call lalg_scal(NP, M_z1/hm(n+1, n), v(:, idim, n+1)) 
        end do
        call zgpadm(6, n, timestep, -M_zI*hm(1:n, 1:n), n, wsp, lwsp, ipiv(1:n), iexph, ns, iflag)
        k = 0
        do i = 1, n
          do j = 1, n
            expo(j, i) = wsp(iexph + k); k = k + 1
          end do
        end do

        res = abs(hm(n+1, n)*abs(expo(1, n)))
        if(abs(hm(n+1, n)) < CNST(1.0e-12)) exit ! (very unlikely) happy breakdown!!! Yuppi!!!
        if(n>2 .and. res<tol) exit
      end do

      if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
        write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
        call write_warning(1)
      end if

      call zgpadm(6, korder+1, timestep, -M_zI*hm(1:korder+1, 1:korder+1), korder+1, wsp, &
        lwsp, ipiv(1:korder+1), iexph, ns, iflag)
      k = 0
      do i = 1, korder+1
        do j = 1, korder+1
          expo(j, i) = wsp( iexph + k); k = k + 1
        end do
      end do
      do i = 1, korder + 1
        v(NP+1:NP_PART, 1:h%d%dim, i) = M_Z0
      end do
      ! zpsi = nrm * V * expo(1:korder, 1) = nrm * V * expo * V^(T) * zpsi
      call lalg_gemv(NP_PART, h%d%dim, korder+1, M_z1*beta, v, expo(1:korder+1, 1), M_z0, zpsi)

      if(present(order)) order = korder
      deallocate(v, hm, expo, ipiv, wsp)
      call pop_sub()
    end subroutine lanczos


#if defined(HAVE_FFT)
    ! ---------------------------------------------------------
    subroutine split
      call push_sub('td_exp.split')

      if(h%gauge == VELOCITY) then
        message(1) = 'Split operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      end if

      call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*timestep/M_TWO)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*timestep/M_TWO, .true.)
      call zexp_kinetic(gr, h, zpsi, te%cf, -M_zI*timestep)
      if(h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, zpsi, -M_zI*timestep/M_TWO, .false.)
      call zexp_vlpsi (gr, h, zpsi, ik, t, -M_zI*timestep/M_TWO)

      if(present(order)) order = 0
      call pop_sub()
    end subroutine split


    ! ---------------------------------------------------------
    subroutine suzuki
      FLOAT :: dt(5), p, pp(5)
      integer :: k

      call push_sub('td_exp.suzuki')

      if(h%gauge == 2) then
        message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      end if

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

end module td_exp_m
