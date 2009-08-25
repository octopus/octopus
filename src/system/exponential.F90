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
!! $Id$

#include "global.h"

module exponential_m
  use batch_m
  use blas_m
  use cube_function_m
  use datasets_m
  use hardware_m
  use fourier_space_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_math_m
  use loct_parser_m
  use mesh_function_m
  use profiling_m
  use states_calc_m
  use exponential_split_m
  use varinfo_m
  use states_m

  implicit none

  private
  public ::                      &
    exponential_t,               &
    exponential_init,            &
    exponential_copy,            &
    exponential_end,             &
    exponential_apply_batch,     &
    exponential_apply,           &
    exponential_apply_all

  integer, public, parameter :: &
    SPLIT_OPERATOR     = 0,     &
    SUZUKI_TROTTER     = 1,     &
    LANCZOS_EXPANSION  = 2,     &
    TAYLOR             = 3,     &
    CHEBYSHEV          = 4

  type exponential_t
    integer     :: exp_method  ! which method is used to apply the exponential
    FLOAT       :: lanczos_tol ! tolerance for the Lanczos method
    integer     :: exp_order   ! order to which the propagator is expanded
    type(zcf_t) :: cf          ! auxiliary cube for split operator methods
  end type exponential_t

contains

  ! ---------------------------------------------------------
  subroutine exponential_init(te, gr)
    type(grid_t),        intent(in)  :: gr
    type(exponential_t), intent(out) :: te

    !%Variable TDExponentialMethod
    !%Type integer
    !%Default taylor
    !%Section Time-Dependent::Propagation
    !%Description
    !% Method used to numerically calculate the exponential of the Hamiltonian,
    !% a core part of the full algorithm used to approximate the evolution
    !% operator, specified through the variable <tt>TDEvolutionMethod</tt>.
    !% In the case of using the Magnus method, described below, the action of the exponential
    !% of the Magnus operator is also calculated through the algorithm specified
    !% by this variable.
    !%Option split 0
    !% It is important to distinguish between applying the split-operator method
    !% to calculate the exponential of the Hamiltonian at a given time -- which
    !% is what this variable is referring to -- from the split-operator method
    !% as an algorithm to approximate the full evolution operator <math>U(t+\delta t, t)</math>,
    !% and which will be described below as one of the possibilities
    !% of the variable <tt>TDEvolutionMethod</tt>.
    !% The equation that describes the split-operator scheme is well known:
    !%
    !% <MATH>\exp_{\rm SO} (-i \delta t H) = \exp (-i \delta t/2 V) \exp (-i \delta t T) \exp (-i \delta t/2 V).</MATH>
    !%
    !% Note that this is a "kinetic-referenced SO", since the kinetic term is sandwiched in the
    !% middle. This is so because in <tt>octopus</tt>, the states spend most of their time in real-space; doing
    !% it "potential-referenced" would imply 4 FFTs instead of 2.
    !% This split-operator technique may be used in combination with, for example,
    !% the exponential midpoint rule as a means to approximate the evolution operator.
    !% In that case, the potential operator <i>V</i> that appears in the equation would be
    !% calculated at time <math>t+\delta t/2</math>, that is, in the middle of the time-step.
    !% However, note that if the split-operator method is invoked as a means to approximate
    !% the evolution operator (<tt>TDEvolutionMethod = 0</tt>), a different procedure is taken -- 
    !% described below -- and in fact the variable <tt>TDExponentialMethod</tt> has no
    !% effect at all.
    !%Option suzuki_trotter 1
    !% This is a higher-order SO-based algorithm. See O. Sugino and Y. Miyamoto,
    !% Phys. Rev. B <b>59</b>, 2579 (1999). Allows for larger time-steps,
    !% but requires five times more time than the normal SO.
    !%
    !% The considerations above for the SO algorithm about the distinction
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
    !% the stopping criterion (controlled by variable <tt>TDLanczosTol</tt>).
    !% The smaller the stopping criterion, the more precisely the exponential
    !% is calculated, but also the larger the dimension of the Arnoldi
    !% subspace. If the maximum dimension allowed by <tt>TDExpOrder</tt> is not
    !% enough to meet the criterion, the above-mentioned warning is emitted.
    !%Option taylor 3
    !% This method amounts to a straightforward application of the definition of
    !% the exponential of an operator, in terms of its Taylor expansion.
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
    !% where <math>J_k</math> are the Bessel functions of the first kind, and H has to be previously
    !% scaled to <math>[-1,1]</math>.
    !% See H. Tal-Ezer and R. Kosloff, J. Chem. Phys. <b>81</b>,
    !% 3967 (1984); R. Kosloff, Annu. Rev. Phys. Chem. <b>45</b>, 145 (1994);
    !% C. W. Clenshaw, MTAC <b>9</b>, 118 (1955).
    !%End
    call loct_parse_int(datasets_check('TDExponentialMethod'), TAYLOR, te%exp_method)

    select case(te%exp_method)
    case(TAYLOR)
    case(CHEBYSHEV)
    case(LANCZOS_EXPANSION)
      !%Variable TDLanczosTol
      !%Type float
      !%Default 1e-5
      !%Section Time-Dependent::Propagation
      !%Description
      !% An internal tolerance variable for the Lanczos method. The smaller, the more
      !% precisely the exponential is calculated, and also the bigger the dimension
      !% of the Krylov subspace needed to perform the algorithm. One should carefully
      !% make sure that this value is not too big, or else the evolution will be
      !% wrong.
      !%End
      call loct_parse_float(datasets_check('TDLanczosTol'), CNST(1e-5), te%lanczos_tol)
      if (te%lanczos_tol <= M_ZERO) call input_error('TDLanczosTol')

    case(SPLIT_OPERATOR)
    case(SUZUKI_TROTTER)

    case default
      call input_error('TDExponentialMethod')
    end select
    call messages_print_var_option(stdout, 'TDExponentialMethod', te%exp_method)

    if(te%exp_method==TAYLOR.or.te%exp_method==CHEBYSHEV.or.te%exp_method==LANCZOS_EXPANSION) then
      !%Variable TDExpOrder
      !%Type integer
      !%Default 4
      !%Section Time-Dependent::Propagation
      !%Description
      !% For <tt>TDExponentialMethod</tt> = <tt>standard</tt> or <tt>chebyshev</tt>, 
      !% the order to which the exponential is expanded. For the Lanczos approximation, 
      !% it is the Lanczos-subspace dimension.
      !%End
      call loct_parse_int(datasets_check('TDExpOrder'), 4, te%exp_order)
      if (te%exp_order < 2) call input_error('TDExpOrder')

    else if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) then
      call zcf_new(gr%mesh%idx%ll, te%cf)
      call zcf_fft_init(te%cf, gr%sb)
    end if

  end subroutine exponential_init

  ! ---------------------------------------------------------
  subroutine exponential_end(te)
    type(exponential_t), intent(inout) :: te

    if(te%exp_method==SPLIT_OPERATOR.or.te%exp_method==SUZUKI_TROTTER) call zcf_free(te%cf)

  end subroutine exponential_end

  ! ---------------------------------------------------------
  subroutine exponential_copy(teo, tei)
    type(exponential_t), intent(inout) :: teo
    type(exponential_t), intent(in)    :: tei

    teo%exp_method  = tei%exp_method
    teo%lanczos_tol = tei%lanczos_tol
    teo%exp_order   = tei%exp_order
    if(teo%exp_method == SPLIT_OPERATOR .or. teo%exp_method == SUZUKI_TROTTER) call zcf_new_from(teo%cf, tei%cf)

  end subroutine exponential_copy

  ! ---------------------------------------------------------
  ! This routine performs the operation:
  !
  ! exp{-i*deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! If imag_time is present and is set to true, it performa instead:
  !
  ! exp{ deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! If the hamiltonian contains an inhomogeneous term, the operation is:
  !
  ! exp{-i*deltat*hm(t)}|zpsi> + deltat*phi{-i*deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! where:
  !
  ! phi(x) = (e^x - 1)/x
  ! ---------------------------------------------------------
  subroutine exponential_apply(te, gr, hm, zpsi, ist, ik, deltat, t, order, vmagnus, imag_time)
    type(exponential_t), intent(inout) :: te
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    FLOAT,               intent(in)    :: deltat, t
    integer, optional,   intent(inout) :: order
    FLOAT,   optional,   intent(in)    :: vmagnus(gr%mesh%np, hm%d%nspin, 2)
    logical, optional,   intent(in)    :: imag_time

    CMPLX   :: timestep
    logical :: apply_magnus
    type(profile_t), save :: exp_prof

    call push_sub('exponential.exponential_td')
    call profiling_in(exp_prof, "EXPONENTIAL")

    ! The only method that is currently taking care of the presence of an inhomogeneous
    ! term is the Lanczos expansion.
    ! However, I disconnect this check, because this routine is sometimes called in an
    ! auxiliary way, in order to compute (1-hm(t)*deltat)|zpsi>.
    ! This should be cleaned up.
    !_ASSERT(.not.(hamiltonian_inh_term(hm) .and. (te%exp_method .ne. LANCZOS_EXPANSION)))

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    ! If we want to use imaginary time, timestep = i*deltat
    ! Otherwise, timestep is simply equal to deltat.
    timestep = cmplx(deltat, M_ZERO)
    if(present(imag_time)) then
      if(imag_time) then
        select case(te%exp_method)
          case(TAYLOR, LANCZOS_EXPANSION)
            timestep = M_zI * deltat
          case default
            write(message(1), '(a)') &
              'Imaginary time evolution can only be performed with the Lanczos'
            write(message(2), '(a)') &
              'exponentiation scheme ("TDExponentialMethod = lanczos") or with the'
            write(message(3), '(a)') &
              'Taylor expansion ("TDExponentialMethod = taylor") method.'
            call write_fatal(3)
        end select
      end if
    end if

    select case(te%exp_method)
    case(TAYLOR)
      call taylor_series
    case(LANCZOS_EXPANSION)
      call lanczos
    case(SPLIT_OPERATOR)
      call split
    case(SUZUKI_TROTTER)
      call suzuki
    case(CHEBYSHEV)
      call cheby
    end select

    call profiling_out(exp_prof)
    call pop_sub()
  contains

    ! ---------------------------------------------------------
    subroutine operate(psi, oppsi)
      CMPLX, intent(inout) :: psi(:, :)
      CMPLX, intent(inout) :: oppsi(:, :)

      if(apply_magnus) then
        call zmagnus(hm, gr, psi, oppsi, ik, vmagnus)
      else
        call zhamiltonian_apply(hm, gr, psi, oppsi, ist, ik, t)
      end if

    end subroutine operate
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine taylor_series()
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:), hzpsi1(:,:)
      integer :: i, idim
      logical :: zfact_is_real

      call push_sub('exponential.taylor_series')

      SAFE_ALLOCATE(zpsi1 (1:gr%mesh%np_part, 1:hm%d%dim))
      SAFE_ALLOCATE(hzpsi1(1:gr%mesh%np,      1:hm%d%dim))

      zfact = M_z1
      zfact_is_real = .true.

      do idim = 1, hm%d%dim
        call lalg_copy(gr%mesh%np, zpsi(:, idim), zpsi1(:, idim))
      end do

      do i = 1, te%exp_order
        zfact = zfact*(-M_zI*timestep)/i
        zfact_is_real = .not. zfact_is_real
        
        call operate(zpsi1, hzpsi1)

        if(zfact_is_real) then
          do idim = 1, hm%d%dim
            call lalg_axpy(gr%mesh%np, real(zfact), hzpsi1(:, idim), zpsi(:, idim))
          end do
        else
          do idim = 1, hm%d%dim
            call lalg_axpy(gr%mesh%np, zfact, hzpsi1(:, idim), zpsi(:, idim))
          end do
        end if

        if(i .ne. te%exp_order) then
          do idim = 1, hm%d%dim
            call lalg_copy(gr%mesh%np, hzpsi1(:, idim), zpsi1(:, idim))
          end do
        end if

      end do
      SAFE_DEALLOCATE_A(zpsi1)
      SAFE_DEALLOCATE_A(hzpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine taylor_series
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine cheby()
      ! Calculates the exponential of the hamiltonian through a expansion in 
      ! Chebyshev polynomials.
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

      call push_sub('exponential.cheby')

      SAFE_ALLOCATE(zpsi1(1:gr%mesh%np_part, 1:hm%d%dim, 0:2))
      zpsi1 = M_z0
      do j = te%exp_order-1, 0, -1
        do idim = 1, hm%d%dim
          call lalg_copy(gr%mesh%np, zpsi1(:, idim, 1), zpsi1(:, idim, 2))
          call lalg_copy(gr%mesh%np, zpsi1(:, idim, 0), zpsi1(:, idim, 1))
        end do

        call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
        zfact = 2*(-M_zI)**j*loct_bessel(j, hm%spectral_half_span*deltat)

        do idim = 1, hm%d%dim
          call lalg_axpy(gr%mesh%np, -hm%spectral_middle_point, zpsi1(:, idim, 1), &
            zpsi1(:, idim, 0))
          call lalg_scal(gr%mesh%np, M_TWO/hm%spectral_half_span, zpsi1(:, idim, 0))
          call lalg_axpy(gr%mesh%np, zfact, zpsi(:, idim), zpsi1(:, idim, 0))
          call lalg_axpy(gr%mesh%np, -M_ONE, zpsi1(:, idim, 2),  zpsi1(:, idim, 0))
        end do
      end do

      zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
      do idim = 1, hm%d%dim
        call lalg_scal(gr%mesh%np, exp(-M_zI*hm%spectral_middle_point*deltat), zpsi(:, idim))
      end do
      SAFE_DEALLOCATE_A(zpsi1)

      if(present(order)) order = te%exp_order
      call pop_sub()
    end subroutine cheby
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine lanczos
      integer ::  iter, l, idim
      CMPLX, allocatable :: hamilt(:,:), v(:,:,:), expo(:,:), tmp(:, :), psi(:, :)
      FLOAT :: beta, res, tol !, nrm
      CMPLX :: pp

      call push_sub('exponential.lanczos')


      SAFE_ALLOCATE(     v(1:gr%mesh%np, 1:hm%d%dim, 1:te%exp_order+1))
      SAFE_ALLOCATE(   tmp(1:gr%mesh%np, 1:hm%d%dim))
      SAFE_ALLOCATE(hamilt(1:te%exp_order+1, 1:te%exp_order+1))
      SAFE_ALLOCATE(  expo(1:te%exp_order+1, 1:te%exp_order+1))
      SAFE_ALLOCATE(   psi(1:gr%mesh%np_part, 1:hm%d%dim))

      tol    = te%lanczos_tol
      pp = deltat
      if(.not. present(imag_time)) pp = -M_zI*pp

      beta = zmf_nrm2(gr%mesh, hm%d%dim, zpsi)
      ! If we have a null vector, no need to compute the action of the exponential.
      if(beta > CNST(1.0e-12)) then

        hamilt = M_z0
        expo = M_z0

        ! Normalize input vector, and put it into v(:, :, 1)
        v(1:gr%mesh%np, 1:hm%d%dim, 1) = zpsi(1:gr%mesh%np, 1:hm%d%dim)/beta

        ! This is the Lanczos loop...
        do iter = 1, te%exp_order

          !copy v(:, :, n) to an array of size 1:gr%mesh%np_part
          do idim = 1, hm%d%dim
            call lalg_copy(gr%mesh%np, v(:, idim, iter), zpsi(:, idim))
          end do

          !to apply the hamiltonian
          call operate(zpsi, v(:, :,  iter + 1))
        
          if(hamiltonian_hermitean(hm)) then
            l = max(1, iter - 1)
          else
            l = 1
          end if

          !orthogonalize against previous vectors
          call zstates_gram_schmidt(gr%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), v(:, :, iter + 1), &
            normalize = .true., overlap = hamilt(l:iter, iter), norm = hamilt(iter + 1, iter))

          call zlalg_exp(iter, pp, hamilt, expo, hamiltonian_hermitean(hm))

          res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

          if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
          if(iter > 2 .and. res < tol) exit
        end do

        if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
          write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
          call write_warning(1)
        end if

        ! zpsi = nrm * V * expo(1:iter, 1) = nrm * V * expo * V^(T) * zpsi
        call lalg_gemv(gr%mesh%np, hm%d%dim, iter, M_z1*beta, v, expo(1:iter, 1), M_z0, tmp)

        do idim = 1, hm%d%dim
          call lalg_copy(gr%mesh%np, tmp(:, idim), zpsi(:, idim))
        end do

      end if

      ! We have an inhomogeneous term.
      if( hamiltonian_inh_term(hm) ) then

        beta = zmf_nrm2(gr%mesh, hm%d%dim, hm%inh_st%zpsi(:, :, ist, ik))
        if(beta > CNST(1.0e-12)) then

          hamilt = M_z0
          expo = M_z0

          v(1:gr%mesh%np, 1:hm%d%dim, 1) = &
            hm%inh_st%zpsi(1:gr%mesh%np, 1:hm%d%dim, ist, ik)/beta

          psi = M_z0
          ! This is the Lanczos loop...
          do iter = 1, te%exp_order
            !copy v(:, :, n) to an array of size 1:gr%mesh%np_part
            do idim = 1, hm%d%dim
              call lalg_copy(gr%mesh%np, v(:, idim, iter), psi(:, idim))
            end do

            !to apply the hamiltonian
            call operate(psi, v(:, :, iter + 1))
  

            if(hamiltonian_hermitean(hm)) then
              l = max(1, iter - 1)
            else
              l = 1
            end if

            !orthogonalize against previous vectors
            call zstates_gram_schmidt(gr%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), &
              v(:, :, iter + 1), normalize = .true., overlap = hamilt(l:iter, iter), &
              norm = hamilt(iter + 1, iter))

            call zlalg_phi(iter, pp, hamilt, expo, hamiltonian_hermitean(hm))
 
            res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

            if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
            if(iter > 2 .and. res < tol) exit
          end do

          if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
            write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
            call write_warning(1)
          end if

          call lalg_gemv(gr%mesh%np, hm%d%dim, iter, M_z1*beta, v, expo(1:iter, 1), M_z0, tmp)

          do idim = 1, hm%d%dim
            call lalg_copy(gr%mesh%np, tmp(:, idim), psi(:, idim))
          end do

          zpsi = zpsi + deltat*psi
        end if
      end if


      if(present(order)) order = te%exp_order
      SAFE_DEALLOCATE_A(v)
      SAFE_DEALLOCATE_A(hamilt)
      SAFE_DEALLOCATE_A(expo)
      SAFE_DEALLOCATE_A(tmp)
      SAFE_DEALLOCATE_A(psi)
      call pop_sub()
    end subroutine lanczos
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine split
      call push_sub('exponential.split')

      if(hm%gauge == VELOCITY) then
        message(1) = 'Split operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      end if

      call zexp_vlpsi (gr, hm, zpsi, ik, t, -M_zI*deltat/M_TWO)
      if(hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, zpsi, -M_zI*deltat/M_TWO, .true.)
      call zexp_kinetic(gr, hm, zpsi, te%cf, -M_zI*deltat)
      if(hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, zpsi, -M_zI*deltat/M_TWO, .false.)
      call zexp_vlpsi (gr, hm, zpsi, ik, t, -M_zI*deltat/M_TWO)

      if(present(order)) order = 0
      call pop_sub()
    end subroutine split
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine suzuki
      FLOAT :: dt(5), p, pp(5)
      integer :: k

      call push_sub('exponential.suzuki')

      if(hm%gauge == 2) then
        message(1) = 'Suzuki-Trotter operator does not work well if velocity gauge is used.'
        call write_fatal(1)
      end if

      p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
      pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
      dt(1:5) = pp(1:5)*deltat

      do k = 1, 5
        call zexp_vlpsi (gr, hm, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
        if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, zpsi, -M_zI*dt(k)/M_TWO, .true.)
        call zexp_kinetic(gr, hm, zpsi, te%cf, -M_zI*dt(k))
        if (hm%ep%non_local) call zexp_vnlpsi (gr%mesh, hm, zpsi, -M_zI*dt(k)/M_TWO, .false.)
        call zexp_vlpsi (gr, hm, zpsi, ik, t, -M_zI*dt(k)/M_TWO)
      end do

      if(present(order)) order = 0
      call pop_sub()
    end subroutine suzuki
    ! ---------------------------------------------------------

  end subroutine exponential_apply

  subroutine exponential_apply_batch(te, gr, hm, psib, ik, deltat, t)
    type(exponential_t), intent(inout) :: te
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    integer,             intent(in)    :: ik
    type(batch_t),       intent(inout) :: psib
    FLOAT,               intent(in)    :: deltat
    FLOAT,               intent(in)    :: t
    
    integer :: ii, ist
    CMPLX, pointer :: psi(:, :)

    if (te%exp_method == TAYLOR) then 
      call taylor_series_batch
    else
      
      do ii = 1, psib%nst
        psi  => psib%states(ii)%zpsi
        ist  =  psib%states(ii)%ist
        
        call exponential_apply(te, gr, hm, psi, ist, ik, deltat, t)
      end do

    end if
    
  contains
    
    subroutine taylor_series_batch()
      CMPLX :: zfact
      CMPLX, allocatable :: psi1(:, :, :), hpsi1(:, :, :)
      integer :: iter, idim
      logical :: zfact_is_real
      integer :: st_start, st_end
      type(batch_t) :: psi1b, hpsi1b
      integer :: bsize, ip
      type(profile_t), save :: prof

      call push_sub('exponential.taylor_series_batch')
      call profiling_in(prof, "EXP_TAYLOR_BATCH")

      SAFE_ALLOCATE(psi1 (1:gr%mesh%np_part, 1:hm%d%dim, 1:psib%nst))
      SAFE_ALLOCATE(hpsi1(1:gr%mesh%np, 1:hm%d%dim, 1:psib%nst))

      st_start = psib%states(1)%ist
      st_end = psib%states(psib%nst)%ist

      zfact = M_z1
      zfact_is_real = .true.

      call batch_init(psi1b, hm%d%dim, st_start, st_end, psi1)
      call batch_init(hpsi1b, hm%d%dim, st_start, st_end, hpsi1)

      !$omp parallel do private(ii, idim)
      do ii = 1, psib%nst
        do idim = 1, psib%dim
          call lalg_copy(gr%mesh%np, psib%states(ii)%zpsi(:, idim), psi1b%states(ii)%zpsi(:, idim))
        end do
      end do

      do iter = 1, te%exp_order
        zfact = zfact*(-M_zI*deltat)/iter
        zfact_is_real = .not. zfact_is_real

        call zhamiltonian_apply_batch(hm, gr, psi1b, hpsi1b, ik, t)

        !$omp parallel do private(ii, idim, ip, bsize)
        do ii = 1, psib%nst
          do idim = 1, hm%d%dim
    
            do ip = 1, gr%mesh%np, hardware%zblock_size
              bsize = min(hardware%zblock_size, gr%mesh%np - ip + 1)
              if(zfact_is_real) then
                call blas_axpy(bsize, real(zfact, REAL_PRECISION), hpsi1(ip, idim, ii), psib%states(ii)%zpsi(ip, idim))
              else
                call blas_axpy(bsize, zfact, hpsi1(ip, idim, ii), 1, psib%states(ii)%zpsi(ip, idim), 1)
              end if
              if(iter /= te%exp_order) call blas_copy(bsize, hpsi1(ip, idim, ii), 1, psi1(ip, idim, ii), 1)
            end do

          end do
        end do
        !$omp end parallel do
      end do
      
      call batch_end(hpsi1b)
      call batch_end(psi1b)

      call profiling_count_operations(psib%nst*hm%d%dim*dble(gr%mesh%np)*te%exp_order*CNST(6.0))

      SAFE_DEALLOCATE_A(psi1)
      SAFE_DEALLOCATE_A(hpsi1)
      
      call profiling_out(prof)
      call pop_sub()
    end subroutine taylor_series_batch
    
  end subroutine exponential_apply_batch


  ! ---------------------------------------------------------
  ! This routine performs the operation:
  !
  ! exp{-i*deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! If imag_time is present and is set to true, it performa instead:
  !
  ! exp{ deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! If the hamiltonian contains an inhomogeneous term, the operation is:
  !
  ! exp{-i*deltat*hm(t)}|zpsi> + deltat*phi{-i*deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! where:
  !
  ! phi(x) = (e^x - 1)/x
  ! ---------------------------------------------------------
  subroutine exponential_apply_all(te, gr, hm, psi, deltat, t, order, vmagnus, imag_time)
    type(exponential_t), intent(inout) :: te
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: psi
    FLOAT,               intent(in)    :: deltat, t
    integer, optional,   intent(inout) :: order
    FLOAT,   optional,   intent(in)    :: vmagnus(gr%mesh%np, hm%d%nspin, 2)
    logical, optional,   intent(in)    :: imag_time

    integer :: ik, ist
    CMPLX   :: timestep
    logical :: apply_magnus
    CMPLX :: zfact
    integer :: i, idim
    logical :: zfact_is_real

    type(states_t) :: psi1, hpsi1

    call push_sub('exponential.exponential_apply_all')

    ASSERT(te%exp_method .eq. TAYLOR)

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    timestep = cmplx(deltat, M_ZERO)
    if(present(imag_time)) then
      if(imag_time) timestep = M_zI * deltat
    end if

    call states_copy(psi1, psi)
    call states_copy(hpsi1, psi)

    forall(ik = psi%d%kpt%start:psi%d%kpt%end, ist = psi%st_start:psi%st_end)
      psi1%zpsi(:, :, ist, ik)  = M_z0
      hpsi1%zpsi(:, :, ist, ik) = M_z0
    end forall

    zfact = M_z1
    zfact_is_real = .true.

    do i = 1, te%exp_order
      zfact = zfact*(-M_zI*timestep)/i
      zfact_is_real = .not. zfact_is_real
      
      if (i == 1) then
        call zhamiltonian_apply_all (hm, gr, psi, hpsi1, t)
      else
        call zhamiltonian_apply_all(hm, gr, psi1, hpsi1, t)
      end if

      if(zfact_is_real) then
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_axpy(gr%mesh%np, real(zfact), hpsi1%zpsi(:, idim, ist, ik), psi%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      else
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_axpy(gr%mesh%np, zfact, hpsi1%zpsi(:, idim, ist, ik), psi%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      end if

      if(i .ne. te%exp_order) then
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_copy(gr%mesh%np, hpsi1%zpsi(:, idim, ist, ik), psi1%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      end if

    end do
    ! End of Taylor expansion loop.

    call states_end(psi1)
    call states_end(hpsi1)

    if(present(order)) order = te%exp_order * psi%d%nik * psi%nst ! This should be the correct number
    call pop_sub()
  end subroutine exponential_apply_all

end module exponential_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
