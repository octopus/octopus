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
  use derivatives_m
  use global_m
  use hardware_m
  use hamiltonian_m
  use fourier_space_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_math_m
  use parser_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use states_m
  use states_calc_m
  use types_m
  use varinfo_m

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

  integer, public, parameter ::  &
    EXP_LANCZOS            = 2,  &
    EXP_TAYLOR             = 3,  &
    EXP_CHEBYSHEV          = 4

  type exponential_t
    integer     :: exp_method  ! which method is used to apply the exponential
    FLOAT       :: lanczos_tol ! tolerance for the Lanczos method
    integer     :: exp_order   ! order to which the propagator is expanded
  end type exponential_t

contains

  ! ---------------------------------------------------------
  subroutine exponential_init(te, der)
    type(exponential_t), intent(out) :: te
    type(derivatives_t), intent(in)  :: der

    PUSH_SUB(exponential_init)

    !%Variable TDExponentialMethod
    !%Type integer
    !%Default taylor
    !%Section Time-Dependent::Propagation
    !%Description
    !% Method used to numerically calculate the exponential of the Hamiltonian,
    !% a core part of the full algorithm used to approximate the evolution
    !% operator, specified through the variable <tt>TDPropagator</tt>.
    !% In the case of using the Magnus method, described below, the action of the exponential
    !% of the Magnus operator is also calculated through the algorithm specified
    !% by this variable.
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
    !% There exists a closed analytic form for the coefficients of the exponential in terms
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
    call parse_integer(datasets_check('TDExponentialMethod'), EXP_TAYLOR, te%exp_method)

    select case(te%exp_method)
    case(EXP_TAYLOR)
    case(EXP_CHEBYSHEV)
    case(EXP_LANCZOS)
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
      call parse_float(datasets_check('TDLanczosTol'), CNST(1e-5), te%lanczos_tol)
      if (te%lanczos_tol <= M_ZERO) call input_error('TDLanczosTol')

    case default
      call input_error('TDExponentialMethod')
    end select
    call messages_print_var_option(stdout, 'TDExponentialMethod', te%exp_method)

    if(te%exp_method==EXP_TAYLOR.or.te%exp_method==EXP_CHEBYSHEV.or.te%exp_method==EXP_LANCZOS) then
      !%Variable TDExpOrder
      !%Type integer
      !%Default 4
      !%Section Time-Dependent::Propagation
      !%Description
      !% For <tt>TDExponentialMethod</tt> = <tt>standard</tt> or <tt>chebyshev</tt>, 
      !% the order to which the exponential is expanded. For the Lanczos approximation, 
      !% it is the Lanczos-subspace dimension.
      !%End
      call parse_integer(datasets_check('TDExpOrder'), 4, te%exp_order)
      if (te%exp_order < 2) call input_error('TDExpOrder')

    end if

    POP_SUB(exponential_init)
  end subroutine exponential_init

  ! ---------------------------------------------------------
  subroutine exponential_end(te)
    type(exponential_t), intent(inout) :: te

    PUSH_SUB(exponential_end)

    POP_SUB(exponential_end)
  end subroutine exponential_end

  ! ---------------------------------------------------------
  subroutine exponential_copy(teo, tei)
    type(exponential_t), intent(inout) :: teo
    type(exponential_t), intent(in)    :: tei

    PUSH_SUB(exponential_copy)

    teo%exp_method  = tei%exp_method
    teo%lanczos_tol = tei%lanczos_tol
    teo%exp_order   = tei%exp_order

    POP_SUB(exponential_copy)
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
  ! If the Hamiltonian contains an inhomogeneous term, the operation is:
  !
  ! exp{-i*deltat*hm(t)}|zpsi> + deltat*phi{-i*deltat*hm(t)}|zpsi>  <-- |zpsi>
  !
  ! where:
  !
  ! phi(x) = (e^x - 1)/x
  ! ---------------------------------------------------------
  subroutine exponential_apply(te, der, hm, zpsi, ist, ik, deltat, time, order, vmagnus, imag_time)
    type(exponential_t), intent(inout) :: te
    type(derivatives_t), intent(in)    :: der
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: ist
    integer,             intent(in)    :: ik
    CMPLX,               intent(inout) :: zpsi(:, :)
    FLOAT,               intent(in)    :: deltat
    FLOAT,               intent(in)    :: time
    integer, optional,   intent(inout) :: order
    FLOAT,   optional,   intent(in)    :: vmagnus(der%mesh%np, hm%d%nspin, 2)
    logical, optional,   intent(in)    :: imag_time

    CMPLX   :: timestep
    logical :: apply_magnus
    type(profile_t), save :: exp_prof

    PUSH_SUB(exponential_apply)
    call profiling_in(exp_prof, "EXPONENTIAL")

    ! The only method that is currently taking care of the presence of an inhomogeneous
    ! term is the Lanczos expansion.
    ! However, I disconnect this check, because this routine is sometimes called in an
    ! auxiliary way, in order to compute (1-hm(t)*deltat)|zpsi>.
    ! This should be cleaned up.
    !_ASSERT(.not.(hamiltonian_inh_term(hm) .and. (te%exp_method .ne. EXP_LANCZOS)))

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    ! If we want to use imaginary time, timestep = i*deltat
    ! Otherwise, timestep is simply equal to deltat.
    timestep = cmplx(deltat, M_ZERO)
    if(present(imag_time)) then
      if(imag_time) then
        select case(te%exp_method)
          case(EXP_TAYLOR, EXP_LANCZOS)
            timestep = M_zI*deltat
          case default
            write(message(1), '(a)') &
              'Imaginary time evolution can only be performed with the Lanczos'
            write(message(2), '(a)') &
              'exponentiation scheme ("TDExponentialMethod = lanczos") or with the'
            write(message(3), '(a)') &
              'Taylor expansion ("TDExponentialMethod = taylor") method.'
            call messages_fatal(3)
        end select
      end if
    end if

    select case(te%exp_method)
    case(EXP_TAYLOR)
      call taylor_series
    case(EXP_LANCZOS)
      call lanczos
    case(EXP_CHEBYSHEV)
      call cheby
    end select

    call profiling_out(exp_prof)
    POP_SUB(exponential_apply)

  contains

    ! ---------------------------------------------------------
    subroutine operate(psi, oppsi)
      CMPLX, intent(inout) :: psi(:, :)
      CMPLX, intent(inout) :: oppsi(:, :)

      PUSH_SUB(exponential_apply.operate)

      if(apply_magnus) then
        call zmagnus(hm, der, psi, oppsi, ik, vmagnus)
      else
        call zhamiltonian_apply(hm, der, psi, oppsi, ist, ik, time)
      end if

      POP_SUB(exponential_apply.operate)
    end subroutine operate
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine taylor_series()
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:), hzpsi1(:,:)
      integer :: i, idim
      logical :: zfact_is_real

      PUSH_SUB(exponential_apply.taylor_series)

      SAFE_ALLOCATE(zpsi1 (1:der%mesh%np_part, 1:hm%d%dim))
      SAFE_ALLOCATE(hzpsi1(1:der%mesh%np,      1:hm%d%dim))

      zfact = M_z1
      zfact_is_real = .true.

      do idim = 1, hm%d%dim
        call lalg_copy(der%mesh%np, zpsi(:, idim), zpsi1(:, idim))
      end do

      do i = 1, te%exp_order
        zfact = zfact*(-M_zI*timestep)/i
        zfact_is_real = .not. zfact_is_real
        
        call operate(zpsi1, hzpsi1)

        if(zfact_is_real) then
          do idim = 1, hm%d%dim
            call lalg_axpy(der%mesh%np, real(zfact), hzpsi1(:, idim), zpsi(:, idim))
          end do
        else
          do idim = 1, hm%d%dim
            call lalg_axpy(der%mesh%np, zfact, hzpsi1(:, idim), zpsi(:, idim))
          end do
        end if

        if(i .ne. te%exp_order) then
          do idim = 1, hm%d%dim
            call lalg_copy(der%mesh%np, hzpsi1(:, idim), zpsi1(:, idim))
          end do
        end if

      end do
      SAFE_DEALLOCATE_A(zpsi1)
      SAFE_DEALLOCATE_A(hzpsi1)

      if(present(order)) order = te%exp_order

      POP_SUB(exponential_apply.taylor_series)
    end subroutine taylor_series
    ! ---------------------------------------------------------


    ! ---------------------------------------------------------
    subroutine cheby()
      ! Calculates the exponential of the Hamiltonian through an expansion in 
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

      PUSH_SUB(exponential_apply.cheby)

      SAFE_ALLOCATE(zpsi1(1:der%mesh%np_part, 1:hm%d%dim, 0:2))
      zpsi1 = M_z0
      do j = te%exp_order - 1, 0, -1
        do idim = 1, hm%d%dim
          call lalg_copy(der%mesh%np, zpsi1(:, idim, 1), zpsi1(:, idim, 2))
          call lalg_copy(der%mesh%np, zpsi1(:, idim, 0), zpsi1(:, idim, 1))
        end do

        call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
        zfact = 2*(-M_zI)**j*loct_bessel(j, hm%spectral_half_span*deltat)

        do idim = 1, hm%d%dim
          call lalg_axpy(der%mesh%np, -hm%spectral_middle_point, zpsi1(:, idim, 1), &
            zpsi1(:, idim, 0))
          call lalg_scal(der%mesh%np, M_TWO/hm%spectral_half_span, zpsi1(:, idim, 0))
          call lalg_axpy(der%mesh%np, zfact, zpsi(:, idim), zpsi1(:, idim, 0))
          call lalg_axpy(der%mesh%np, -M_ONE, zpsi1(:, idim, 2),  zpsi1(:, idim, 0))
        end do
      end do

      zpsi(:, :) = M_HALF*(zpsi1(:, :, 0) - zpsi1(:, :, 2))
      do idim = 1, hm%d%dim
        call lalg_scal(der%mesh%np, exp(-M_zI*hm%spectral_middle_point*deltat), zpsi(:, idim))
      end do
      SAFE_DEALLOCATE_A(zpsi1)

      if(present(order)) order = te%exp_order

      POP_SUB(exponential_apply.cheby)
    end subroutine cheby
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    subroutine lanczos
      integer ::  iter, l, idim
      CMPLX, allocatable :: hamilt(:,:), v(:,:,:), expo(:,:), tmp(:, :), psi(:, :)
      FLOAT :: beta, res, tol !, nrm
      CMPLX :: pp

      PUSH_SUB(exponential_apply.lanczos)

      SAFE_ALLOCATE(     v(1:der%mesh%np, 1:hm%d%dim, 1:te%exp_order+1))
      SAFE_ALLOCATE(   tmp(1:der%mesh%np, 1:hm%d%dim))
      SAFE_ALLOCATE(hamilt(1:te%exp_order+1, 1:te%exp_order+1))
      SAFE_ALLOCATE(  expo(1:te%exp_order+1, 1:te%exp_order+1))
      SAFE_ALLOCATE(   psi(1:der%mesh%np_part, 1:hm%d%dim))

      tol    = te%lanczos_tol
      pp = deltat
      if(.not. present(imag_time)) pp = -M_zI*pp

      beta = zmf_nrm2(der%mesh, hm%d%dim, zpsi)
      ! If we have a null vector, no need to compute the action of the exponential.
      if(beta > CNST(1.0e-12)) then

        hamilt = M_z0
        expo = M_z0

        ! Normalize input vector, and put it into v(:, :, 1)
        v(1:der%mesh%np, 1:hm%d%dim, 1) = zpsi(1:der%mesh%np, 1:hm%d%dim)/beta

        ! This is the Lanczos loop...
        do iter = 1, te%exp_order

          !copy v(:, :, n) to an array of size 1:der%mesh%np_part
          do idim = 1, hm%d%dim
            call lalg_copy(der%mesh%np, v(:, idim, iter), zpsi(:, idim))
          end do

          !to apply the Hamiltonian
          call operate(zpsi, v(:, :,  iter + 1))
        
          if(hamiltonian_hermitian(hm)) then
            l = max(1, iter - 1)
          else
            l = 1
          end if

          !orthogonalize against previous vectors
          call zstates_orthogonalization(der%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), v(:, :, iter + 1), &
            normalize = .true., overlap = hamilt(l:iter, iter), norm = hamilt(iter + 1, iter))

          call zlalg_exp(iter, pp, hamilt, expo, hamiltonian_hermitian(hm))

          res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

          if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
          if(iter > 2 .and. res < tol) exit
        end do

        if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
          write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
          call messages_warning(1)
        end if

        ! zpsi = nrm * V * expo(1:iter, 1) = nrm * V * expo * V^(T) * zpsi
        call lalg_gemv(der%mesh%np, hm%d%dim, iter, M_z1*beta, v, expo(1:iter, 1), M_z0, tmp)

        do idim = 1, hm%d%dim
          call lalg_copy(der%mesh%np, tmp(:, idim), zpsi(:, idim))
        end do

      end if

      ! We have an inhomogeneous term.
      if( hamiltonian_inh_term(hm) ) then

        call states_get_state(hm%inh_st, der%mesh, ist, ik, v(:, :, 1))
        beta = zmf_nrm2(der%mesh, hm%d%dim, v(:, :, 1))

        if(beta > CNST(1.0e-12)) then

          hamilt = M_z0
          expo = M_z0

          v(1:der%mesh%np, 1:hm%d%dim, 1) = v(1:der%mesh%np, 1:hm%d%dim, 1)/beta

          psi = M_z0
          ! This is the Lanczos loop...
          do iter = 1, te%exp_order
            !copy v(:, :, n) to an array of size 1:der%mesh%np_part
            do idim = 1, hm%d%dim
              call lalg_copy(der%mesh%np, v(:, idim, iter), psi(:, idim))
            end do

            !to apply the Hamiltonian
            call operate(psi, v(:, :, iter + 1))
  

            if(hamiltonian_hermitian(hm)) then
              l = max(1, iter - 1)
            else
              l = 1
            end if

            !orthogonalize against previous vectors
            call zstates_orthogonalization(der%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), &
              v(:, :, iter + 1), normalize = .true., overlap = hamilt(l:iter, iter), &
              norm = hamilt(iter + 1, iter))

            call zlalg_phi(iter, pp, hamilt, expo, hamiltonian_hermitian(hm))
 
            res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

            if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
            if(iter > 2 .and. res < tol) exit
          end do

          if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
            write(message(1),'(a,es8.2)') 'Lanczos exponential expansion did not converge: ', res
            call messages_warning(1)
          end if

          call lalg_gemv(der%mesh%np, hm%d%dim, iter, M_z1*beta, v, expo(1:iter, 1), M_z0, tmp)

          do idim = 1, hm%d%dim
            call lalg_copy(der%mesh%np, tmp(:, idim), psi(:, idim))
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

      POP_SUB(exponential_apply.lanczos)
    end subroutine lanczos
    ! ---------------------------------------------------------

  end subroutine exponential_apply

  subroutine exponential_apply_batch(te, der, hm, psib, ik, deltat, time)
    type(exponential_t), intent(inout) :: te
    type(derivatives_t), intent(inout) :: der
    type(hamiltonian_t), intent(inout) :: hm
    integer,             intent(in)    :: ik
    type(batch_t),       intent(inout) :: psib
    FLOAT,               intent(in)    :: deltat
    FLOAT,               intent(in)    :: time
    
    integer :: ii, ist
    CMPLX, pointer :: psi(:, :)

    PUSH_SUB(exponential_apply_batch)

    ASSERT(batch_type(psib) == TYPE_CMPLX)

    if (te%exp_method == EXP_TAYLOR) then 
      call taylor_series_batch
    else
      
      do ii = 1, psib%nst
        psi  => psib%states(ii)%zpsi
        ist  =  psib%states(ii)%ist
        
        call exponential_apply(te, der, hm, psi, ist, ik, deltat, time)
      end do

    end if
    
    POP_SUB(exponential_apply_batch)

  contains
    
    subroutine taylor_series_batch()
      CMPLX :: zfact
      CMPLX, allocatable :: psi1(:, :, :), hpsi1(:, :, :)
      integer :: iter
      logical :: zfact_is_real
      integer :: st_start, st_end
      type(batch_t) :: psi1b, hpsi1b
      type(profile_t), save :: prof

      PUSH_SUB(exponential_apply_batch.taylor_series_batch)
      call profiling_in(prof, "EXP_TAYLOR_BATCH")

      SAFE_ALLOCATE(psi1 (1:der%mesh%np_part, 1:hm%d%dim, 1:psib%nst))
      SAFE_ALLOCATE(hpsi1(1:der%mesh%np, 1:hm%d%dim, 1:psib%nst))

      st_start = psib%states(1)%ist
      st_end = psib%states(psib%nst)%ist

      zfact = M_z1
      zfact_is_real = .true.

      call batch_init(psi1b, hm%d%dim, st_start, st_end, psi1)
      call batch_init(hpsi1b, hm%d%dim, st_start, st_end, hpsi1)

      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_pack(psib)
        call batch_pack(psi1b, copy = .false.)
        call batch_pack(hpsi1b, copy = .false.)
      end if
      
      call batch_copy_data(der%mesh%np, psib, psi1b)

      do iter = 1, te%exp_order
        zfact = zfact*(-M_zI*deltat)/iter
        zfact_is_real = .not. zfact_is_real

        ! FIXME: need a test here for runaway exponential, e.g. for too large dt.
        !  in runaway case the problem is really hard to trace back: the positions
        !  go haywire on the first step of dynamics (often NaN) and with debugging options
        !  the code stops in ZAXPY below without saying why.

        call zhamiltonian_apply_batch(hm, der, psi1b, hpsi1b, ik, time)

        if(zfact_is_real) then
          call batch_axpy(der%mesh%np, real(zfact, REAL_PRECISION), hpsi1b, psib)
        else
          call batch_axpy(der%mesh%np, zfact, hpsi1b, psib)
        end if

        if(iter /= te%exp_order) call batch_copy_data(der%mesh%np, hpsi1b, psi1b)

      end do

      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_unpack(psi1b, copy = .false.)
        call batch_unpack(hpsi1b, copy = .false.)
        call batch_unpack(psib)
      end if

      call batch_end(hpsi1b)
      call batch_end(psi1b)

      call profiling_count_operations(psib%nst*hm%d%dim*dble(der%mesh%np)*te%exp_order*CNST(6.0))

      SAFE_DEALLOCATE_A(psi1)
      SAFE_DEALLOCATE_A(hpsi1)
      
      call profiling_out(prof)
      POP_SUB(exponential_apply_batch.taylor_series_batch)
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
  subroutine exponential_apply_all(te, der, hm, psi, deltat, t, order, vmagnus, imag_time)
    type(exponential_t), intent(inout) :: te
    type(derivatives_t), intent(inout) :: der
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: psi
    FLOAT,               intent(in)    :: deltat, t
    integer, optional,   intent(inout) :: order
    FLOAT,   optional,   intent(in)    :: vmagnus(der%mesh%np, hm%d%nspin, 2)
    logical, optional,   intent(in)    :: imag_time

    integer :: ik, ist
    CMPLX   :: timestep
    logical :: apply_magnus
    CMPLX :: zfact
    integer :: i, idim
    logical :: zfact_is_real

    type(states_t) :: psi1, hpsi1

    PUSH_SUB(exponential_apply_all)

    ASSERT(te%exp_method .eq. EXP_TAYLOR)

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
        call zhamiltonian_apply_all(hm, der, psi, hpsi1, t)
      else
        call zhamiltonian_apply_all(hm, der, psi1, hpsi1, t)
      end if

      if(zfact_is_real) then
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_axpy(der%mesh%np, real(zfact), hpsi1%zpsi(:, idim, ist, ik), psi%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      else
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_axpy(der%mesh%np, zfact, hpsi1%zpsi(:, idim, ist, ik), psi%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      end if

      if(i .ne. te%exp_order) then
        do ik = psi%d%kpt%start, psi%d%kpt%end
          do ist = psi%st_start, psi%st_end
            do idim = 1, hm%d%dim
              call lalg_copy(der%mesh%np, hpsi1%zpsi(:, idim, ist, ik), psi1%zpsi(:, idim, ist, ik))
            end do
          end do
        end do
      end if

    end do
    ! End of Taylor expansion loop.

    call states_end(psi1)
    call states_end(hpsi1)

    if(present(order)) order = te%exp_order*psi%d%nik*psi%nst ! This should be the correct number

    POP_SUB(exponential_apply_all)
  end subroutine exponential_apply_all

end module exponential_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
