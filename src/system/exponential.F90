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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module exponential_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use derivatives_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use loct_math_oct_m
  use mesh_function_oct_m
  use mesh_batch_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use types_oct_m
  use xc_oct_m

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
    private
    integer, public :: exp_method  !< which method is used to apply the exponential
    FLOAT           :: lanczos_tol !< tolerance for the Lanczos method
    integer, public :: exp_order   !< order to which the propagator is expanded
    integer         :: arnoldi_gs  !< Orthogonalization scheme used for Arnoldi
    integer         :: tmp_nst, tmp_nst_linear
  end type exponential_t

contains

  ! ---------------------------------------------------------
  subroutine exponential_init(te, namespace)
    type(exponential_t), intent(out) :: te
    type(namespace_t),   intent(in)  :: namespace
    
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
    !% (see M. Hochbruck and C. Lubich, <i>SIAM J. Numer. Anal.</i> <b>34</b>, 1911 (1997) for details). 
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
    !% <math>\exp_{\rm STD} (-i\delta t H) = \sum_{i=0}^{k} {(-i\delta t)^i\over{i!}} H^i.</math>
    !%
    !% The order <i>k</i> is determined by variable <tt>TDExpOrder</tt>.
    !% Some numerical considerations from <a href=http://www.phys.washington.edu/~bertsch/num3.ps>
    !% Jeff Giansiracusa and George F. Bertsch</a>
    !% suggest the 4th order as especially suitable and stable.
    !%Option chebyshev 4
    !% In principle, the Chebyshev expansion
    !% of the exponential represents it more accurately than the canonical or standard expansion. 
    !% As in the latter case, <tt>TDExpOrder</tt> determines the order of the expansion.
    !%
    !% There exists a closed analytic form for the coefficients of the exponential in terms
    !% of Chebyshev polynomials:
    !%
    !% <math>\exp_{\rm CHEB} \left( -i\delta t H \right) = \sum_{k=0}^{\infty} (2-\delta_{k0})(-i)^{k}J_k(\delta t) T_k(H),</math>
    !%
    !% where <math>J_k</math> are the Bessel functions of the first kind, and H has to be previously
    !% scaled to <math>[-1,1]</math>.
    !% See H. Tal-Ezer and R. Kosloff, <i>J. Chem. Phys.</i> <b>81</b>,
    !% 3967 (1984); R. Kosloff, <i>Annu. Rev. Phys. Chem.</i> <b>45</b>, 145 (1994);
    !% C. W. Clenshaw, <i>MTAC</i> <b>9</b>, 118 (1955).
    !%End
    call parse_variable(namespace, 'TDExponentialMethod', EXP_TAYLOR, te%exp_method)

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
      call parse_variable(namespace, 'TDLanczosTol', CNST(1e-5), te%lanczos_tol)
      if (te%lanczos_tol <= M_ZERO) call messages_input_error('TDLanczosTol')

    case default
      call messages_input_error('TDExponentialMethod')
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
      call parse_variable(namespace, 'TDExpOrder', DEFAULT__TDEXPORDER, te%exp_order)
      if (te%exp_order < 2) call messages_input_error('TDExpOrder')

    end if

    te%arnoldi_gs = OPTION__ARNOLDIORTHOGONALIZATION__CGS
    if(te%exp_method == EXP_LANCZOS) then
      !%Variable ArnoldiOrthogonalization
      !%Type integer
      !%Section Time-Dependent::Propagation
      !%Description
      !% The orthogonalization method used for the Arnoldi procedure.
      !% Only for TDExponentialMethod = lanczos. 
      !%Option cgs 3
      !% Classical Gram-Schmidt (CGS) orthogonalization.
      !% The algorithm is defined in Giraud et al., Computers and Mathematics with Applications 50, 1069 (2005).
      !%Option drcgs 5
      !% Classical Gram-Schmidt orthogonalization with double-step reorthogonalization.
      !% The algorithm is taken from Giraud et al., Computers and Mathematics with Applications 50, 1069 (2005). 
      !% According to this reference, this is much more precise than CGS or MGS algorithms.
      !%End
      call parse_variable(namespace, 'ArnoldiOrthogonalization', OPTION__ARNOLDIORTHOGONALIZATION__CGS, &
                              te%arnoldi_gs)
    end if

    te%tmp_nst = -1
    te%tmp_nst_linear = -1

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
    teo%arnoldi_gs  = tei%arnoldi_gs 

    POP_SUB(exponential_copy)
  end subroutine exponential_copy

  ! ---------------------------------------------------------
  !> This routine performs the operation:
  !! \f[
  !! \exp{-i*\Delta t*hm(t)}|\psi_z>  <-- |\psi_z>
  !! \f]
  !! If imag_time is present and is set to true, it performs instead:
  !! \f[
  !! \exp{ \Delta t*hm(t)}|\psi_z>  <-- |\psi_z>
  !! \f]
  !! If the Hamiltonian contains an inhomogeneous term, the operation is:
  !! \f[
  !! \exp{-i*\Delta t*hm(t)}|\psi_z> + \Delta t*\phi{-i*\Delta t*hm(t)}|\psi_z>  <-- |\psi_z>
  !! \f]
  !! where:
  !! \f[
  !! \phi(x) = (e^x - 1)/x
  !! \f]
  ! ---------------------------------------------------------
  subroutine exponential_apply(te, der, hm, psolver, zpsi, ist, ik, deltat, order, vmagnus, imag_time)
    type(exponential_t),      intent(inout) :: te
    type(derivatives_t),      intent(in)    :: der
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(poisson_t),          intent(in)    :: psolver
    integer,                  intent(in)    :: ist
    integer,                  intent(in)    :: ik
    CMPLX,                    intent(inout) :: zpsi(:, :)
    FLOAT,                    intent(in)    :: deltat
    integer, optional,        intent(inout) :: order
    FLOAT,   optional,        intent(in)    :: vmagnus(der%mesh%np, hm%d%nspin, 2)
    logical, optional,        intent(in)    :: imag_time

    CMPLX   :: timestep
    logical :: apply_magnus, phase_correction
    type(profile_t), save :: exp_prof

    PUSH_SUB(exponential_apply)
    call profiling_in(exp_prof, "EXPONENTIAL")

    ! The only method that is currently taking care of the presence of an inhomogeneous
    ! term is the Lanczos expansion.
    ! However, I disconnect this check, because this routine is sometimes called in an
    ! auxiliary way, in order to compute (1-hm(t)*deltat)|zpsi>.
    ! This should be cleaned up.
    !_ASSERT(.not.(hamiltonian_elec_inh_term(hm) .and. (te%exp_method /= EXP_LANCZOS)))

    apply_magnus = .false.
    if(present(vmagnus)) apply_magnus = .true.

    phase_correction = .false.
    if(associated(hm%hm_base%phase)) phase_correction = .true.
    if(accel_is_enabled()) phase_correction = .false.

    ! If we want to use imaginary time, timestep = i*deltat
    ! Otherwise, timestep is simply equal to deltat.
    timestep = TOCMPLX(deltat, M_ZERO)

    if(present(imag_time)) then
      if(imag_time) then
        select case(te%exp_method)
          case(EXP_TAYLOR, EXP_LANCZOS)
            timestep = M_zI*deltat
          case default
            write(message(1), '(a)') &
              'Imaginary  time evolution can only be performed with the Lanczos'
            write(message(2), '(a)') &
              'exponentiation scheme ("TDExponentialMethod = lanczos") or with the'
            write(message(3), '(a)') &
              'Taylor expansion ("TDExponentialMethod = taylor") method.'
            call messages_fatal(3)
        end select
      end if
    end if

   !We apply the phase only to np points, and the phase for the np+1 to np_part points
   !will be treated as a phase correction in the Hamiltonian
    if(phase_correction) then
      call states_elec_set_phase(hm%d, zpsi, hm%hm_base%phase(1:der%mesh%np, ik), der%mesh%np, .false.)
    end if

    select case(te%exp_method)
    case(EXP_TAYLOR)
      call taylor_series()
    case(EXP_LANCZOS)
      call lanczos()
    case(EXP_CHEBYSHEV)
      call cheby()
    end select

    if(phase_correction) then
      call states_elec_set_phase(hm%d, zpsi, hm%hm_base%phase(1:der%mesh%np, ik), der%mesh%np, .true.)
    end if


    call profiling_out(exp_prof)
    POP_SUB(exponential_apply)

  contains

    ! ---------------------------------------------------------
    subroutine operate(psi, oppsi)
      CMPLX,   intent(inout) :: psi(:, :)
      CMPLX,   intent(inout) :: oppsi(:, :)

      PUSH_SUB(exponential_apply.operate)

      if(apply_magnus) then
        call zmagnus(hm, der, psolver, psi, oppsi, ik, vmagnus, set_phase = .not.phase_correction)
        else
        call zhamiltonian_elec_apply(hm, der, psolver, psi, oppsi, ist, ik, set_phase = .not.phase_correction)
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
            call lalg_axpy(der%mesh%np, real(zfact, REAL_PRECISION), hzpsi1(:, idim), zpsi(:, idim))
          end do
        else
          do idim = 1, hm%d%dim
            call lalg_axpy(der%mesh%np, zfact, hzpsi1(:, idim), zpsi(:, idim))
          end do
        end if

        if(i /= te%exp_order) then
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
    !> Calculates the exponential of the Hamiltonian through an expansion in 
    !! Chebyshev polynomials.
    !!
    !! For that purposes it uses the closed form of the coefficients[1] and Clenshaw-Gordons[2]
    !! recursive algorithm.
    !! [1] H. Tal-Ezer and R. Kosloff, J. Chem. Phys 81, 3967 (1984).
    !! [2] C. W. Clenshaw, MTAC 9, 118 (1955).
    !! Since I don't have access to MTAC, I copied next Maple algorithm from Dr. F. G. Lether's
    !! (University of Georgia) homepage: (www.math.uga.edu/~fglether):
    !! \verbatim
    !! twot := t + t; u0 := 0; u1 := 0;
    !!  for k from n to 0 by -1 do
    !!    u2 := u1; u1 := u0;
    !!    u0 := twot*u1 - u2 + c[k];
    !!  od;
    !!  ChebySum := 0.5*(u0 - u2);
    !! \endverbatim
    subroutine cheby()
      integer :: j, idim
      CMPLX :: zfact
      CMPLX, allocatable :: zpsi1(:,:,:)

      integer :: np

      PUSH_SUB(exponential_apply.cheby)

      np = der%mesh%np

      !TODO: We can save memory here as we only need one array of size np_part and not 4
      SAFE_ALLOCATE(zpsi1(1:der%mesh%np_part, 1:hm%d%dim, 0:2))
      zpsi1 = M_z0
      do j = te%exp_order - 1, 0, -1
        do idim = 1, hm%d%dim
          call lalg_copy(der%mesh%np, zpsi1(1:np, idim, 1), zpsi1(1:np, idim, 2))
          call lalg_copy(der%mesh%np, zpsi1(1:np, idim, 0), zpsi1(1:np, idim, 1))
        end do

        call operate(zpsi1(:, :, 1), zpsi1(:, :, 0))
        zfact = 2*(-M_zI)**j*loct_bessel(j, hm%spectral_half_span*deltat)

        do idim = 1, hm%d%dim
          call lalg_axpy(np, -hm%spectral_middle_point, zpsi1(1:np, idim, 1), &
            zpsi1(1:np, idim, 0))
          call lalg_scal(np, M_TWO/hm%spectral_half_span, zpsi1(1:np, idim, 0))
          call lalg_axpy(np, zfact, zpsi(:, idim), zpsi1(1:np, idim, 0))
          call lalg_axpy(der%mesh%np, -M_ONE, zpsi1(1:np, idim, 2),  zpsi1(1:np, idim, 0))
        end do
      end do

      zpsi(1:np, 1:hm%d%dim) = M_HALF*(zpsi1(1:np, 1:hm%d%dim, 0) - zpsi1(1:np, 1:hm%d%dim, 2))
      do idim = 1, hm%d%dim
        call lalg_scal(np, exp(-M_zI*hm%spectral_middle_point*deltat), zpsi(1:np, idim))
      end do
      SAFE_DEALLOCATE_A(zpsi1)

      if(present(order)) order = te%exp_order

      POP_SUB(exponential_apply.cheby)
    end subroutine cheby
    ! ---------------------------------------------------------

    ! ---------------------------------------------------------
    !TODO: Add a reference
    subroutine lanczos()
      integer ::  iter, l, idim
      CMPLX, allocatable :: hamilt(:,:), v(:,:,:), expo(:,:), psi(:, :)
      FLOAT :: beta, res, tol !, nrm
      CMPLX :: pp

      PUSH_SUB(exponential_apply.lanczos)

      SAFE_ALLOCATE(     v(1:der%mesh%np, 1:hm%d%dim, 1:te%exp_order+1))
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
        
          if(hm%is_hermitian()) then
            l = max(1, iter - 1)
          else
            l = 1
          end if

          !orthogonalize against previous vectors
          call zstates_elec_orthogonalization(der%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), v(:, :, iter + 1), &
            normalize = .false., overlap = hamilt(l:iter, iter), norm = hamilt(iter + 1, iter), &
            gs_scheme = te%arnoldi_gs)

          call zlalg_exp(iter, pp, hamilt, expo, hm%is_hermitian())

          res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

          if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
          !We normalize only if the norm is non-zero
          ! see http://www.netlib.org/utk/people/JackDongarra/etemplates/node216.html#alg:arn0
          do idim = 1, hm%d%dim
            call lalg_scal(der%mesh%np, M_ONE / hamilt(iter + 1, iter), v(:, idim, iter+1))
          end do
           
          if(iter > 3 .and. res < tol) exit
        end do

        if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
          write(message(1),'(a,es9.2)') 'Lanczos exponential expansion did not converge: ', res
          call messages_warning(1)
        end if

        ! zpsi = nrm * V * expo(1:iter, 1) = nrm * V * expo * V^(T) * zpsi
        do idim = 1, hm%d%dim
          call blas_gemv('N', der%mesh%np, iter, M_z1*beta, v(1,idim,1), der%mesh%np*hm%d%dim, expo(1,1), 1, M_z0, zpsi(1,idim), 1)
        end do

      end if

      ! We have an inhomogeneous term.
      if( hamiltonian_elec_inh_term(hm) ) then

        call states_elec_get_state(hm%inh_st, der%mesh, ist, ik, v(:, :, 1))
        beta = zmf_nrm2(der%mesh, hm%d%dim, v(:, :, 1))

        if(beta > CNST(1.0e-12)) then

          hamilt = M_z0
          expo = M_z0

          v(1:der%mesh%np, 1:hm%d%dim, 1) = v(1:der%mesh%np, 1:hm%d%dim, 1)/beta

          ! This is the Lanczos loop...
          do iter = 1, te%exp_order
            !copy v(:, :, n) to an array of size 1:der%mesh%np_part
            do idim = 1, hm%d%dim
              call lalg_copy(der%mesh%np, v(:, idim, iter), psi(:, idim))
            end do

            !to apply the Hamiltonian
            call operate(psi, v(:, :, iter + 1))
  

            if(hm%is_hermitian()) then
              l = max(1, iter - 1)
            else
              l = 1
            end if

            !orthogonalize against previous vectors
            call zstates_elec_orthogonalization(der%mesh, iter - l + 1, hm%d%dim, v(:, :, l:iter), &
              v(:, :, iter + 1), normalize = .true., overlap = hamilt(l:iter, iter), &
              norm = hamilt(iter + 1, iter), gs_scheme = te%arnoldi_gs)

            call zlalg_phi(iter, pp, hamilt, expo, hm%is_hermitian())
 
            res = abs(hamilt(iter + 1, iter)*abs(expo(iter, 1)))

            if(abs(hamilt(iter + 1, iter)) < CNST(1.0e4)*M_EPSILON) exit ! "Happy breakdown"
            if(iter > 3 .and. res < tol) exit
          end do

          if(res > tol) then ! Here one should consider the possibility of the happy breakdown.
            write(message(1),'(a,es9.2)') 'Lanczos exponential expansion did not converge: ', res
            call messages_warning(1)
          end if

          do idim = 1, hm%d%dim
            call blas_gemv('N', der%mesh%np, iter, deltat*M_z1*beta, v(1,idim,1), &
                           der%mesh%np*hm%d%dim, expo(1,1), 1, M_z1, zpsi(1,idim), 1)
          end do

        end if
      end if


      if(present(order)) order = te%exp_order
      SAFE_DEALLOCATE_A(v)
      SAFE_DEALLOCATE_A(hamilt)
      SAFE_DEALLOCATE_A(expo)
      SAFE_DEALLOCATE_A(psi)

      POP_SUB(exponential_apply.lanczos)
    end subroutine lanczos
    ! ---------------------------------------------------------

  end subroutine exponential_apply

  subroutine exponential_apply_batch(te, der, hm, psolver, psib, ik, deltat, psib2, deltat2)
    type(exponential_t),             intent(inout) :: te
    type(derivatives_t),             intent(inout) :: der
    type(hamiltonian_elec_t),        intent(inout) :: hm
    type(poisson_t),                 intent(in)    :: psolver
    integer,                         intent(in)    :: ik
    type(batch_t), target,           intent(inout) :: psib
    FLOAT,                           intent(in)    :: deltat
    type(batch_t), target, optional, intent(inout) :: psib2
    FLOAT, optional,                 intent(in)    :: deltat2
    
    integer :: ii, ist
    CMPLX, pointer :: psi(:, :), psi2(:, :)
    logical :: phase_correction

    PUSH_SUB(exponential_apply_batch)

    ASSERT(batch_type(psib) == TYPE_CMPLX)
    ASSERT(present(psib2) .eqv. present(deltat2))

    
    ! check if we only want a phase correction for the boundary points
    phase_correction = .false.
    if(associated(hm%hm_base%phase)) phase_correction = .true.
    if(accel_is_enabled()) phase_correction = .false.

    if (te%exp_method == EXP_TAYLOR .or. &
          (te%exp_method == EXP_LANCZOS.and..not. hamiltonian_elec_inh_term(hm)&
           .and..not. present(psib2))) then 
     !We apply the phase only to np points, and the phase for the np+1 to np_part points
     !will be treated as a phase correction in the Hamiltonian
      if(phase_correction) then
        call zhamiltonian_elec_base_phase(hm%hm_base, der, der%mesh%np, ik, .false., psib)
      end if

      select case(te%exp_method)
      case(EXP_TAYLOR)
        call taylor_series_batch()
      case(EXP_LANCZOS) 
        call lanczos_batch()
      end select

      if(phase_correction) then
        call zhamiltonian_elec_base_phase(hm%hm_base, der, der%mesh%np, ik, .true., psib)
        if(present(psib2)) then
          call zhamiltonian_elec_base_phase(hm%hm_base, der, der%mesh%np, ik, .true., psib2)
        end if
      end if
    else

      if(present(psib2)) call batch_copy_data(der%mesh%np, psib, psib2)

      ! only allocate for packed cases
      if (batch_status(psib) /= BATCH_NOT_PACKED) then
        SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:hm%d%dim))
      end if
      if(present(psib2)) then
        if (batch_status(psib2) /= BATCH_NOT_PACKED) then
          SAFE_ALLOCATE(psi2(1:der%mesh%np_part, 1:hm%d%dim))
        end if
      end if

      do ii = 1, psib%nst
        ist = psib%states(ii)%ist

        ! avoid copy for the unpacked case, simply set pointer
        ! -> this should be removed by having batched versions of
        ! exponential_apply
        if (batch_status(psib) /= BATCH_NOT_PACKED) then
          call batch_get_state(psib, ii, der%mesh%np, psi)
        else
          psi => psib%states(ii)%zpsi
        end if
        call exponential_apply(te, der, hm, psolver, psi, ist, ik, deltat)
        if (batch_status(psib) /= BATCH_NOT_PACKED) then
          call batch_set_state(psib, ii, der%mesh%np, psi)
        end if
        
        if(present(psib2)) then
          ! also avoid copying unpacked batches here
          if (batch_status(psib2) /= BATCH_NOT_PACKED) then
            call batch_get_state(psib2, ii, der%mesh%np, psi2)
          else
            psi2 => psib2%states(ii)%zpsi
          end if
          call exponential_apply(te, der, hm, psolver, psi2, ist, ik, deltat2)
          if (batch_status(psib2) /= BATCH_NOT_PACKED) then
            call batch_set_state(psib2, ii, der%mesh%np, psi2)
          end if
        end if
      end do

      if (batch_status(psib) /= BATCH_NOT_PACKED) then
        SAFE_DEALLOCATE_P(psi)
      end if
      if(present(psib2)) then
        if (batch_status(psib2) /= BATCH_NOT_PACKED) then
           SAFE_DEALLOCATE_P(psi2)
        end if
      end if

   end if
    
   POP_SUB(exponential_apply_batch)

  contains
    
    subroutine taylor_series_batch()
      CMPLX :: zfact, zfact2
      integer :: iter
      logical :: zfact_is_real
      type(profile_t), save :: prof
      type(batch_t) :: psi1b, hpsi1b

      PUSH_SUB(exponential_apply_batch.taylor_series_batch)
      call profiling_in(prof, "EXP_TAYLOR_BATCH")

      if(hamiltonian_elec_apply_packed(hm, der%mesh)) then
        call batch_pack(psib)
        if(present(psib2)) call batch_pack(psib2, copy = .false.)
      end if

      call batch_copy(psib, psi1b)
      call batch_copy(psib, hpsi1b)

      zfact = M_z1
      zfact2 = M_z1
      zfact_is_real = .true.
      
      if(present(psib2)) call batch_copy_data(der%mesh%np, psib, psib2)

      do iter = 1, te%exp_order
        zfact = zfact*(-M_zI*deltat)/iter
        if(present(deltat2)) zfact2 = zfact2*(-M_zI*deltat2)/iter
        zfact_is_real = .not. zfact_is_real
        ! FIXME: need a test here for runaway exponential, e.g. for too large dt.
        !  in runaway case the problem is really hard to trace back: the positions
        !  go haywire on the first step of dynamics (often NaN) and with debugging options
        !  the code stops in ZAXPY below without saying why.

        if(iter /= 1) then
          call zhamiltonian_elec_apply_batch(hm, der, psolver, psi1b, hpsi1b, ik, set_phase = .not.phase_correction)
        else
          call zhamiltonian_elec_apply_batch(hm, der, psolver, psib, hpsi1b, ik, set_phase = .not.phase_correction)
        end if
        
        if(zfact_is_real) then
          call batch_axpy(der%mesh%np, real(zfact, REAL_PRECISION), hpsi1b, psib)
          if(present(psib2)) call batch_axpy(der%mesh%np, real(zfact2, REAL_PRECISION), hpsi1b, psib2)
        else
          call batch_axpy(der%mesh%np, zfact, hpsi1b, psib)
          if(present(psib2)) call batch_axpy(der%mesh%np, zfact2, hpsi1b, psib2)
        end if

        if(iter /= te%exp_order) call batch_copy_data(der%mesh%np, hpsi1b, psi1b)

      end do

      call batch_end(psi1b)
      call batch_end(hpsi1b)
      
      if(hamiltonian_elec_apply_packed(hm, der%mesh)) then
        if(present(psib2)) call batch_unpack(psib2)
        call batch_unpack(psib)
      end if

      call profiling_count_operations(psib%nst*hm%d%dim*dble(der%mesh%np)*te%exp_order*CNST(6.0))
      
      call profiling_out(prof)
      POP_SUB(exponential_apply_batch.taylor_series_batch)

    end subroutine taylor_series_batch


    ! ---------------------------------------------------------
    !TODO: Add a reference
    subroutine lanczos_batch()

      integer ::  iter, l, idim, bind
      CMPLX, allocatable :: hamilt(:,:,:), expo(:,:,:)
      FLOAT, allocatable :: beta(:), res(:), norm(:)
      type(batch_t), allocatable :: vb(:)
      type(profile_t), save :: prof

      PUSH_SUB(exponential_apply_batch.lanczos_batch)
      call profiling_in(prof, "EXP_LANCZOS_BATCH")

      SAFE_ALLOCATE(beta(1:psib%nst))
      SAFE_ALLOCATE(res(1:psib%nst))
      SAFE_ALLOCATE(norm(1:psib%nst))
      call mesh_batch_nrm2(der%mesh, psib, beta) 

      ! If we have a null vector, no need to compute the action of the exponential.
      if(all(abs(beta) <= CNST(1.0e-12))) then
        call profiling_out(prof)
        POP_SUB(exponential_apply_batch.lanczos_batch)
        return
      end if

      if(hamiltonian_elec_apply_packed(hm, der%mesh)) then
        call batch_pack(psib)
      end if

      SAFE_ALLOCATE(vb(1:te%exp_order+1))
      do iter = 1, te%exp_order+1
        call batch_copy(psib, vb(iter))
      end do
      call batch_copy_data(der%mesh%np, psib, vb(1))
      call batch_scal(der%mesh%np, M_ONE/beta, vb(1), a_full = .false.)

      SAFE_ALLOCATE(hamilt(1:te%exp_order+1, 1:te%exp_order+1, 1:psib%nst))
      SAFE_ALLOCATE(  expo(1:te%exp_order+1, 1:te%exp_order+1, 1:psib%nst))
      hamilt = M_z0
      expo = M_z0

      ! This is the Lanczos loop...
      do iter = 1, te%exp_order

        !to apply the Hamiltonian
        call zhamiltonian_elec_apply_batch(hm, der, psolver, vb(iter), vb(iter+1), ik, set_phase = .not.phase_correction)

        if(hm%is_hermitian()) then
          l = max(1, iter - 1)
        else
          l = 1
        end if

         !orthogonalize against previous vectors
        call zmesh_batch_orthogonalization(der%mesh, iter - l + 1, vb(l:iter), vb(iter+1), &
            normalize = .false., overlap = hamilt(l:iter, iter, 1:psib%nst), norm = hamilt(iter + 1, iter, 1:psib%nst), &
            gs_scheme = te%arnoldi_gs)

        do ii = 1, psib%nst
          call zlalg_exp(iter, -M_zI*deltat, hamilt(:,:,ii), expo(:,:,ii), hm%is_hermitian())

          res(ii) = abs(hamilt(iter + 1, iter, ii)*abs(expo(iter, 1, ii)))
        end do !ii

        if(all(abs(hamilt(iter + 1, iter, :)) < CNST(1.0e4)*M_EPSILON)) exit ! "Happy breakdown"
        !We normalize only if the norm is non-zero
        ! see http://www.netlib.org/utk/people/JackDongarra/etemplates/node216.html#alg:arn0 
        norm = M_ONE
        do ist = 1, psib%nst
          if( abs(hamilt(iter + 1, iter, ist)) >= CNST(1.0e4)*M_EPSILON ) then
            norm(ist) = M_ONE / abs(hamilt(iter + 1, iter, ist))
          end if
        end do
        call batch_scal(der%mesh%np, norm, vb(iter+1), a_full = .false.)

        if(iter > 3 .and. all(res < te%lanczos_tol)) exit

      end do !iter 

      if(any(res > te%lanczos_tol)) then ! Here one should consider the possibility of the happy breakdown.
        write(message(1),'(a,es9.2)') 'Lanczos exponential expansion did not converge: ', maxval(res)
        call messages_warning(1)
      end if

      ! zpsi = nrm * V * expo(1:iter, 1) = nrm * V * expo * V^(T) * zpsi
      call batch_scal(der%mesh%np, expo(1,1,1:psib%nst), psib, a_full = .false.)
      !TODO: We should have a routine batch_gemv fro improve performances
      do ii = 2, iter
        call batch_axpy(der%mesh%np, beta(1:psib%nst)*expo(ii,1,1:psib%nst), vb(ii), psib, a_full = .false.)
        !In order to apply the two exponentials, we mush store the eigenvales and eigenvectors given by zlalg_exp
        !And to recontruct here the exp(i*dt*H) for deltat2
      end do

      do iter = 1, te%exp_order+1
        call batch_end(vb(iter))
      end do

      SAFE_DEALLOCATE_A(hamilt)
      SAFE_DEALLOCATE_A(expo)
      SAFE_DEALLOCATE_A(beta)
      SAFE_DEALLOCATE_A(res)
      SAFE_DEALLOCATE_A(norm)

      if(hamiltonian_elec_apply_packed(hm, der%mesh)) then
        call batch_unpack(psib)
      end if

      call profiling_out(prof)

      POP_SUB(exponential_apply_batch.lanczos_batch)
    end subroutine lanczos_batch

  end subroutine exponential_apply_batch

  ! ---------------------------------------------------------
  !> Note that this routine not only computes the exponential, but
  !! also an extra term if there is a inhomogeneous term in the
  !! Hamiltonian hm.
  subroutine exponential_apply_all(te, der, hm, psolver, xc, st, deltat, order)
    type(exponential_t),      intent(inout) :: te
    type(derivatives_t),      intent(inout) :: der
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(poisson_t),          intent(in)    :: psolver
    type(xc_t),               intent(in)    :: xc
    type(states_elec_t),      intent(inout) :: st
    FLOAT,                    intent(in)    :: deltat
    integer, optional,        intent(inout) :: order

    integer :: ik, ib, i
    FLOAT :: zfact

    type(states_elec_t) :: st1, hst1

    PUSH_SUB(exponential_apply_all)

    ASSERT(te%exp_method  ==  EXP_TAYLOR)

    call states_elec_copy(st1, st)
    call states_elec_copy(hst1, st)

    zfact = M_ONE
    do i = 1, te%exp_order
      zfact = zfact * deltat / i
      
      if (i == 1) then
        call zhamiltonian_elec_apply_all(hm, xc, der, psolver, st, hst1)
      else
        call zhamiltonian_elec_apply_all(hm, xc, der, psolver, st1, hst1)
      end if

      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
            call batch_set_zero(st1%group%psib(ib, ik))
            call batch_axpy(der%mesh%np, -M_zI, hst1%group%psib(ib, ik), st1%group%psib(ib, ik))
            call batch_axpy(der%mesh%np, zfact, st1%group%psib(ib, ik), st%group%psib(ib, ik))
        end do
      end do

    end do
    ! End of Taylor expansion loop.

    call states_elec_end(st1)
    call states_elec_end(hst1)

    ! We now add the inhomogeneous part, if present.
    if(hamiltonian_elec_inh_term(hm)) then
      !write(*, *) 'Now we apply the inhomogeneous term...'

      call states_elec_copy(st1, hm%inh_st)
      call states_elec_copy(hst1, hm%inh_st)


      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call batch_axpy(der%mesh%np, deltat, st1%group%psib(ib, ik), st%group%psib(ib, ik))
        end do
      end do

      zfact = M_ONE
      do i = 1, te%exp_order
        zfact = zfact * deltat / (i+1)
      
        if (i == 1) then
          call zhamiltonian_elec_apply_all(hm, xc, der, psolver, hm%inh_st, hst1)
        else
          call zhamiltonian_elec_apply_all(hm, xc, der, psolver, st1, hst1)
        end if

        do ik = st%d%kpt%start, st%d%kpt%end
          do ib = st%group%block_start, st%group%block_end
            call batch_set_zero(st1%group%psib(ib, ik))
            call batch_axpy(der%mesh%np, -M_zI, hst1%group%psib(ib, ik), st1%group%psib(ib, ik))
            call batch_axpy(der%mesh%np, deltat * zfact, st1%group%psib(ib, ik), st%group%psib(ib, ik))
          end do
        end do

      end do

      call states_elec_end(st1)
      call states_elec_end(hst1)

    end if

    if(present(order)) order = te%exp_order*st%d%nik*st%nst ! This should be the correct number

    POP_SUB(exponential_apply_all)
  end subroutine exponential_apply_all

end module exponential_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
