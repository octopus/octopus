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

module td_rti_m
  use global_m
  use messages_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use math_m
  use mesh_m
  use cube_function_m
  use functions_m
  use mesh_function_m
  use states_m
  use v_ks_m
  use hamiltonian_m
  use external_pot_m
  use td_exp_m
  use td_exp_split_m
  use grid_m
  use varinfo_m

  implicit none

  private
  public ::                   &
    td_rti_t,                 &
    td_rti_init,              &
    td_rti_end,               &
    td_rti_run_zero_iter,     &
    td_rti_dt

  integer, parameter ::       &
    SPLIT_OPERATOR       = 0, &
    SUZUKI_TROTTER       = 1, &
    REVERSAL             = 2, &
    APP_REVERSAL         = 3, &
    EXPONENTIAL_MIDPOINT = 4, &
    MAGNUS               = 5

  FLOAT, parameter :: scf_threshold = CNST(1.0e-3)

  type td_rti_t
    integer           :: method  ! which evolution method to use
    type(td_exp_t) :: te      ! how to apply the propagator (e^{-i H \Delta t})

    FLOAT, pointer :: v_old(:, :, :) ! storage of the KS potential of previous iterations
    FLOAT, pointer :: vmagnus(:, :, :) ! auxiliary function to store the Magnus potentials.

    type(zcf_t) :: cf              ! auxiliary cube for split operator methods
  end type td_rti_t


contains

  ! ---------------------------------------------------------
  subroutine td_rti_init(gr, st, tr)
    type(grid_t),   intent(in)    :: gr
    type(states_t), intent(in)    :: st
    type(td_rti_t), intent(inout) :: tr

    !%Variable TDEvolutionMethod
    !%Type integer
    !%Default etrs
    !%Section Time Dependent::Propagation
    !%Description
    !% This variable determines which algorithm will be used to approximate
    !% the evolution operator <math>U(t+\delta t, t)</math>. That is, known
    !% <math>\psi(\tau)</math> and <math>H(\tau)</math> for <math>tau \le t</math>,
    !% calculate <math>t+\delta t</math>. Note that in general the Hamiltonian
    !% is not known at times in the interior of the interval <math>[t,t+\delta t]</math>.
    !% This is due to the self-consistent nature of the time-dependent Kohn-Sham problem:
    !% the Hamiltonian at a given time <math>\tau</math> is built from the
    !% "solution" wavefunctions at that time.
    !%
    !% Some methods, however, do require the knowledge of the Hamiltonian at some
    !% point of the interval <math>[t,t+\delta t]</math>. This problem is solved by making
    !% use of extrapolation: given a number <math>l</math> of time steps previous to time
    !% <math>t</math>, this information is used to build the Hamiltonian at arbitrary times
    !% within <math>[t,t+\delta t]</math>. To be fully precise, one should then proceed
    !% <i>self-consistently</i>: the obtained Hamiltonian at time <math>t+\delta t</math>
    !% may then be used to interpolate the Hamiltonian, and repeat the evolution
    !% algorithm with this new information. Whenever iterating the procedure does
    !% not change the solution wave-functions, the cycle is stopped. In practice,
    !% in <tt>octopus</tt> we perform a second-order extrapolation without
    !% self-consistente check, except for the first two iterations, where obviously
    !% the extrapolation is not reliable.
    !%
    !% The proliferation of methods is certainly excessive; The reason for it is that 
    !% the propagation algorithm is currently a topic of active development. We
    !% hope that in the future the optimal schemes are clearly identified. In the
    !% mean time, if you do not feel like testing, use the default choices and
    !% make sure the time step is small enough.
    !%Option split 0
    !% Split Operator (SO).
    !% This is one of the most traditional methods. It splits the full Hamiltonian
    !% into a kinetic and a potential part, performing the first in Fourier-space,
    !% and the latter in real space. The necessary transformations are performed
    !% with the FFT algorithm.
    !%
    !% <MATH>
    !%    U_{\rm SO}(t+\delta t, t) = \exp \lbrace - {i \over 2}\delta t T \rbrace
    !%       \exp \lbrace -i\delta t V^* \rbrace
    !%       \exp \lbrace - {i \over 2}\delta t T \rbrace
    !% </MATH>
    !%
    !% Since those three exponentials may be calculated exactly, one does not
    !% need to use any of the methods specified by variable <tt>TDExponentialMethod</tt>
    !% to perform them. 
    !%
    !% The key point is the specification of <math>V^*</math>. Let <math>V(t)</math> be divided into
    !% <math>V_{\rm int}(t)</math>, the "internal" potential which depends self-consistently
    !% on the density, and <math>V_{\rm ext}(t)</math>, the external potential that we know
    !% at all times since it is imposed to the system by us (e.g. a laser field):
    !% <math>V(t)=V_{\rm int}(t)+V_{\rm ext}(t)</math>. Then we define to be <math>V^*</math> to
    !% be the sum of <math>V_{\rm ext}(t+\delta t/2)</math> and the internal potential built
    !% from the wavefunctions <i>after</i> applying the right-most kinetic term
    !% of the equation, <math>\exp \lbrace -i\delta t/2 T \rbrace</math>.
    !%
    !% It may the be demonstrated that the order of the error of the algorithm is the
    !% same that the one that we would have by making use of the Exponential Midpoint Rule
    !% (EM, described below), the SO algorithm to calculate the action of the 
    !% exponential of the Hamiltonian, and a full self-consistent procedure.
    !%Option suzuki_trotter 1
    !% This is the generalization of the Suzuki-Trotter algorithm, described
    !% as one of the choices of the <tt>TDExponentialMethod</tt>,
    !% to time-dependent problem. Consult the paper by O. Sugino and M. Miyamoto,
    !% Phys. Rev. B <b>59</b>, 2579 (1999), for details.
    !%
    !% It requires of Hamiltonian extrapolations.
    !%Option etrs 2
    !% The idea is to make use of the time-reversal symmetry from the beginning:
    !%
    !% <MATH>
    !%   \exp \left(-i\delta t/2 H_{n}\right)\psi_n = exp \left(i\delta t/2 H_{n+1}\right)\psi_{n+1},
    !% </MATH>
    !%
    !% and then invert to obtain:
    !%
    !% <MATH>
    !%   \psi_{n+1} = \exp \left(-i\delta t/2 H_{n+1}\right) exp \left(-i\delta t/2 H_{n}\right)\psi_{n}.
    !% </MATH>
    !%
    !% But we need to know <math>H_{n+1}</math>, which can only be known exactly through the solution
    !% <math>\psi_{n+1}</math>. What we do is to estimate it by performing a single exponential:
    !% <math>\psi^{*}_{n+1}=\exp \left( -i\delta t H_{n} \right) \psi_n</math>, and then
    !% <math>H_{n+1} = H[\psi^{*}_{n+1}]</math>. Thus no extrapolation is performed in this case.
    !%Option aetrs 3
    !% Approximated Enforced Time-Reversal Symmetry (AETRS).
    !% A modification of previous method to make it faster.
    !% It is based on extrapolation of the time-dependent potentials. It is faster
    !% by about 40%.
    !%
    !% The only difference is the procedure to estimate @math{H_{n+1}}: in this case
    !% it is extrapolated trough a second-order polynomial by making use of the
    !% Hamiltonian at time @math{t-2\delta t}, @math{t-\delta t} and @math{t}.
    !%Option exp_mid 4
    !% Exponential Midpoint Rule (EM).
    !% This is maybe the simplest method, but is is very well grounded theretically:
    !% it is unitary (if the exponential is performed correctly) and preserves
    !% time symmetry (if the self-consistency problem is dealt with correctly).
    !% It is defined as:
    !%
    !% <MATH>
    !%   U_{\rm EM}(t+\delta t, t) = \exp \left( -i\delta t H_{t+\delta t/2}\right)\,.
    !% </MATH>
    !%Option magnus 5
    !% Magnus Expansion (M4).
    !% This is the most sophisticated approach. It is a fourth order scheme (feature
    !% that shares with the ST scheme; the other schemes are in principle second order).
    !% It is tailored for making use of very large time steps, or equivalently,
    !% dealing with problem with very high-frequency time dependence.
    !% It is still in a experimental state; we are not yet sure of when it is
    !% advantageous.
    !%End
    call loct_parse_int(check_inp('TDEvolutionMethod'), REVERSAL, tr%method)

    select case(tr%method)
    case(SPLIT_OPERATOR)
      call zcf_new(gr%m%l, tr%cf)
      call zcf_fft_init(tr%cf, gr%sb)
    case(SUZUKI_TROTTER)
      call zcf_new(gr%m%l, tr%cf)
      call zcf_fft_init(tr%cf, gr%sb)
    case(REVERSAL)
    case(APP_REVERSAL)
    case(EXPONENTIAL_MIDPOINT)
    case(MAGNUS)
      ALLOCATE(tr%vmagnus(NP, st%d%nspin, 2), NP*st%d%nspin*2)
    case default
      call input_error('TDEvolutionMethod')
    end select
    call messages_print_var_option(stdout, 'TDEvolutionMethod', tr%method)

    ! allocate memory to store the old KS potentials
    ALLOCATE(tr%v_old(NP, st%d%nspin, 0:3), NP*st%d%nspin*(3+1))
    call td_exp_init(gr, tr%te)             ! initialize propagator

  end subroutine td_rti_init


  ! ---------------------------------------------------------
  subroutine td_rti_end(tr)
    type(td_rti_t), intent(inout) :: tr

    ASSERT(associated(tr%v_old)) ! sanity check
    deallocate(tr%v_old)         ! clean ols KS potentials
    nullify(tr%v_old)

    if(tr%method == MAGNUS) then
      ASSERT(associated(tr%vmagnus))
      deallocate(tr%vmagnus); nullify(tr%vmagnus)
    end if

    if(tr%method == SUZUKI_TROTTER .or. tr%method == SPLIT_OPERATOR) then
      call zcf_free(tr%cf)
    end if

    call td_exp_end(tr%te)       ! clean propagator method
  end subroutine td_rti_end


  ! ---------------------------------------------------------
  subroutine td_rti_run_zero_iter(h, tr)
    type(hamiltonian_t), intent(in)    :: h
    type(td_rti_t),      intent(inout) :: tr

    tr%v_old(:, :, 2) = h%vhxc(:, :)
    tr%v_old(:, :, 3) = h%vhxc(:, :)
    tr%v_old(:, :, 1) = h%vhxc(:, :)

  end subroutine td_rti_run_zero_iter


  ! ---------------------------------------------------------
  subroutine td_rti_dt(ks, h, gr, st, tr, t, dt)
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(td_rti_t),      intent(inout) :: tr
    FLOAT,                  intent(in) :: t, dt

    integer :: is, iter
    FLOAT   :: d, d_max
    logical :: self_consistent
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:)

    call push_sub('td_rti.td_rti')

    self_consistent = .false.
    if(t<3*dt .and. (.not.h%ip_app)) then
      self_consistent = .true.
      ALLOCATE(zpsi1(NP_PART, st%d%dim, st%st_start:st%st_end, st%d%nik), NP_PART*st%d%dim*(st%st_end-st%st_start+1)*st%d%nik)
      zpsi1 = st%zpsi
    end if

    call lalg_copy(NP, st%d%nspin, tr%v_old(:, :, 2), tr%v_old(:, :, 3))
    call lalg_copy(NP, st%d%nspin, tr%v_old(:, :, 1), tr%v_old(:, :, 2))
    call lalg_copy(NP, st%d%nspin, h%vhxc(:, :),      tr%v_old(:, :, 1))
    call dextrapolate(2, NP, st%d%nspin, tr%v_old(:, :, 1:3), tr%v_old(:, :, 0), dt, dt)

    select case(tr%method)
    case(SPLIT_OPERATOR);       call td_rti0
    case(SUZUKI_TROTTER);       call td_rti1
    case(REVERSAL);             call td_rti2
    case(APP_REVERSAL);         call td_rti3
    case(EXPONENTIAL_MIDPOINT); call td_rti4
    case(MAGNUS);               call td_rti5
    end select

    if(self_consistent) then
      do iter = 1, 3
        call lalg_copy(NP, st%d%nspin, tr%v_old(:, :, 0), tr%v_old(:, :, 3))

        call zstates_calc_dens(st, NP, st%rho(1:NP,:))
        call zv_ks_calc(gr, ks, h, st)
        tr%v_old(1:NP, :, 0) = h%vhxc  (1:NP, :)
        h%vhxc  (1:NP, :)    = tr%v_old(1:NP, :, 1)

        d_max = M_ZERO
        ALLOCATE(dtmp(NP), NP)
        do is = 1, st%d%nspin
          dtmp(1:NP) = tr%v_old(1:NP, is, 3) - tr%v_old(1:NP, is, 0)
          d          = dmf_nrm2(gr%m, dtmp)
          if(d > d_max) d_max = d
        end do
        deallocate(dtmp)

        if(d_max < scf_threshold) exit

        st%zpsi = zpsi1
        select case(tr%method)
        case(SPLIT_OPERATOR);       call td_rti0
        case(SUZUKI_TROTTER);       call td_rti1
        case(REVERSAL);             call td_rti2
        case(APP_REVERSAL);         call td_rti3
        case(EXPONENTIAL_MIDPOINT); call td_rti4
        case(MAGNUS);               call td_rti5
        end select
      end do

      deallocate(zpsi1)
    end if

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    ! Split operator.
    subroutine td_rti0
      integer :: ik, ist
      call push_sub('td_rti.td_rti0')

      do ik = 1, st%d%nik
        do ist = 1, st%nst
          call zexp_kinetic(gr, h, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
        end do
      end do
      call zstates_calc_dens(st, NP, st%rho(1:NP,:))
      call zv_ks_calc(gr, ks, h, st)
      do ik = 1, st%d%nik
        do ist = 1, st%nst
          if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, st%zpsi(:, :, ist, ik), -M_zI*dt, .true.)
          call zexp_vlpsi (gr, h, st%zpsi(:, :, ist, ik), ik, t-dt*M_HALF, -M_zI*dt)
          if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, st%zpsi(:, :, ist, ik), -M_zI*dt, .false.)
        end do
      end do
      do ik = 1, st%d%nik
        do ist = 1, st%nst
          call zexp_kinetic(gr, h, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
        end do
      end do

      call pop_sub()
    end subroutine td_rti0


    ! ---------------------------------------------------------
    ! Suzuki-Trotter.
    subroutine td_rti1
      FLOAT :: p, pp(5), time(5), dtime(5)
      integer :: ik, ist, k

      call push_sub('td_rti.td_rti1')

      p = M_ONE/(M_FOUR - M_FOUR**(M_THIRD))
      pp = (/ p, p, M_ONE-M_FOUR*p, p, p /)
      dtime = pp*dt
      time(1) = t-dt+pp(1)/M_TWO*dt
      time(2) = t-dt+(pp(1)+pp(2)/M_TWO)*dt
      time(3) = t-dt+(pp(1)+pp(2)+pp(3)/M_TWO)*dt
      time(4) = t-dt+(pp(1)+pp(2)+pp(3)+pp(4)/M_TWO)*dt
      time(5) = t-dt+(pp(1)+pp(2)+pp(3)+pp(4)+pp(5)/M_TWO)*dt

      do k = 1, 5
        call dextrapolate(2, NP, st%d%nspin, tr%v_old(:, :, 0:2), h%vhxc, dt, time(k))
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            call zexp_vlpsi (gr, h, st%zpsi(:, :, ist, ik), ik, time(k), -M_zI*dtime(k)/M_TWO)
            if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, &
              st%zpsi(:, :, ist, ik), -M_zI*dtime(k)/M_TWO, .true.)

            call zexp_kinetic(gr, h, st%zpsi(:, :, ist, ik), tr%cf, -M_zI*dtime(k))
            if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, &
              st%zpsi(:, :, ist, ik), -M_zI*dtime(k)/M_TWO, .false.)
            call zexp_vlpsi (gr, h, st%zpsi(:, :, ist, ik), ik, time(k), -M_zI*dtime(k)/M_TWO)
          end do
        end do
      end do

      call pop_sub()
    end subroutine td_rti1


    ! ---------------------------------------------------------
    subroutine td_rti2
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      CMPLX, allocatable :: zpsi1(:,:,:,:)
      integer :: ik, ist

      call push_sub('td_rti.td_rti2')

      if(.not.h%ip_app) then
        ALLOCATE(zpsi1(NP_PART, st%d%dim, st%st_start:st%st_end, st%d%nik), NP_PART*(st%st_end-st%st_start+1)*st%d%nik)
        zpsi1 = st%zpsi ! store zpsi

        ALLOCATE(vhxc_t1(NP, st%d%nspin), NP*st%d%nspin)
        ALLOCATE(vhxc_t2(NP, st%d%nspin), NP*st%d%nspin)
        vhxc_t1 = h%vhxc

        ! propagate dt with H(t-dt)
        do ik = 1, st%d%nik
          do ist = st%st_start, st%st_end
            call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt, t-dt)
          end do
        end do

        call zstates_calc_dens(st, NP, st%rho(1:NP,:))
        call zv_ks_calc(gr, ks, h, st)

        st%zpsi = zpsi1
        deallocate(zpsi1)

        vhxc_t2 = h%vhxc
        h%vhxc = vhxc_t1
      end if

      ! propagate dt/2 with H(t-dt)
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt/M_TWO, t-dt)
        end do
      end do

      ! propagate dt/2 with H(t)
      if(.not.h%ip_app) h%vhxc = vhxc_t2
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt/M_TWO, t)
        end do
      end do

      if(.not.h%ip_app) deallocate(vhxc_t1, vhxc_t2)

      call pop_sub()
    end subroutine td_rti2


    ! ---------------------------------------------------------
    subroutine td_rti3
      integer ik, ist

      call push_sub('td_rti.td_rti3')

      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt/M_TWO, t-dt)
        end do
      end do

      h%vhxc = tr%v_old(:, :, 0)
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt/M_TWO, t)
        end do
      end do

      call pop_sub()
    end subroutine td_rti3


    ! ---------------------------------------------------------
    subroutine td_rti4
      integer :: ist, ik

      call push_sub('td_rti.td_rti4')

      if(.not.h%ip_app) then
        call dextrapolate(2, NP, st%d%nspin, tr%v_old(:, :, 0:2), h%vhxc, dt, -dt/M_TWO)
      end if

      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt, t - dt/M_TWO)
        end do
      end do

      call pop_sub()
    end subroutine td_rti4


    ! ---------------------------------------------------------
    subroutine td_rti5
      integer :: j, is, ist, k
      FLOAT :: time(2)
      FLOAT, allocatable :: vaux(:, :, :)

      call push_sub('td_rti.td_rti5')

      ALLOCATE(vaux(NP, st%d%nspin, 2), NP*st%d%nspin*2)

      time(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
      time(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

      if(.not.h%ip_app) then
        do j = 1, 2
          call dextrapolate(2, NP, st%d%nspin, tr%v_old(:,:, 0:2), vaux(:,:, j), dt, time(j)-dt)
        end do
      else
        vaux = M_ZERO
      end if

      do j = 1, 2
        if(h%ep%no_lasers > 0) then
          select case(h%gauge)
          case(1) ! length gauge
             do is = 1, st%d%nspin
                vaux(:, is, j) = vaux(:, is, j) + epot_laser_scalar_pot(gr%m%np, gr, h%ep, t-dt+time(j))
             end do
          case(2) ! velocity gauge
            write(message(1),'(a)') 'Inconsistency in td_rti5'
            call write_fatal(1)
          end select
        end if
      end do

      tr%vmagnus(:, :, 2)  = M_HALF*(vaux(:, :, 1) + vaux(:, :, 2))
      tr%vmagnus(:, :, 1) = (sqrt(M_THREE)/CNST(12.0))*dt*(vaux(:, :, 2) - vaux(:, :, 1))

      do k = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, k), k, dt, M_ZERO, &
            vmagnus = tr%vmagnus)
        end do
      end do

      deallocate(vaux)
      call pop_sub()
    end subroutine td_rti5

  end subroutine td_rti_dt

end module td_rti_m
