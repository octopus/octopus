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

module td_rti
  use global
  use messages
  use lib_oct_parser
  use lib_basic_alg
  use math
  use mesh
  use cube_function
  use functions
  use mesh_function
  use states
  use v_ks
  use hamiltonian
  use external_pot
  use td_exp
  use td_exp_split
  use grid

  implicit none

  type td_rti_type
     integer           :: method         ! which evolution method to use
     type(td_exp_type) :: te             ! how to apply the propagator (e^{-i H \Delta t})


     FLOAT, pointer :: v_old(:, :, :) ! storage of the KS potential of previous iterations
     FLOAT, pointer :: vmagnus(:, :, :) ! auxiliary function to store the Magnus potentials.

     type(zcf) :: cf             ! auxiliary cube for split operator methods
  end type td_rti_type

  integer, parameter ::            &
       SPLIT_OPERATOR       = 0,   &
       SUZUKI_TROTTER       = 1,   &
       REVERSAL             = 2,   &
       APP_REVERSAL         = 3,   &
       EXPONENTIAL_MIDPOINT = 4,   &
       MAGNUS               = 5

  FLOAT, parameter :: scf_threshold = CNST(1.0e-3)

  private
  public :: td_rti_type,     &
       td_rti_init,          &
       td_rti_end,           &
       td_rti_run_zero_iter, &
       td_rti_dt


contains

  !-------------------------------------------------------------------
  subroutine td_rti_init(gr, st, tr)
    type(grid_type),   intent(in)    :: gr
    type(states_type), intent(in)    :: st
    type(td_rti_type), intent(inout) :: tr

    call loct_parse_int(check_inp('TDEvolutionMethod'), REVERSAL, tr%method)
    select case(tr%method)
    case(SPLIT_OPERATOR)
       call zcf_new(gr%m%l, tr%cf)
       call zcf_fft_init(tr%cf, gr%sb)
       message(1) = 'Info: Evolution method:  Split-Operator'
    case(SUZUKI_TROTTER)
       call zcf_new(gr%m%l, tr%cf)
       call zcf_fft_init(tr%cf, gr%sb)
       message(1) = 'Info: Evolution method:  Suzuki-Trotter'
    case(REVERSAL);             message(1) = 'Info: Evolution method:  Enforced Time-Reversal Symmetry'
    case(APP_REVERSAL);         message(1) = 'Info: Evolution method:  Approx.Enforced Time-Reversal Symmetry'
    case(EXPONENTIAL_MIDPOINT); message(1) = 'Info: Evolution method:  Exponential Midpoint Rule.'
    case(MAGNUS);               message(1) = 'Info: Evolution method:  Magnus expansion.'
       allocate(tr%vmagnus(NP, st%d%nspin, 2))
    case default
       write(message(1), '(a,i6,a)') "Input: '", tr%method, "' is not a valid TDEvolutionMethod"
       message(2) = '(0 <= TDEvolutionMethod <= 5)'
       call write_fatal(2)
    end select
    call write_info(1)

    allocate(tr%v_old(NP, st%d%nspin, 0:3)) ! allocate memory to store the old KS potentials
    call td_exp_init(gr, tr%te)            ! initialize propagator

  end subroutine td_rti_init


  !-------------------------------------------------------------------
  subroutine td_rti_end(tr)
    type(td_rti_type), intent(inout) :: tr

    ASSERT(associated(tr%v_old)) ! sanity check
    deallocate(tr%v_old)         ! clean ols KS potentials
    nullify(tr%v_old)

    if(tr%method == MAGNUS) then
       ASSERT(associated(tr%vmagnus))
       deallocate(tr%vmagnus); nullify(tr%vmagnus)
    endif

    if(tr%method == SUZUKI_TROTTER .or. tr%method == SPLIT_OPERATOR) then
       call zcf_free(tr%cf)
    endif

    call td_exp_end(tr%te)       ! clean propagator method
  end subroutine td_rti_end


  !-------------------------------------------------------------------
  subroutine td_rti_run_zero_iter(h, tr)
    type(hamiltonian_type), intent(in)    :: h
    type(td_rti_type),      intent(inout) :: tr
    tr%v_old(:, :, 2) = h%vhxc(:, :)
    tr%v_old(:, :, 3) = h%vhxc(:, :)
    tr%v_old(:, :, 1) = h%vhxc(:, :)
  end subroutine td_rti_run_zero_iter


  !-------------------------------------------------------------------
  subroutine td_rti_dt(ks, h, gr, st, tr, t, dt)
    type(v_ks_type),        intent(inout) :: ks
    type(hamiltonian_type), intent(inout) :: h
    type(grid_type),        intent(inout) :: gr
    type(states_type),      intent(inout) :: st
    type(td_rti_type),      intent(inout) :: tr
    FLOAT,                  intent(in)    :: t, dt

    integer :: is, iter
    FLOAT   :: d, d_max
    logical :: self_consistent
    CMPLX, allocatable :: zpsi1(:, :, :, :)
    FLOAT, allocatable :: dtmp(:)

    call push_sub('td_rti.td_rti')

    self_consistent = .false.
    if(t<3*dt .and. (.not.h%ip_app)) then
       self_consistent = .true.
       allocate(zpsi1(NP, st%d%dim, st%st_start:st%st_end, st%d%nik))
       zpsi1 = st%zpsi
    endif

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

          call zstates_calc_dens(st, NP, st%rho, .true.)
          call zv_ks_calc(gr, ks, h, st)
          tr%v_old(:, :, 0) = h%vhxc
          h%vhxc = tr%v_old(:, :, 1)

          d_max = M_ZERO
          allocate(dtmp(NP))
          do is = 1, st%d%nspin
             dtmp = tr%v_old(:, is, 3) - tr%v_old(:, is, 0)
             d    = dmf_nrm2(gr%m, dtmp)
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
       enddo
       deallocate(zpsi1)
    endif

    call pop_sub()
  contains

    ! Split operator.
    subroutine td_rti0
      integer :: ik, ist
      call push_sub('td_rti.td_rti0')

      do ik = 1, st%d%nik
         do ist = 1, st%nst
            call zexp_kinetic(gr, h, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
         enddo
      enddo
      call zstates_calc_dens(st, NP, st%rho, .true.)
      call zv_ks_calc(gr, ks, h, st)
      do ik = 1, st%d%nik
         do ist = 1, st%nst
            if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, st%zpsi(:, :, ist, ik), -M_zI*dt, .true.)
            call zexp_vlpsi (gr, h, st%zpsi(:, :, ist, ik), ik, t-dt*M_HALF, -M_zI*dt)
            if (h%ep%nvnl > 0) call zexp_vnlpsi (gr%m, h, st%zpsi(:, :, ist, ik), -M_zI*dt, .false.)
         enddo
      enddo
      do ik = 1, st%d%nik
         do ist = 1, st%nst
            call zexp_kinetic(gr, h, st%zpsi(:, :, ist, ik), tr%cf, -M_HALF*M_zI*dt)
         enddo
      enddo

      call pop_sub()
    end subroutine td_rti0

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
            enddo
         enddo
      enddo

      call pop_sub()
    end subroutine td_rti1

    subroutine td_rti2
      FLOAT, allocatable :: vhxc_t1(:,:), vhxc_t2(:,:)
      CMPLX, allocatable :: zpsi1(:,:,:,:)
      integer :: ik, ist

      call push_sub('td_rti.td_rti2')

      if(.not.h%ip_app) then
         allocate(zpsi1(NP, st%d%dim, st%st_start:st%st_end, st%d%nik))
         zpsi1 = st%zpsi ! store zpsi

         allocate(vhxc_t1(NP, st%d%nspin), vhxc_t2(NP, st%d%nspin))
         vhxc_t1 = h%vhxc

         ! propagate dt with H(t-dt)
         do ik = 1, st%d%nik
            do ist = st%st_start, st%st_end
               call td_exp_dt(tr%te, gr, h, st%zpsi(:,:, ist, ik), ik, dt, t-dt)
            end do
         end do

         call zstates_calc_dens(st, NP, st%rho, .true.)
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

    subroutine td_rti5
      integer :: j, is, ist, k
      FLOAT :: time(2)
      FLOAT, allocatable :: vaux(:, :, :)

      call push_sub('td_rti.td_rti5')

      allocate(vaux(NP, st%d%nspin, 2))

      time(1) = (M_HALF-sqrt(M_THREE)/M_SIX)*dt
      time(2) = (M_HALF+sqrt(M_THREE)/M_SIX)*dt

      if(.not.h%ip_app) then
         do j = 1, 2
            call dextrapolate(2, NP, st%d%nspin, tr%v_old(:,:, 0:2), vaux(:,:, j), dt, time(j)-dt)
         enddo
      else
         vaux = M_ZERO
      endif

      do j = 1, 2
         if(h%ep%no_lasers > 0) then
            select case(h%gauge)
            case(1) ! length gauge
               vaux(:, is, j) = vaux(:, is, j) + epot_laser_scalar_pot(gr%m%np, gr, h%ep, t-dt+time(j))
            case(2) ! velocity gauge
               write(message(1),'(a)') 'Inconsistency in td_rti5'
               call write_fatal(1)
            end select
         endif
      enddo

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

end module td_rti
