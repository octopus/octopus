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
!! $Id: opt_control.F90 2868 2007-04-26 17:11:47Z acastro $


  ! ---------------------------------------------------------
  ! Calculates the J1 functional, i.e.:
  ! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  ! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  ! the time-dependent case.
  ! ---------------------------------------------------------
  FLOAT function j1_functional(gr, geo, ep, psi, target) result(j1)
    type(grid_t),   intent(inout)   :: gr
    type(geometry_t), intent(inout) :: geo
    type(epot_t), intent(inout)     :: ep
    type(states_t), intent(inout)   :: psi
    type(target_t), intent(in)      :: target

    integer :: i, p, j
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: opsi(:, :)

    call push_sub('aux.j1_functional')

    select case(target%type)
    case(oct_tg_density)

      ALLOCATE(local_function(NP), NP)
      do i = 1, NP
        local_function(i) = - ( sqrt(psi%rho(i, 1)) - sqrt(target%rho(i)) )**2
      end do
      j1 = dmf_integrate(gr%m, local_function)
      deallocate(local_function)

    case(oct_tg_local)

      select case(psi%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi%d%nik.eq.1)
        ALLOCATE(opsi(NP_PART, 1), NP_PART)
        opsi = M_z0
        j1 = M_ZERO
        do p  = psi%st_start, psi%st_end
          do j = 1, NP
            opsi(j, 1) = target%rho(j) * psi%zpsi(j, 1, p, 1)
          end do
          j1 = j1 + zstates_dotp(gr%m, psi%d%dim, psi%zpsi(:, :, p, 1), opsi(:, :))
        end do
        deallocate(opsi)
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_td_local)
      j1 = sum(target%td_fitness) 

    case(oct_tg_excited)
      j1 = abs(zstates_mpdotp(gr%m, target%est, psi))**2

    case(oct_tg_exclude_state)
      j1 = M_ONE
      do i = 1, target%excluded_states
        j1 = j1 - abs(zstates_dotp(gr%m, psi%d%dim, target%st%zpsi(:, :, i, 1), psi%zpsi(:, :, 1, 1)))**2
      end do

    case default
      j1 = abs(zstates_mpdotp(gr%m, psi, target%st))**2

    end select

    call pop_sub()
  end function j1_functional


  ! ---------------------------------------------------------
  ! calculate |chi(T)> = \hat{O}(T) |psi(T)>
  ! ---------------------------------------------------------
  subroutine calc_chi(oct, gr, geo, ep, target, psi_in, chi_out)
    type(oct_t),       intent(in)  :: oct
    type(grid_t),      intent(inout)  :: gr
    type(geometry_t),  intent(inout) :: geo
    type(epot_t),  intent(in) :: ep
    type(states_t),    intent(inout)  :: psi_in
    type(target_t),    intent(in)  :: target
    type(states_t),    intent(inout) :: chi_out
    
    CMPLX   :: olap, zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, p, dim, k, j, no_electrons, ia, ib, n_pairs, nst, kpoints, jj

    call push_sub('opt_control.calc_chi')

    no_electrons = -nint(psi_in%val_charge)

    select case(target%type)

    case(oct_tg_density)

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)

        ASSERT(psi_in%d%nik.eq.1)

        if(no_electrons .eq. 1) then
          do j = 1, NP
            chi_out%zpsi(j, 1, 1, 1) = sqrt(target%rho(j)) * &
              exp( M_z1 * atan2(aimag(psi_in%zpsi(j, 1, 1, 1)), &
                                real(psi_in%zpsi(j, 1, 1, 1)  )) )
          end do
        else
          do p  = psi_in%st_start, psi_in%st_end
            do j = 1, NP_PART
              if(psi_in%rho(j, 1) > CNST(1.0e-8)) then
                chi_out%zpsi(j, 1, p, 1) = sqrt(target%rho(j)/psi_in%rho(j, 1)) * &
                  psi_in%zpsi(j, 1, p, 1)
              else
                chi_out%zpsi(j, 1, p, 1) = M_ZERO!sqrt(target%rho(j))
              end if
            end do
          end do
        end if

      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_local)

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)
        ASSERT(psi_in%d%nik.eq.1)
        do p  = psi_in%st_start, psi_in%st_end
          do j = 1, NP
            chi_out%zpsi(j, 1, p, 1) = target%rho(j) * psi_in%zpsi(j, 1, p, 1)
          end do
        end do
      case(SPIN_POLARIZED); stop 'Error'
      case(SPINORS);        stop 'Error'
      end select

    case(oct_tg_td_local)
      !We assume that there is no time-independent operator.
      do ik = 1, chi_out%d%nik
        do p = chi_out%st_start, chi_out%st_end
          do dim = 1, chi_out%d%dim
            do j = 1, NP
              chi_out%zpsi(j, dim, p, ik) = M_z0
            end do
          end do
        end do
      end do

    case(oct_tg_excited) 

      n_pairs = target%est%n_pairs
      kpoints = psi_in%d%nik
      nst = psi_in%nst

      ALLOCATE(cI(n_pairs), n_pairs)
      ALLOCATE(dI(n_pairs), n_pairs)
      ALLOCATE(mat(target%est%st%nst, nst, psi_in%d%nik), target%est%st%nst*nst*psi_in%d%nik)
      ALLOCATE(mm(nst, nst, kpoints, n_pairs), nst*nst*kpoints*n_pairs)
      ALLOCATE(mk(NP_PART, psi_in%d%dim), NP_PART * psi_in%d%dim)
      ALLOCATE(lambda(n_pairs, n_pairs), n_pairs*n_pairs)

      call zstates_matrix(gr%m, target%est%st, psi_in, mat)

      do ia = 1, n_pairs
        cI(ia) = target%est%weight(ia)
        call zstates_matrix_swap(mat, target%est%pair(ia))
        mm(1:nst, 1:nst, 1:kpoints, ia) = mat(1:nst, 1:kpoints, 1:kpoints)
        dI(ia) = zstates_mpdotp(gr%m, target%est%st, psi_in, mat)
        if(abs(dI(ia)) > CNST(1.0e-12)) then
          do ik = 1, kpoints
            zdet = lalg_inverter(nst, mm(1:nst, 1:nst, ik, ia))
          end do
        end if
        call zstates_matrix_swap(mat, target%est%pair(ia))
      end do

      do ia = 1, n_pairs
        do ib = 1, n_pairs
          lambda(ia, ib) = conjg(cI(ib)) * cI(ia) * conjg(dI(ia)) * dI(ib)
        end do
      end do

      select case(psi_in%d%ispin)
      case(UNPOLARIZED)
        write(message(1), '(a)') 'Internal error in aux.calc_chi'
        call write_fatal(2)

      case(SPIN_POLARIZED)
        ASSERT(chi_out%d%nik .eq. 2)

        do ik = 1, kpoints
          do k = chi_out%st_start, chi_out%st_end
            chi_out%zpsi(:, :, k, ik) = M_z0
            do ia = 1, n_pairs
              if(ik .ne. target%est%pair(ia)%sigma) cycle
              if(abs(dI(ia)) < CNST(1.0e-12)) cycle
              do ib = 1, n_pairs
                if(abs(dI(ib)) < CNST(1.0e-12)) cycle
                mk = M_z0
                do j = 1, nst
                  if(j .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                  mk(:, :) = mk(:, :) + conjg(mm(k, j, ik, ib)) * target%est%st%zpsi(:, :, jj, ik)
                end do
                call lalg_axpy(NP_PART, psi_in%d%dim, M_z1, lambda(ib, ia)*mk(:, :), chi_out%zpsi(:, :, k, ik))
              end do
            end do
          end do
        end do
        
      case(SPINORS)
        ASSERT(chi_out%d%nik .eq. 1)

        do k = chi_out%st_start, chi_out%st_end
          chi_out%zpsi(:, :, k, 1) = M_z0

          do ia = 1, n_pairs
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle

            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle

              mk = M_z0
              do j = 1, nst
                if(j .eq. target%est%pair(ib)%i) jj = target%est%pair(ia)%a
                mk(:, :) = mk(:, :) + conjg(mm(k, j, 1, ib)) * target%est%st%zpsi(:, :, jj, 1)
              end do

              call lalg_axpy(NP_PART, 2, M_z1, lambda(ib, ia)*mk(:, :), chi_out%zpsi(:, :, k, 1))
            end do
          end do
        end do

      end select

      deallocate(cI, dI, mat, mm, mk, lambda)

    case(oct_tg_exclude_state)

      chi_out%zpsi(:, :, 1, 1) = psi_in%zpsi(:, :, 1, 1)
      do p = 1, target%excluded_states
        olap = zstates_dotp(gr%m, psi_in%d%dim, target%st%zpsi(:, :, p, 1), psi_in%zpsi(:, :, 1, 1))
        chi_out%zpsi(:, :, 1, 1) = chi_out%zpsi(:, :, 1, 1) - olap*target%st%zpsi(:, :, p, 1)
      end do

    case default

      olap = zstates_mpdotp(gr%m, target%st, psi_in)
      do ik = 1, psi_in%d%nik
        do p  = psi_in%st_start, psi_in%st_end
          select case(oct%algorithm_type)
            case(oct_algorithm_zbr98)
              chi_out%zpsi(:, :, p, ik) = target%st%zpsi(:, :, p, ik)
            case default
              chi_out%zpsi(:, :, p, ik) = olap*target%st%zpsi(:, :, p, ik)
          end select
        end do
      end do

    end select

    call pop_sub()
  end subroutine calc_chi


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
