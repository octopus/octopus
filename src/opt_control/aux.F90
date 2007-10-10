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
  FLOAT function j1(m, max_iter, psi, target)
    type(mesh_t), intent(in)   :: m
    integer, intent(in)        :: max_iter
    type(states_t), intent(in) :: psi
    type(target_t), intent(in) :: target

    integer :: i, p, j
    FLOAT, allocatable :: local_function(:)
    CMPLX, allocatable :: opsi(:, :)

    call push_sub('opt_control.overlap_function')

    if(target%targetmode==oct_targetmode_td) then
       ! 1/T * int(<Psi| O | Psi>)
       j1 = sum(target%td_fitness) / real(max_iter, REAL_PRECISION) 
    else
      select case(target%totype)
      case(oct_tg_density)

        ALLOCATE(local_function(m%np), m%np)
        do i = 1, m%np
          local_function(i) = - ( sqrt(psi%rho(i, 1)) - sqrt(target%rho(i)) )**2
        end do
        j1 = dmf_integrate(m, local_function)
        deallocate(local_function)

      case(oct_tg_local)

        select case(psi%d%ispin)
        case(UNPOLARIZED)
          ASSERT(psi%d%nik.eq.1)
          ALLOCATE(opsi(m%np_part, 1), m%np_part)
          opsi = M_z0
          j1 = M_ZERO
          do p  = psi%st_start, psi%st_end
            do j = 1, m%np
              opsi(j, 1) = target%rho(j) * psi%zpsi(j, 1, p, 1)
            end do
            j1 = j1 + zstates_dotp(m, psi%d%dim, psi%zpsi(:, :, p, 1), opsi(:, :))
          end do
          deallocate(opsi)
        case(SPIN_POLARIZED); stop 'Error'
        case(SPINORS);        stop 'Error'
        end select

      case default

        j1 = abs(zstates_mpdotp(m, psi, target%st))**2

      end select
    end if

    call pop_sub()
  end function j1


  ! ---------------------------------------------------------
  ! calculate |chi> = \hat{O} |psi>
  ! ---------------------------------------------------------
  subroutine calc_chi(oct, gr, target, psi_in, chi_out)
    type(oct_t),       intent(in)  :: oct
    type(grid_t),      intent(in)  :: gr
    type(states_t),    intent(in)  :: psi_in
    type(target_t),    intent(in)  :: target

    type(states_t),    intent(inout) :: chi_out
    
    CMPLX   :: olap
    integer :: ik, p, dim, j, no_electrons
  
    call push_sub('opt_control.calc_chi')

    no_electrons = -nint(psi_in%val_charge)

    if(target%targetmode == oct_targetmode_td) then
      message(1) = 'Info: Time-dependent Target selected'
      call write_info(1)
      chi_out%zpsi = M_z0
      call pop_sub()
      return
    end if


    select case(target%totype)

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

    case default

      olap = zstates_mpdotp(gr%m, target%st, psi_in)
      do ik = 1, psi_in%d%nik
        do p  = psi_in%st_start, psi_in%st_end
          if(oct%algorithm_type == oct_algorithm_zr98) &
            chi_out%zpsi(:,:,p,ik) = olap*target%st%zpsi(:, :, p, ik)
          if(oct%algorithm_type == oct_algorithm_zbr98) then
            chi_out%zpsi(:,:,p,ik) = target%st%zpsi(:, :, p, ik)
          end if
          if(oct%algorithm_type == oct_algorithm_wg05) &
            chi_out%zpsi(:,:,p,ik) = olap*target%st%zpsi(:, :, p, ik)
        end do
      end do

    end select

      
    call pop_sub()
  end subroutine calc_chi


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
