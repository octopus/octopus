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

module states
use global
use lib_oct_parser
use io
use lib_alg
use math
use mesh
use functions
use output

implicit none

type states_type
  integer :: dim           ! Dimension of the state (one or two for spinors)
  integer :: nst           ! Number of states in each irreducible subspace
  integer :: nik           ! Number of irreducible subspaces
  integer :: nik_axis(3)   ! Number of kpoints per axis
  integer :: ispin         ! spin mode (unpolarized, spin polarized, spinors)
  integer :: nspin         ! dimension of rho (1, 2 or 4)
  integer :: spin_channels ! 1 or 2, wether spin is or not considered.
  logical :: select_axis(3)! which axes are used fo k points

  ! pointers to the wavefunctions
  FLOAT, pointer :: dpsi(:,:,:,:)
  CMPLX, pointer :: zpsi(:,:,:,:)

  ! the densities (after all we are doing DFT :)
  FLOAT, pointer :: rho(:,:)

  FLOAT, pointer :: rho_core(:)! core charge for nl core corrections

  FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
  logical           :: fixed_occ ! should the occupation numbers be fixed?
  FLOAT, pointer :: occ(:,:)  ! the occupation numbers
  FLOAT, pointer :: mag(:, :, :)

  FLOAT :: qtot    ! (-) The total charge in the system (used in Fermi)

  FLOAT :: el_temp ! electronic temperature for the Fermi function
  FLOAT :: ef      ! the fermi energy

  integer :: st_start, st_end ! needed for some parallel parts

  FLOAT, pointer :: kpoints(:,:) ! obviously the kpoints
  FLOAT, pointer :: kweights(:)  ! weights for the kpoint integrations

end type states_type

! Parameters...
integer, parameter :: UNPOLARIZED    = 1, &
                      SPIN_POLARIZED = 2, &
                      SPINORS        = 3

interface assignment (=)
  module procedure states_copy
end interface

contains

subroutine states_null(st)
  type(states_type) :: st
  nullify(st%dpsi, st%zpsi, st%rho, st%rho_core, st%eigenval, st%occ, st%mag, st%kpoints, st%kweights)
end subroutine states_null

subroutine states_init(st, m, val_charge)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  FLOAT, intent(in) :: val_charge

  FLOAT :: excess_charge, r
  integer :: nempty, i, j
  character(len=80) :: str

  call push_sub('states_init')

  call states_null(st)

  call loct_parse_int('SpinComponents', UNPOLARIZED, st%ispin)
  if (st%ispin < UNPOLARIZED .or. st%ispin > SPINORS) then
    write(message(1),'(a,i4,a)') "Input: '", st%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if

  if (conf%periodic_dim>0) then
    if(loct_parse_block_n('NumberKPoints')<1) then
      message(1) = 'Block "NumberKPoints" not found in input file.'
      call write_fatal(1)
    end if
    st%nik_axis = 1
    do i = 1, conf%periodic_dim
      call loct_parse_block_int('NumberKPoints', 0, i-1, st%nik_axis(i))
    end do    
    if (any(st%nik_axis < 1)) then
      message(1) = 'Input: NumberKPoints is not valid'
      message(2) = '(NumberKPoints >= 1)'
      call write_fatal(2)
    end if
    
    select case(loct_parse_block_n('SelectKAxis'))
    case(1)
      do i = 1, conf%periodic_dim
         call loct_parse_block_logical('SelectKAxis', 0, i-1, st%select_axis(i))
      end do
      do i=1, conf%periodic_dim
         if (.not.st%select_axis(i)) then
          write(message(1),'(a,i1,a)')'Info: K points in axis ',i,' are neglected'
          call write_info(1) 
          st%nik_axis(i)=1
         end if
      end do
    case default
      st%select_axis=.true.
    end select
    st%nik=PRODUCT(st%nik_axis)
  else
    st%nik = 1
  end if

  
  call loct_parse_float('ExcessCharge', M_ZERO, excess_charge)

  call loct_parse_int('ExtraStates', 0, nempty)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
    message(2) = '(0 <= ExtraStates)'
    call write_fatal(2)
  end if
  
!!$  st%qtot = -(val_charge + excess_charge)
!!$  st%nst  = int(st%qtot/M_TWO)
!!$  if(st%nst*M_TWO < st%qtot) &
!!$       st%nst = st%nst + 1
!!$
!!$  select case(st%ispin)
!!$  case(UNPOLARIZED)
!!$    st%dim = 1
!!$    st%nst = st%nst + nempty
!!$  case(SPIN_POLARIZED)
!!$    st%dim = 1
!!$    st%nst = st%nst + nempty
!!$    st%nik = st%nik*2
!!$  case(SPINORS)
!!$    st%dim = 2
!!$    st%nst = st%nst*2 + nempty
!!$  end select

  st%qtot = -(val_charge + excess_charge)

  select case(st%ispin)
  case(UNPOLARIZED)
    st%dim = 1
    st%nst = int(st%qtot/2)
    if(st%nst*2 < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%nspin = 1
    st%spin_channels = 1
  case(SPIN_POLARIZED)
    st%dim = 1
    st%nst = int(st%qtot/2)
    if(st%nst*2 < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%nik = st%nik*2
    st%nspin = 2
    st%spin_channels = 2
  case(SPINORS)
    st%dim = 2
    st%nst = int(st%qtot)
    if(st%nst < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%nspin = 4
    st%spin_channels = 2
  end select

  ! For non-periodic systems this should just return the Gamma point
  call states_choose_kpoints(st, m)

  ! we now allocate some arrays
!!$  st%nspin = min(2, st%ispin)
  allocate(st%rho(m%np, st%nspin), &
           st%occ(st%nst, st%nik), &
           st%eigenval(st%nst, st%nik))
  if(st%ispin == SPINORS) then
    allocate(st%mag(st%nst, st%nik, 2))
  end if

  str = "Occupations"
  occ_fix: if(loct_parse_isdef(str) .ne. 0) then
    ! read in occupations
    st%fixed_occ = .true.

    do i = 1, st%nik
      do j = 1, st%nst
        call loct_parse_block_float(str, i-1, j-1, st%occ(j, i))
      end do
    end do
  else
    st%fixed_occ = .false.

    ! first guest for occupation...paramagnetic configuration
    if(st%ispin == UNPOLARIZED) then
      r = M_TWO
    else
      r = M_ONE
    endif
!!$    r = M_TWO / st%nspin
    st%occ  = M_ZERO
    st%qtot = M_ZERO
    do j = 1, st%nst
      do i = 1, st%nik
        st%occ(j, i) = min(r, -(val_charge + excess_charge) - st%qtot)
        st%qtot = st%qtot + st%occ(j, i)
      end do
    end do

    ! read in fermi distribution temperature
    call loct_parse_float('ElectronicTemperature', M_ZERO, st%el_temp)
  end if occ_fix

  st%st_start = 1; st%st_end = st%nst
  nullify(st%dpsi, st%zpsi, st%rho_core)

  call pop_sub()
end subroutine states_init

subroutine states_copy(stout, stin)
  type(states_type), intent(in)   :: stin
  type(states_type), intent(out)  :: stout

  call states_null(stout)

  stout%dim           = stin%dim
  stout%nst           = stin%nst
  stout%nik           = stin%nik
  stout%ispin         = stin%ispin
  stout%nspin         = stin%nspin
  stout%spin_channels = stin%spin_channels
  stout%qtot = stin%qtot
  stout%el_temp = stin%el_temp
  stout%ef = stin%ef
  stout%st_start = stin%st_start
  stout%st_end = stin%st_end
  if(associated(stin%dpsi)) then
    allocate(stout%dpsi(size(stin%dpsi, 1), stin%dim, stin%st_start:stin%st_end, stin%nik))
    stout%dpsi = stin%dpsi
  endif
  if(associated(stin%zpsi)) then
    allocate(stout%zpsi(size(stin%zpsi, 1), stin%dim, stin%st_start:stin%st_end, stin%nik))
    stout%zpsi = stin%zpsi
  endif
  if(associated(stin%rho)) then
    allocate(stout%rho(size(stin%rho, 1), size(stin%rho, 2)))
    stout%rho = stin%rho
  endif
  if(associated(stin%rho_core)) then
    allocate(stout%rho_core(size(stin%rho_core, 1)))
    stout%rho_core = stin%rho_core
  endif
  if(associated(stin%eigenval)) then
    allocate(stout%eigenval(stin%st_start:stin%st_end, stin%nik))
    stout%eigenval = stin%eigenval
  endif
  stout%fixed_occ = stin%fixed_occ
  if(associated(stin%occ)) then
    allocate(stout%occ(size(stin%occ, 1), size(stin%occ, 2)))
    stout%occ = stin%occ
  endif
  if(associated(stin%mag)) then
    allocate(stout%mag(size(stin%mag, 1), size(stin%mag, 2), size(stin%mag, 3)))
    stout%mag = stin%mag
  endif
  if(associated(stin%kpoints)) then
    allocate(stout%kpoints(size(stin%kpoints, 1), size(stin%kpoints, 2)))
    stout%kpoints = stin%kpoints
  endif
  if(associated(stin%kweights)) then
    allocate(stout%kweights(size(stin%kweights, 1)))
    stout%kweights = stin%kweights
  endif

end subroutine states_copy

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  call push_sub('states_end')

  if(associated(st%rho)) then
    deallocate(st%rho, st%occ, st%eigenval)
    nullify   (st%rho, st%occ, st%eigenval)
  end if

  if(associated(st%rho_core)) then
    deallocate(st%rho_core)
    nullify(st%rho_core)
  end if

  if(st%ispin==3 .and. associated(st%mag)) then
    deallocate(st%mag); nullify(st%mag)
  end if

  if(associated(st%dpsi)) then
    deallocate(st%dpsi); nullify(st%dpsi)
  end if

  if(associated(st%zpsi)) then
    deallocate(st%zpsi); nullify(st%zpsi)
  end if

  if(associated(st%kpoints)) then
    nullify(st%kpoints, st%kweights)
  end if

  call pop_sub()
end subroutine states_end

! generate a hydrogen s-wavefunction around a random point
subroutine states_generate_random(st, m, ist_start)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  integer, intent(in), optional :: ist_start

  integer, save :: iseed = 123
  integer :: ist, ik, id, i, ist_s
  FLOAT :: a(3), rnd, r

  call push_sub('states_generate_random')

  ist_s = 1
  if(present(ist_start)) ist_s = ist_start

  do ik = 1, st%nik
    do ist = ist_s, st%nst
      do id = 1, st%dim
         call X(states_random)(m, st%X(psi)(1:m%np, id, ist, ik))
      end do
    end do
    call X(states_gram_schmidt)(st%nst, m, st%dim, st%X(psi)(1:,:,:,ik))
  end do
  st%eigenval = M_ZERO

  call pop_sub()
end subroutine states_generate_random

subroutine states_fermi(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(in) :: m

! Local variables
  integer :: ie, ik, iter
  integer, parameter :: nitmax = 200
  FLOAT :: drange, t, emin, emax, sumq
  FLOAT, parameter :: tol = CNST(1.0e-10)
  logical :: conv

  call push_sub('fermi')

  if(st%fixed_occ) then ! nothing to do
     ! Calculate magnetizations...
     if(st%ispin == 3) then
       do ik = 1, st%nik
         do ie = 1, st%nst
            st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
         enddo
       enddo
     endif
    call pop_sub()
    return
  end if

! Initializations
  emin = minval(st%eigenval)
  emax = maxval(st%eigenval)

  if(st%ispin == 3) then
     sumq = real(st%nst, PRECISION)
  else
     sumq = M_TWO*st%nst
  endif

  t = max(st%el_temp, CNST(1.0e-6))
  st%ef = emax

  conv = .true.
  if (abs(sumq - st%qtot) > tol) conv = .false.
  if (conv) then ! all orbitals are full; nothing to be done
     st%occ = M_TWO/st%spin_channels!st%nspin
     ! Calculate magnetizations...
     if(st%ispin == 3) then
       do ik = 1, st%nik
         do ie = 1, st%nst
            st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
         enddo
       enddo
     endif
     call pop_sub()
     return
  endif

  if (sumq < st%qtot) then ! not enough states
    message(1) = 'Fermi: Not enough states'
    write(message(2),'(6x,a,f12.6,a,f12.6)')'(total charge = ', st%qtot, &
        ' max charge = ', sumq
    call write_fatal(2)
  endif

  drange = t*sqrt(-log(tol*CNST(.01)))

  emin = emin - drange
  emax = emax + drange

  do iter = 1, nitmax
    st%ef = M_HALF*(emin + emax)
    sumq  = M_ZERO

    do ik = 1, st%nik
      do ie =1, st%nst
        sumq = sumq + st%kweights(ik)/st%spin_channels * & !st%nspin * &
             stepf((st%eigenval(ie, ik) - st%ef)/t)
      end do
    end do

    conv = .true.
    if(abs(sumq - st%qtot) > tol) conv = .false.
    if(conv) exit
    
    if(sumq <= st%qtot ) emin = st%ef
    if(sumq >= st%qtot ) emax = st%ef
  end do
  
  if(iter == nitmax) then
    message(1) = 'Fermi: did not converge'
    call write_fatal(1)
  end if

  do ik = 1, st%nik
    do ie = 1, st%nst
      st%occ(ie, ik) = stepf((st%eigenval(ie, ik) - st%ef)/t)/st%spin_channels!st%nspin
    end do
  end do

  ! Calculate magnetizations...
  if(st%ispin == 3) then
    do ik = 1, st%nik
       do ie = 1, st%nst
          st%mag(ie, ik, 1) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
          st%mag(ie, ik, 2) = X(mf_nrm2) (m, st%X(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
       enddo
    enddo
  endif
  
  call pop_sub(); return
end subroutine states_fermi

subroutine states_calculate_multipoles(m, st, pol, dipole, lmax, multipole)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  FLOAT, intent(in) :: pol(3)
  FLOAT, intent(out) :: dipole(st%nspin)
  integer, intent(in), optional :: lmax
  !FLOAT, intent(out), optional :: multipole((lmax + 1)**2, st%nspin)
  FLOAT, intent(out), optional :: multipole(:, :)

  integer :: i, is, l, lm, add_lm
  FLOAT :: x(3), r, ylm, mult

  call push_sub('states_calculate_multipoles')

  dipole = M_ZERO
  do is = 1, st%nspin
    do i = 1, m%np
      call mesh_xyz(m, i, x)
      dipole(is) = dipole(is) + st%rho(i, is)*sum(x*pol)
    end do
    dipole(is) = dipole(is) * m%vol_pp

    if(present(lmax).and.present(multipole)) then
      add_lm = 1
      do l = 0, lmax
         do lm = -l, l
            mult = M_ZERO
            do i = 1, m%np
               call mesh_r(m, i, r, x=x)
               ylm = loct_ylm(x(1), x(2), x(3), l, lm)
               if(l == 0) then
                 mult = mult + st%rho(i, is) * ylm
               else
                 mult = mult + st%rho(i, is) * ylm * r**l
               end if
            end do
            multipole(add_lm, is) = mult * m%vol_pp
            add_lm = add_lm + 1
         end do
      end do
    endif

  end do

  call pop_sub(); return
end subroutine states_calculate_multipoles

! function to calculate the eigenvalues sum using occupations as weights
function states_eigenvalues_sum(st)
  FLOAT :: states_eigenvalues_sum
  type(states_type), intent(in) :: st

  states_eigenvalues_sum = sum(st%eigenval * st%occ)

end function states_eigenvalues_sum

subroutine states_write_eigenvalues(iunit, nst, st, error)
  integer, intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st
  FLOAT, intent(in), optional :: error(nst, st%nik)

  integer ik, j, ns, is
  FLOAT :: o, oplus, ominus

  if(iunit==stdout.and.conf%verbose<=20) return

  ns = 1
  if(st%nspin == 2) ns = 2

  message(1) = 'Eigenvalues ['+trim(units_out%energy%abbrev)+']'
  call write_info(1, iunit)
  if (conf%periodic_dim>0) then 
  end if
  if (st%nik > ns) then
    message(1) = 'Kpoints ['+trim(units_out%length%abbrev)+'^-1]'
    call write_info(1, iunit)
  end if

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    do ik = 1, st%nik, ns
      if(st%nik > ns) then
        write(iunit, '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
        st%kpoints(1, ik)*units_out%length%factor, ',',           &
        st%kpoints(2, ik)*units_out%length%factor, ',',           &
        st%kpoints(3, ik)*units_out%length%factor, ')'
      end if
      
      do is = 1, ns
        write(iunit, '(a4)', advance='no') '#st'
        if(present(error)) then
          write(iunit, '(1x,a12,3x,a12,2x,a10,i3,a1)', advance='no') &
             ' Eigenvalue', 'Occupation ', 'Error (', is, ')'
        else
          write(iunit, '(1x,a12,3x,a12,2x)', advance='no') &
             ' Eigenvalue', 'Occupation '
        endif
      end do
      write(iunit, '(1x)', advance='yes')

      do j = 1, nst
        do is = 0, ns-1
          if(j > st%nst) then
            o = M_ZERO
            if(st%ispin == 3) oplus = M_ZERO; ominus = M_ZERO
          else
            o = st%occ(j, ik+is)
            if(st%ispin == 3) then 
              oplus  = st%mag(j, ik+is, 1)
              ominus = st%mag(j, ik+is, 2)
            endif
          end if
      
          write(iunit, '(i4)', advance='no') j
          if (conf%periodic_dim>0) then 
            if(st%ispin == 3) then
              write(iunit, '(1x,f12.6,3x,f5.3,a1,f5.3)', advance='no') &
                   (st%eigenval(j, ik)-st%ef)/units_out%energy%factor, oplus, '/', ominus
              if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
            else
              write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                    (st%eigenval(j, ik+is))/units_out%energy%factor, o
             if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
            endif
          else
            if(st%ispin == 3) then
              write(iunit, '(1x,f12.6,3x,f5.3,a1,f5.3)', advance='no') &
                   st%eigenval(j, ik)/units_out%energy%factor, oplus, '/', ominus
              if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
            else
              write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                   st%eigenval(j, ik+is)/units_out%energy%factor, o
              if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
            endif
          end if
        end do
        write(iunit, '(1x)', advance='yes')
      end do
    end do

#ifdef HAVE_MPI
  end if
#endif
  
end subroutine states_write_eigenvalues

subroutine states_write_bands(iunit, nst, st)
  integer, intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st

  integer ik, j, ns, is, pd
!  type(mesh_type):: m

  if(iunit==stdout.and.conf%verbose<=20) return
  
  ! shortcuts
  pd=conf%periodic_dim
  ns = 1
  if(st%nspin == 2) ns = 2

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif
      do j = 1, nst
        do ik = 1, st%nik, ns
            write(iunit, '(1x,3f12.4,3x,f12.6))', advance='yes') &
            st%kpoints(:,ik)*units_out%length%factor,            &
            (st%eigenval(j, ik))/units_out%energy%factor
        end do
        write(iunit, '(a)')' '
      end do                 
                 
#ifdef HAVE_MPI
  end if
#endif

end subroutine states_write_bands


integer function states_spin_channel(ispin, ik, dim)
  integer, intent(in) :: ispin, ik, dim

  if( ispin < 1 .or. ispin>3 .or. &
      ik<0 .or.                   &
      dim<1 .or. dim>2 .or.       &
      (ispin.ne.3 .and. dim==2)) then
      message(1) = 'Internal bug: Unacceptable entry in states_spin_channel function'
      call write_fatal(1)
  else
    select case(ispin)
      case(1); states_spin_channel = 1
      case(2); states_spin_channel = mod(ik+1, 2)+1
      case(3); states_spin_channel = dim
      case default; states_spin_channel = -1
    end select
  endif

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates:
! p(uist, ist, ik) = < phi0(uist, k) | phi(ist, ik) (t) >
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_projection(u_st, st, m, p)
  type(states_type), intent(in) :: u_st, st
  type(mesh_type),   intent(in) :: m
  CMPLX, intent(out)      :: p(u_st%nst, st%st_start:st%st_end, st%nik)

  integer :: uist, uik, ist, ik

  do ik = 1, st%nik
     do ist = st%st_start, st%st_end
        do uist = 1, u_st%nst
          p(uist, ist, ik) = zstates_dotp( m, st%dim,          &
              cmplx(u_st%X(psi)(:, :, uist, ik), kind=PRECISION) , &
              st%zpsi(:, :, ist, ik) )
        end do
     end do
  end do
end subroutine calc_projection

! This routine (obviously) assumes complex wave-functions
subroutine calc_current(m, st, j)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  FLOAT, intent(out) :: j(3, m%np, st%nspin)
  
  integer :: ik, p, sp, ierr, k
  CMPLX, allocatable :: aux(:,:)
#if defined(HAVE_MPI) && defined(MPI_TD)
  FLOAT, allocatable :: red(:,:,:)
#endif

  call push_sub('calc_current')
  
  if(st%ispin == 2) then
    sp = 2
  else
    sp = 1
  end if
  
  j = M_ZERO
  allocate(aux(conf%dim, m%np))
  
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      call zf_gradient(m, st%zpsi(:, 1, p, ik), aux)
      
      ! spin-up density
      do k = 1, m%np
        j(1:conf%dim, k, 1) = j(1:conf%dim, k, 1) + st%kweights(ik)*st%occ(p, ik)  &
             * aimag(st%zpsi(k, 1, p, ik) * conjg(aux(1:conf%dim, k)))
      end do
        
      ! spin-down density
      if(st%ispin == 2) then
        call zf_gradient(m, st%zpsi(:, 1, p, ik+1), aux)

        do k = 1, m%np
          j(1:conf%dim, k, 2) = j(1:conf%dim, k, 2) + st%kweights(ik+1)*st%occ(p, ik+1) &
             * aimag(st%zpsi(k, 1, p, ik+1) * conjg(aux(1:conf%dim, k)))
        end do

        ! WARNING: the next lines DO NOT work properly
      else if(st%ispin == 3) then ! off-diagonal densities
        call zf_gradient(m, st%zpsi(:, 2, p, ik), aux)

        do k = 1, m%np
          j(1:conf%dim, k, 2) = j(1:conf%dim, k, 2) + st%kweights(ik)*st%occ(p, ik) &
               * aimag(st%zpsi(k, 2, p, ik) * conjg(aux(1:conf%dim, k)))
        end do
      end if

    end do
  end do
  deallocate(aux)

#if defined(HAVE_MPI) && defined(MPI_TD)
  allocate(red(3, m%np, st%nspin))
  call MPI_ALLREDUCE(j(1, 1, 1), red(1, 1, 1), 3*m%np*st%nspin, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  j = red
  deallocate(red)
#endif

  call pop_sub()
end subroutine calc_current

subroutine zstates_project_gs(st, m, p)
  type(states_type), intent(in) :: st
  type(mesh_type), intent(in)   :: m
  CMPLX, intent(out) :: p

  type(states_type) :: stgs

  call push_sub('zstates_project_gs')

  stgs = st
  if(.not.zstates_load_restart("tmp/restart.static", m, stgs)) then
    message(1) = 'Error loading GS in zstates_project_gs'
    call write_fatal(1)
  endif
  p = zstates_mpdotp(m, 1, stgs, st)

  call states_end(stgs)
  call pop_sub()
end subroutine zstates_project_gs

#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
