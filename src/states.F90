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
use lib_basic_alg
use lib_adv_alg
use math
use mesh
use functions
use mesh_function
use output
use geometry
use crystal

implicit none

private
public :: states_type, &
          states_dim_type, &
          states_init, &
          states_null, &
          states_end,  &
          states_copy, &
          states_generate_random, &
          states_fermi, &
          states_calculate_multipoles, &
          states_eigenvalues_sum, &
          states_write_eigenvalues, &
          states_write_bands, &
          states_spin_channel, &
          calc_projection, &
          calc_current_physical, &
          kpoints_write_info, &
          dcalcdens, zcalcdens, &
          dstates_gram_schmidt, zstates_gram_schmidt, &
          dstates_dotp, zstates_dotp, &
          dstates_nrm2, zstates_nrm2, &
          dstates_residue, zstates_residue, &
          dstates_output, zstates_output, &
          dstates_mpdotp, zstates_mpdotp, &
          dstates_calculate_magnetization, zstates_calculate_magnetization, &
          dstates_calculate_angular, zstates_calculate_angular


type states_type

  type(states_dim_type), pointer :: d

  ! This is (by now!) replicated from the states_dim_type
  integer :: dim
  integer :: nik           ! Number of irreducible subspaces
  integer :: nspin         ! dimension of rho (1, 2 or 4)
  
  integer :: nst           ! Number of states in each irreducible subspace

  ! pointers to the wavefunctions
  FLOAT, pointer :: dpsi(:,:,:,:)
  CMPLX, pointer :: zpsi(:,:,:,:)

  ! the densities and currents (after all we are doing DFT :)
  FLOAT, pointer :: rho(:,:)
  FLOAT, pointer :: j(:,:,:)

  FLOAT, pointer :: rho_core(:)! core charge for nl core corrections

  FLOAT, pointer :: eigenval(:,:) ! obviously the eigenvalues
  logical           :: fixed_occ ! should the occupation numbers be fixed?
  FLOAT, pointer :: occ(:,:)  ! the occupation numbers
  FLOAT, pointer :: mag(:, :, :)

  FLOAT :: qtot    ! (-) The total charge in the system (used in Fermi)

  FLOAT :: el_temp ! electronic temperature for the Fermi function
  FLOAT :: ef      ! the fermi energy

  integer :: st_start, st_end ! needed for some parallel parts

end type states_type

type states_dim_type
  integer :: dim           ! * Dimension of the state (one or two for spinors)
  integer :: nik           ! * Number of irreducible subspaces
  integer :: nik_axis(3)   ! * Number of kpoints per axis
  integer :: ispin         ! * spin mode (unpolarized, spin polarized, spinors)
  integer :: nspin         ! * dimension of rho (1, 2 or 4)
  integer :: spin_channels ! * 1 or 2, wether spin is or not considered.
  logical :: cdft          ! * Are we using Current-DFT or not?
  FLOAT, pointer :: kpoints(:,:) ! * obviously the kpoints
  FLOAT, pointer :: kweights(:)  ! * weights for the kpoint integrations
end type states_dim_type

! Parameters...
integer, public, parameter :: UNPOLARIZED    = 1, &
                              SPIN_POLARIZED = 2, &
                              SPINORS        = 3

interface assignment (=)
  module procedure states_copy
end interface

contains

subroutine states_null(st)
  type(states_type), intent(out) :: st

  nullify(st%dpsi, st%zpsi, st%rho, st%j, st%rho_core, st%eigenval, st%occ, st%mag)
  nullify(st%d); allocate(st%d)
  nullify(st%d%kpoints, st%d%kweights)
end subroutine states_null


subroutine states_init(st, m, geo, val_charge, nlcc)
  type(states_type),   intent(inout) :: st
  type(mesh_type),     intent(IN)    :: m
  type(geometry_type), intent(IN)    :: geo ! this is needed to generate the k points
  FLOAT,               intent(in)    :: val_charge
  logical, optional,   intent(in)    :: nlcc

  FLOAT :: excess_charge, r
  integer :: nempty, i, j
  integer(POINTER_SIZE) :: blk

  call push_sub('states_init')

  call states_null(st)

  call loct_parse_int('SpinComponents', UNPOLARIZED, st%d%ispin)
  if (st%d%ispin < UNPOLARIZED .or. st%d%ispin > SPINORS) then
    write(message(1),'(a,i4,a)') "Input: '", st%d%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if

  call loct_parse_float('ExcessCharge', M_ZERO, excess_charge)

  call loct_parse_int('ExtraStates', 0, nempty)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
    message(2) = '(0 <= ExtraStates)'
    call write_fatal(2)
  end if
  
  st%qtot = -(val_charge + excess_charge)

  select case(st%d%ispin)
  case(UNPOLARIZED)
    st%d%dim = 1
    st%nst = int(st%qtot/2)
    if(st%nst*2 < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%d%nspin = 1
    st%d%spin_channels = 1
  case(SPIN_POLARIZED)
    st%d%dim = 1
    st%nst = int(st%qtot/2)
    if(st%nst*2 < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%d%nik = st%d%nik*2
    st%d%nspin = 2
    st%d%spin_channels = 2
  case(SPINORS)
    st%d%dim = 2
    st%nst = int(st%qtot)
    if(st%nst < st%qtot) st%nst = st%nst + 1
    st%nst = st%nst + nempty
    st%d%nspin = 4
    st%d%spin_channels = 2
  end select

  ! current
  call loct_parse_logical("CurrentDFT", .false., st%d%cdft)
#ifdef COMPLEX_WFNS
  if (st%d%cdft .and. st%d%ispin == SPINORS) then
    message(1) = "Sorry, Current DFT not working yet for spinors"
    call write_fatal(1)
  elseif (st%d%cdft) then
    message(1) = "Info: Using Current DFT"
    call write_info(1)
  end if
#else
  if (st%d%cdft) then
    message(1) = "Cannot use Current DFT with an executable compiled"
    message(2) = "for real wavefunctions."
    call write_fatal(2)
  end if
#endif

  ! For non-periodic systems this should just return the Gamma point
  call states_choose_kpoints(st%d, m, geo)

  ! we now allocate some arrays
  allocate(st%rho(m%np, st%d%nspin), &
           st%occ(st%nst, st%d%nik), &
           st%eigenval(st%nst, st%d%nik))
  if(st%d%ispin == SPINORS) then
    allocate(st%mag(st%nst, st%d%nik, 2))
  end if
  if (st%d%cdft) then
    allocate(st%j(m%np, conf%dim, st%d%nspin))
    st%j = M_ZERO
  end if
  if (present(nlcc)) then
    if(nlcc) allocate(st%rho_core(m%np))
  end if

  occ_fix: if(loct_parse_block("Occupations", blk)==0) then
    ! read in occupations
    st%fixed_occ = .true.

    do i = 1, st%d%nik
      do j = 1, st%nst
        call loct_parse_block_float(blk, i-1, j-1, st%occ(j, i))
      end do
    end do
    call loct_parse_block_end(blk)
  else
    st%fixed_occ = .false.

    ! first guest for occupation...paramagnetic configuration
    if(st%d%ispin == UNPOLARIZED) then
      r = M_TWO
    else
      r = M_ONE
    endif
    st%occ  = M_ZERO
    st%qtot = M_ZERO

    do j = 1, st%nst
      do i = 1, st%d%nik
        st%occ(j, i) = min(r, -(val_charge + excess_charge) - st%qtot)
        st%qtot = st%qtot + st%occ(j, i)

      end do
    end do

    ! read in fermi distribution temperature
    call loct_parse_float('ElectronicTemperature', M_ZERO, st%el_temp)
    st%el_temp = st%el_temp * units_inp%energy%factor
  end if occ_fix

  st%st_start = 1; st%st_end = st%nst
  nullify(st%dpsi, st%zpsi)

  ! Duplicate this
  st%dim   = st%d%dim
  st%nik   = st%d%nik
  st%nspin = st%d%nspin

  call pop_sub()
end subroutine states_init

subroutine states_copy(stout, stin)
  type(states_type), intent(in)  :: stin
  type(states_type), intent(out) :: stout

  call states_null(stout)

  stout%dim           = stin%dim
  stout%nst           = stin%nst
  stout%nik           = stin%nik
  stout%nspin         = stin%nspin
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
  if(associated(stin%j)) then
    allocate(stout%j(size(stin%j, 1), size(stin%j, 2), size(stin%j, 3)))
    stout%j = stin%j
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
  stout%d%dim = stin%d%dim
  stout%d%nik = stin%d%nik
  stout%d%nik_axis(:) = stout%d%nik_axis(:)
  stout%d%ispin = stin%d%ispin
  stout%d%nspin = stin%d%nspin
  stout%d%spin_channels = stin%d%spin_channels
  stout%d%cdft = stin%d%cdft
  if(associated(stin%d%kpoints)) then
    allocate(stout%d%kpoints(size(stin%d%kpoints, 1), size(stin%d%kpoints, 2)))
    stout%d%kpoints = stin%d%kpoints
  endif
  if(associated(stin%d%kweights)) then
    allocate(stout%d%kweights(size(stin%d%kpoints, 1)))
    stout%d%kweights = stin%d%kweights
  endif
end subroutine states_copy

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  call push_sub('states_end')

  if(associated(st%rho)) then
    deallocate(st%rho, st%occ, st%eigenval)
    nullify   (st%rho, st%occ, st%eigenval)
  end if

  if(associated(st%j)) then
    deallocate(st%j)
    nullify(st%j)
  end if

  if(associated(st%rho_core)) then
    deallocate(st%rho_core)
    nullify(st%rho_core)
  end if

  if(st%d%ispin==SPINORS .and. associated(st%mag)) then
    deallocate(st%mag); nullify(st%mag)
  end if

  if(associated(st%dpsi)) then
    deallocate(st%dpsi); nullify(st%dpsi)
  end if

  if(associated(st%zpsi)) then
    deallocate(st%zpsi); nullify(st%zpsi)
  end if

  if(associated(st%d%kpoints)) then
    deallocate(st%d%kpoints); nullify(st%d%kpoints)
  end if

  if(associated(st%d%kweights)) then
    deallocate(st%d%kweights); nullify(st%d%kweights)
  end if
  
  call pop_sub()
end subroutine states_end

! generate a hydrogen s-wavefunction around a random point
subroutine states_generate_random(st, m, ist_start_, ist_end_)
  type(states_type), intent(inout) :: st
  type(mesh_type),   intent(IN)    :: m
  integer, optional, intent(in)    :: ist_start_, ist_end_

  integer :: ist, ik, id, ist_start, ist_end

  call push_sub('states_generate_random')

  ist_start = 1
  if(present(ist_start_)) ist_start = ist_start_
  ist_end = st%nst
  if(present(ist_end_)) ist_end = ist_end_

  do ik = 1, st%nik
    do ist = ist_start, ist_end
      do id = 1, st%dim
        call X(mf_random)(m, st%X(psi)(:, id, ist, ik))
      end do
      st%eigenval(ist, ik) = M_ZERO
    end do
    call X(states_gram_schmidt)(st%nst, m, st%dim, st%X(psi)(:,:,:,ik))
  end do

  call pop_sub()
end subroutine states_generate_random

subroutine states_fermi(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type),   intent(IN)    :: m

! Local variables
  integer :: ie, ik, iter
  integer, parameter :: nitmax = 200
  FLOAT :: drange, t, emin, emax, sumq
  FLOAT, parameter :: tol = CNST(1.0e-10)
  logical :: conv

  call push_sub('fermi')

  if(st%fixed_occ) then ! nothing to do
     ! Calculate magnetizations...
     if(st%d%ispin == SPINORS) then
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

  if(st%d%ispin == SPINORS) then
     sumq = real(st%nst, PRECISION)
  else
     sumq = M_TWO*st%nst
  endif

  t = max(st%el_temp, CNST(1.0e-6))
  st%ef = emax

  conv = .true.
  if (abs(sumq - st%qtot) > tol) conv = .false.
  if (conv) then ! all orbitals are full; nothing to be done
     st%occ = M_TWO/st%d%spin_channels!st%nspin
     ! Calculate magnetizations...
     if(st%d%ispin == SPINORS) then
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
        sumq = sumq + st%d%kweights(ik)/st%d%spin_channels * & !st%nspin * &
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
      st%occ(ie, ik) = stepf((st%eigenval(ie, ik) - st%ef)/t)/st%d%spin_channels!st%nspin
    end do
  end do

  ! Calculate magnetizations...
  if(st%d%ispin == SPINORS) then
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
  type(mesh_type),   intent(IN)  :: m
  type(states_type), intent(IN)  :: st
  FLOAT,             intent(in)  :: pol(:)          ! pol(3)
  FLOAT,             intent(out) :: dipole(:)       ! dipole(st%d%nspin)
  integer, optional, intent(in)  :: lmax
  FLOAT,   optional, intent(out) :: multipole(:, :) ! multipole((lmax + 1)**2, st%nspin)

  integer :: i, is, l, lm, add_lm
  FLOAT :: x(3), r, ylm, mult

  call push_sub('states_calculate_multipoles')

  dipole(1:st%d%nspin) = M_ZERO
  do is = 1, st%d%nspin
    do i = 1, m%np
      call mesh_xyz(m, i, x)
      dipole(is) = dipole(is) + st%rho(i, is)*sum(x(1:conf%dim)*pol(1:conf%dim))*m%vol_pp(i)
    end do
    dipole(is) = dipole(is)

    if(present(lmax).and.present(multipole)) then
      add_lm = 1
      do l = 0, lmax
         do lm = -l, l
            mult = M_ZERO
            do i = 1, m%np
               call mesh_r(m, i, r, x=x)
               ylm = loct_ylm(x(1), x(2), x(3), l, lm) * m%vol_pp(i)
               if(l == 0) then
                 mult = mult + st%rho(i, is) * ylm
               else
                 mult = mult + st%rho(i, is) * ylm * r**l
               end if
            end do
            multipole(add_lm, is) = mult
            add_lm = add_lm + 1
         end do
      end do
    endif

  end do

  call pop_sub(); return
end subroutine states_calculate_multipoles

! function to calculate the eigenvalues sum using occupations as weights
function states_eigenvalues_sum(st) result(e)
  type(states_type), intent(IN) :: st
  FLOAT                         :: e

  integer :: ik

  e = M_ZERO
  do ik = 1, st%d%nik
    e = e + st%d%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
         st%eigenval(st%st_start:st%st_end, ik))
  end do

end function states_eigenvalues_sum

subroutine states_write_eigenvalues(iunit, nst, st, error)
  integer,           intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st
  FLOAT,             intent(IN), optional :: error(nst, st%nik)

  integer ik, j, ns, is
  FLOAT :: o, oplus, ominus

  if(iunit==stdout.and.conf%verbose<=20) return

  ns = 1
  if(st%d%nspin == 2) ns = 2

  message(1) = 'Eigenvalues [' // trim(units_out%energy%abbrev) // ']'
  call write_info(1, iunit)
  if (conf%periodic_dim>0) then 
  end if
  if (st%nik > ns) then
    message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
    call write_info(1, iunit)
  end if

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif
    
    do ik = 1, st%nik, ns
      if(st%nik > ns) then
        write(iunit, '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
           st%d%kpoints(1, ik)*units_out%length%factor, ',',           &
           st%d%kpoints(2, ik)*units_out%length%factor, ',',           &
           st%d%kpoints(3, ik)*units_out%length%factor, ')'
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
            if(st%d%ispin == SPINORS) oplus = M_ZERO; ominus = M_ZERO
          else
            o = st%occ(j, ik+is)
            if(st%d%ispin == SPINORS) then 
              oplus  = st%mag(j, ik+is, 1)
              ominus = st%mag(j, ik+is, 2)
            endif
          end if
      
          write(iunit, '(i4)', advance='no') j
          if (conf%periodic_dim>0) then 
            if(st%d%ispin == SPINORS) then
              write(iunit, '(1x,f12.6,3x,f5.2,a1,f5.2)', advance='no') &
                   (st%eigenval(j, ik)-st%ef)/units_out%energy%factor, oplus, '/', ominus
              if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
            else
              write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                    (st%eigenval(j, ik+is))/units_out%energy%factor, o
             if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
            endif
          else
            if(st%d%ispin == SPINORS) then
              write(iunit, '(1x,f12.6,3x,f5.2,a1,f5.2)', advance='no') &
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
  integer,           intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st

  integer ik, j, ns

  if(iunit==stdout.and.conf%verbose<=20) return
  
  ! shortcuts
  ns = 1
  if(st%d%nspin == 2) ns = 2

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    do j = 1, nst
      do ik = 1, st%nik, ns
        write(iunit, '(1x,3f12.4,3x,f12.6))', advance='yes') &
           st%d%kpoints(:,ik)*units_out%length%factor,            &
           (st%eigenval(j, ik))/units_out%energy%factor
      end do
      write(iunit, '(a)')' '
    end do
    
#ifdef HAVE_MPI
  end if
#endif

end subroutine states_write_bands


! ---------------------------------------------------------
integer function states_spin_channel(ispin, ik, dim)
  integer, intent(in) :: ispin, ik, dim

  ASSERT(ispin >= UNPOLARIZED .or. ispin <= SPINORS)
  ASSERT(ik > 0)
  ASSERT(dim==1 .or. dim==2)
  ASSERT(.not.(ispin.ne.3 .and. dim==2))

  select case(ispin)
  case(1); states_spin_channel = 1
  case(2); states_spin_channel = mod(ik+1, 2)+1
  case(3); states_spin_channel = dim
  case default; states_spin_channel = -1
  end select

end function states_spin_channel

! ---------------------------------------------------------
! This subroutine calculates:
! p(uist, ist, ik) = < phi0(uist, k) | phi(ist, ik) (t) >
! ---------------------------------------------------------
subroutine calc_projection(u_st, st, m, p)
  type(states_type), intent(IN)  :: u_st, st
  type(mesh_type),   intent(IN)  :: m
  CMPLX,             intent(out) :: p(u_st%nst, st%st_start:st%st_end, st%nik)

  integer :: uist, ist, ik
  CMPLX, allocatable :: tmp(:,:)

  allocate(tmp(m%np, st%dim))
  do ik = 1, st%nik
     do ist = st%st_start, st%st_end
        do uist = 1, u_st%nst
          tmp = cmplx(u_st%X(psi)(:, :, uist, ik), kind=PRECISION)
          p(uist, ist, ik) = zstates_dotp(m, st%dim, tmp, st%zpsi(:, :, ist, ik) )
        end do
     end do
  end do
  deallocate(tmp)

end subroutine calc_projection

! This routine (obviously) assumes complex wave-functions
subroutine calc_current_paramagnetic(m, f_der, st, jp)
  type(mesh_type),   intent(in)    :: m
  type(f_der_type),  intent(inout) :: f_der
  type(states_type), intent(in)    :: st
  FLOAT,             intent(out)   :: jp(:,:,:)  ! (m%np, conf%dim, st%d%nspin)

  integer :: ik, p, sp, k
  CMPLX, allocatable :: grad(:,:)
#if defined(HAVE_MPI)
  integer :: ierr
  FLOAT, allocatable :: red(:,:,:)
#endif

  call push_sub('calc_current_paramagnetic')
  
  if(st%d%ispin == SPIN_POLARIZED) then
    sp = 2
  else
    sp = 1
  end if
  
  jp = M_ZERO
  allocate(grad(m%np, conf%dim))
  
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      call zf_gradient(f_der, st%zpsi(:, 1, p, ik), grad)
      
      ! spin-up density
      do k = 1, conf%dim
        jp(:, k, 1) = jp(:, k, 1) + st%d%kweights(ik)*st%occ(p, ik)  &
           * aimag(conjg(st%zpsi(:, 1, p, ik)) * grad(:, k))
      end do

      ! spin-down density
      if(st%d%ispin == SPIN_POLARIZED) then
        call zf_gradient(f_der, st%zpsi(:, 1, p, ik+1), grad)
        
        do k = 1, conf%dim
          jp(:, k, 2) = jp(:, k, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1) &
             * aimag(conjg(st%zpsi(:, 1, p, ik+1)) * grad(:, k))
        end do

        ! WARNING: the next lines DO NOT work properly
      else if(st%d%ispin == SPINORS) then ! off-diagonal densities
        call zf_gradient(f_der, st%zpsi(:, 2, p, ik), grad)

        do k = 1, conf%dim
          jp(:, k, 2) = jp(:, k, 2) + st%d%kweights(ik)*st%occ(p, ik) &
             * aimag(conjg(st%zpsi(:, 2, p, ik)) * grad(:, k))
        end do
      end if

    end do
  end do
  deallocate(grad)

#if defined(HAVE_MPI)
  allocate(red(m%np, conf%dim, st%d%nspin))
  call MPI_ALLREDUCE(jp(1, 1, 1), red(1, 1, 1), conf%dim*m%np*st%d%nspin, &
       MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
  jp = red
  deallocate(red)
#endif

  call pop_sub()
end subroutine calc_current_paramagnetic

subroutine calc_current_physical(m, f_der, st, j)
  type(mesh_type),     intent(in)    :: m
  type(f_der_type),    intent(inout) :: f_der
  type(states_type),   intent(in)    :: st
  FLOAT,               intent(out)   :: j(:,:,:)   ! j(m%np, conf%dim, st%d%nspin)

  call push_sub('calc_current_physical')
  
  ! Paramagnetic contribution to the physical current
  call calc_current_paramagnetic(m, f_der, st, j)

  ! TODO
  ! Diamagnetic contribution to the physical current 

  call pop_sub()
end subroutine calc_current_physical

#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
