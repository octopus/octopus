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

#include "config_F90.h"

module states
use global
use liboct
use io
use math
use output
use mesh

implicit none

type states_type
  integer :: dim           ! Dimension of the state (one or two for spinors)
  integer :: nst           ! Number of states in each irreducible subspace
  integer :: nik           ! Number of irreducible subspaces
  integer :: ispin         ! spin mode (unpolarized, spin polarized, spinors)
  integer :: nspin         ! dimension of rho (1, 2 or 4)
  integer :: spin_channels ! 1 or 2, wether spin is or not considered.

  ! pointers to the wavefunctions
  real(r8), pointer :: dpsi(:,:,:,:)
  complex(r8), pointer :: zpsi(:,:,:,:)

  ! the densities (after all we are doing DFT :)
  real(r8), pointer :: rho(:,:)

  real(r8), pointer :: rho_core(:)! core charge for nl core corrections

  real(r8), pointer :: eigenval(:,:) ! obviously the eigenvalues
  logical           :: fixed_occ ! should the occupation numbers be fixed?
  real(r8), pointer :: occ(:,:)  ! the occupation numbers
  real(r8), pointer :: mag(:, :, :)

  real(r8) :: qtot    ! (-) The total charge in the system (used in Fermi)

  real(r8) :: el_temp ! electronic temperature for the Fermi function
  real(r8) :: ef      ! the fermi energy

  integer :: st_start, st_end ! needed for some parallel parts

  real(r8), pointer :: kpoints(:,:) ! obviously the kpoints
  real(r8), pointer :: kweights(:)  ! weights for the kpoint integrations
end type states_type

! Parameters...
integer, parameter :: UNPOLARIZED    = 1, &
                      SPIN_POLARIZED = 2, &
                      SPINORS        = 3

contains

subroutine states_init(st, m, val_charge)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  real(r8), intent(in) :: val_charge

  real(r8) :: excess_charge, r
  integer :: nempty, i, j
  character(len=80) :: str

  sub_name = 'states_init'; call push_sub()

  call oct_parse_int(C_string('SpinComponents'), UNPOLARIZED, st%ispin)
  if (st%ispin < UNPOLARIZED .or. st%ispin > SPINORS) then
    write(message(1),'(a,i4,a)') "Input: '", st%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if

  if (conf%periodic_dim>0) then
    call oct_parse_int(C_string('NumberKPoints'), 4, st%nik)
    if (st%nik < 1) then
      write(message(1),'(a,i4,a)') "Input: '", st%nik,"' is not a valid NumberKPoints"
      message(2) = '(NumberKPoints >= 1)'
      call write_fatal(2)
    else if (st%nik == 1) then
      message(1) = "The system is periodic, but the number of k points is set to 1"
      message(2) = "output will be generated for the Gamma point only"
      call write_warning(2)
    end if
  else
    st%nik = 1
  end if
  
  call oct_parse_double(C_string('ExcessCharge'), 0.0_r8, excess_charge)

  if(conf%periodic_dim>0) then
    i = 2 ! we take two extra states, just in case we have a metal
  else
    i = 0 ! no extra states
  end if

  call oct_parse_int(C_string('ExtraStates'), i, nempty)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid ExtraStates"
    message(2) = '(0 <= ExtraStates)'
    call write_fatal(2)
  end if
  
!!$  st%qtot = -(val_charge + excess_charge)
!!$  st%nst  = int(st%qtot/2._r8)
!!$  if(st%nst*2._r8 < st%qtot) &
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

  str = C_string("Occupations")
  occ_fix: if(oct_parse_isdef(str) .ne. 0) then
    ! read in occupations
    st%fixed_occ = .true.

    do i = 1, st%nik
      do j = 1, st%nst
        call oct_parse_block_double(str, i-1, j-1, st%occ(j, i))
      end do
    end do
  else
    st%fixed_occ = .false.

    ! first guest for occupation...paramagnetic configuration
    if(st%ispin == UNPOLARIZED) then
      r = 2._r8
    else
      r = 1._r8
    endif
!!$    r = 2.0_r8 / st%nspin
    st%occ  = 0.0_r8
    st%qtot = 0.0_r8
    do j = 1, st%nst
      do i = 1, st%nik
        st%occ(j, i) = min(r, -(val_charge + excess_charge) - st%qtot)
        st%qtot = st%qtot + st%occ(j, i)
      end do
    end do

    ! read in fermi distribution temperature
    call oct_parse_double(C_string('ElectronicTemperature'), 0.0_r8, st%el_temp)
  end if occ_fix

  st%st_start = 1; st%st_end = st%nst
  nullify(st%dpsi, st%zpsi, st%rho_core)

  call pop_sub()
end subroutine states_init

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  sub_name = 'states_end'; call push_sub()

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
  real(r8) :: a(3), rnd, r

  sub_name = 'states_generate_random'; call push_sub()

  ist_s = 1
  if(present(ist_start)) ist_s = ist_start

  st%R_FUNC(psi) (0, :, :, :) = M_ZERO
  do ik = 1, st%nik
    do ist = ist_s, st%nst
      do id = 1, st%dim
         call R_FUNC(mesh_random)(m, st%R_FUNC(psi)(1:m%np, id, ist, ik))
      end do
    end do
    call R_FUNC(states_gram_schmidt)(st%nst, m, st%dim, st%R_FUNC(psi)(:,:,:,ik))
  end do
  st%eigenval = M_ZERO

  call pop_sub(); return
end subroutine states_generate_random

subroutine states_fermi(st, m)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(in) :: m

! Local variables
  integer :: ie, ik, iter
  integer, parameter :: nitmax = 200
  real(r8) :: drange, t, emin, emax, sumq
  real(r8), parameter :: tol = 1.0e-10_r8
  logical :: conv

  sub_name = 'fermi'; call push_sub()

  if(st%fixed_occ) then ! nothing to do
     ! Calculate magnetizations...
     if(st%ispin == 3) then
       do ik = 1, st%nik
         do ie = 1, st%nst
            st%mag(ie, ik, 1) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
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
     sumq = real(st%nst, r8)
  else
     sumq = 2.0_r8*st%nst
  endif

  t = max(st%el_temp, 1.0e-6_r8)
  st%ef = emax

  conv = .true.
  if (abs(sumq - st%qtot) > tol) conv = .false.
  if (conv) then ! all orbitals are full; nothing to be done
     st%occ = 2.0_r8/st%spin_channels!st%nspin
     ! Calculate magnetizations...
     if(st%ispin == 3) then
       do ik = 1, st%nik
         do ie = 1, st%nst
            st%mag(ie, ik, 1) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
            st%mag(ie, ik, 2) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
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

  drange = t*sqrt(-log(tol*.01_r8))

  emin = emin - drange
  emax = emax + drange

  do iter = 1, nitmax
    st%ef = 0.5_r8*(emin + emax)
    sumq  = 0.0_r8

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
          st%mag(ie, ik, 1) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 1, ie, ik))**2 * st%occ(ie, ik)
          st%mag(ie, ik, 2) = R_FUNC(mesh_nrm2) (m, st%R_FUNC(psi)(1:m%np, 2, ie, ik))**2 * st%occ(ie, ik)
       enddo
    enddo
  endif
  
  call pop_sub()
  return
end subroutine states_fermi

subroutine states_calculate_multipoles(m, st, pol, lmax, dipole, multipole)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  integer, intent(in) :: lmax
  real(r8), intent(in) :: pol(3)
  real(r8), intent(out) :: dipole(st%nspin), multipole((lmax + 1)**2, st%nspin)

  integer :: i, is, l, lm, add_lm
  real(r8) :: x(3), r, ylm, mult

  dipole = 0._r8
  do is = 1, st%nspin
    do i = 1, m%np
      call mesh_xyz(m, i, x)
      dipole(is) = dipole(is) + st%rho(i, is)*sum(x*pol)
    end do
    dipole(is) = dipole(is) * m%vol_pp


    add_lm = 1
    do l = 0, lmax
      do lm = -l, l
        mult = 0._r8

        do i = 1, m%np
          call mesh_r(m, i, r, x=x)
          ylm = oct_ylm(x(1), x(2), x(3), l, lm)
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
  end do
end subroutine states_calculate_multipoles

subroutine states_write_eigenvalues(iunit, nst, st, error)
  integer, intent(in) :: iunit, nst
  type(states_type), intent(IN) :: st
  real(r8), intent(in), optional :: error(nst, st%nik)

  integer ik, j, ns, is
  real(r8) :: o, oplus, ominus

  if(iunit==stdout.and.conf%verbose<=20) return

  message(1) = 'Eigenvalues ['//trim(units_out%energy%abbrev)//']'
  call write_info(1, iunit)

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    ns = 1
    if(st%nspin == 2) ns = 2

    do ik = 1, st%nik, ns
      if(st%nik > ns) then
        write(iunit, '(3(a,f12.6),a)') 'k = (', st%kpoints(1, ik), ',', st%kpoints(2, ik), ',', st%kpoints(3, ik), ')'
      end if
      
      do is = 1, ns
        write(iunit, '(a4)', advance='no') '#'
        write(iunit, '(1x,a12,3x,a12,2x,a10,i3,a1)', advance='no') &
             ' Eigenvalue', 'Occupation ', 'Error (', is, ')'
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
          if(st%ispin == 3) then
            write(iunit, '(1x,f12.6,3x,f5.3,a1,f5.3)', advance='no') &
                 st%eigenval(j, ik+is)/units_out%energy%factor, oplus, '/', ominus
            if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik+is), ')'
          else
            write(iunit, '(1x,f12.6,3x,f12.6)', advance='no') &
                 st%eigenval(j, ik)/units_out%energy%factor, o
            if(present(error)) write(iunit, '(a7,es7.1,a1)', advance='no')'      (', error(j, ik), ')'
          endif
        end do
        write(iunit, '(1x)', advance='yes')
      end do
    end do

#ifdef HAVE_MPI
  end if
#endif
  
end subroutine states_write_eigenvalues

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
  complex(r8), intent(out)      :: p(u_st%nst, st%st_start:st%st_end, st%nik)

  integer :: uist, uik, ist, ik

  do ik = 1, st%nik
     do ist = st%st_start, st%st_end
        do uist = 1, u_st%nst
          p(uist, ist, ik) = zstates_dotp( m, st%dim,                                          &
                                           cmplx(u_st%R_FUNC(psi)(1:, :, uist, ik), kind=r8) , &
                                           st%zpsi(1:, :, ist, ik) )
        end do
     end do
  end do
end subroutine calc_projection

! This routine (obviously) assumes complex wave-functions
subroutine calc_current(m, st, j)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: j(3, m%np, st%nspin)
  
  integer :: ik, p, sp, ierr, k
  complex(r8), allocatable :: aux(:,:)
  
  sub_name = 'calc_current'; call push_sub()
  
  if(st%ispin == 2) then
    sp = 2
  else
    sp = 1
  end if
  
  j = M_ZERO
  allocate(aux(3, m%np))
  
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      call zmesh_derivatives(m, st%zpsi(0:m%np, 1, p, ik), grad=aux(:,:))
      
      ! spin-up density
      do k = 1, m%np
        j(1:conf%dim, k, 1) = j(1:conf%dim, k, 1) + st%kweights(ik)*st%occ(p, ik)  &
             * aimag(st%zpsi(k, 1, p, ik) * conjg(aux(1:conf%dim, k)))
      end do
        
      ! spin-down density
      if(st%ispin == 2) then
        call zmesh_derivatives(m, st%zpsi(0:m%np, 1, p, ik+1), aux(:,:))

        do k = 1, m%np
          j(1:conf%dim, k, 2) = j(1:conf%dim, k, 2) + st%kweights(ik+1)*st%occ(p, ik+1) &
             * aimag(st%zpsi(k, 1, p, ik+1) * conjg(aux(1:conf%dim, k)))
        end do

        ! WARNING: the next lines DO NOT work properly
      else if(st%ispin == 3) then ! off-diagonal densities
        call zmesh_derivatives(m, st%zpsi(0:m%np, 2, p, ik), aux(:,:))

        do k = 1, m%np
          j(1:conf%dim, k, 2) = j(1:conf%dim, k, 2) + st%kweights(ik)*st%occ(p, ik) &
               * aimag(st%zpsi(k, 2, p, ik) * conjg(aux(1:conf%dim, k)))
        end do
      end if

    end do
  end do

#if defined(HAVE_MPI) && defined(MPI_TD)
  ! reduce current (assumes memory is contiguous)
  call MPI_ALLREDUCE(j(1, 1), aux(1, 1), m%np*st%nspin, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  j = aux
#endif

  deallocate(aux)
  call pop_sub()
end subroutine calc_current

#include "states_kpoints.F90"

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
