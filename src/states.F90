#include "config.h"

module states
use global
use liboct
use io
use math
use mesh

implicit none

! If you change everything in this structure, do not forget to
! change states_dup
type states_type
  integer :: dim ! Dimension of the state
  integer :: nst ! Number of states in each irreducible subspace
  integer :: nik ! Number of irreducible subspaces
  integer :: ispin ! spin mode
  integer :: nspin ! dimension of rho

  ! pointers to the wavefunctions
  real(r8), pointer :: dpsi(:,:,:,:)
  complex(r8), pointer :: zpsi(:,:,:,:)

  ! the density (after all we are doing DFT :)
  real(r8), pointer :: rho(:,:)

  real(r8), pointer :: eigenval(:,:) ! obviously the eigenvalues

  logical           :: fixed_occ ! should the occupation numbers be fixed?
  real(r8), pointer :: occ(:,:)  ! the occupation numbers

  real(r8) :: qtot    ! (-) The total charge in the system (used in Fermi)

  real(r8) :: el_temp ! electronic temperature for the Fermi function
  real(r8) :: ef      ! the fermi energy

  integer :: st_start, st_end ! needed for some parallel parts
end type states_type

contains

subroutine states_init(st, m, val_charge)
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  real(r8), intent(in) :: val_charge

  real(r8) :: excess_charge, r
  integer :: nempty, i, j
  character(len=80) :: str

  sub_name = 'states_init'; call push_sub()

  call oct_parse_int(C_string('SpinComponents'), 1, st%ispin)
  if (st%ispin /= 1 .and. st%ispin /= 2 .and. st%ispin /= 4) then
    write(message(1),'(a,i4,a)') "Input: '", st%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 4)'
    call write_fatal(2)
  end if

#ifdef PERIODIC_1D
  call oct_parse_int(C_string('NumberKPoints'), 1, st%nik)
  if (st%nik < 1) then
    write(message(1),'(a,i4,a)') "Input: '", st%nik,"' is not a valid NumberKPoints"
    message(2) = '(NumberKPoints >= 1)'
    call write_fatal(2)
  end if
#else
  st%nik = 1
#endif
  
  call oct_parse_double(C_string('ExcessCharge'), 0.0_r8, excess_charge)

  call oct_parse_int(C_string('EmptyStates'), 0, nempty)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid EmptyStates"
    message(2) = '(0 <= EmptyStates)'
    call write_fatal(2)
  end if
  
  st%qtot = -(val_charge + excess_charge)
  st%nst  = int(st%qtot/2._r8)
  if(st%nst*2._r8 < st%qtot) &
       st%nst = st%nst + 1
  st%nst = st%nst + nempty

  select case(st%ispin)
  case(1)
    st%dim = 1
  case(2)
    st%dim = 1
    st%nik = st%nik*2
  case(4)
    st%dim = 2
    st%nik = st%nik*2
  end select

  ! we now allocate some arrays
  st%nspin = min(2, st%ispin)
  allocate(st%rho(m%np, st%nspin), &
       st%occ(st%nst, st%nik), st%eigenval(st%nst, st%nik))

  str = C_string("Occupations")
  occ_fix: if(oct_parse_isdef(str) .ne. 0) then
    print *, "HHHELLOO"
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
    ! TODO: for polimers this has to be changed
    r = 2.0_r8
    if (st%ispin > 1) r = 1.0_r8
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

  call pop_sub()
end subroutine states_init

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  sub_name = 'states_end'; call push_sub()

  if(associated(st%rho)) then
    deallocate(st%rho); nullify(st%rho)
    deallocate(st%occ); nullify(st%occ)
    deallocate(st%eigenval); nullify(st%eigenval)
  end if

  if(associated(st%dpsi)) then
    deallocate(st%dpsi); nullify(st%dpsi)
  end if

  if(associated(st%zpsi)) then
    deallocate(st%zpsi); nullify(st%zpsi)
  end if

  call pop_sub()
end subroutine states_end

subroutine states_dup(i, o)
  type(states_type), intent(IN) :: i
  type(states_type), intent(out) :: o

  o%dim   = i%dim
  o%nst   = i%nst
  o%nik   = i%nik
  o%ispin = i%ispin
  o%nspin = i%nspin

  o%dpsi     => i%dpsi
  o%zpsi     => i%zpsi
  o%rho      => i%rho
  o%eigenval => i%eigenval
  o%occ      => i%occ

  o%qtot     = i%qtot
  o%el_temp  = i%el_temp
  o%ef       = i%ef
  o%st_start = i%st_start
  o%st_end   = i%st_end
end subroutine states_dup

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

  st%dpsi(0, :, :, :) = 0.0_r8
  do ik = 1, st%nik
    do ist = ist_s, st%nst
      do id = 1, st%dim
        call quickrnd(iseed, rnd)
        a(1) = 2.0_r8*(2*rnd - 1)
        call quickrnd(iseed, rnd)
        a(2) = 2.0_r8*(2*rnd - 1)
        call quickrnd(iseed, rnd)
        a(3) = 2.0_r8*(2*rnd - 1)

        do i = 1, m%np
          call mesh_r(m, i, r, a=a)
          st%R_FUNC(psi)(i, id, ist, ik) = exp(-0.5_r8*r*r)
        end do
      end do
    end do
    call R_FUNC(states_gram_schmidt)(st%nst, m, st%dim, st%R_FUNC(psi)(:,:,:,ik))
  end do
  st%eigenval = 0._r8

  call pop_sub()
end subroutine states_generate_random

subroutine states_fermi(st)
  type(states_type), intent(inout) :: st

! Local variables
  integer :: ie, ik, iter
  integer, parameter :: nitmax = 200
  real(r8) :: drange, t, emin, emax, sumq
  real(r8), parameter :: tol = 1.0e-10_r8
  logical :: conv

  sub_name = 'fermi'; call push_sub()

  if(st%fixed_occ) then ! nothing to do
    call pop_sub()
    return
  end if

! Initializations
  emin = minval(st%eigenval)
  emax = maxval(st%eigenval)

  sumq = 2.0_r8*st%nst
  t = max(st%el_temp, 1.0e-6_r8)
  st%ef = emax

  conv = .true.
  if (abs(sumq - st%qtot) > tol) conv = .false.
  if (conv) then ! all orbitals are full; nothing to be done
     st%occ = 2.0_r8/min(2, st%ispin)

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
        sumq = sumq + stepf((st%eigenval(ie, ik) - st%ef)/t)/min(2, st%ispin)
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
      st%occ(ie, ik) = stepf((st%eigenval(ie, ik) - st%ef)/t)/min(2, st%ispin)
    end do
  end do
 
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
  real(r8) :: x(3), r, ylm, mult, lixo(3)

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

  integer ik, j
  real(r8) :: o(st%nik)

  if(iunit==stdout.and.conf%verbose<=20) return

  message(1) = 'Eigenvalues ['//trim(units_out%energy%abbrev)//']'
  call write_info(1, iunit)

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    write(iunit, '(a4)', advance='no') '#'
    do ik = 1, st%nik
      write(iunit, '(1x,a12,1x,a12,2x,a10,i1,a1)', advance='no') &
           'Eigenvalue ', 'Occupation ', 'Error (', ik, ')'
    end do
    write(iunit, '(1x)', advance='yes')
    
    do j = 1, nst
      if(j > st%nst) then
        o = 0._r8
      else
        o = st%occ(j, :)
      end if
      
      write(iunit, '(i4)', advance='no') j
      do ik = 1, st%nik
        write(iunit, '(1x,f12.6,1x,f12.6)', advance='no') &
             st%eigenval(j, ik)/units_out%energy%factor, o(ik)
        if(present(error)) then
          write(iunit, '(a2,f12.8,a1)', advance='no')' (', error(j, ik), ')'
        end if
      end do
      write(iunit, '(1x)', advance='yes')
    end do

#ifdef HAVE_MPI
  end if
#endif
  
end subroutine states_write_eigenvalues

#include "undef.F90"
#include "real.F90"
#include "states_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "states_inc.F90"
#include "undef.F90"

end module states
