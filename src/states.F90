#include "config.h"

module states
use global
use fdf

implicit none

type states_type
  integer :: dim ! Dimension of the state
  integer :: nst ! Number of states in each irreducible subspace
  integer :: nik ! Number of irreducible subspaces
  integer :: ispin ! spin mode

  real(r8), pointer :: rpsi(:,:,:,:)
  complex(r8), pointer :: zpsi(:,:,:,:)

  real(r8), pointer :: eigenval(:,:)
  real(r8), pointer :: occ(:,:)

  integer :: st_start, st_end ! needed for some parallel parts
end type states_type

contains

subroutine states_init(st, val_charge)
  type(states_type), intent(inout) :: st
  real(r8), intent(in) :: val_charge

  real(r8) :: excess_charge
  integer :: nempty

  st%ispin = fdf_integer('SpinComponents', 1)
  if (st%ispin /= 1 .and. st%ispin /= 2 .and. st%ispin /= 4) then
    write(message(1),'(a,i4,a)') "Input: '", st%ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 4)'
    call write_fatal(2)
  end if

#ifdef PERIODIC_1D
  st%nik = fdf_integer('NumberKPoints', 1)
  if (st%nik < 1) then
    write(message(1),'(a,i4,a)') "Input: '", st%nik,"' is not a valid NumberKPoints"
    message(2) = '(NumberKPoints >= 1)'
    call write_fatal(2)
  end if
#else
  st%nik = 1
#endif
  
  excess_charge = fdf_double('ExcessCharge', 0.0_r8)

  nempty = fdf_integer('EmptyStates', 0)
  if (nempty < 0) then
    write(message(1), '(a,i5,a)') "Input: '", nempty, "' is not a valid EmptyStates"
    message(2) = '(0 <= EmptyStates)'
    call write_fatal(2)
  end if
  
  st%nst = -int((val_charge + excess_charge)/2._r8)
  if(st%nst*2._r8 < -(val_charge + excess_charge)) &
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

  ! we now define the occupation
  allocate(st%occ(st%nst, st%nik), st%eigenval(st%nst, st%nik))

  st%st_start = 1; st%st_end = st%nst

end subroutine states_init

subroutine states_end(st)
  type(states_type), intent(inout) :: st
  
  sub_name = 'states_end'; call push_sub()

  if(associated(st%occ)) then
    deallocate(st%occ); nullify(st%occ)
    deallocate(st%eigenval); nullify(st%eigenval)
  end if

  if(associated(st%rpsi)) then
    deallocate(st%rpsi); nullify(st%rpsi)
  end if

  if(associated(st%zpsi)) then
    deallocate(st%zpsi); nullify(st%zpsi)
  end if

  call pop_sub()
end subroutine states_end

end module states
