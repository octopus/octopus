#include "global.h"

module iatom_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t, atom_classical_t, atom_same_species, atom_get_label, atom_end
  use json_m,    only: JSON_OK, json_object_t, json_init, json_get, json_set
  use species_m, only: LABEL_LEN, species_t, species_label, species_type

  implicit none

  private
  public ::                            &
    atom_t,                            &
    operator(==),                      &
    operator(/=),                      &
    atom_classical_t,                  &
    atom_init,                         &
    atom_hash,                         &
    atom_create_data_object,           &
    atom_classical_create_data_object, &
    atom_end

  interface operator(==)
    module procedure atom_equal
  end interface operator(==)

  interface operator(/=)
    module procedure atom_equal
  end interface operator(/=)

  interface atom_init
    module procedure atom_init_from_data_object
    module procedure atom_classical_init_from_data_object
  end interface atom_init

contains

  ! ---------------------------------------------------------
  subroutine atom_init_from_data_object(this, spec, json)
    type(atom_t),            intent(out) :: this
    type(species_t), target, intent(in)  :: spec
    type(json_object_t),     intent(in)  :: json
    !
    integer :: ierr
    !
    PUSH_SUB(atom_init_from_data_object)
    call json_get(json, "label", this%label, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "label" from atom data object.'
      call messages_fatal(1)
      return
    end if
    ASSERT(atom_get_label(this)==species_label(spec))
    this%spec=>spec
    call json_get(json, "x", this%x, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "x" (coordinates) from atom data object.'
      call messages_fatal(1)
      return
    end if
    this%v=M_ZERO
    this%f=M_ZERO
    this%move=.true.
    POP_SUB(atom_init_from_data_object)
    return
  end subroutine atom_init_from_data_object

  ! ---------------------------------------------------------
  elemental function atom_equal(this, that) result(is)
    type(atom_t), intent(in) :: this
    type(atom_t), intent(in) :: that
    !
    logical :: is
    !
    is=(atom_get_label(this)==atom_get_label(that))
    is=is.and.(species_type(this%spec)==species_type(that%spec))
    return
  end function atom_equal

  ! ---------------------------------------------------------
  elemental function atom_not_equal(this, that) result(is)
    type(atom_t), intent(in) :: this
    type(atom_t), intent(in) :: that
    !
    logical :: is
    !
    is=(.not.atom_equal(this, that))
    return
  end function atom_not_equal

  ! ---------------------------------------------------------
  ! Daniel J. Bernstein Hash Function
  elemental function atom_hash(this, size) result(hash)
    type(atom_t), intent(in) :: this
    integer,      intent(in) :: size
    !
    integer :: hash
    !
    character(len=LABEL_LEN) :: label
    integer                  :: i
    !
    hash=5381
    label=atom_get_label(this)
    do i = 1, len_trim(label)
      hash = ieor(33*hash, iachar(label(i:i)))
    end do
    hash = ieor(33*hash, species_type(this%spec))
    hash=modulo(hash, size)+1
    return
  end function atom_hash

  ! ---------------------------------------------------------
  subroutine atom_create_data_object(this, json)
    type(atom_t),        intent(in)  :: this
    type(json_object_t), intent(out) :: json
    !
    PUSH_SUB(atom_create_data_object)
    call json_init(json)
    call json_set(json, "label", trim(adjustl(this%label)))
    call json_set(json, "x", this%x)
    POP_SUB(atom_create_data_object)
    return
  end subroutine atom_create_data_object

  ! ---------------------------------------------------------
  subroutine atom_classical_init_from_data_object(this, json)
    type(atom_classical_t), intent(out) :: this
    type(json_object_t),    intent(in)  :: json
    !
    integer :: ierr
    !
    PUSH_SUB(atom_classical_init_from_data_object)
    call json_get(json, "label", this%label, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "label" from atom classical data object.'
      call messages_fatal(1)
      return
    end if
    call json_get(json, "x", this%x, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "x" (coordinates) from atom classical data object.'
      call messages_fatal(1)
      return
    end if
    this%v=M_ZERO
    this%f=M_ZERO
    call json_get(json, "charge", this%charge, ierr)
    if(ierr/=JSON_OK)then
      message(1) = 'Could not read "charge" from atom classical data object.'
      call messages_fatal(1)
      return
    end if
    POP_SUB(atom_classical_init_from_data_object)
    return
  end subroutine atom_classical_init_from_data_object

  ! ---------------------------------------------------------
  subroutine atom_classical_create_data_object(this, json)
    type(atom_classical_t), intent(in)  :: this
    type(json_object_t),    intent(out) :: json
    !
    PUSH_SUB(atom_classical_create_data_object)
    call json_init(json)
    call json_set(json, "label", trim(adjustl(this%label)))
    call json_set(json, "x", this%x)
    call json_set(json, "charge", this%charge)
    POP_SUB(atom_classical_create_data_object)
    return
  end subroutine atom_classical_create_data_object

end module iatom_m

#define TEMPLATE_NAME atom
#define MODULE_TYPE atom_m
#define TYPE atom_t
#include "tlist.F90"
#undef TYPE
#undef MODULE_TYPE
#undef TEMPLATE_NAME
#define TEMPLATE_NAME atomspec
#define MODULE_TYPE_KEY iatom_m
#define TYPE_KEY atom_t
#define HASH_FUNCTION atom_hash
#define MODULE_TYPE_VAL species_m
#define TYPE_VAL species_t
#include "thash.F90"
#undef TYPE_VAL
#undef MODULE_TYPE_VAL
#undef HASH_FUNCTION
#undef TYPE_KEY
#undef MODULE_TYPE_KEY
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

