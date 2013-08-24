#include "global.h"

module iatom_m

  use global_m
  use messages_m
  use profiling_m

  use species_m, only: LABEL_LEN, species_type

  use atom_m, only:             &
    atom_t,                     &
    atom_get_label,             &
    atom_init_from_data_object

  use atom_m, only:  &
    end => atom_end

  use atom_m, only:                       &
    atom_classical_t,                     &
    atom_classical_get_label,             &
    atom_classical_init_from_data_object, &
    atom_classical_end

  implicit none

  private
  public ::       &
    operator(==), &
    operator(/=)

  public ::    &
    atom_t,    &
    atom_init, &
    atom_hash, &
    atom_end

  interface operator(==)
    module procedure atom_equal
    module procedure atom_classical_equal
  end interface operator(==)

  interface operator(/=)
    module procedure atom_not_equal
    module procedure atom_classical_not_equal
  end interface operator(/=)

  interface atom_init
     module procedure atom_init_from_data_object
     module procedure atom_classical_init_from_data_object
  end interface atom_init

  interface atom_end
     module procedure end
     module procedure atom_classical_end
  end interface atom_end

contains

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
  elemental function atom_classical_equal(this, that) result(is)
    type(atom_classical_t), intent(in) :: this
    type(atom_classical_t), intent(in) :: that
    !
    logical :: is
    !
    is=(atom_classical_get_label(this)==atom_classical_get_label(that))
    return
  end function atom_classical_equal

  ! ---------------------------------------------------------
  elemental function atom_classical_not_equal(this, that) result(is)
    type(atom_classical_t), intent(in) :: this
    type(atom_classical_t), intent(in) :: that
    !
    logical :: is
    !
    is=(.not.atom_classical_equal(this, that))
    return
  end function atom_classical_not_equal

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

