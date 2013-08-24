#include "global.h"

module iatom_m

  use global_m
  use messages_m
  use profiling_m

  use species_m, only: LABEL_LEN, species_type

  use atom_m, only: &
    operator(==),   &
    operator(/=)

  use atom_m, only: &
    atom_t,         &
    atom_get_label, &
    atom_end

  use atom_m, only:                       &
    atom_init => atom_create_data_object

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

contains

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

