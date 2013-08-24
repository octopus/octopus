#include "global.h"
#include "template.h"

#ifdef MODULE_TYPE_KEY
#define MODULE_INVOCATION_KEY use MODULE_TYPE_KEY, only: operator(==), key_t=>TYPE_KEY, HASH_FUNCTION
#define ITYPE_KEY type(key_t)
#else
#define MODULE_INVOCATION_KEY
#define ITYPE_KEY TYPE_KEY
#endif

#ifdef MODULE_TYPE_VAL
#define MODULE_INVOCATION_VAL use MODULE_TYPE_VAL, only: val_t=>TYPE_VAL
#define ITYPE_VAL type(val_t)
#else
#define MODULE_INVOCATION_VAL
#define ITYPE_VAL TYPE_VAL
#endif

#define MODULE_TYPE_1 MODULE_TYPE_KEY
#define TYPE_1 TYPE_KEY
#define MODULE_TYPE_2 MODULE_TYPE_VAL
#define TYPE_2 TYPE_VAL
#include "tpair.F90"
#undef MODULE_TYPE_1
#undef TYPE_1
#undef MODULE_TYPE_2
#undef TYPE_2

#define MODULE_TYPE TEMPLATE(pair_m)
#define TYPE TEMPLATE(pair_t)
#include "tlist.F90"
#undef MODULE_TYPE
#undef TYPE

module TEMPLATE(table_m)

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  MODULE_INVOCATION_KEY
  MODULE_INVOCATION_VAL

  use TEMPLATE(pair_m), only:                     &
    pair_t          => TEMPLATE(pair_t),          &
    pair_init       => TEMPLATE(pair_init),       &
    pair_get_first  => TEMPLATE(pair_get_first),  &
    pair_get_second => TEMPLATE(pair_get_second), &
    pair_set_second => TEMPLATE(pair_set_second), &
    pair_end        => TEMPLATE(pair_end)

  use TEMPLATE(list_m), only:         &
    list_t    => TEMPLATE(list_t),    &
    list_init => TEMPLATE(list_init), &
    list_push => TEMPLATE(list_push), &
    list_pop  => TEMPLATE(list_pop),  &
    list_copy => TEMPLATE(list_copy), &
    list_end  => TEMPLATE(list_end)

  use TEMPLATE(list_m), only:                           &
    list_iterator_t    => TEMPLATE(list_iterator_t),    &
    list_iterator_init => TEMPLATE(list_iterator_init), &
    list_iterator_next => TEMPLATE(list_iterator_next), &
    list_iterator_end  => TEMPLATE(list_iterator_end)

  implicit none

  private
  public ::               &
    TEMPLATE(table_init), &
    TEMPLATE(table_len),  &
    TEMPLATE(table_pop),  &
    TEMPLATE(table_set),  &
    TEMPLATE(table_get),  &
    TEMPLATE(table_copy), &
    TEMPLATE(table_end)
    

  real(kind=wp), parameter :: TABLE_GROWTH_FACTOR = 1.5_wp
  integer,       parameter :: TABLE_INIT_LEN      = 127

  integer, public, parameter :: TABLE_KEY_OK    = 0
  integer, public, parameter :: TABLE_KEY_ERROR =-1

  type, public :: TEMPLATE(table_t)
    private
    integer                             :: size  = 0
    integer                             :: used  = 0
    type(list_t), dimension(:), pointer :: table =>null()
  end type TEMPLATE(table_t)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_init)(this, size)
    type(TEMPLATE(table_t)), intent(out) :: this
    integer,       optional, intent(in)  :: size
    !
    integer :: i
    !
    PUSH_SUB(TEMPLATE(table_init))
    this%used=0
    this%size=TABLE_INIT_LEN
    if(present(size))this%size=size
    SAFE_ALLOCATE(this%table(this%size))
    do i=1, this%size
      call list_init(this%table(i))
    end do
    POP_SUB(TEMPLATE(table_init))
    return
  end subroutine TEMPLATE(table_init)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(table_len)(this) result(len)
    type(TEMPLATE(table_t)), intent(in) :: this
    !
    integer :: len
    !
    len=this%used
    return
  end function TEMPLATE(table_len)

  ! ---------------------------------------------------------
  subroutine table_reallocate(this)
    type(TEMPLATE(table_t)), intent(inout) :: this
    !
    type(TEMPLATE(table_t)) :: buff
    type(pair_t),   pointer :: pair
    ITYPE_KEY,      pointer :: key
    real(kind=wp)           :: need
    integer                 :: i, n
    !
    PUSH_SUB(TEMPLATE(table_reallocate))
    need=real(this%used+1, kind=wp)
    if(this%size<int(TABLE_GROWTH_FACTOR*need))then
      n=max(ceiling((log(need)-log(real(this%size,kind=wp)))/log(TABLE_GROWTH_FACTOR)),1)
      call TEMPLATE(table_init)(buff, ceiling((TABLE_GROWTH_FACTOR**n)*real(this%size,kind=wp)))
      pair=>null()
      do i = 1, this%size
        do
          call list_pop(this%table(i), pair)
          if(.not.associated(pair))exit
          call pair_get_first(pair, key)
          n=HASH_FUNCTION(key, buff%size)
          call list_push(buff%table(n), pair)
          buff%used=buff%used+1
          this%used=this%used-1
          pair=>null()
        end do
        call list_end(this%table(i))
      end do
      ASSERT(this%used==0)
      call TEMPLATE(table_end)(this)
      this%size=buff%size
      this%used=buff%used
      this%table=>buff%table
    end if
    POP_SUB(TEMPLATE(table_reallocate))
    return
  end subroutine table_reallocate

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_pop)(this, key, val)
    type(TEMPLATE(table_t)), intent(inout) :: this
    ITYPE_KEY,              pointer        :: key
    ITYPE_VAL,              pointer        :: val
    !
    type(pair_t), pointer :: pair
    integer               :: i
    !
    PUSH_SUB(TEMPLATE(table_pop))
    nullify(key, val)
    if(this%used>0)then
      pair=>null()
      do i = 1, this%size
        call list_pop(this%table(i), pair)
        if(associated(pair))then
          this%used=this%used-1
          call pair_get_first(pair, key)
          call pair_get_second(pair, val)
          call pair_end(pair)
          SAFE_DEALLOCATE_P(pair)
          exit
        end if
        pair=>null()
      end do
    end if
    POP_SUB(TEMPLATE(table_pop))
    return
  end subroutine TEMPLATE(table_pop)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_set)(this, key, val)
    type(TEMPLATE(table_t)), intent(inout) :: this
    ITYPE_KEY,               intent(in)    :: key
    ITYPE_VAL,               intent(in)    :: val
    !
    type(pair_t), pointer :: pair
    integer               :: n
    !
    PUSH_SUB(TEMPLATE(table_set))
    call table_reallocate(this)
    n=HASH_FUNCTION(key, this%size)
    call table_list_get_pair(this%table(n), key, pair)
    if(associated(pair))then
      call pair_set_second(pair, val)
    else
      SAFE_ALLOCATE(pair)
      call pair_init(pair, key, val)
      call list_push(this%table(n), pair)
      this%used=this%used+1
    end if
    POP_SUB(TEMPLATE(table_set))
    return
  end subroutine TEMPLATE(table_set)

  ! ---------------------------------------------------------
  subroutine table_list_get_pair(this, key, pair)
    type(list_t),          intent(in) :: this
    ITYPE_KEY,     target, intent(in) :: key
    type(pair_t), pointer             :: pair
    !
    type(list_iterator_t) :: iter
    ITYPE_KEY,    pointer :: k
    !
    PUSH_SUB(TEMPLATE(table_list_get_pair))
    nullify(k, pair)
    call list_iterator_init(iter, this)
    do
      call list_iterator_next(iter, pair)
      if(.not.associated(pair))exit
      call pair_get_first(pair, k)
      if(k==key)exit
      nullify(k, pair)
    end do
    call list_iterator_end(iter)
    POP_SUB(TEMPLATE(table_list_get_pair))
    return
  end subroutine table_list_get_pair

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_get)(this, key, val, ierr)
    type(TEMPLATE(table_t)), intent(in)  :: this
    ITYPE_KEY ,              intent(in)  :: key
    ITYPE_VAL,              pointer      :: val
    integer,       optional, intent(out) :: ierr
    !
    type(pair_t), pointer :: pair
    integer               :: n
    !
    PUSH_SUB(TEMPLATE(table_get))
    val=>null()
    if(present(ierr))ierr=TABLE_KEY_ERROR
    n=HASH_FUNCTION(key, this%size)
    call table_list_get_pair(this%table(n), key, pair)
    if(associated(pair))then
      if(present(ierr))ierr=TABLE_KEY_OK
      call pair_get_second(pair, val)
    end if
    POP_SUB(TEMPLATE(table_get))
    return
  end subroutine TEMPLATE(table_get)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_copy)(this_out, this_in)
    type(TEMPLATE(table_t)), intent(out) :: this_out
    type(TEMPLATE(table_t)), intent(in)  :: this_in
    !
    integer :: i
    !
    PUSH_SUB(TEMPLATE(table_copy))
    call TEMPLATE(table_init)(this_out, this_in%size)
    do i = 1, this_out%size
      call list_copy(this_out%table(i), this_in%table(i))
    end do
    this_out%used=this_in%used
    POP_SUB(TEMPLATE(table_copy))
    return
  end subroutine TEMPLATE(table_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(table_end)(this)
    type(TEMPLATE(table_t)), intent(inout) :: this
    !
    type(pair_t), pointer :: pair
    integer               :: i
    !
    PUSH_SUB(TEMPLATE(table_end))
    pair=>null()
    do i = 1, this%size
      do
        call list_pop(this%table(i), pair)
        if(.not.associated(pair))exit
        call pair_end(pair)
        SAFE_DEALLOCATE_P(pair)
        this%used=this%used-1
        pair=>null()
      end do
      call list_end(this%table(i))
    end do
    SAFE_DEALLOCATE_P(this%table)
    this%table=>null()
    ASSERT(this%used==0)
    this%size=0
    POP_SUB(TEMPLATE(table_end))
    return
  end subroutine TEMPLATE(table_end)

!!$  ! ---------------------------------------------------------
!!$  subroutine table_iterator_init(this, table)
!!$    type(table_iterator_t), intent(out) :: this
!!$    type(table_t), target,  intent(in)  :: table
!!$    !
!!$    integer :: i
!!$    !
!!$    call json_object_iterator_nullify(this)
!!$    if(object%used>0)then
!!$      this%size=object%size
!!$      this%table=>object%table
!!$      do i = 1, object%size
!!$        if(associated(object%table(i)%head))then
!!$          this%pos=i
!!$          this%node=>object%table(i)%head
!!$          exit
!!$        end if
!!$      end do
!!$    end if
!!$    return
!!$  end subroutine json_object_iterator_init
!!$
!!$  ! ---------------------------------------------------------
!!$  subroutine table_iterator_next(this, val)
!!$    type(json_object_iterator_t), target, intent(inout) :: this
!!$    type(json_member_t),         pointer, intent(out)   :: member
!!$    !
!!$    integer :: i
!!$    !
!!$    member=>null()
!!$    if(associated(this%node))then
!!$      member=>this%node%member
!!$      if(associated(this%node%next))then
!!$        this%node=>this%node%next
!!$      else
!!$        this%node=>null()
!!$        do i = this%pos+1, this%size
!!$          if(associated(this%table(i)%head))then
!!$            this%pos=i
!!$            this%node=>this%table(i)%head
!!$            exit
!!$          end if
!!$        end do
!!$      end if
!!$    end if
!!$    return
!!$  end subroutine table_iterator_next
!!$
!!$  ! ---------------------------------------------------------
!!$  elemental subroutine table_iterator_copy(this)
!!$    type(json_object_iterator_t), intent(inout) :: this
!!$    !
!!$    this%pos=0
!!$    this%node=>null()
!!$    this%table=>null()
!!$    return
!!$  end subroutine table_iterator_copy
!!$
!!$  ! ---------------------------------------------------------
!!$  elemental subroutine table_iterator_end(this)
!!$    type(json_object_iterator_t), intent(inout) :: this
!!$    !
!!$    this%pos=0
!!$    this%node=>null()
!!$    this%table=>null()
!!$    return
!!$  end subroutine table_iterator_end

end module TEMPLATE(table_m)

#undef ITYPE_VAL
#undef MODULE_INVOCATION_VAL
#undef ITYPE_KEY
#undef MODULE_INVOCATION_KEY

!! Local Variables:
!! mode: f90
!! End:
