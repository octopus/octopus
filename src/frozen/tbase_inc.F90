#include "global.h"
#include "util.h"

!BASE: BASE_TEMPLATE_NAME
!BASE: BASE_TYPE_NAME
!BASE: BASE_TYPE_MODULE_NAME

#if defined(BASE_TEMPLATE_NAME)
#if !defined(BASE_TYPE_NAME)
#define BASE_TYPE_NAME DECORATE(BASE_TEMPLATE_NAME,t)
#endif
#if !defined(BASE_TYPE_MODULE_NAME)
#define BASE_TYPE_MODULE_NAME DECORATE(BASE_TEMPLATE_NAME,oct_m)
#endif
#else
#error "'BASE_TEMPLATE_NAME' must be defined"
#endif

#undef BASE_INCLUDE_MODULE
#if !defined(BASE_INCLUDE_PREFIX) && !defined(BASE_INCLUDE_HEADER) && !defined(BASE_INCLUDE_BODY)
#define BASE_INCLUDE_MODULE
#endif

#if defined(BASE_INCLUDE_PREFIX) && defined(BASE_INCLUDE_HEADER)
#error "Only one off 'BASE_INCLUDE_PREFIX' or 'BASE_INCLUDE_HEADER' can be defined."
#endif

#if defined(BASE_INCLUDE_PREFIX) && defined(BASE_INCLUDE_BODY)
#error "Only one off 'BASE_INCLUDE_PREFIX' or 'BASE_INCLUDE_BODY' can be defined."
#endif

#if defined(BASE_INCLUDE_HEADER) && defined(BASE_INCLUDE_BODY)
#error "Only one off 'BASE_INCLUDE_HEADER' or 'BASE_INCLUDE_BODY' can be defined."
#endif

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

#if defined(BASE_INCLUDE_MODULE)

module TEMPLATE(base_oct_m)

  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m

#endif
#if defined(BASE_INCLUDE_PREFIX) || defined(BASE_INCLUDE_MODULE)

#define LIST_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

#endif
#if defined(BASE_INCLUDE_MODULE)

  implicit none

  private

#endif
#if defined(BASE_INCLUDE_HEADER) || defined(BASE_INCLUDE_MODULE)

#define LIST_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

  public ::             &
    TEMPLATE(NAME_LEN)

  public ::                &
    TEMPLATE(OK),          &
    TEMPLATE(KEY_ERROR),   &
    TEMPLATE(EMPTY_ERROR)

  public ::               &
    TEMPLATE(iterator_t)

  public ::          &
    SPECIAL(build),  &
    SPECIAL(apply),  &
    SPECIAL(reduce)
    
  public ::         &
    TEMPLATE(sets), &
    TEMPLATE(gets), &
    TEMPLATE(dels), &
    TEMPLATE(next)

  integer, parameter :: TEMPLATE(NAME_LEN) = TEMPLATE(DICT_NAME_LEN)

  integer, parameter :: TEMPLATE(OK)          = TEMPLATE(DICT_OK)
  integer, parameter :: TEMPLATE(KEY_ERROR)   = TEMPLATE(DICT_KEY_ERROR)
  integer, parameter :: TEMPLATE(EMPTY_ERROR) = TEMPLATE(DICT_EMPTY_ERROR)

#if !defined(BASE_EXTENDED_ITERATOR_TYPE)

  type :: TEMPLATE(iterator_t)
    private
    type(BASE_TYPE_NAME),   pointer :: self =>null()
    type(TEMPLATE(dict_iterator_t)) :: iter
  end type TEMPLATE(iterator_t)

#endif

  interface TEMPLATE(del)
    module procedure TEMPLATE(del_type)
    module procedure TEMPLATE(del_pass)
  end interface TEMPLATE(del)

  interface TEMPLATE(init)
    module procedure TEMPLATE(iterator_init_type)
    module procedure TEMPLATE(iterator_init_copy)
  end interface TEMPLATE(init)

  interface TEMPLATE(gets)
    module procedure TEMPLATE(gets_type)
    module procedure TEMPLATE(gets_config)
  end interface TEMPLATE(gets)

  interface TEMPLATE(next)
    module procedure TEMPLATE(iterator_next_name)
    module procedure TEMPLATE(iterator_next_type)
    module procedure TEMPLATE(iterator_next_name_type)
  end interface TEMPLATE(next)

  interface TEMPLATE(copy)
    module procedure TEMPLATE(iterator_copy)
  end interface TEMPLATE(copy)

  interface TEMPLATE(end)
    module procedure TEMPLATE(end_type)
    module procedure TEMPLATE(end_pass)
    module procedure TEMPLATE(iterator_end)
  end interface TEMPLATE(end)

#endif
#if defined(BASE_INCLUDE_MODULE)

contains

#endif
#if defined(BASE_INCLUDE_BODY) || defined(BASE_INCLUDE_MODULE)

#define LIST_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(del_type)(this)
    type(BASE_TYPE_NAME), pointer :: this

    PUSH_SUB(TEMPLATE(del_type))

    call TEMPLATE(del)(this, TEMPLATE(end_type))
    
    POP_SUB(TEMPLATE(del_type))
  end subroutine TEMPLATE(del_type)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(del_pass)(this, finis)
    type(BASE_TYPE_NAME), pointer :: this

    interface
      subroutine finis(this)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
      end subroutine finis
    end interface
    
    PUSH_SUB(TEMPLATE(del_pass))

    if(associated(this))then
      if(associated(this%prnt))then
        call TEMPLATE(list_del)(this%prnt%list, this)
        call finis(this)
        call SPECIAL(del)(this)
      end if
    end if

    POP_SUB(TEMPLATE(del_pass))
  end subroutine TEMPLATE(del_pass)

  ! ---------------------------------------------------------
  subroutine SPECIAL(build)(this, build)
    type(BASE_TYPE_NAME), intent(inout) :: this

    interface
      subroutine build(this, name, that)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
        character(len=*),     intent(in)    :: name
        type(BASE_TYPE_NAME), intent(in)    :: that
      end subroutine build
    end interface

    type(TEMPLATE(dict_iterator_t))   :: iter
    character(len=TEMPLATE(NAME_LEN)) :: name
    type(BASE_TYPE_NAME),     pointer :: subs

    PUSH_SUB(SPECIAL(build))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_init)(iter, this%dict)
    do
      nullify(subs)
      call TEMPLATE(dict_next)(iter, name, subs)
      if(.not.associated(subs))exit
      call build(this, trim(adjustl(name)), subs)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(subs)

    POP_SUB(SPECIAL(build))
  end subroutine SPECIAL(build)

  ! ---------------------------------------------------------
  subroutine SPECIAL(apply)(this, operation, parent)
    type(BASE_TYPE_NAME), intent(inout) :: this
    logical,    optional, intent(in)    :: parent

    interface
      subroutine operation(this)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
      end subroutine operation
    end interface

    type(TEMPLATE(dict_iterator_t)) :: iter
    type(BASE_TYPE_NAME),   pointer :: subs
    logical                         :: prnt

    PUSH_SUB(SPECIAL(apply))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_init)(iter, this%dict)
    do
      nullify(subs)
      call TEMPLATE(dict_next)(iter, subs)
      if(.not.associated(subs))exit
      call operation(subs)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(subs)
    prnt = .true.
    if(present(parent)) prnt = parent
    if(prnt) call operation(this)

    POP_SUB(SPECIAL(apply))
  end subroutine SPECIAL(apply)

  ! ---------------------------------------------------------
  subroutine SPECIAL(reduce)(this, operation)
    type(BASE_TYPE_NAME), intent(inout) :: this

    interface
      subroutine operation(this, that)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
        type(BASE_TYPE_NAME), intent(in)    :: that
      end subroutine operation
    end interface

    type(TEMPLATE(dict_iterator_t)) :: iter
    type(BASE_TYPE_NAME),   pointer :: subs

    PUSH_SUB(SPECIAL(reduce))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_init)(iter, this%dict)
    do
      nullify(subs)
      call TEMPLATE(dict_next)(iter, subs)
      if(.not.associated(subs))exit
      call operation(this, subs)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(subs)

    POP_SUB(SPECIAL(reduce))
  end subroutine SPECIAL(reduce)

#if defined(BASE_LEAF_TYPE)

  ! ---------------------------------------------------------
  subroutine SPECIAL(sets)(this, name, that, config, lock, active)
    type(BASE_TYPE_NAME),          intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(BASE_TYPE_NAME),          intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    PUSH_SUB(SPECIAL(sets))

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    if(present(config)) continue
    if(present(lock)) continue
    if(present(active)) continue
    
    POP_SUB(SPECIAL(sets))
  end subroutine SPECIAL(sets)
    
  ! ---------------------------------------------------------
  subroutine SPECIAL(dels)(this, name, ierr)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    integer,              intent(out)   :: ierr

    PUSH_SUB(SPECIAL(dels))

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ierr = TEMPLATE(OK)

    POP_SUB(SPECIAL(dels))
  end subroutine SPECIAL(dels)

#endif

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(sets)(this, name, that)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(BASE_TYPE_NAME), intent(in)    :: that

    PUSH_SUB(TEMPLATE(sets))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_set)(this%dict, trim(adjustl(name)), that)
    call SPECIAL(sets)(this, trim(adjustl(name)), that)

    POP_SUB(TEMPLATE(sets))
  end subroutine TEMPLATE(sets)
    
  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(gets_type)(this, name, that)
    type(BASE_TYPE_NAME),  intent(in) :: this
    character(len=*),      intent(in) :: name
    type(BASE_TYPE_NAME), pointer     :: that

    type(BASE_TYPE_NAME), pointer :: subs
    integer                       :: ipos
    
    PUSH_SUB(TEMPLATE(gets_type))

    ASSERT(associated(this%config))
    nullify(that, subs)
    ipos = index(name, "/")
    if(ipos>0)then
      call TEMPLATE(dict_get)(this%dict, trim(adjustl(name(1:ipos-1))), subs)
      if(associated(subs)) call TEMPLATE(gets)(subs, trim(adjustl(name(ipos+1:))), that)
      nullify(subs)
    else
      call TEMPLATE(dict_get)(this%dict, trim(adjustl(name)), that)
    end if

    POP_SUB(TEMPLATE(gets_type))
  end subroutine TEMPLATE(gets_type)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(dels)(this, name, ierr)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    integer,    optional, intent(out)   :: ierr

    type(BASE_TYPE_NAME), pointer :: subs
    integer                       :: jerr

    PUSH_SUB(TEMPLATE(dels))

    ASSERT(associated(this%config))
    nullify(subs)
    call TEMPLATE(dict_del)(this%dict, trim(adjustl(name)), subs, jerr)
    if(associated(subs))then
      call SPECIAL(dels)(this, trim(adjustl(name)), jerr)
      call TEMPLATE(del)(subs)
      nullify(subs)
    end if
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(dels))
  end subroutine TEMPLATE(dels)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(gets_config)(this, name, that)
    type(BASE_TYPE_NAME), intent(in) :: this
    character(len=*),     intent(in) :: name
    type(json_object_t), pointer     :: that
    
    type(BASE_TYPE_NAME), pointer :: subs

    PUSH_SUB(TEMPLATE(gets_config))

    nullify(that, subs)
    call TEMPLATE(gets)(this, trim(adjustl(name)), subs)
    if(associated(subs)) call TEMPLATE(get)(subs, that)
    nullify(subs)

    POP_SUB(TEMPLATE(gets_config))
  end subroutine TEMPLATE(gets_config)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(end_type)(this)
    type(BASE_TYPE_NAME), intent(inout) :: this

    PUSH_SUB(TEMPLATE(end_type))

    call TEMPLATE(end)(this, TEMPLATE(end_type))
    
    POP_SUB(TEMPLATE(end_type))
  end subroutine TEMPLATE(end_type)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(end_pass)(this, finis)
    type(BASE_TYPE_NAME), intent(inout) :: this

    interface
      subroutine finis(this)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
      end subroutine finis
    end interface
    
    type(BASE_TYPE_NAME), pointer :: subs

    PUSH_SUB(TEMPLATE(end_pass))

    do
      nullify(subs)
      call TEMPLATE(list_pop)(this%list, subs)
      if(.not.associated(subs))exit
      call TEMPLATE(del)(subs, finis)
    end do
    nullify(subs)
    call SPECIAL(end)(this)

    POP_SUB(TEMPLATE(end_pass))
  end subroutine TEMPLATE(end_pass)

#if !defined(BASE_EXTENDED_ITERATOR_TYPE)

  ! ---------------------------------------------------------
  subroutine SPECIAL(iterator,init)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(BASE_TYPE_NAME),       intent(in)  :: that

    PUSH_SUB(SPECIAL(iterator,init))

    ASSERT(.not.associated(this%self))
    ASSERT(associated(that%config))
    
    POP_SUB(SPECIAL(iterator,init))
  end subroutine SPECIAL(iterator,init)

  ! ---------------------------------------------------------
  subroutine SPECIAL(iterator,copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(TEMPLATE(iterator_t)), intent(in)    :: that

    PUSH_SUB(SPECIAL(iterator,copy))

    if(associated(this%self))then
      ASSERT(associated(this%self, that%self))
    end if
    
    POP_SUB(SPECIAL(iterator,copy))
  end subroutine SPECIAL(iterator,copy)

  ! ---------------------------------------------------------
  subroutine SPECIAL(iterator,end)(this)
    type(TEMPLATE(iterator_t)), intent(inout) :: this

    PUSH_SUB(SPECIAL(iterator,end))

    ASSERT(.not.associated(this%self))
    
    POP_SUB(SPECIAL(iterator,end))
  end subroutine SPECIAL(iterator,end)

#endif
  
  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_init_type)(this, that)
    type(TEMPLATE(iterator_t)),   intent(out) :: this
    type(BASE_TYPE_NAME), target, intent(in)  :: that

    PUSH_SUB(TEMPLATE(iterator_init_type))

    call SPECIAL(iterator,init)(this, that)
    this%self => that
    call TEMPLATE(dict_init)(this%iter, this%self%dict)

    POP_SUB(TEMPLATE(iterator_init_type))
  end subroutine TEMPLATE(iterator_init_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_init_copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(TEMPLATE(iterator_t)), intent(in)  :: that

    PUSH_SUB(TEMPLATE(iterator_init_copy))

    ASSERT(associated(that%self))
    call TEMPLATE(iterator_copy)(this, that)

    POP_SUB(TEMPLATE(iterator_init_copy))
  end subroutine TEMPLATE(iterator_init_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name)(this, name, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_name))

    call TEMPLATE(dict_next)(this%iter, name, ierr)

    POP_SUB(TEMPLATE(iterator_next_name))
  end subroutine TEMPLATE(iterator_next_name)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_type)(this, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(BASE_TYPE_NAME),      pointer        :: that
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_type))

    nullify(that)
    call TEMPLATE(dict_next)(this%iter, that, ierr)

    POP_SUB(TEMPLATE(iterator_next_type))
  end subroutine TEMPLATE(iterator_next_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name_type)(this, name, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    type(TEMPLATE(t)),         pointer        :: that
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_name_type))

    nullify(that)
    call TEMPLATE(dict_next)(this%iter, name, that, ierr)

    POP_SUB(TEMPLATE(iterator_next_name_type))
  end subroutine TEMPLATE(iterator_next_name_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(TEMPLATE(iterator_t)), intent(in)    :: that

    PUSH_SUB(TEMPLATE(iterator_copy))

    call TEMPLATE(iterator_end)(this)
    this%self => that%self
    call SPECIAL(iterator,copy)(this, that)
    call TEMPLATE(dict_copy)(this%iter, that%iter)

    POP_SUB(TEMPLATE(iterator_copy))
  end subroutine TEMPLATE(iterator_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_end)(this)
    type(TEMPLATE(iterator_t)), intent(inout) :: this

    PUSH_SUB(TEMPLATE(iterator_end))

    nullify(this%self)
    call SPECIAL(iterator,end)(this)
    call TEMPLATE(dict_end)(this%iter)

    POP_SUB(TEMPLATE(iterator_end))
  end subroutine TEMPLATE(iterator_end)

#endif
#if defined(BASE_INCLUDE_MODULE)

end module TEMPLATE(base_oct_m)

#endif

#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
