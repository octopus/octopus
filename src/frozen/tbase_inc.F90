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

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME BASE_TEMPLATE_NAME
#define DICT_TYPE_NAME DECORATE(BASE_TEMPLATE_NAME,husk_t)

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

#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

  use refcount_oct_m

#endif
#if defined(BASE_INCLUDE_MODULE)

  implicit none

  private

#endif
#if defined(BASE_INCLUDE_HEADER) || defined(BASE_INCLUDE_MODULE)

#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER

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

  public ::            &
    SPECIAL(register)

  public ::          &
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

  integer, parameter :: TEMPLATE(HUSK_DISA) = 0
  integer, parameter :: TEMPLATE(HUSK_NULL) = 1
  integer, parameter :: TEMPLATE(HUSK_ASSC) = 2

  type :: TEMPLATE(husk_t)
    private
    type(BASE_TYPE_NAME), pointer :: self =>null()
    type(refcount_t),     pointer :: rcnt =>null()
    type(json_object_t),  pointer :: cnfg =>null()
    logical                       :: lock = .false.
    logical                       :: actv = .true.
    integer                       :: type = TEMPLATE(HUSK_DISA)
  end type TEMPLATE(husk_t)

#if !defined(BASE_EXTENDED_ITERATOR_TYPE)

  type :: TEMPLATE(iterator_t)
    private
    type(BASE_TYPE_NAME),   pointer :: self =>null()
    type(TEMPLATE(dict_iterator_t)) :: iter
  end type TEMPLATE(iterator_t)

#endif

  interface TEMPLATE(husk_assoc)
    module procedure TEMPLATE(husk_assoc_type)
    module procedure TEMPLATE(husk_assoc_self)
    module procedure TEMPLATE(husk_assoc_husk)
  end interface TEMPLATE(husk_assoc)

  interface TEMPLATE(husk_set)
    module procedure TEMPLATE(husk_set_info)
    module procedure TEMPLATE(husk_set_type)
  end interface TEMPLATE(husk_set)

  interface TEMPLATE(husk_get)
    module procedure TEMPLATE(husk_get_info)
    module procedure TEMPLATE(husk_get_type)
  end interface TEMPLATE(husk_get)

  interface TEMPLATE(new)
    module procedure TEMPLATE(new_copy)
  end interface TEMPLATE(new)

  interface TEMPLATE(del)
    module procedure TEMPLATE(del_type)
    module procedure TEMPLATE(del_pass)
  end interface TEMPLATE(del)

  interface TEMPLATE(init)
    module procedure TEMPLATE(init_copy)
    module procedure TEMPLATE(iterator_init_type)
    module procedure TEMPLATE(iterator_init_copy)
  end interface TEMPLATE(init)

  interface TEMPLATE(sets)
    module procedure TEMPLATE(sets_info)
    module procedure TEMPLATE(sets_type)
  end interface TEMPLATE(sets)

  interface TEMPLATE(get)
    module procedure TEMPLATE(get_sub_config)
  end interface TEMPLATE(get)

  interface TEMPLATE(next)
    module procedure TEMPLATE(iterator_next_name)
    module procedure TEMPLATE(iterator_next_type)
    module procedure TEMPLATE(iterator_next_name_type)
  end interface TEMPLATE(next)

  interface TEMPLATE(copy)
    module procedure TEMPLATE(copy_type)
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

#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY

#undef TEMPLATE_NAME
#define TEMPLATE_NAME BASE_TEMPLATE_NAME
#include "template.h"

  ! ---------------------------------------------------------
  function TEMPLATE(husk_new)(that, config, lock, active) result(this)
    type(BASE_TYPE_NAME), intent(in) :: that
    type(json_object_t),  intent(in) :: config
    logical,    optional, intent(in) :: lock
    logical,    optional, intent(in) :: active

    type(TEMPLATE(husk_t)), pointer :: this

    PUSH_SUB(TEMPLATE(husk_new))

    nullify(this)
    SAFE_ALLOCATE(this)
    call TEMPLATE(husk_init)(this, that, config, lock, active)

    POP_SUB(TEMPLATE(husk_new))
  end function TEMPLATE(husk_new)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_del)(this)
    type(TEMPLATE(husk_t)), pointer, intent(inout) :: this

    PUSH_SUB(TEMPLATE(husk_del))

    if(associated(this))then
      call TEMPLATE(husk_end)(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(TEMPLATE(husk_del))
  end subroutine TEMPLATE(husk_del)

  ! ---------------------------------------------------------
  function TEMPLATE(husk_assoc_type)(this) result(that)
    type(TEMPLATE(husk_t)), intent(in) :: this

    logical :: that

    PUSH_SUB(TEMPLATE(husk_assoc_type))

    select case(this%type)
    case(TEMPLATE(HUSK_DISA))
      that = .false.
    case(TEMPLATE(HUSK_NULL))
      ASSERT(.not.associated(this%self))
      ASSERT(.not.associated(this%rcnt))
      ASSERT(.not.associated(this%cnfg))
      that = .false.
    case(TEMPLATE(HUSK_ASSC))
      ASSERT(associated(this%self))
      ASSERT(associated(this%rcnt))
      ASSERT(associated(this%cnfg))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(TEMPLATE(husk_assoc_type))
  end function TEMPLATE(husk_assoc_type)

  ! ---------------------------------------------------------
  function TEMPLATE(husk_assoc_self)(this, that) result(assc)
    type(TEMPLATE(husk_t)),       intent(in) :: this
    type(BASE_TYPE_NAME), target, intent(in) :: that

    logical :: assc

    PUSH_SUB(TEMPLATE(husk_assoc_self))

    assc = .false.
    if(TEMPLATE(husk_assoc)(this))then
      ASSERT(this%type==TEMPLATE(HUSK_ASSC))
      assc = associated(this%self, that)
    end if

    POP_SUB(TEMPLATE(husk_assoc_self))
  end function TEMPLATE(husk_assoc_self)

  ! ---------------------------------------------------------
  function TEMPLATE(husk_assoc_husk)(this, that) result(assc)
    type(TEMPLATE(husk_t)), intent(in) :: this
    type(TEMPLATE(husk_t)), intent(in) :: that

    logical :: assc

    PUSH_SUB(TEMPLATE(husk_assoc_husk))

    assc = .false.
    if(TEMPLATE(husk_assoc)(this))then
      ASSERT(this%type==TEMPLATE(HUSK_ASSC))
      assc = associated(this%self, that%self)
    end if

    POP_SUB(TEMPLATE(husk_assoc_husk))
  end function TEMPLATE(husk_assoc_husk)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_init)(this, that, config, lock, active)
    type(TEMPLATE(husk_t)), intent(out) :: this
    type(BASE_TYPE_NAME),   intent(in)  :: that
    type(json_object_t),    intent(in)  :: config
    logical,      optional, intent(in)  :: lock
    logical,      optional, intent(in)  :: active

    PUSH_SUB(TEMPLATE(husk_init))

    nullify(this%self, this%rcnt, this%cnfg)
    this%lock = .false.
    this%actv = .true.
    this%type = TEMPLATE(HUSK_NULL)
    call TEMPLATE(husk_set)(this, type=that, config=config, lock=lock, active=active)

    POP_SUB(TEMPLATE(husk_init))
  end subroutine TEMPLATE(husk_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_set_info)(this, lock, active)
    type(TEMPLATE(husk_t)), intent(inout) :: this
    logical,      optional, intent(in)    :: lock
    logical,      optional, intent(in)    :: active

    PUSH_SUB(TEMPLATE(husk_set_info))

    ASSERT(TEMPLATE(husk_assoc)(this))
    if(present(lock))   this%lock = lock
    if(present(active)) this%actv = active

    POP_SUB(TEMPLATE(husk_set_info))
  end subroutine TEMPLATE(husk_set_info)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_set_type)(this, type, config, lock, active)
    type(TEMPLATE(husk_t)),       intent(inout) :: this
    type(BASE_TYPE_NAME), target, intent(in)    :: type
    type(json_object_t),  target, intent(in)    :: config
    logical,            optional, intent(in)    :: lock
    logical,            optional, intent(in)    :: active

    PUSH_SUB(TEMPLATE(husk_set_type))

    ASSERT(this%type==TEMPLATE(HUSK_NULL))
    ASSERT(.not.TEMPLATE(husk_assoc)(this))
    this%self => type
    this%cnfg => config
    call SPECIAL(register)(type, this%rcnt)
    ASSERT(associated(this%rcnt))
    call refcount_inc(this%rcnt)
    this%type = TEMPLATE(HUSK_ASSC)
    call TEMPLATE(husk_set)(this, lock=lock, active=active)

    POP_SUB(TEMPLATE(husk_set_type))
  end subroutine TEMPLATE(husk_set_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_get_info)(this, config, lock, active)
    type(TEMPLATE(husk_t)),                  intent(in)  :: this
    type(json_object_t),  pointer, optional, intent(out) :: config
    logical,                       optional, intent(out) :: lock
    logical,                       optional, intent(out) :: active

    PUSH_SUB(TEMPLATE(husk_get_info))

    ASSERT(TEMPLATE(husk_assoc)(this))
    if(present(config)) config =>this%cnfg
    if(present(lock))     lock = this%lock
    if(present(active)) active = this%actv

    POP_SUB(TEMPLATE(husk_get_info))
  end subroutine TEMPLATE(husk_get_info)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_get_type)(this, type, config, lock, active)
    type(TEMPLATE(husk_t)),                  intent(in)  :: this
    type(BASE_TYPE_NAME), pointer,           intent(out) :: type
    type(json_object_t),  pointer, optional, intent(out) :: config
    logical,                       optional, intent(out) :: lock
    logical,                       optional, intent(out) :: active

    PUSH_SUB(TEMPLATE(husk_get_type))

    ASSERT(this%type>TEMPLATE(HUSK_DISA))
    nullify(type)
    if(TEMPLATE(husk_assoc)(this)) type => this%self
    call TEMPLATE(husk_get)(this, config=config, lock=lock, active=active)

    POP_SUB(TEMPLATE(husk_get_type))
  end subroutine TEMPLATE(husk_get_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_copy)(this, that)
    type(TEMPLATE(husk_t)), intent(inout) :: this
    type(TEMPLATE(husk_t)), intent(in)    :: that

    PUSH_SUB(TEMPLATE(husk_copy))

    call TEMPLATE(husk_end)(this)
    this%type = that%type
    if(TEMPLATE(husk_assoc)(that))&
      call TEMPLATE(husk_set)(this, that%self, that%cnfg, lock=that%lock, active=that%actv)

    POP_SUB(TEMPLATE(husk_copy))
  end subroutine TEMPLATE(husk_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(husk_end)(this)
    type(TEMPLATE(husk_t)), intent(inout) :: this

    PUSH_SUB(TEMPLATE(husk_end))

    if(TEMPLATE(husk_assoc)(this)) call refcount_dec(this%rcnt)
    nullify(this%self, this%rcnt, this%cnfg)
    this%lock = .false.
    this%actv = .true.
    this%type = TEMPLATE(HUSK_DISA)

    POP_SUB(TEMPLATE(husk_end))
  end subroutine TEMPLATE(husk_end)

  ! ---------------------------------------------------------
  recursive function TEMPLATE(new_copy)(source, mold) result(this)
    type(BASE_TYPE_NAME), optional, intent(in) :: source
    type(BASE_TYPE_NAME), optional, intent(in) :: mold

    type(BASE_TYPE_NAME), pointer :: this

    PUSH_SUB(TEMPLATE(new_copy))

    ASSERT(.not.(present(source).eqv.present(mold)))
    nullify(this)
    SAFE_ALLOCATE(this)
    if(present(source))then
      call TEMPLATE(copy)(this, source)
    elseif(present(mold))then
      call TEMPLATE(init)(this, mold)
    else
      ASSERT(.FALSE.)
    end if
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(TEMPLATE(new_copy))
  end function TEMPLATE(new_copy)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(del_type)(this)
    type(BASE_TYPE_NAME), pointer, intent(inout) :: this

    PUSH_SUB(TEMPLATE(del_type))

    call TEMPLATE(del)(this, TEMPLATE(end_type))
    
    POP_SUB(TEMPLATE(del_type))
  end subroutine TEMPLATE(del_type)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(del_pass)(this, finis)
    type(BASE_TYPE_NAME), pointer, intent(inout) :: this

    interface
      subroutine finis(this)
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME), intent(inout) :: this
      end subroutine finis
    end interface
    
    logical :: free

    PUSH_SUB(TEMPLATE(del_pass))

    if(associated(this))then
      ASSERT(associated(this%rcnt))
      call refcount_get(this%rcnt, free=free)
      if(free)then
        call finis(this)
        SAFE_DEALLOCATE_P(this)
        nullify(this)
      end if
    end if

    POP_SUB(TEMPLATE(del_pass))
  end subroutine TEMPLATE(del_pass)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(init_copy)(this, that)
    type(BASE_TYPE_NAME), intent(out) :: this
    type(BASE_TYPE_NAME), intent(in)  :: that

    type(TEMPLATE(dict_iterator_t))   :: iter
    character(len=TEMPLATE(NAME_LEN)) :: name
    type(TEMPLATE(husk_t)),   pointer :: husk
    type(BASE_TYPE_NAME),     pointer :: subs
    type(json_object_t),      pointer :: cnfg
    logical                           :: lock, actv
    integer                           :: ierr

    PUSH_SUB(TEMPLATE(init_copy))

    ASSERT(associated(that%config))
    call SPECIAL(init)(this, that)
    call TEMPLATE(dict_init)(iter, that%dict)
    do
      nullify(husk, subs, cnfg)
      call TEMPLATE(dict_next)(iter, name, husk, ierr)
      if(ierr/=TEMPLATE(OK))exit
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=subs, config=cnfg, lock=lock, active=actv)
      call TEMPLATE(sets)(this, trim(adjustl(name)), TEMPLATE(new)(mold=subs), config=cnfg, lock=lock, active=actv)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(husk, subs, cnfg)

    POP_SUB(TEMPLATE(init_copy))
  end subroutine TEMPLATE(init_copy)

  ! ---------------------------------------------------------
  subroutine SPECIAL(register)(this, that)
    type(BASE_TYPE_NAME), intent(in) :: this
    type(refcount_t),    pointer     :: that

    PUSH_SUB(SPECIAL(register))

    nullify(that)
    if(associated(this%rcnt)) that => this%rcnt

    POP_SUB(SPECIAL(register))
  end subroutine SPECIAL(register)

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
    type(TEMPLATE(husk_t)), pointer :: husk
    type(BASE_TYPE_NAME),   pointer :: subs
    logical                         :: prnt, lock

    PUSH_SUB(SPECIAL(apply))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_init)(iter, this%dict)
    do
      nullify(husk, subs)
      call TEMPLATE(dict_next)(iter, husk)
      if(.not.associated(husk))exit
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=subs, lock=lock)
      if(.not.lock) call operation(subs)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(husk, subs)
    prnt = .true.
    if(present(parent)) prnt = parent
    if(prnt) call operation(this)

    POP_SUB(SPECIAL(apply))
  end subroutine SPECIAL(apply)

  ! ---------------------------------------------------------
  subroutine SPECIAL(reduce)(this, operation)
    type(BASE_TYPE_NAME), intent(inout) :: this

    interface
      subroutine operation(this, that, config)
        use json_oct_m
        import :: BASE_TYPE_NAME
        type(BASE_TYPE_NAME),          intent(inout) :: this
        type(BASE_TYPE_NAME),          intent(in)    :: that
        type(json_object_t), optional, intent(in)    :: config
      end subroutine operation
    end interface

    type(TEMPLATE(dict_iterator_t)) :: iter
    type(TEMPLATE(husk_t)), pointer :: husk
    type(BASE_TYPE_NAME),   pointer :: subs
    type(json_object_t),    pointer :: cnfg
    logical                         :: actv

    PUSH_SUB(SPECIAL(reduce))

    ASSERT(associated(this%config))
    call TEMPLATE(dict_init)(iter, this%dict)
    do
      nullify(husk, subs)
      call TEMPLATE(dict_next)(iter, husk)
      if(.not.associated(husk))exit
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=subs, active=actv, config=cnfg)
      if(actv) call operation(this, subs, config=cnfg)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(husk, subs)

    POP_SUB(SPECIAL(reduce))
  end subroutine SPECIAL(reduce)

#if defined(BASE_LEAF_TYPE)

  ! ---------------------------------------------------------
  subroutine SPECIAL(sets)(this, name, that, config, lock, active)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(BASE_TYPE_NAME), optional, intent(in)    :: that
    type(json_object_t),  optional, intent(in)    :: config
    logical,    optional, intent(in)    :: lock
    logical,    optional, intent(in)    :: active

    PUSH_SUB(SPECIAL(sets))

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    if(present(that))   continue
    if(present(config)) continue
    if(present(lock))   continue
    if(present(active)) continue
    
    POP_SUB(SPECIAL(sets))
  end subroutine SPECIAL(sets)
    
  ! ---------------------------------------------------------
  subroutine SPECIAL(dels)(this, name, that)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(BASE_TYPE_NAME), intent(in)    :: that

    PUSH_SUB(SPECIAL(dels))

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))

    POP_SUB(SPECIAL(dels))
  end subroutine SPECIAL(dels)

#endif

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(sets_info)(this, name, lock, active)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    logical,    optional, intent(in)    :: lock
    logical,    optional, intent(in)    :: active

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: ierr

    PUSH_SUB(TEMPLATE(sets_info))

    ASSERT(associated(this%config))
    nullify(husk)
    call TEMPLATE(dict_get)(this%dict, trim(adjustl(name)), husk, ierr)
    ASSERT(ierr==TEMPLATE(OK))
    ASSERT(associated(husk))
    ASSERT(TEMPLATE(husk_assoc)(husk))
    call TEMPLATE(husk_set)(husk, lock=lock, active=active)
    call SPECIAL(sets)(this, trim(adjustl(name)), lock=lock, active=active)
    nullify(husk)

    POP_SUB(TEMPLATE(sets_info))
  end subroutine TEMPLATE(sets_info)
    
  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(sets_type)(this, name, that, config, lock, active)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(BASE_TYPE_NAME), intent(in)    :: that
    type(json_object_t),  intent(in)    :: config
    logical,    optional, intent(in)    :: lock
    logical,    optional, intent(in)    :: active

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: ierr

    PUSH_SUB(TEMPLATE(sets_type))

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    nullify(husk)
    call TEMPLATE(dict_del)(this%dict, trim(adjustl(name)), husk, ierr)
    if(ierr==TEMPLATE(OK))then
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_del)(husk)
      nullify(husk)
    end if
    call TEMPLATE(dict_set)(this%dict, trim(adjustl(name)), TEMPLATE(husk_new)(that, config, lock=lock, active=active))
    call SPECIAL(sets)(this, trim(adjustl(name)), that, config, lock=lock, active=active)

    POP_SUB(TEMPLATE(sets_type))
  end subroutine TEMPLATE(sets_type)
    
  ! ---------------------------------------------------------
  recursive subroutine INTERNAL(gets)(this, name, that, ierr)
    type(BASE_TYPE_NAME),            intent(in)  :: this
    character(len=*),                intent(in)  :: name
    type(TEMPLATE(husk_t)), pointer, intent(out) :: that
    integer,                         intent(out) :: ierr

    type(TEMPLATE(husk_t)), pointer :: husk
    type(BASE_TYPE_NAME),   pointer :: subs
    integer                         :: ipos
    
    PUSH_SUB(INTERNAL(gets))

    ASSERT(associated(this%config))
    nullify(that, husk, subs)
    ipos = index(name, "/")
    if(ipos>0)then
      call TEMPLATE(dict_get)(this%dict, trim(adjustl(name(1:ipos-1))), husk, ierr)
      if(ierr==TEMPLATE(OK))then
        ASSERT(associated(husk))
        ASSERT(TEMPLATE(husk_assoc)(husk))
        call TEMPLATE(husk_get)(husk, type=subs)
        call INTERNAL(gets)(subs, trim(adjustl(name(ipos+1:))), that, ierr)
        nullify(subs)
      end if
      nullify(husk)
    else
      call TEMPLATE(dict_get)(this%dict, trim(adjustl(name)), that, ierr)
    end if

    POP_SUB(INTERNAL(gets))
  end subroutine INTERNAL(gets)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(gets)(this, name, type, config, lock, active)
    type(BASE_TYPE_NAME),                    intent(in)  :: this
    character(len=*),                        intent(in)  :: name
    type(BASE_TYPE_NAME), pointer, optional, intent(out) :: type
    type(json_object_t),  pointer, optional, intent(out) :: config
    logical,                       optional, intent(out) :: lock
    logical,                       optional, intent(out) :: active

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: ierr
    
    PUSH_SUB(TEMPLATE(gets_type))

    ASSERT(associated(this%config))
    nullify(type, husk)
    call INTERNAL(gets)(this, name, husk, ierr)
    if(ierr==TEMPLATE(OK))then
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      if(present(type))then
        call TEMPLATE(husk_get)(husk, type=type, config=config, lock=lock, active=active)
      else
        call TEMPLATE(husk_get)(husk, config=config, lock=lock, active=active)
      end if
    end if
    nullify(husk)

    POP_SUB(TEMPLATE(gets))
  end subroutine TEMPLATE(gets)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(dels)(this, name, that)
    type(BASE_TYPE_NAME), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(BASE_TYPE_NAME), intent(in)    :: that

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: ierr

    PUSH_SUB(TEMPLATE(dels))

    ASSERT(associated(this%config))
    nullify(husk)
    call TEMPLATE(dict_del)(this%dict, trim(adjustl(name)), husk, ierr)
    ASSERT(ierr==TEMPLATE(OK))
    ASSERT(TEMPLATE(husk_assoc)(husk, that))
    call TEMPLATE(husk_del)(husk)
    nullify(husk)
    call SPECIAL(dels)(this, trim(adjustl(name)), that)

    POP_SUB(TEMPLATE(dels))
  end subroutine TEMPLATE(dels)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(get_sub_config)(this, name, that)
    type(BASE_TYPE_NAME),         intent(in)  :: this
    character(len=*),             intent(in)  :: name
    type(json_object_t), pointer, intent(out) :: that
    
    type(BASE_TYPE_NAME), pointer :: subs

    PUSH_SUB(TEMPLATE(get_sub_config))

    nullify(that, subs)
    call TEMPLATE(gets)(this, trim(adjustl(name)), subs)
    if(associated(subs)) call TEMPLATE(get)(subs, that)
    nullify(subs)

    POP_SUB(TEMPLATE(get_sub_config))
  end subroutine TEMPLATE(get_sub_config)

  ! ---------------------------------------------------------
  recursive subroutine TEMPLATE(copy_type)(this, that)
    type(BASE_TYPE_NAME), intent(inout) :: this
    type(BASE_TYPE_NAME), intent(in)    :: that

    type(TEMPLATE(dict_iterator_t))   :: iter
    character(len=TEMPLATE(NAME_LEN)) :: name
    type(TEMPLATE(husk_t)),   pointer :: husk
    type(BASE_TYPE_NAME),     pointer :: subs
    type(json_object_t),      pointer :: cnfg
    type(refcount_t),         pointer :: rcnt
    logical                           :: lock, actv
    integer                           :: ierr

    PUSH_SUB(TEMPLATE(copy_type))

    rcnt => this%rcnt
    nullify(this%rcnt)
    call TEMPLATE(end)(this)
    call SPECIAL(copy)(this, that)
    call TEMPLATE(dict_init)(iter, that%dict)
    do
      nullify(husk, subs, cnfg)
      call TEMPLATE(dict_next)(iter, name, husk, ierr)
      if(ierr/=TEMPLATE(OK))exit
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=subs, config=cnfg, lock=lock, active=actv)
      call TEMPLATE(sets)(this, trim(adjustl(name)), TEMPLATE(new)(source=subs), config=cnfg, lock=lock, active=actv)
    end do
    call TEMPLATE(dict_end)(iter)
    nullify(husk, subs, cnfg)
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(TEMPLATE(copy_type))
  end subroutine TEMPLATE(copy_type)

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
    
    character(len=TEMPLATE(NAME_LEN)) :: name
    type(TEMPLATE(husk_t)),   pointer :: husk
    type(BASE_TYPE_NAME),     pointer :: subs
    integer                           :: ierr

    PUSH_SUB(TEMPLATE(end_pass))

    do
      nullify(husk, subs)
      call TEMPLATE(dict_pop)(this%dict, name, husk, ierr)
      if(ierr/=TEMPLATE(OK))exit
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=subs)
      call TEMPLATE(husk_del)(husk)
      call SPECIAL(dels)(this, trim(adjustl(name)), subs)
      call TEMPLATE(del)(subs, finis)
    end do
    nullify(husk, subs)
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
    call TEMPLATE(dict_init)(this%iter, that%dict)

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
    type(TEMPLATE(iterator_t)),    intent(inout) :: this
    type(BASE_TYPE_NAME), pointer, intent(out)   :: that
    integer,             optional, intent(out)   :: ierr

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_type))

    nullify(that, husk)
    call TEMPLATE(dict_next)(this%iter, husk, jerr)
    if(jerr==TEMPLATE(OK))then
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=that)
    end if
    nullify(husk)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(iterator_next_type))
  end subroutine TEMPLATE(iterator_next_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name_type)(this, name, that, ierr)
    type(TEMPLATE(iterator_t)),    intent(inout) :: this
    character(len=*),              intent(out)   :: name
    type(BASE_TYPE_NAME), pointer, intent(out)   :: that
    integer,             optional, intent(out)   :: ierr

    type(TEMPLATE(husk_t)), pointer :: husk
    integer                         :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_name_type))

    nullify(that, husk)
    call TEMPLATE(dict_next)(this%iter, name, husk, jerr)
    if(jerr==TEMPLATE(OK))then
      ASSERT(associated(husk))
      ASSERT(TEMPLATE(husk_assoc)(husk))
      call TEMPLATE(husk_get)(husk, type=that)
    end if
    nullify(husk)
    if(present(ierr)) ierr = jerr

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

#undef DICT_TYPE_NAME
#undef DICT_TEMPLATE_NAME

#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
