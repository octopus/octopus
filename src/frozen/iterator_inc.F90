#include "template.h"

#if defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)
  use config_dict_m
#endif

#if !defined(INCLUDE_PREFIX) && defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

  public ::             &
    TEMPLATE(NAME_LEN)

  public ::               &
    TEMPLATE(iterator_t)

  public ::         &
    TEMPLATE(next)

  integer, parameter :: TEMPLATE(NAME_LEN) = CONFIG_DICT_NAME_LEN

#if !defined(EXTENDED_TYPE)

  type :: TEMPLATE(iterator_t)
    private
    type(TEMPLATE(t)),   pointer :: self =>null()
    type(config_dict_iterator_t) :: iter
  end type TEMPLATE(iterator_t)

#endif

  interface TEMPLATE(init)
    module procedure TEMPLATE(iterator_init_type)
    module procedure TEMPLATE(iterator_init_copy)
  end interface TEMPLATE(init)

  interface TEMPLATE(next)
    module procedure TEMPLATE(iterator_next_name_config)
    module procedure TEMPLATE(iterator_next_name)
    module procedure TEMPLATE(iterator_next_config)
    module procedure TEMPLATE(iterator_next_name_config_that)
    module procedure TEMPLATE(iterator_next_name_that)
    module procedure TEMPLATE(iterator_next_config_that)
    module procedure TEMPLATE(iterator_next_that)
  end interface TEMPLATE(next)

  interface TEMPLATE(copy)
    module procedure TEMPLATE(iterator_copy)
  end interface TEMPLATE(copy)

  interface TEMPLATE(end)
    module procedure TEMPLATE(iterator_end)
  end interface TEMPLATE(end)

#endif
#if !defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && defined(INCLUDE_BODY)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_init_type)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(TEMPLATE(t)),  target, intent(in)  :: that

    PUSH_SUB(TEMPLATE(iterator_init_type))

#if defined(EXTENDED_TYPE)
    call TEMPLATE(iterator__init__)(this, that)
#endif
    this%self => that
    call config_dict_init(this%iter, this%self%dict)

    POP_SUB(TEMPLATE(iterator_init_type))
  end subroutine TEMPLATE(iterator_init_type)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_init_copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(out) :: this
    type(TEMPLATE(iterator_t)), intent(in)  :: that

    PUSH_SUB(TEMPLATE(iterator_init_copy))

    ASSERT(associated(that%self))
    call TEMPLATE(iterator_copy)(this, that)

    POP_SUB(TEMPLATE(iterator_init_type))
  end subroutine TEMPLATE(iterator_init_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name_config)(this, name, config, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    type(json_object_t),       pointer        :: config
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_name_config))

    call config_dict_next(this%iter, name, config, ierr)

    POP_SUB(TEMPLATE(iterator_next_name_config))
  end subroutine TEMPLATE(iterator_next_name_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name)(this, name, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_name))

    call config_dict_next(this%iter, name, ierr)

    POP_SUB(TEMPLATE(iterator_next_name))
  end subroutine TEMPLATE(iterator_next_name)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_config)(this, config, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(json_object_t),       pointer        :: config
    integer,          optional, intent(out)   :: ierr

    PUSH_SUB(TEMPLATE(iterator_next_config))

    call config_dict_next(this%iter, config, ierr)

    POP_SUB(TEMPLATE(iterator_next_config))
  end subroutine TEMPLATE(iterator_next_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name_config_that)(this, name, config, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    type(json_object_t),       pointer        :: config
    type(TEMPLATE(t)),         pointer        :: that
    integer,          optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_name_config_that))

    nullify(config, that)
    call config_dict_next(this%iter, name, config, jerr)
    if(jerr==CONFIG_DICT_OK) call TEMPLATE(gets)(this%self, config, that)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(iterator_next_name_config_that))
  end subroutine TEMPLATE(iterator_next_name_config_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_config_that)(this, config, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(json_object_t),       pointer        :: config
    type(TEMPLATE(t)),         pointer        :: that
    integer,          optional, intent(out)   :: ierr

    integer :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_name_config_that))

    nullify(config, that)
    call config_dict_next(this%iter, config, jerr)
    if(jerr==CONFIG_DICT_OK) call TEMPLATE(gets)(this%self, config, that)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(iterator_next_config_that))
  end subroutine TEMPLATE(iterator_next_config_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_name_that)(this, name, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    character(len=*),           intent(out)   :: name
    type(TEMPLATE(t)),         pointer        :: that
    integer,          optional, intent(out)   :: ierr

    type(json_object_t), pointer :: config
    integer                      :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_name_that))

    nullify(config, that)
    call config_dict_next(this%iter, name, config, jerr)
    if(jerr==CONFIG_DICT_OK) call TEMPLATE(gets)(this%self, config, that)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(iterator_next_name_that))
  end subroutine TEMPLATE(iterator_next_name_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_next_that)(this, that, ierr)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(TEMPLATE(t)),         pointer        :: that
    integer,          optional, intent(out)   :: ierr

    type(json_object_t), pointer :: config
    integer                      :: jerr

    PUSH_SUB(TEMPLATE(iterator_next_that))

    nullify(config, that)
    call config_dict_next(this%iter, config, jerr)
    if(jerr==CONFIG_DICT_OK) call TEMPLATE(gets)(this%self, config, that)
    if(present(ierr)) ierr = jerr

    POP_SUB(TEMPLATE(iterator_next_that))
  end subroutine TEMPLATE(iterator_next_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_copy)(this, that)
    type(TEMPLATE(iterator_t)), intent(inout) :: this
    type(TEMPLATE(iterator_t)), intent(in)    :: that

    PUSH_SUB(TEMPLATE(iterator_copy))

    this%self => that%self
    call config_dict_copy(this%iter, that%iter)
#if defined(EXTENDED_TYPE)
    call TEMPLATE(iterator__copy__)(this, that)
#endif

    POP_SUB(TEMPLATE(iterator_copy))
  end subroutine TEMPLATE(iterator_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(iterator_end)(this)
    type(TEMPLATE(iterator_t)), intent(inout) :: this

    PUSH_SUB(TEMPLATE(iterator_end))

    nullify(this%self)
    call config_dict_end(this%iter)
#if defined(EXTENDED_TYPE)
    call TEMPLATE(iterator__end__)(this)
#endif

    POP_SUB(TEMPLATE(iterator_end))
  end subroutine TEMPLATE(iterator_end)

#endif

!! Local Variables:
!! mode: f90
!! End:
