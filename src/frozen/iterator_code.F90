
#include "template.h"

#if defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

#if defined(EXCLUDE_TYPE)

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_next

#else

  use config_dict_m, only:  &
    config_dict_iterator_t

  use config_dict_m, only:  &
    config_dict_init,       &
    config_dict_next,       &
    config_dict_copy,       &
    config_dict_end

#endif

  use config_dict_m, only: &
    config_dict_t

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

#if defined(INTERNAL_ACCESS)

  use TEMPLATE(root,m), only: &
    TEMPLATE(root,_aget__)

#endif

#endif

#if !defined(INCLUDE_PREFIX) && defined(INCLUDE_HEADER) && !defined(INCLUDE_BODY)

#if !defined(EXCLUDE_TYPE)

  type, public :: TEMPLATE(base,iterator_t)
    private
    type(TEMPLATE(base,t)), pointer :: self =>null()
    type(config_dict_iterator_t)    :: iter
  end type TEMPLATE(base,iterator_t)

  interface TEMPLATE(base,init)
    module procedure TEMPLATE(base,iterator_init)
  end interface TEMPLATE(base,init)

#endif 

  interface TEMPLATE(base,next)
    module procedure TEMPLATE(base,iterator_next_name_config)
    module procedure TEMPLATE(base,iterator_next_name)
    module procedure TEMPLATE(base,iterator_next_config)
    module procedure TEMPLATE(base,iterator_next_name_config_that)
    module procedure TEMPLATE(base,iterator_next_name_that)
    module procedure TEMPLATE(base,iterator_next_config_that)
    module procedure TEMPLATE(base,iterator_next_that)
  end interface TEMPLATE(base,next)

#if !defined(EXCLUDE_TYPE)

  interface TEMPLATE(base,copy)
    module procedure TEMPLATE(base,iterator_copy)
  end interface TEMPLATE(base,copy)

  interface TEMPLATE(base,end)
    module procedure TEMPLATE(base,iterator_end)
  end interface TEMPLATE(base,end)

#endif

  integer, public, parameter :: TEMPLATE(base,NAME_LEN) = CONFIG_DICT_NAME_LEN

#endif
#if !defined(INCLUDE_PREFIX) && !defined(INCLUDE_HEADER) && defined(INCLUDE_BODY)

#if !defined(EXCLUDE_TYPE)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_init)(this, that)
    type(TEMPLATE(base,iterator_t)), intent(out) :: this
    type(TEMPLATE(base,t)),  target, intent(in)  :: that
    !
    type(config_dict_t), pointer :: dict
    !
    PUSH_SUB(TEMPLATE(base,iterator_init))
    this%self=>that
    call TEMPLATE(base,iterator__iget__)(this, dict)
    ASSERT(associated(dict))
    call config_dict_init(this%iter, dict)
    nullify(dict)
    POP_SUB(TEMPLATE(base,iterator_init))
    return
  end subroutine TEMPLATE(base,iterator_init)

#endif

#if defined(INTERNAL_ACCESS)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator__iget__)(this, that)
    type(TEMPLATE(base,iterator_t)), intent(in) :: this
    type(config_dict_t),            pointer     :: that
    !
    PUSH_SUB(TEMPLATE(base,iterator__iget__))
    call TEMPLATE(root,_aget__)(this%self, that)
    POP_SUB(TEMPLATE(base,iterator__iget__))
    return
  end subroutine TEMPLATE(base,iterator__iget__)

#else

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator__iget__)(this, that)
    type(TEMPLATE(base,iterator_t)), target, intent(in) :: this
    type(config_dict_t),            pointer             :: that
    !
    PUSH_SUB(TEMPLATE(base,iterator__iget__))
    that=>this%self%dict
    POP_SUB(TEMPLATE(base,iterator__iget__))
    return
  end subroutine TEMPLATE(base,iterator__iget__)

#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_name_config)(this, name, config, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: name
    type(json_object_t),            pointer        :: config
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_name_config))
    call config_dict_next(this%iter, name, config, ierr)
    POP_SUB(TEMPLATE(base,iterator_next_name_config))
    return
  end subroutine TEMPLATE(base,iterator_next_name_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_name)(this, name, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: name
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_name))
    call config_dict_next(this%iter, name, ierr)
    POP_SUB(TEMPLATE(base,iterator_next_name))
    return
  end subroutine TEMPLATE(base,iterator_next_name)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_config)(this, config, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    type(json_object_t),            pointer        :: config
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_config))
    call config_dict_next(this%iter, config, ierr)
    POP_SUB(TEMPLATE(base,iterator_next_config))
    return
  end subroutine TEMPLATE(base,iterator_next_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_name_config_that)(this, name, config, that, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: name
    type(json_object_t),            pointer        :: config
    type(TEMPLATE(base,t)),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    integer :: jerr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_name_config_that))
    nullify(config, that)
    call config_dict_next(this%iter, name, config, jerr)
    if(jerr==CONFIG_DICT_OK)&
      call TEMPLATE(base,get)(this%self, config, that)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(base,iterator_next_name_config_that))
    return
  end subroutine TEMPLATE(base,iterator_next_name_config_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_config_that)(this, config, that, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    type(json_object_t),            pointer        :: config
    type(TEMPLATE(base,t)),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    integer :: jerr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_name_config_that))
    nullify(config, that)
    call config_dict_next(this%iter, config, jerr)
    if(jerr==CONFIG_DICT_OK)&
      call TEMPLATE(base,get)(this%self, config, that)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(base,iterator_next_config_that))
    return
  end subroutine TEMPLATE(base,iterator_next_config_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_name_that)(this, name, that, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    character(len=*),                intent(out)   :: name
    type(TEMPLATE(base,t)),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(json_object_t), pointer :: config
    integer                      :: jerr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_name_that))
    nullify(config, that)
    call config_dict_next(this%iter, name, config, jerr)
    if(jerr==CONFIG_DICT_OK)&
      call TEMPLATE(base,get)(this%self, config, that)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(base,iterator_next_name_that))
    return
  end subroutine TEMPLATE(base,iterator_next_name_that)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_next_that)(this, that, ierr)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    type(TEMPLATE(base,t)),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    type(json_object_t), pointer :: config
    integer                      :: jerr
    !
    PUSH_SUB(TEMPLATE(base,iterator_next_that))
    nullify(config, that)
    call config_dict_next(this%iter, config, jerr)
    if(jerr==CONFIG_DICT_OK)&
      call TEMPLATE(base,get)(this%self, config, that)
    if(present(ierr))ierr=jerr
    POP_SUB(TEMPLATE(base,iterator_next_that))
    return
  end subroutine TEMPLATE(base,iterator_next_that)

#if !defined(EXCLUDE_TYPE)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_copy)(this, that)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    type(TEMPLATE(base,iterator_t)), intent(in)    :: that
    !
    PUSH_SUB(TEMPLATE(base,iterator_copy))
    this%self=>that%self
    call config_dict_copy(this%iter, that%iter)
    POP_SUB(TEMPLATE(base,iterator_copy))
    return
  end subroutine TEMPLATE(base,iterator_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base,iterator_end)(this)
    type(TEMPLATE(base,iterator_t)), intent(inout) :: this
    !
    PUSH_SUB(TEMPLATE(base,iterator_end))
    nullify(this%self)
    call config_dict_end(this%iter)
    POP_SUB(TEMPLATE(base,iterator_end))
    return
  end subroutine TEMPLATE(base,iterator_end)

#endif

#endif

!! Local Variables:
!! mode: f90
!! End:
