#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

module base_density_oct_m

  use dnst_oct_m
  use dnst_intrf_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_density
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

#define TEMPLATE_PREFIX base_density
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                   &
    BASE_DENSITY_OK,          &
    BASE_DENSITY_KEY_ERROR,   &
    BASE_DENSITY_EMPTY_ERROR

  public ::         &
    base_density_t

  public ::                 &
    base_density__init__,   &
    base_density__build__,  &
    base_density__start__,  &
    base_density__acc__,    &
    base_density__update__, &
    base_density__reset__,  &
    base_density__stop__,   &
    base_density__copy__,   &
    base_density__end__

  public ::              &
    base_density_new,    &
    base_density_del,    &
    base_density_init,   &
    base_density_start,  &
    base_density_acc,    &
    base_density_update, &
    base_density_reset,  &
    base_density_stop,   &
    base_density_set,    &
    base_density_get,    &
    base_density_copy,   &
    base_density_end

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_density
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

  integer, parameter :: BASE_DENSITY_OK          = BASE_DENSITY_DICT_OK
  integer, parameter :: BASE_DENSITY_KEY_ERROR   = BASE_DENSITY_DICT_KEY_ERROR
  integer, parameter :: BASE_DENSITY_EMPTY_ERROR = BASE_DENSITY_DICT_EMPTY_ERROR

  type :: base_density_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(base_density_t), pointer :: prnt   =>null()
    logical                       :: accu   = .false.
    type(dnst_intrf_t)            :: dnst
    type(base_density_dict_t)     :: dict
    type(base_density_list_t)     :: list
  end type base_density_t

  interface base_density__init__
    module procedure base_density__init__type
    module procedure base_density__init__copy
  end interface base_density__init__

  interface base_density_init
    module procedure base_density_init_type
    module procedure base_density_init_copy
  end interface base_density_init

  interface base_density_set
    module procedure base_density_set_info
    module procedure base_density_set_charge
    module procedure base_density_set_sub
  end interface base_density_set

  interface base_density_get
    module procedure base_density_get_info
    module procedure base_density_get_charge
    module procedure base_density_get_config
    module procedure base_density_get_simulation
    module procedure base_density_get_storage
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
    module procedure base_density_get_dnst_intrf
    module procedure base_density_get_sub
    module procedure base_density_get_sub_density_1d
    module procedure base_density_get_sub_density_2d
  end interface base_density_get

  interface base_density_copy
    module procedure base_density_copy_type
  end interface base_density_copy

  interface base_density_end
    module procedure base_density_end_type
  end interface base_density_end

#define TEMPLATE_PREFIX base_density
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_density
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_density
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_density__new__(this)
    type(base_density_t), pointer :: this

    PUSH_SUB(base_density__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_density__new__)
  end subroutine base_density__new__

  ! ---------------------------------------------------------
  subroutine base_density__del__(this)
    type(base_density_t), pointer :: this

    PUSH_SUB(base_density__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_density__del__)
  end subroutine base_density__del__

  ! ---------------------------------------------------------
  subroutine base_density_new(this, that)
    type(base_density_t),  target, intent(inout) :: this
    type(base_density_t), pointer                :: that

    PUSH_SUB(base_density_new)

    nullify(that)
    call base_density__new__(that)
    that%prnt => this
    call base_density_list_push(this%list, that)

    POP_SUB(base_density_new)
  end subroutine base_density_new

  ! ---------------------------------------------------------
  subroutine base_density_del(this)
    type(base_density_t), pointer :: this

    PUSH_SUB(base_density_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_density_list_del(this%prnt%list, this)
        call base_density_end(this)
        call base_density__del__(this)
      end if
    end if

    POP_SUB(base_density_del)
  end subroutine base_density_del

  ! ---------------------------------------------------------
  subroutine base_density__iinit__(this, config)
    type(base_density_t), target, intent(out) :: this
    type(json_object_t),  target, intent(in)  :: config

    integer :: ierr

    PUSH_SUB(base_density__iinit__)

    this%config => config
    call json_get(this%config, "reduce", this%accu, ierr)
    if(ierr/=JSON_OK) this%accu = .false.
    call dnst_intrf_init(this%dnst, config)
    call base_density_dict_init(this%dict)
    call base_density_list_init(this%list)

    POP_SUB(base_density__iinit__)
  end subroutine base_density__iinit__

  ! ---------------------------------------------------------
  subroutine base_density__init__type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    logical :: dflt
    integer :: ierr

    PUSH_SUB(base_density__init__type)

    call base_density__iinit__(this, config)
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    if(dflt) call dnst_intrf_new(this%dnst, init)

    POP_SUB(base_density__init__type)

  contains
    
    subroutine init(this, config)
      type(dnst_t),        intent(out) :: this
      type(json_object_t), intent(in)  :: config

      PUSH_SUB(base_density__init__type.init)

      call dnst_init(this, config)
      
      POP_SUB(base_density__init__type.init)
    end subroutine init

  end subroutine base_density__init__type

  ! ---------------------------------------------------------
  subroutine base_density__init__copy(this, that)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that

    type(dnst_t), pointer :: pdns
    
    PUSH_SUB(base_density__init__copy)

    ASSERT(associated(that%config))
    nullify(pdns)
    call base_density__iinit__(this, that%config)
    if(dnst_intrf_assoc(that%dnst))then
      if(dnst_intrf_alloc(that%dnst))then
        call dnst_intrf_new(this%dnst, init)
      else
        call dnst_intrf_get(that%dnst, pdns)
        ASSERT(associated(pdns))
        call dnst_intrf_set(this%dnst, pdns)
        nullify(pdns)
      end if
      if(associated(that%sim)) call base_density__start__(this, that%sim)
    end if

    POP_SUB(base_density__init__copy)
    
  contains
    
    subroutine init(this, config)
      type(dnst_t),        intent(out) :: this
      type(json_object_t), intent(in)  :: config

      PUSH_SUB(base_density__init__copy.init)

      call dnst_init(this, config)
      
      POP_SUB(base_density__init__copy.init)
    end subroutine init

  end subroutine base_density__init__copy

  ! ---------------------------------------------------------
  subroutine base_density__build__(this, build)
    type(base_density_t), intent(inout) :: this

    interface
      subroutine build(this, that, name)
        import :: base_density_t
        type(base_density_t), intent(inout) :: this
        type(base_density_t), intent(in)    :: that
        character(len=*),     intent(in)    :: name
      end subroutine build
    end interface

    type(base_density_iterator_t)        :: iter
    character(len=BASE_DENSITY_NAME_LEN) :: name
    type(base_density_t),        pointer :: subs
    integer                              :: ierr

    PUSH_SUB(base_density__build__)

    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, name, subs, ierr)
      if(ierr/=BASE_DENSITY_OK)exit
      call build(this, subs, name)
    end do
    call base_density_end(iter)
    nullify(subs)

    POP_SUB(base_density__build__)
  end subroutine base_density__build__

  ! ---------------------------------------------------------
  subroutine base_density_init_type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(base_density_init_type)

    call base_density__init__(this, config)

    POP_SUB(base_density_init_type)
  end subroutine base_density_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_density_init_copy(this, that)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that

    PUSH_SUB(base_density_init_copy)

    call base_density__init__(this, that)
    call base_density__build__(this, init)
    
    POP_SUB(base_density_init_copy)

  contains

    recursive subroutine init(this, isub, name)
      type(base_density_t), intent(inout) :: this
      type(base_density_t), intent(in)    :: isub
      character(len=*),     intent(in)    :: name

      type(base_density_t), pointer :: osub
      
      POP_SUB(base_density_init_copy.init)

      nullify(osub)
      if(base_density_list_index(that%list, isub)>0)then
        call base_density_new(this, osub)
        call base_density_init(osub, isub)
        call base_density_set(this, name, osub)
        nullify(osub)
      else
        call base_density_set(this, name, isub)
      end if

      PUSH_SUB(base_density_init_copy.init)
    end subroutine init

  end subroutine base_density_init_copy

  ! ---------------------------------------------------------
  subroutine base_density__start__(this, sim)
    type(base_density_t),       intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density__start__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(.not.associated(this%sim))
    nullify(dnst)
    this%sim => sim
    if(dnst_intrf_alloc(this%dnst))then
      call dnst_intrf_get(this%dnst, dnst)
      ASSERT(associated(dnst))
      call dnst_start(dnst, sim)
      nullify(dnst)
    end if

    POP_SUB(base_density__start__)
  end subroutine base_density__start__

  ! ---------------------------------------------------------
  subroutine base_density__acc__(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    type(dnst_t), pointer :: odns, idns

    PUSH_SUB(base_density__acc__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    nullify(odns, idns)
    if(dnst_intrf_alloc(this%dnst))then
      ASSERT(associated(that%config))
      ASSERT(dnst_intrf_assoc(that%dnst))
      ASSERT(associated(that%sim))
      call dnst_intrf_get(this%dnst, odns)
      ASSERT(associated(odns))
      call dnst_intrf_get(that%dnst, idns)
      ASSERT(associated(idns))
      call dnst_acc(odns, idns)
      nullify(odns, idns)
    end if

    POP_SUB(base_density__acc__)
  end subroutine base_density__acc__

  ! ---------------------------------------------------------
  subroutine base_density__update__(this)
    type(base_density_t), intent(inout) :: this

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density__update__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    nullify(dnst)
    if(dnst_intrf_alloc(this%dnst))then
      call dnst_intrf_get(this%dnst, dnst)
      ASSERT(associated(dnst))
      call dnst_update(dnst)
      nullify(dnst)
    end if

    POP_SUB(base_density__update__)
  end subroutine base_density__update__

  ! ---------------------------------------------------------
  subroutine base_density__reset__(this)
    type(base_density_t), intent(inout) :: this

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density__reset__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    nullify(dnst)
    if(dnst_intrf_alloc(this%dnst))then
      call dnst_intrf_get(this%dnst, dnst)
      ASSERT(associated(dnst))
      call dnst_reset(dnst)
      nullify(dnst)
    end if

    POP_SUB(base_density__reset__)
  end subroutine base_density__reset__

  ! ---------------------------------------------------------
  subroutine base_density__stop__(this)
    type(base_density_t), intent(inout) :: this

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density__stop__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    nullify(this%sim, dnst)
    if(dnst_intrf_alloc(this%dnst))then
      call dnst_intrf_get(this%dnst, dnst)
      ASSERT(associated(dnst))
      call dnst_stop(dnst)
      nullify(dnst)
    end if

    POP_SUB(base_density__stop__)
  end subroutine base_density__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_density__apply__(this, operation)
    type(base_density_t), intent(inout) :: this

    interface
      subroutine operation(this)
        import :: base_density_t
        type(base_density_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(base_density_list_iterator_t) :: iter
    type(base_density_t),      pointer :: subs

    PUSH_SUB(base_density__apply__)

    call base_density_list_init(iter, this%list)
    do
      nullify(subs)
      call base_density_list_next(iter, subs)
      if(.not.associated(subs))exit
      call base_density__apply__(subs, operation)
    end do
    call base_density_list_end(iter)
    nullify(subs)
    call operation(this)

    POP_SUB(base_density__apply__)
  end subroutine base_density__apply__

  ! ---------------------------------------------------------
  subroutine base_density_start(this, sim)
    type(base_density_t), intent(inout) :: this
    type(simulation_t),   intent(in)    :: sim

    PUSH_SUB(base_density_start)

    call base_density__apply__(this, start)
    
    POP_SUB(base_density_start)
    
  contains

    subroutine start(this)
      type(base_density_t), intent(inout) :: this

      PUSH_SUB(base_density_start.start)
      
      call base_density__start__(this, sim)

      POP_SUB(base_density_start.start)
    end subroutine start

  end subroutine base_density_start

  ! ---------------------------------------------------------
  recursive subroutine base_density_acc(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_acc)

    ASSERT(associated(this%config))
    if(this%accu) call base_density__reset__(this)
    call base_density_init(iter, this)
    do
      nullify(subs)
      call base_density_next(iter, subs)
      if(.not.associated(subs))exit
      call base_density_acc(subs)
      if(this%accu) call base_density__acc__(this, subs)
    end do
    call base_density_end(iter)
    nullify(subs)
    if(this%accu) call base_density__update__(this)
    
    POP_SUB(base_density_acc)
  end subroutine base_density_acc

  ! ---------------------------------------------------------
  subroutine base_density_update(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density_update)

    call base_density__apply__(this, base_density__update__)
    
    POP_SUB(base_density_update)
  end subroutine base_density_update

  ! ---------------------------------------------------------
  subroutine base_density_reset(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density_reset)

    call base_density__apply__(this, base_density__reset__)
    
    POP_SUB(base_density_reset)
  end subroutine base_density_reset

  ! ---------------------------------------------------------
  subroutine base_density_stop(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density_stop)

    call base_density__apply__(this, base_density__stop__)
    
    POP_SUB(base_density_stop)
  end subroutine base_density_stop

  ! ---------------------------------------------------------
  subroutine base_density_set_sub(this, name, that)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density_set_sub)

    ASSERT(associated(this%config))
    call base_density_dict_set(this%dict, trim(adjustl(name)), that)

    POP_SUB(base_density_set_sub)
  end subroutine base_density_set_sub

  ! ---------------------------------------------------------
  subroutine base_density_get_sub(this, name, that)
    type(base_density_t),  intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    PUSH_SUB(base_density_get_sub)

    ASSERT(associated(this%config))
    nullify(that)
    call base_density_dict_get(this%dict, trim(adjustl(name)), that)

    POP_SUB(base_density_get_sub)
  end subroutine base_density_get_sub

  ! ---------------------------------------------------------
  subroutine base_density_get_sub_density_1d(this, name, that, spin, total)
    type(base_density_t),         intent(in) :: this
    character(len=*),             intent(in) :: name
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_get_sub_density_1d)

    nullify(that, subs)
    call base_density_get(this, name, subs)
    if(associated(subs)) call base_density_get(subs, that, spin, total)

    POP_SUB(base_density_get_sub_density_1d)
  end subroutine base_density_get_sub_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_sub_density_2d(this, name, that, total)
    type(base_density_t),           intent(in) :: this
    character(len=*),               intent(in) :: name
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_get_sub_density_2d)

    nullify(that, subs)
    call base_density_get(this, name, subs)
    if(associated(subs)) call base_density_get(subs, that, total)

    POP_SUB(base_density_get_sub_density_2d)
  end subroutine base_density_get_sub_density_2d

  ! ---------------------------------------------------------
  subroutine base_density_set_info(this, static)
    type(base_density_t), intent(inout) :: this
    logical,              intent(in)    :: static

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_set_info)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_alloc(this%dnst))
    nullify(dnst)
    call dnst_intrf_get(this%dnst, dnst)
    ASSERT(associated(dnst))
    call dnst_set(dnst, static=static)
    nullify(dnst)

    POP_SUB(base_density_set_info)
  end subroutine base_density_set_info

  ! ---------------------------------------------------------
  subroutine base_density_set_charge(this, charge, spin, total)
    type(base_density_t), intent(inout) :: this
    real(kind=wp),        intent(in)    :: charge
    integer,    optional, intent(in)    :: spin
    logical,    optional, intent(in)    :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_set_charge)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_alloc(this%dnst))
    nullify(dnst)
    call dnst_intrf_get(this%dnst, dnst)
    ASSERT(associated(dnst))
    call dnst_set(dnst, charge, spin, total)
    nullify(dnst)

    POP_SUB(base_density_set_charge)
  end subroutine base_density_set_charge

  ! ---------------------------------------------------------
  subroutine base_density_get_info(this, size, nspin, reduce, fine, use)
    type(base_density_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin
    logical,    optional, intent(out) :: reduce
    logical,    optional, intent(out) :: fine
    logical,    optional, intent(out) :: use

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_info)

    nullify(dnst)
    if(present(reduce)) reduce = this%accu
    call dnst_intrf_get(this%dnst, dnst)
    ASSERT(associated(dnst))
    call dnst_get(dnst, size=size, nspin=nspin, fine=fine, use=use)
    nullify(dnst)

    POP_SUB(base_density_get_info)
  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_charge(this, charge, spin, total)
    type(base_density_t), intent(in)  :: this
    real(kind=wp),        intent(out) :: charge
    integer,    optional, intent(in)  :: spin
    logical,    optional, intent(in)  :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_charge)

    nullify(dnst)
    call dnst_intrf_get(this%dnst, dnst)
    ASSERT(associated(dnst))
    call dnst_get(dnst, charge, spin, total)
    nullify(dnst)

    POP_SUB(base_density_get_charge)
  end subroutine base_density_get_charge

  ! ---------------------------------------------------------
  subroutine base_density_get_config(this, that)
    type(base_density_t), target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_density_get_config)

    nullify(that)
    if(associated(this%config)) that=>this%config

    POP_SUB(base_density_get_config)
  end subroutine base_density_get_config

  ! ---------------------------------------------------------
  subroutine base_density_get_simulation(this, that)
    type(base_density_t), target, intent(in) :: this
    type(simulation_t),  pointer             :: that

    PUSH_SUB(base_density_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_density_get_simulation)
  end subroutine base_density_get_simulation

  ! ---------------------------------------------------------
  subroutine base_density_get_storage(this, that, total)
    type(base_density_t), target, intent(in) :: this
    type(storage_t),     pointer             :: that
    logical,            optional, intent(in) :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_storage)

    nullify(that, dnst)
    call dnst_intrf_get(this%dnst, dnst)
    if(associated(dnst))then
      call dnst_get(dnst, that, total)
      nullify(dnst)
    end if

    POP_SUB(base_density_get_storage)
  end subroutine base_density_get_storage

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, spin, total)
    type(base_density_t),         intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_density_1d)

    nullify(that, dnst)
    call dnst_intrf_get(this%dnst, dnst)
    if(associated(dnst))then
      call dnst_get(dnst, that, spin, total)
      nullify(dnst)
    end if

    POP_SUB(base_density_get_density_1d)
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that, total)
    type(base_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_density_2d)

    nullify(that, dnst)
    call dnst_intrf_get(this%dnst, dnst)
    if(associated(dnst))then
      call dnst_get(dnst, that, total)
      nullify(dnst)
    end if
    
    POP_SUB(base_density_get_density_2d)
  end subroutine base_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_density_get_dnst_intrf(this, that)
    type(base_density_t), target, intent(in) :: this
    type(dnst_intrf_t),  pointer             :: that

    PUSH_SUB(base_density_get_dnst_intrf)

    nullify(that)
    if(associated(this%config)) that => this%dnst

    POP_SUB(base_density_get_dnst_intrf)
  end subroutine base_density_get_dnst_intrf

  ! ---------------------------------------------------------
  subroutine base_density__copy__(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density__copy__)

    call base_density__end__(this)
    if(associated(that%config))then
      call base_density__iinit__(this, that%config)
      this%sim => that%sim
      call dnst_intrf_copy(this%dnst, that%dnst, dnst_copy)
    end if

    POP_SUB(base_density__copy__)
  end subroutine base_density__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_density_copy_type(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density_copy_type)

    call base_density_end(this)
    call base_density__copy__(this, that)
    call base_density__build__(this, copy)
    
    POP_SUB(base_density_copy_type)

  contains

    recursive subroutine copy(this, isub, name)
      type(base_density_t), intent(inout) :: this
      type(base_density_t), intent(in)    :: isub
      character(len=*),     intent(in)    :: name

      type(base_density_t), pointer :: osub
      
      POP_SUB(base_density_copy_type.copy)

      nullify(osub)
      if(base_density_list_index(that%list, isub)>0)then
        call base_density_new(this, osub)
        call base_density_copy(osub, isub)
        call base_density_set(this, name, osub)
        nullify(osub)
      else
        call base_density_set(this, name, isub)
      end if

      PUSH_SUB(base_density_copy_type.copy)
    end subroutine copy

  end subroutine base_density_copy_type

  ! ---------------------------------------------------------
  subroutine base_density__end__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__end__)

    nullify(this%config, this%sim, this%prnt)
    call dnst_intrf_end(this%dnst, dnst_end)
    call base_density_dict_end(this%dict)
    call base_density_list_end(this%list)

    POP_SUB(base_density__end__)
  end subroutine base_density__end__

  ! ---------------------------------------------------------
  recursive subroutine base_density_end_type(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_end_type)

    do
      nullify(subs)
      call base_density_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_density_end(subs)
      call base_density__del__(subs)
    end do
    nullify(subs)
    call base_density__end__(this)

    POP_SUB(base_density_end_type)
  end subroutine base_density_end_type

#define TEMPLATE_PREFIX base_density
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_density_oct_m

!! Local Variables:
!! mode: f90
!! End:

