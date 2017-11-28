#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_density_oct_m

  use dnst_oct_m
  use dnst_intrf_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use message_oct_m
  use messages_oct_m
  use msgbus_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

#define BASE_TEMPLATE_NAME base_density
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::         &
    base_density_t

  public ::                 &
    base_density__init__,   &
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
    base_density_notify, &
    base_density_update, &
    base_density_reset,  &
    base_density_stop,   &
    base_density_set,    &
    base_density_get,    &
    base_density_copy,   &
    base_density_end

#define BASE_TEMPLATE_NAME base_density
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_density_t
    private
    type(json_object_t),  pointer :: config =>null()
    type(simulation_t),   pointer :: sim    =>null()
    type(base_density_t), pointer :: prnt   =>null()
    type(dnst_intrf_t)            :: dnst
    type(msgbus_t)                :: msgb
    type(base_density_dict_t)     :: dict
    type(base_density_list_t)     :: list
  end type base_density_t

  interface base_density__init__
    module procedure base_density__init__type
    module procedure base_density__init__pass
    module procedure base_density__init__copy
  end interface base_density__init__

  interface base_density_new
    module procedure base_density_new_type
  end interface base_density_new

  interface base_density_init
    module procedure base_density_init_type
  end interface base_density_init

  interface base_density_notify
    module procedure base_density_notify_type
    module procedure base_density_notify_subs
  end interface base_density_notify

  interface base_density_set
    module procedure base_density_set_info
    module procedure base_density_set_charge
    module procedure base_density_set_dnst
  end interface base_density_set

  interface base_density_get
    module procedure base_density_get_info
    module procedure base_density_get_charge
    module procedure base_density_get_config
    module procedure base_density_get_simulation
    module procedure base_density_get_msgbus
    module procedure base_density_get_dnst
    module procedure base_density_get_storage
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
  end interface base_density_get

  interface base_density_gets
    module procedure base_density_gets_dnst
    module procedure base_density_gets_storage
    module procedure base_density_gets_density_1d
    module procedure base_density_gets_density_2d
  end interface base_density_gets

contains

#define BASE_TEMPLATE_NAME base_density
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_density_new_type(this, that)
    type(base_density_t),  target, intent(inout) :: this
    type(base_density_t), pointer, intent(out)   :: that

    PUSH_SUB(base_density_new_type)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt => this
    call base_density_list_push(this%list, that)

    POP_SUB(base_density_new_type)
  end subroutine base_density_new_type

  ! ---------------------------------------------------------
  subroutine base_density__iinit__(this, config)
    type(base_density_t),        intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_density__iinit__)

    this%config => config
    call msgbus_init(this%msgb)
    call base_density_dict_init(this%dict)
    call base_density_list_init(this%list)

    POP_SUB(base_density__iinit__)
  end subroutine base_density__iinit__

  ! ---------------------------------------------------------
  subroutine base_density__init__type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    character(len=BASE_DENSITY_NAME_LEN) :: name
    integer                              :: ierr
    logical                              :: dflt, refs, accu

    PUSH_SUB(base_density__init__type)

    call base_density__iinit__(this, config)
    call dnst_intrf_init(this%dnst)
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    if(dflt)then
      refs = .false.
      call json_get(this%config, "reference", name, ierr)
      if(ierr==JSON_OK) refs = .true.
      call json_get(this%config, "reduce", accu, ierr)
      if(ierr/=JSON_OK) accu = .false.
      ASSERT(.not.(refs.and.accu))
      if(.not.refs) call dnst_intrf_new(this%dnst, this%config)
    end if

    POP_SUB(base_density__init__type)
  end subroutine base_density__init__type

  ! ---------------------------------------------------------
  subroutine base_density__init__pass(this, init)
    type(base_density_t), intent(inout) :: this

    interface
      subroutine init(this, config)
        use dnst_oct_m
        use json_oct_m
        type(dnst_t),        intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    integer :: ierr
    logical :: dflt

    PUSH_SUB(base_density__init__pass)

    ASSERT(associated(this%config))
    ASSERT(.not.dnst_intrf_assoc(this%dnst))
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    ASSERT(.not.dflt)
    call dnst_intrf_new(this%dnst, iinit)
    ASSERT(dnst_intrf_alloc(this%dnst))

    POP_SUB(base_density__init__pass)

  contains

    subroutine iinit(dnst)
      type(dnst_t), intent(out) :: dnst

      PUSH_SUB(base_density__init__pass.iinit)

      call init(dnst, this%config)

      POP_SUB(base_density__init__pass.iinit)
    end subroutine iinit

  end subroutine base_density__init__pass

  ! ---------------------------------------------------------
  subroutine base_density__init__copy(this, that, start)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that
    logical,    optional, intent(in)  :: start

    logical :: istr

    PUSH_SUB(base_density__init__copy)

    ASSERT(associated(that%config))
    call base_density__iinit__(this, that%config)
    call dnst_intrf_init(this%dnst, that%dnst, start=.false.)
    istr = .true.
    if(present(start)) istr = start
    if(istr)then
      if(present(start))then
        ASSERT(associated(that%sim))
      end if
      if(associated(that%sim)) call base_density__start__(this, that%sim)
    end if

    POP_SUB(base_density__init__copy)
  end subroutine base_density__init__copy

  ! ---------------------------------------------------------
  subroutine base_density_init_type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(base_density_init_type)

    call base_density__init__(this, config)

    POP_SUB(base_density_init_type)
  end subroutine base_density_init_type

  ! ---------------------------------------------------------
  subroutine base_density__start__(this, sim)
    type(base_density_t),       intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    ASSERT(dnst_intrf_assoc(this%dnst))
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
    logical               :: sttc

    PUSH_SUB(base_density__update__)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    nullify(dnst)
    call base_density_get(this, static=sttc)
    if(.not.sttc.and.dnst_intrf_alloc(this%dnst))then
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
  subroutine base_density__recv__(this, update)
    type(base_density_t), intent(in)  :: this
    logical,              intent(out) :: update

    type(msgbus_iterator_t)      :: iter
    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg
    logical                      :: updt
    integer                      :: ierr

    PUSH_SUB(base_density__recv__)

    update = .false.
    call msgbus_init(iter, this%msgb)
    do
      nullify(data, mssg)
      call msgbus_next(iter, mssg)
      if(.not.associated(mssg))exit
      call message_get(mssg, data)
      ASSERT(associated(data))
      call json_get(data, "update", updt, ierr)
      if(ierr==JSON_OK)then
        update = (update .or. updt)
        call msgbus_remove(iter, ierr=ierr)
        ASSERT(ierr==MSGBUS_OK)
      end if
    end do
    nullify(data, mssg)
    
    POP_SUB(base_density__recv__)
  end subroutine base_density__recv__
    
  ! ---------------------------------------------------------
  subroutine base_density_notify_type(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg

    PUSH_SUB(base_density_notify_type)

    nullify(data, mssg)
    mssg => message_new()
    call message_get(mssg, data)
    ASSERT(associated(data))
    call json_set(data, "update", .true.)
    nullify(data)
    call msgbus_notify(this%msgb, mssg)
    call msgbus_publish(this%msgb, mssg)
    nullify(mssg)
    
    POP_SUB(base_density_notify_type)
  end subroutine base_density_notify_type
    
  ! ---------------------------------------------------------
  subroutine base_density_notify_subs(this, name)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: name

    type(base_density_t), pointer :: subs
    
    PUSH_SUB(base_density_notify_subs)

    nullify(subs)
    call base_density_gets(this, trim(adjustl(name)), subs)
    ASSERT(associated(subs))
    call base_density_notify(subs)
    nullify(subs)
    
    POP_SUB(base_density_notify_subs)
  end subroutine base_density_notify_subs
    
  ! ---------------------------------------------------------
  subroutine base_density_start(this, sim)
    type(base_density_t), intent(inout) :: this
    type(simulation_t),   intent(in)    :: sim

    PUSH_SUB(base_density_start)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(.not.associated(this%sim))
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

    PUSH_SUB(base_density_acc)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    ASSERT(base_density_dict_len(this%dict)>0)
    call base_density__reset__(this)
    call base_density__reduce__(this, base_density__acc__)
    call base_density__update__(this)
    
    POP_SUB(base_density_acc)
  end subroutine base_density_acc

  ! ---------------------------------------------------------
  recursive subroutine base_density_update(this)
    type(base_density_t), intent(inout) :: this

    type(base_density_iterator_t) :: iter
    type(base_density_t), pointer :: subs
    logical                       :: sttc, updt, accu
    integer                       :: ierr

    PUSH_SUB(base_density_update)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    ASSERT(associated(this%sim))
    call base_density_get(this, static=sttc)
    call base_density__recv__(this, updt)
    ASSERT(.not.(sttc.and.updt))
    if(.not.sttc.and.updt)then
      call json_get(this%config, "reduce", accu, ierr)
      if(ierr/=JSON_OK) accu = .false.
      accu = (accu .and. dnst_intrf_alloc(this%dnst))
      if(accu)then
        ASSERT(base_density_dict_len(this%dict)>0)
        call base_density__reset__(this)
      end if
      call base_density_init(iter, this)
      do
        nullify(subs)
        call base_density_next(iter, subs)
        if(.not.associated(subs))exit
        call base_density_update(subs)
        if(accu) call base_density__acc__(this, subs)
      end do
      call base_density_end(iter)
      nullify(subs)
      call base_density__update__(this)
    end if

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
  subroutine base_density__sets__(this, name, that)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(base_density_t), intent(in)    :: that

    character(len=BASE_DENSITY_NAME_LEN) :: rnam
    type(dnst_t),                pointer :: dnst
    integer                              :: indx, ierr
    logical                              :: dflt, accu

    PUSH_SUB(base_density__sets__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    nullify(dnst)
    call msgbus_attach(this%msgb, that%msgb)
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    if(dflt)then
      call json_get(this%config, "reduce", accu, ierr)
      if(ierr/=JSON_OK) accu = .false.
      call json_get(this%config, "reference", rnam, ierr)
      if(ierr==JSON_OK)then
        ASSERT(.not.accu)
        indx = index(trim(adjustl(rnam)), trim(adjustl(name)))
        if(indx==1)then
          ASSERT(.not.dnst_intrf_assoc(this%dnst))
          call base_density_gets(this, trim(adjustl(rnam)), dnst)
          ASSERT(associated(dnst))
          call dnst_intrf_set(this%dnst, dnst)
          nullify(dnst)
        end if
      end if
    end if

    POP_SUB(base_density__sets__)
  end subroutine base_density__sets__
  
  ! ---------------------------------------------------------
  subroutine base_density__dels__(this, name, that)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: name
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density__dels__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    call msgbus_detach(this%msgb, that%msgb)

    POP_SUB(base_density__dels__)
  end subroutine base_density__dels__

  ! ---------------------------------------------------------
  subroutine base_density_gets_dnst(this, name, that)
    type(base_density_t),  intent(in)  :: this
    character(len=*),      intent(in)  :: name
    type(dnst_t), pointer, intent(out) :: that

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_gets_dnst)

    nullify(that, subs)
    call base_density_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_density_get(subs, that)
    nullify(subs)

    POP_SUB(base_density_gets_dnst)
  end subroutine base_density_gets_dnst

  ! ---------------------------------------------------------
  subroutine base_density_gets_storage(this, name, that, total)
    type(base_density_t),     intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(storage_t), pointer, intent(out) :: that
    logical,        optional, intent(in)  :: total

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_gets_storage)

    nullify(that, subs)
    call base_density_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_density_get(subs, that, total)
    nullify(subs)

    POP_SUB(base_density_gets_storage)
  end subroutine base_density_gets_storage

  ! ---------------------------------------------------------
  subroutine base_density_gets_density_1d(this, name, that, spin, total)
    type(base_density_t),                 intent(in)  :: this
    character(len=*),                     intent(in)  :: name
    real(kind=wp), dimension(:), pointer, intent(out) :: that
    integer,                    optional, intent(in)  :: spin
    logical,                    optional, intent(in)  :: total

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_gets_density_1d)

    nullify(that, subs)
    call base_density_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_density_get(subs, that, spin, total)
    nullify(subs)

    POP_SUB(base_density_gets_density_1d)
  end subroutine base_density_gets_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_gets_density_2d(this, name, that, total)
    type(base_density_t),                   intent(in)  :: this
    character(len=*),                       intent(in)  :: name
    real(kind=wp), dimension(:,:), pointer, intent(out) :: that
    logical,                      optional, intent(in)  :: total

    type(base_density_t), pointer :: subs

    PUSH_SUB(base_density_gets_density_2d)

    nullify(that, subs)
    call base_density_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_density_get(subs, that, total)
    nullify(subs)

    POP_SUB(base_density_gets_density_2d)
  end subroutine base_density_gets_density_2d

  ! ---------------------------------------------------------
  subroutine base_density_set_info(this, static)
    type(base_density_t), intent(inout) :: this
    logical,              intent(in)    :: static

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_set_info)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_alloc(this%dnst))
    nullify(dnst)
    call base_density_get(this, dnst)
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
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dnst_set(dnst, charge, spin, total)
    nullify(dnst)

    POP_SUB(base_density_set_charge)
  end subroutine base_density_set_charge

  ! ---------------------------------------------------------
  subroutine base_density_set_dnst(this, that)
    type(base_density_t), intent(inout) :: this
    type(dnst_t),         intent(in)    :: that

    logical :: dflt
    integer :: ierr

    PUSH_SUB(base_density_set_dnst)

    ASSERT(associated(this%config))
    ASSERT(.not.dnst_intrf_assoc(this%dnst))
    call json_get(this%config, "default", dflt, ierr)
    if(ierr/=JSON_OK) dflt = .true.
    ASSERT(.not.dflt)
    call dnst_intrf_set(this%dnst, that)
    ASSERT(dnst_intrf_assoc(this%dnst))

    POP_SUB(base_density_set_dnst)
  end subroutine base_density_set_dnst

  ! ---------------------------------------------------------
  subroutine base_density_get_info(this, size, nspin, static, fine, use)
    type(base_density_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin
    logical,    optional, intent(out) :: static
    logical,    optional, intent(out) :: fine
    logical,    optional, intent(out) :: use

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_info)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    nullify(dnst)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dnst_get(dnst, size=size, nspin=nspin, static=static, fine=fine, use=use)
    nullify(dnst)

    POP_SUB(base_density_get_info)
  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_config(this, that)
    type(base_density_t), target, intent(in)  :: this
    type(json_object_t), pointer, intent(out) :: that

    PUSH_SUB(base_density_get_config)

    nullify(that)
    if(associated(this%config)) that=>this%config

    POP_SUB(base_density_get_config)
  end subroutine base_density_get_config

  ! ---------------------------------------------------------
  subroutine base_density_get_simulation(this, that)
    type(base_density_t), target, intent(in)  :: this
    type(simulation_t),  pointer, intent(out) :: that

    PUSH_SUB(base_density_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_density_get_simulation)
  end subroutine base_density_get_simulation

  ! ---------------------------------------------------------
  subroutine base_density_get_msgbus(this, that)
    type(base_density_t), target, intent(in)  :: this
    type(msgbus_t),      pointer, intent(out) :: that

    PUSH_SUB(base_density_get_msgbus)

    nullify(that)
    if(associated(this%config)) that => this%msgb

    POP_SUB(base_density_get_msgbus)
  end subroutine base_density_get_msgbus

  ! ---------------------------------------------------------
  subroutine base_density_get_charge(this, charge, spin, total)
    type(base_density_t), intent(in)  :: this
    real(kind=wp),        intent(out) :: charge
    integer,    optional, intent(in)  :: spin
    logical,    optional, intent(in)  :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_charge)

    ASSERT(associated(this%config))
    ASSERT(dnst_intrf_assoc(this%dnst))
    nullify(dnst)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dnst_get(dnst, charge, spin, total)
    nullify(dnst)

    POP_SUB(base_density_get_charge)
  end subroutine base_density_get_charge

  ! ---------------------------------------------------------
  subroutine base_density_get_dnst(this, that)
    type(base_density_t),  intent(in)  :: this
    type(dnst_t), pointer, intent(out) :: that

    PUSH_SUB(base_density_get_dnst)

    nullify(that)
    if(associated(this%config))&
      call dnst_intrf_get(this%dnst, that)

    POP_SUB(base_density_get_dnst)
  end subroutine base_density_get_dnst

  ! ---------------------------------------------------------
  subroutine base_density_get_storage(this, that, total)
    type(base_density_t), target, intent(in)  :: this
    type(storage_t),     pointer, intent(out) :: that
    logical,            optional, intent(in)  :: total

    type(dnst_t), pointer :: dnst

    PUSH_SUB(base_density_get_storage)

    nullify(that, dnst)
    if(dnst_intrf_assoc(this%dnst))then
      call base_density_get(this, dnst)
      ASSERT(associated(dnst))
      call dnst_get(dnst, that, total)
      nullify(dnst)
    end if

    POP_SUB(base_density_get_storage)
  end subroutine base_density_get_storage

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, spin, total)
    type(base_density_t),                 intent(inout) :: this
    real(kind=wp), dimension(:), pointer, intent(out)   :: that
    integer,                    optional, intent(in)    :: spin
    logical,                    optional, intent(in)    :: total

    type(dnst_t), pointer :: dnst
    logical               :: sttc

    PUSH_SUB(base_density_get_density_1d)

    nullify(that, dnst)
    if(dnst_intrf_assoc(this%dnst))then
      call base_density_update(this)
      call base_density_get(this, static=sttc)
      call base_density_set(this, static=.false.)
      call base_density_get(this, dnst)
      ASSERT(associated(dnst))
      call dnst_get(dnst, that, spin, total)
      nullify(dnst)
      call base_density_set(this, static=sttc)
    end if

    POP_SUB(base_density_get_density_1d)
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that, total)
    type(base_density_t),                   intent(inout) :: this
    real(kind=wp), dimension(:,:), pointer, intent(out)   :: that
    logical,                      optional, intent(in)    :: total

    type(dnst_t), pointer :: dnst
    logical               :: sttc

    PUSH_SUB(base_density_get_density_2d)

    nullify(that, dnst)
    if(dnst_intrf_assoc(this%dnst))then
      call base_density_update(this)
      call base_density_get(this, static=sttc)
      call base_density_set(this, static=.false.)
      call base_density_get(this, dnst)
      ASSERT(associated(dnst))
      call dnst_get(dnst, that, total)
      nullify(dnst)
      call base_density_set(this, static=sttc)
    end if
    
    POP_SUB(base_density_get_density_2d)
  end subroutine base_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_density__copy__(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(base_density__copy__)

    call base_density__end__(this)
    if(associated(that%config))then
      call base_density__iinit__(this, that%config)
      this%sim => that%sim
      call dnst_intrf_copy(this%dnst, that%dnst)
    end if

    POP_SUB(base_density__copy__)
  end subroutine base_density__copy__

  ! ---------------------------------------------------------
  subroutine base_density__end__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(base_density__end__)

    nullify(this%config, this%sim, this%prnt)
    call dnst_intrf_end(this%dnst)
    call msgbus_end(this%msgb)
    call base_density_dict_end(this%dict)
    ASSERT(base_density_list_len(this%list)==0)
    call base_density_list_end(this%list)

    POP_SUB(base_density__end__)
  end subroutine base_density__end__

end module base_density_oct_m

!! Local Variables:
!! mode: f90
!! End:

