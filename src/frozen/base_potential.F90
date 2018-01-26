#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_potential_oct_m

  use atom_oct_m
  use base_density_oct_m
  use base_geometry_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use global_oct_m
  use intrpl_oct_m
  use json_oct_m
  use kinds_oct_m
  use memo_oct_m
  use message_oct_m
  use messages_oct_m
  use msgbus_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use species_oct_m
  use storage_oct_m

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private
  
  public ::         &
    POTN_KIND_NONE, &
    POTN_KIND_XTRN

  public ::           &
    base_potential_t

  public ::                   &
    base_potential__init__,   &
    base_potential__start__,  &
    base_potential__acc__,    &
    base_potential__update__, &
    base_potential__reset__,  &
    base_potential__stop__,   &
    base_potential__copy__,   &
    base_potential__end__

  public ::                &
    base_potential_new,    &
    base_potential_del,    &
    base_potential_init,   &
    base_potential_start,  &
    base_potential_calc,   &
    base_potential_notify, &
    base_potential_update, &
    base_potential_stop,   &
    base_potential_set,    &
    base_potential_get,    &
    base_potential_copy,   &
    base_potential_end

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  integer, parameter :: default_nspin = 1

  integer, parameter :: POTN_KIND_NONE = 0
  integer, parameter :: POTN_KIND_XTRN = 1

  type :: base_potential_t
    private
    type(json_object_t),  pointer :: config  =>null()
    type(base_system_t),  pointer :: sys     =>null()
    type(simulation_t),   pointer :: sim     =>null()
    type(refcount_t),     pointer :: rcnt    =>null()
    integer                       :: kind    = POTN_KIND_NONE
    integer                       :: nspin   = default_nspin
    type(memo_t)                  :: memo
    type(storage_t)               :: data
    type(msgbus_t)                :: msgb
    type(base_potential_dict_t)   :: dict
  end type base_potential_t

  interface base_potential__init__
    module procedure base_potential__init__type
    module procedure base_potential__init__copy
  end interface base_potential__init__

  interface base_potential__load__
    module procedure base_potential__load__r1
    module procedure base_potential__load__r2
  end interface base_potential__load__
  
  interface base_potential__sets__
    module procedure base_potential__sets__info
    module procedure base_potential__sets__type
  end interface base_potential__sets__

  interface base_potential_new
    module procedure base_potential_new_type
    module procedure base_potential_new_pass
  end interface base_potential_new

  interface base_potential_init
    module procedure base_potential_init_type
  end interface base_potential_init

  interface base_potential_notify
    module procedure base_potential_notify_type
    module procedure base_potential_notify_subs
  end interface base_potential_notify

  interface base_potential_set
    module procedure base_potential_set_info
  end interface base_potential_set

  interface base_potential_get
    module procedure base_potential_get_info
    module procedure base_potential_get_energy
    module procedure base_potential_get_config
    module procedure base_potential_get_system
    module procedure base_potential_get_density
    module procedure base_potential_get_simulation
    module procedure base_potential_get_msgbus
    module procedure base_potential_get_storage
    module procedure base_potential_get_potential_r1
    module procedure base_potential_get_potential_r2
    module procedure base_potential_get_sub_storage
    module procedure base_potential_get_sub_potential_r1
    module procedure base_potential_get_sub_potential_r2
  end interface base_potential_get

contains

#define BASE_TEMPLATE_NAME base_potential
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_potential_new_type(sys, config) result(this)
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    type(base_potential_t), pointer :: this

    PUSH_SUB(base_potential_new_type)

    this => base_potential_new(sys, config, base_potential_init_type)

    POP_SUB(base_potential_new_type)
  end function base_potential_new_type

  ! ---------------------------------------------------------
  function base_potential_new_pass(sys, config, init) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_potential_t), pointer :: this
    
    interface
      subroutine init(this, sys, config)
        use json_oct_m
        use base_system_oct_m
        import :: base_potential_t
        type(base_potential_t), intent(out) :: this
        type(base_system_t),    intent(in)  :: sys
        type(json_object_t),    intent(in)  :: config
      end subroutine init
    end interface
    
    PUSH_SUB(base_potential_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, sys, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_potential_new_pass)
  end function base_potential_new_pass

  ! ---------------------------------------------------------
  subroutine base_potential__init__type(this, sys, config)
    type(base_potential_t),      intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    logical                      :: uspn
    integer                      :: ierr

    PUSH_SUB(base_potential__init__type)

    nullify(cnfg)
    this%config => config
    this%sys => sys
    this%rcnt => refcount_new()
    call json_get(this%config, "kind", this%kind, ierr)
    if(ierr/=JSON_OK) this%kind = POTN_KIND_NONE
    this%nspin = default_nspin
    call json_get(this%config, "spin", uspn, ierr)
    if(ierr/=JSON_OK) uspn = .false.
    if(uspn) call base_system_get(this%sys, nspin=this%nspin)
    ASSERT(this%nspin>0)
    ASSERT(this%nspin<3)
    select case(this%kind)
    case(POTN_KIND_XTRN)
      ASSERT(this%nspin==1)
    case default
      ASSERT(.false.)
    end select
    call memo_init(this%memo)
    call json_get(this%config, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, "full", .false.)
    call json_set(cnfg, "dimensions", this%nspin)
    call storage_init(this%data, cnfg)
    nullify(cnfg)
    call msgbus_init(this%msgb, number=2)
    call base_potential_dict_init(this%dict)

    POP_SUB(base_potential__init__type)
  end subroutine base_potential__init__type

  ! ---------------------------------------------------------
  subroutine base_potential__init__copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(base_potential__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_potential__init__(this, that%sys, that%config)
    if(associated(that%sim)) call base_potential__start__(this, that%sim)

    POP_SUB(base_potential__init__copy)
  end subroutine base_potential__init__copy

  ! ---------------------------------------------------------
  subroutine base_potential_init_type(this, sys, config)
    type(base_potential_t), intent(out) :: this
    type(base_system_t),    intent(in)  :: sys
    type(json_object_t),    intent(in)  :: config

    PUSH_SUB(base_potential_init_type)

    call base_potential__init__(this, sys, config)

    POP_SUB(base_potential_init_type)
  end subroutine base_potential_init_type

  ! ---------------------------------------------------------
  subroutine base_potential__start__(this, sim)
    type(base_potential_t),     intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    type(base_density_t), pointer :: dnst
    type(msgbus_t),       pointer :: msgb

    PUSH_SUB(base_potential__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    nullify(dnst, msgb)
    this%sim => sim
    call storage_start(this%data, this%sim)
    call base_potential_get(this, dnst)
    ASSERT(associated(dnst))
    call base_density_get(dnst, msgb)
    ASSERT(associated(msgb))
    call msgbus_attach(this%msgb, msgb, id=2)
    nullify(dnst, msgb)

    POP_SUB(base_potential__start__)
  end subroutine base_potential__start__

  ! ---------------------------------------------------------
  subroutine base_potential__xtrn__(this, x, v)
    type(base_potential_t),      intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v

    type(base_geometry_t), pointer :: geom
    type(atom_t),          pointer :: atom
    type(base_geometry_iterator_t) :: iter
    integer                        :: nd

    PUSH_SUB(base_potential__xtrn__)

    v = 0.0_wp
    call simulation_get(this%sim, ndim=nd)
    call base_system_get(this%sys, geom)
    ASSERT(associated(geom))
    call base_geometry_init(iter, geom)
    do
      nullify(atom)
      call base_geometry_next(iter, atom)
      if(.not.associated(atom))exit
      ASSERT(associated(atom%species))
      v = v + calc(x, atom%x(1:nd), species_zval(atom%species))
    end do
    nullify(atom)
    call base_geometry_end(iter)
    nullify(geom)

    POP_SUB(base_potential__xtrn__)

  contains

    pure function calc(x, y, c) result(v)
      real(kind=wp), dimension(:), intent(in)  :: x
      real(kind=wp), dimension(:), intent(in)  :: y
      real(kind=wp),               intent(in)  :: c

      real(kind=wp) :: v

      real(kind=wp) :: r

      r = sqrt(sum((x-y)**2))
      if(r<r_small) r = r_small
      v = -c / r

    end function calc

  end subroutine base_potential__xtrn__

  ! ---------------------------------------------------------
  subroutine base_potential__interpolate__(this, that, type, list, operation)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    integer,                intent(in)    :: type
    type(json_array_t),     intent(in)    :: list

    interface
      subroutine operation(this, that)
        import :: base_potential_t
        type(base_potential_t), intent(inout) :: this
        type(base_potential_t), intent(in)    :: that
      end subroutine operation
    end interface

    type(json_object_t), pointer :: cnfg
    type(json_array_iterator_t)  :: iter
    type(base_potential_t)       :: potn
    type(intrpl_t)               :: intp
    integer                      :: ierr

    PUSH_SUB(base_potential__interpolate__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%kind==that%kind)
    ASSERT(this%nspin==that%nspin)
    nullify(cnfg)
    call base_potential__init__(potn, this)
    call base_potential_set(potn, static=.false.)
    call intrpl_init(intp, that%data, type)
    ASSERT(json_len(list)>0)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call intrpl_attach(intp, cnfg)
      select case(that%kind)
      case(POTN_KIND_XTRN)
        call base_potential__load__(potn, sfnc)
      case default
        ASSERT(.false.)
      end select
      call intrpl_detach(intp)
      call operation(this, potn)
    end do
    call json_end(iter)
    nullify(cnfg)
    call intrpl_end(intp)
    call base_potential_end(potn)
    
    POP_SUB(base_potential__interpolate__)

  contains

    function sfnc(x) result(val)
      real(kind=wp), dimension(:), intent(in) :: x

      real(kind=wp) :: val

      PUSH_SUB(base_potential__interpolate__.sfnc)

      call intrpl_eval(intp, x, val, ifnc)

      POP_SUB(base_potential__interpolate__.sfnc)
    end function sfnc

    function ifnc(x) result(val)
      real(kind=wp), dimension(:), intent(in) :: x

      real(kind=wp) :: val

      PUSH_SUB(base_potential__interpolate__.ifnc)

      call base_potential__xtrn__(that, x, val)

      POP_SUB(base_potential__interpolate__.ifnc)
    end function ifnc

  end subroutine base_potential__interpolate__

  ! ---------------------------------------------------------
  subroutine base_potential__op__(this, that, config, operation)
    type(base_potential_t),        intent(inout) :: this
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    interface
      subroutine operation(this, that)
        import :: base_potential_t
        type(base_potential_t), intent(inout) :: this
        type(base_potential_t), intent(in)    :: that
      end subroutine operation
    end interface

    type(json_array_t), pointer :: list
    integer                     :: type, ierr
    
    PUSH_SUB(base_potential__op__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    nullify(list)
    if(present(config))then
      call json_get(config, "interpolation", type, ierr)
      if(ierr/=JSON_OK) type = INTRPL_DEFAULT
      call json_get(config, "positions", list, ierr)
      if(ierr/=JSON_OK) nullify(list)
    end if
    if(associated(list))then
      call base_potential__interpolate__(this, that, type, list, operation)
    else
      call operation(this, that)
    end if

    POP_SUB(base_potential__op__)
  end subroutine base_potential__op__

  ! ---------------------------------------------------------
  subroutine base_potential__iacc__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(base_potential__iacc__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    call storage_add(this%data, that%data)

    POP_SUB(base_potential__iacc__)
  end subroutine base_potential__iacc__

  ! ---------------------------------------------------------
  subroutine base_potential__acc__(this, that, config)
    type(base_potential_t),        intent(inout) :: this
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_potential__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    call base_potential__op__(this, that, config, base_potential__iacc__)

    POP_SUB(base_potential__acc__)
  end subroutine base_potential__acc__

  ! ---------------------------------------------------------
  subroutine base_potential__isub__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(base_potential__isub__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    call storage_sub(this%data, that%data)

    POP_SUB(base_potential__isub__)
  end subroutine base_potential__isub__

  ! ---------------------------------------------------------
  subroutine base_potential__sub__(this, that, config)
    type(base_potential_t),        intent(inout) :: this
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_potential__sub__)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(this%nspin==that%nspin)
    call base_potential__op__(this, that, config, base_potential__isub__)

    POP_SUB(base_potential__sub__)
  end subroutine base_potential__sub__

  ! ---------------------------------------------------------
  function base_potential__calc__(this, that) result(energy)
    type(base_potential_t), intent(in) :: this
    type(base_density_t),   intent(in) :: that
    
    real(kind=wp) :: energy

    type(storage_t), pointer :: data
    logical                  :: fuse
    integer                  :: nspn

    PUSH_SUB(base_potential__calc__)

    nullify(data)
    energy = 0.0_wp
    call base_density_get(that, nspin=nspn, use=fuse)
    ASSERT(this%nspin<=nspn)
    if(fuse)then
      call base_density_get(that, data, total=(this%nspin<nspn))
      ASSERT(associated(data))
      call storage_integrate(this%data, data, energy)
      nullify(data)
    end if

    POP_SUB(base_potential__calc__)
  end function base_potential__calc__

  ! ---------------------------------------------------------
  subroutine base_potential__load__r1(this, func)
    type(base_potential_t), intent(inout) :: this

    interface
      function func(x) result(f)
        use kinds_oct_m
        implicit none
        real(kind=wp), dimension(:), intent(in) :: x
        real(kind=wp)                           :: f
      end function func
    end interface

    logical :: accu
    integer :: ierr

    PUSH_SUB(base_potential__load__r1)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(this%nspin==1)
    call json_get(this%config, "reduce", accu, ierr)
    if(ierr/=JSON_OK) accu = .false.
    ASSERT(.not.accu)
    call base_potential__reset__(this)
    call storage_load(this%data, func)
    call base_potential__update__(this)
    
    POP_SUB(base_potential__load__r1)
  end subroutine base_potential__load__r1
  
  ! ---------------------------------------------------------
  subroutine base_potential__load__r2(this, func)
    type(base_potential_t), intent(inout) :: this

    interface
      subroutine func(x, f)
        use kinds_oct_m
        real(kind=wp), dimension(:), intent(in)  :: x
        real(kind=wp), dimension(:), intent(out) :: f
      end subroutine func
    end interface

    logical :: accu
    integer :: ierr

    PUSH_SUB(base_potential__load__r2)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call json_get(this%config, "reduce", accu, ierr)
    if(ierr/=JSON_OK) accu = .false.
    ASSERT(.not.accu)
    call base_potential__reset__(this)
    call storage_load(this%data, func)
    call base_potential__update__(this)

    POP_SUB(base_potential__load__r2)
  end subroutine base_potential__load__r2
  
  ! ---------------------------------------------------------
  subroutine base_potential__update__(this, energy)
    type(base_potential_t), intent(inout) :: this
    logical,      optional, intent(in)    :: energy

    type(base_density_t), pointer :: dnst
    logical                       :: fuse, calc, updt, accu
    integer                       :: ierr

    PUSH_SUB(base_potential__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(dnst)
    call base_potential_get(this, use=fuse)
    if(fuse)then
      calc = .false.
      if(present(energy)) calc = energy
      call base_potential__recv__(this, updt, channel=1)
      if(updt)then
        call json_get(this%config, "reduce", accu, ierr)
        if(ierr/=JSON_OK) accu = .false.
        if(accu)then
          call base_potential__reset__(this)
          call base_potential__reduce__(this, reduce)
        else
          call memo_del(this%memo, "energy")
        end if
        call storage_update(this%data)
        call base_potential__publish__(this)
      end if
      if(calc)then
        call base_potential_get(this, dnst)
        ASSERT(associated(dnst))
        call base_density__update__(dnst)
        call base_potential__recv__(this, updt, channel=2)
        if(.not.memo_in(this%memo, "energy").or.updt)&
          call memo_set(this%memo, "energy", base_potential__calc__(this, dnst))
        nullify(dnst)
      end if
    end if

    POP_SUB(base_potential__update__)

  contains

    subroutine reduce(this, that, config)
      type(base_potential_t),        intent(inout) :: this
      type(base_potential_t),        intent(in)    :: that
      type(json_object_t), optional, intent(in)    :: config

      logical :: fuse

      PUSH_SUB(base_potential__update__.reduce)

      ASSERT(associated(that%config))
      ASSERT(associated(that%sim))
      call base_potential_get(that, use=fuse)
      if(fuse) call base_potential__acc__(this, that, config)

      POP_SUB(base_potential__update__.reduce)
    end subroutine reduce

  end subroutine base_potential__update__

  ! ---------------------------------------------------------
  subroutine base_potential__reset__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call memo_del(this%memo, "energy")
    call storage_reset(this%data)

    POP_SUB(base_potential__reset__)
  end subroutine base_potential__reset__

  ! ---------------------------------------------------------
  subroutine base_potential__stop__(this)
    type(base_potential_t), intent(inout) :: this

    type(base_density_t), pointer :: dnst
    type(msgbus_t),       pointer :: msgb

    PUSH_SUB(base_potential__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim, dnst, msgb)
    call storage_stop(this%data)
    call base_potential_get(this, dnst)
    ASSERT(associated(dnst))
    call base_density_get(dnst, msgb)
    ASSERT(associated(msgb))
    call msgbus_detach(this%msgb, msgb, id=2)
    nullify(dnst, msgb)

    POP_SUB(base_potential__stop__)
  end subroutine base_potential__stop__

  ! ---------------------------------------------------------
  subroutine base_potential__recv__(this, update, channel)
    type(base_potential_t), intent(in)  :: this
    logical,                intent(out) :: update
    integer,      optional, intent(in)  :: channel

    type(msgbus_iterator_t)      :: iter
    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg
    logical                      :: updt
    integer                      :: chid, ierr

    PUSH_SUB(base_potential__recv__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    chid = 1
    if(present(channel)) chid = channel
    update = .false.
    call msgbus_init(iter, this%msgb, id=chid)
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

    POP_SUB(base_potential__recv__)
  end subroutine base_potential__recv__

  ! ---------------------------------------------------------
  subroutine base_potential__publish__(this)
    type(base_potential_t), intent(inout) :: this

    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg

    PUSH_SUB(base_potential__publish__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(data, mssg)
    mssg => message_new()
    call message_get(mssg, data)
    ASSERT(associated(data))
    call json_set(data, "update", .true.)
    nullify(data)
    call msgbus_publish(this%msgb, mssg)
    call message_del(mssg)
    nullify(mssg)

    POP_SUB(base_potential__publish__)
  end subroutine base_potential__publish__

  ! ---------------------------------------------------------
  subroutine base_potential_notify_type(this)
    type(base_potential_t), intent(inout) :: this

    type(json_object_t), pointer :: data
    type(message_t),     pointer :: mssg

    PUSH_SUB(base_potential_notify_type)

    ASSERT(associated(this%config))
    nullify(data, mssg)
    mssg => message_new()
    call message_get(mssg, data)
    ASSERT(associated(data))
    call json_set(data, "update", .true.)
    nullify(data)
    call msgbus_notify(this%msgb, mssg)
    call message_del(mssg)
    nullify(mssg)

    POP_SUB(base_potential_notify_type)
  end subroutine base_potential_notify_type

  ! ---------------------------------------------------------
  subroutine base_potential_notify_subs(this, name)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: name

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_notify_subs)

    ASSERT(associated(this%config))
    nullify(subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    ASSERT(associated(subs))
    call base_potential_notify(subs)
    nullify(subs)

    POP_SUB(base_potential_notify_subs)
  end subroutine base_potential_notify_subs

  ! ---------------------------------------------------------
  subroutine base_potential_start(this, sim)
    type(base_potential_t), intent(inout) :: this
    type(simulation_t),     intent(in)    :: sim

    PUSH_SUB(base_potential_start)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    call base_potential__apply__(this, start)

    POP_SUB(base_potential_start)

  contains

    subroutine start(this)
      type(base_potential_t), intent(inout) :: this

      PUSH_SUB(base_potential_start.start)

      call base_potential__start__(this, sim)

      POP_SUB(base_potential_start.start)
    end subroutine start

  end subroutine base_potential_start

  ! ---------------------------------------------------------
  function base_potential_calc(this, that) result(energy)
    type(base_potential_t),         intent(inout) :: this
    type(base_density_t), optional, intent(in)    :: that
    
    real(kind=wp) :: energy

    logical :: fuse

    PUSH_SUB(base_potential_calc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    energy = 0.0_wp
    call base_potential_get(this, use=fuse)
    if(fuse)then
      if(present(that))then
        call base_potential_update(this)
        energy = base_potential__calc__(this, that)
      else
        call base_potential_get(this, energy)
      end if
    end if

    POP_SUB(base_potential_calc)
  end function base_potential_calc

  ! ---------------------------------------------------------
  subroutine base_potential_update(this, energy)
    type(base_potential_t), intent(inout) :: this
    logical,      optional, intent(in)    :: energy

    PUSH_SUB(base_potential_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential__apply__(this, update, enforce_active=.true.)

    POP_SUB(base_potential_update)

  contains

    subroutine update(this)
      type(base_potential_t), intent(inout) :: this

      PUSH_SUB(base_potential_update.update)

      call base_potential__update__(this, energy)

      POP_SUB(base_potential_update.update)
    end subroutine update

  end subroutine base_potential_update

  ! ---------------------------------------------------------
  subroutine base_potential_reset(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_reset)

    call base_potential__apply__(this, base_potential__reset__)

    POP_SUB(base_potential_reset)
  end subroutine base_potential_reset

  ! ---------------------------------------------------------
  recursive subroutine base_potential_stop(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential_stop)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_potential__apply__(this, base_potential__stop__)

    POP_SUB(base_potential_stop)
  end subroutine base_potential_stop

  ! ---------------------------------------------------------
  subroutine base_potential__sets__info(this, name, lock, active)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: name
    logical,      optional, intent(in)    :: lock
    logical,      optional, intent(in)    :: active

    type(base_potential_t), pointer :: subs
    type(msgbus_t),         pointer :: msgb
    logical                         :: actv, accu
    integer                         :: ierr

    PUSH_SUB(base_potential__sets__info)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    nullify(subs, msgb)
    if(present(lock)) continue
    call json_get(this%config, "reduce", accu, ierr)
    if(ierr/=JSON_OK) accu = .false.
    if(accu.and.present(active))then
      call base_potential_gets(this, trim(adjustl(name)), subs, active=actv)
      ASSERT(associated(subs))
      if(actv.neqv.active)then
        call base_potential_notify(this)
        call base_potential_get(subs, msgb)
        ASSERT(associated(msgb))
        if(active)then
          call msgbus_attach(this%msgb, msgb)
        else
          call msgbus_detach(this%msgb, msgb)
        end if
        nullify(subs, msgb)
      end if
    end if

    POP_SUB(base_potential__sets__info)
  end subroutine base_potential__sets__info

  ! ---------------------------------------------------------
  subroutine base_potential__sets__type(this, name, that, config, lock, active)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: name
    type(base_potential_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config
    logical,      optional, intent(in)    :: lock
    logical,      optional, intent(in)    :: active

    type(msgbus_t), pointer :: msgb
    logical                 :: actv, accu
    integer                 :: ierr

    PUSH_SUB(base_potential__sets__type)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(json_len(config)>0)
    nullify(msgb)
    if(present(lock)) continue
    actv = .true.
    if(present(active)) actv = active
    call json_get(this%config, "reduce", accu, ierr)
    if(ierr/=JSON_OK) accu = .false.
    if(actv.and.accu)then
      call base_potential_notify(this)
      call base_potential_get(that, msgb)
      ASSERT(associated(msgb))
      call msgbus_attach(this%msgb, msgb)
      nullify(msgb)
    end if

    POP_SUB(base_potential__sets__type)
  end subroutine base_potential__sets__type

  ! ---------------------------------------------------------
  subroutine base_potential__dels__(this, name, that, config, lock, active)
    type(base_potential_t),        intent(inout) :: this
    character(len=*),              intent(in)    :: name
    type(base_potential_t),        intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config
    logical,             optional, intent(in)    :: lock
    logical,             optional, intent(in)    :: active

    type(msgbus_t), pointer :: msgb
    logical                 :: actv, accu
    integer                 :: ierr

    PUSH_SUB(base_potential__dels__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    nullify(msgb)
    if(present(config)) continue
    if(present(lock))   continue
    actv = .true.
    if(present(active)) actv = active
    call json_get(this%config, "reduce", accu, ierr)
    if(ierr/=JSON_OK) accu = .false.
    if(actv.and.accu)then
      call base_potential_notify(this)
      call base_potential_get(that, msgb)
      ASSERT(associated(msgb))
      call msgbus_detach(this%msgb, msgb)
      nullify(msgb)
    end if

    POP_SUB(base_potential__dels__)
  end subroutine base_potential__dels__

  ! ---------------------------------------------------------
  subroutine base_potential_get_sub_storage(this, name, that)
    type(base_potential_t),   intent(in)  :: this
    character(len=*),         intent(in)  :: name
    type(storage_t), pointer, intent(out) :: that

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_get_sub_storage)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that)
    nullify(subs)

    POP_SUB(base_potential_get_sub_storage)
  end subroutine base_potential_get_sub_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_sub_potential_r1(this, name, that, spin)
    type(base_potential_t),               intent(in)  :: this
    character(len=*),                     intent(in)  :: name
    real(kind=wp), dimension(:), pointer, intent(out) :: that
    integer,                    optional, intent(in)  :: spin

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_get_sub_potential_r1)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that, spin)

    POP_SUB(base_potential_get_sub_potential_r1)
  end subroutine base_potential_get_sub_potential_r1

  ! ---------------------------------------------------------
  subroutine base_potential_get_sub_potential_r2(this, name, that)
    type(base_potential_t),                 intent(in)  :: this
    character(len=*),                       intent(in)  :: name
    real(kind=wp), dimension(:,:), pointer, intent(out) :: that

    type(base_potential_t), pointer :: subs

    PUSH_SUB(base_potential_get_sub_potential_r2)

    nullify(that, subs)
    call base_potential_gets(this, trim(adjustl(name)), subs)
    if(associated(subs)) call base_potential_get(subs, that)

    POP_SUB(base_potential_get_sub_potential_r2)
  end subroutine base_potential_get_sub_potential_r2

  ! ---------------------------------------------------------
  subroutine base_potential_set_info(this, static)
    type(base_potential_t),  intent(inout) :: this
    logical,       optional, intent(in)    :: static

    PUSH_SUB(base_potential_set_info)

    ASSERT(associated(this%config))
    call storage_set(this%data, lock=static)

    POP_SUB(base_potential_set_info)
  end subroutine base_potential_set_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_info(this, kind, size, nspin, static, fine, use)
    type(base_potential_t), intent(in)  :: this
    integer,      optional, intent(out) :: kind
    integer,      optional, intent(out) :: size
    integer,      optional, intent(out) :: nspin
    logical,      optional, intent(out) :: static
    logical,      optional, intent(out) :: fine
    logical,      optional, intent(out) :: use

    PUSH_SUB(base_potential_get_info)

    ASSERT(associated(this%config))
    if(present(kind)) kind = this%kind
    if(present(nspin)) nspin = this%nspin
    call storage_get(this%data, size=size, lock=static, fine=fine, alloc=use)

    POP_SUB(base_potential_get_info)
  end subroutine base_potential_get_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_energy(this, energy)
    type(base_potential_t), intent(inout) :: this
    real(kind=wp),          intent(out)   :: energy

    logical :: fuse
    integer :: ierr

    PUSH_SUB(base_potential_get_energy)

    ASSERT(associated(this%config))
    energy = 0.0_wp
    call base_potential_get(this, use=fuse)
    if(fuse)then
      call base_potential_update(this, energy=.true.)
      call memo_get(this%memo, "energy", energy, ierr)
      ASSERT(ierr==MEMO_OK)
    end if

    POP_SUB(base_potential_get_energy)
  end subroutine base_potential_get_energy

  ! ---------------------------------------------------------
  subroutine base_potential_get_config(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(json_object_t),   pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_potential_get_config)
  end subroutine base_potential_get_config

  ! ---------------------------------------------------------
  subroutine base_potential_get_system(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(base_system_t),   pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_potential_get_system)
  end subroutine base_potential_get_system

  ! ---------------------------------------------------------
  subroutine base_potential_get_density(this, that)
    type(base_potential_t),        intent(in)  :: this
    type(base_density_t), pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_density)

    nullify(that)
    if(associated(this%sys))&
      call base_system_get(this%sys, that)

    POP_SUB(base_potential_get_density)
  end subroutine base_potential_get_density

  ! ---------------------------------------------------------
  subroutine base_potential_get_simulation(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(simulation_t),    pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_potential_get_simulation)
  end subroutine base_potential_get_simulation

  ! ---------------------------------------------------------
  subroutine base_potential_get_msgbus(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(msgbus_t),        pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_msgbus)

    nullify(that)
    if(associated(this%config)) that => this%msgb

    POP_SUB(base_potential_get_msgbus)
  end subroutine base_potential_get_msgbus

  ! ---------------------------------------------------------
  subroutine base_potential_get_storage(this, that)
    type(base_potential_t), target, intent(in)  :: this
    type(storage_t),       pointer, intent(out) :: that

    PUSH_SUB(base_potential_get_storage)

    nullify(that)
    if(associated(this%config)) that => this%data

    POP_SUB(base_potential_get_storage)
  end subroutine base_potential_get_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_r1(this, that, spin)
    type(base_potential_t),               intent(inout) :: this
    real(kind=wp), dimension(:), pointer, intent(out)   :: that
    integer,                    optional, intent(in)    :: spin

    logical :: fuse
    integer :: ispn

    PUSH_SUB(base_potential_get_potential_r1)

    nullify(that)
    ispn = default_nspin
    if(present(spin)) ispn = spin
    call base_potential_get(this, use=fuse)
    if(fuse)then
      call base_potential_update(this)
      call storage_get(this%data, that, ispn)
    end if

    POP_SUB(base_potential_get_potential_r1)
  end subroutine base_potential_get_potential_r1

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_r2(this, that)
    type(base_potential_t),                 intent(inout) :: this
    real(kind=wp), dimension(:,:), pointer, intent(out)   :: that

    logical :: fuse

    PUSH_SUB(base_potential_get_potential_r2)

    nullify(that)
    call base_potential_get(this, use=fuse)
    if(fuse)then
      call base_potential_update(this)
      call storage_get(this%data, that)
    end if

    POP_SUB(base_potential_get_potential_r2)
  end subroutine base_potential_get_potential_r2

  ! ---------------------------------------------------------
  subroutine base_potential__copy__(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_potential__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_potential__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_potential__init__(this, that)
      call refcount_del(this%rcnt)
      call memo_copy(this%memo, that%memo)
      if(associated(that%sim)) call storage_copy(this%data, that%data)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_potential__copy__)
  end subroutine base_potential__copy__

  ! ---------------------------------------------------------
  subroutine base_potential__end__(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(base_potential__end__)

    nullify(this%config, this%sys, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    nullify(this%rcnt)
    this%kind = POTN_KIND_NONE
    this%nspin = default_nspin
    call memo_end(this%memo)
    call storage_end(this%data)
    call msgbus_end(this%msgb)
    ASSERT(base_potential_dict_len(this%dict)==0)
    call base_potential_dict_end(this%dict)

    POP_SUB(base_potential__end__)
  end subroutine base_potential__end__

end module base_potential_oct_m

!! Local Variables:
!! mode: f90
!! End:
